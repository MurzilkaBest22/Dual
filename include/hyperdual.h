#pragma once // Защита от множественного включения

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <map>
#include <functional>
#include "real_utils.h"
#include "rational.h" // Рациональные числа
#include "multidual.h"

namespace dual {

	class HyperDualCombination {
	private:
		std::vector<long double> cells; // Ячейки гипердуальной комбинации
		std::vector<size_t> orders; // Порядки нильпотентных элементов
		std::vector<size_t> coefficients; // Коэффициенты для преобразований между линейным и многомерным индексами

		// Нахождение коэффициентов для преобразований индексов
		void define_coefficients() {
			coefficients.resize(orders.size());
			if (orders.size() == 0) {
				return;
			}
			coefficients.back() = 1;
			for (size_t j = orders.size(); j-- > 1;) {
				coefficients[j - 1] = coefficients[j] * orders[j];
			}
		}

		// Проверка совместимостей двух гипердуальных комбинаций
		void Check(const HyperDualCombination& other) const {
			if (orders != other.orders) {
				throw std::invalid_argument("HyperDualCombination: order vectors must match");
			}
		}
	public:
		// Конструкторы
		HyperDualCombination() = default;
		explicit HyperDualCombination(const std::vector<size_t>& m) : orders(m) {
			define_coefficients();
			size_t M = 1;
			for (size_t order : orders) M *= order;
			cells.assign(M, 0.0L);
		}
		HyperDualCombination(const std::vector<size_t>& m, long double value) : HyperDualCombination(m) {
			cells[0] = value;
		}
		HyperDualCombination(const std::vector<size_t>& m, const rational::Rational& value) : HyperDualCombination(m) {
			cells[0] = value.toLongDouble();
		}
		HyperDualCombination(const std::vector<size_t>& m, const std::vector<long double>& values) : orders(m), cells(values) {
			define_coefficients();
			size_t M = 1;
			for (size_t order : orders) M *= order;
			if (cells.size() != M) {
				throw std::invalid_argument("HyperDualCombination: cells vector size mismatch");
			}
		}
		
		// Возвращение значений
		size_t number_of_orders() const { return orders.size(); }
		const std::vector<size_t>& get_orders() const { return orders; }
		size_t total_size() const { return cells.size(); }

		// Преобразование многомерного индекса в линейный
		size_t linear_index(const std::vector<size_t>& multi_index) const {
			if (multi_index.size() != orders.size()) {
				throw std::invalid_argument("Linear_index: invalid multi_index size");
			}
			size_t result = 0;
			for (size_t j = 0; j < orders.size(); ++j) {
				if (multi_index[j] >= orders[j]) {
					throw std::out_of_range("Multi_index out of range");
				}
				result += multi_index[j] * coefficients[j];
			}
			return result;
		}

		// Преобразование линейного индекса в многомерный
		std::vector<size_t> multi_index(size_t linear_index) const {
			if (linear_index >= cells.size()) {
				throw std::out_of_range("Linear_index out of range");
			}
			std::vector<size_t> result(orders.size());
			size_t g = linear_index;
			for (size_t j = 0; j < orders.size(); ++j) {
				result[j] = g / coefficients[j];
				g %= coefficients[j];
			}
			return result;
		}

		// Доступ к элементам по линейному индексу
		long double& operator[](size_t i) { return cells[i]; }
		const long double& operator[](size_t i) const { return cells[i]; }

		// Доступ к элементам по многомерному индексу
		long double& at(const std::vector<size_t>& multi_index) {
			return cells[linear_index(multi_index)];
		}
		const long double& at(const std::vector<size_t>& multi_index) const {
			return cells[linear_index(multi_index)];
		}

		// Унарный минус
		HyperDualCombination operator-() const {
			HyperDualCombination result(orders);
			for (size_t i = 0; i < cells.size(); ++i) {
				result.cells[i] = -cells[i];
			}
			return result;
		}

		// Составные операторы присваивания
		HyperDualCombination& operator+= (const HyperDualCombination& other) {
			Check(other);
			for (size_t i = 0; i < cells.size(); ++i) {
				cells[i] += other.cells[i];
			}
			return *this;
		}
		HyperDualCombination& operator+= (long double other) {
			cells[0] += other;
			return *this;
		}
		HyperDualCombination& operator+= (const rational::Rational& other) {
			cells[0] += other.toLongDouble();
			return *this;
		}

		HyperDualCombination& operator-= (const HyperDualCombination& other) {
			Check(other);
			for (size_t i = 0; i < cells.size(); ++i) {
				cells[i] -= other.cells[i];
			}
			return *this;
		}
		HyperDualCombination& operator-= (long double other) {
			cells[0] -= other;
			return *this;
		}
		HyperDualCombination& operator-= (const rational::Rational& other) {
			cells[0] -= other.toLongDouble();
			return *this;
		}

		HyperDualCombination& operator*= (const HyperDualCombination& other) {
			*this = *this * other;
			return *this;
		}
		HyperDualCombination& operator*= (long double other) {
			for (size_t i = 0; i < cells.size(); ++i) {
				cells[i] *= other;
			}
			return *this;
		}
		HyperDualCombination& operator*= (const rational::Rational& other) {
			for (size_t i = 0; i < cells.size(); ++i) {
				cells[i] = cells[i] * other;
			}
			return *this;
		}

		HyperDualCombination& operator/= (const HyperDualCombination& other) {
			if (other[0] == 0.0L) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/=");
			}
			*this = *this / other;
			return *this;
		}
		HyperDualCombination& operator/= (long double other) {
			if (other == 0.0L) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/=");
			}
			for (size_t i = 0; i < cells.size(); ++i) {
				cells[i] /= other;
			}
			return *this;
		}
		HyperDualCombination& operator/= (const rational::Rational& other) {
			if (other == rational::Rational(0,1)) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/=");
			}
			for (size_t i = 0; i < cells.size(); ++i) {
				cells[i] = cells[i] / other;
			}
			return *this;
		}

		// Дружественные операторы
		friend HyperDualCombination operator+(const HyperDualCombination& H1, const HyperDualCombination& H2) {
			H1.Check(H2);
			HyperDualCombination result(H1.orders);
			for (size_t i = 0; i < H1.cells.size(); ++i) {
				result.cells[i] = H1.cells[i] + H2.cells[i];
			}
			return result;
		}
		friend HyperDualCombination operator+(long double p, const HyperDualCombination& H) {
			HyperDualCombination result = H;
			result.cells[0] += p;
			return result;
		}
		friend HyperDualCombination operator+(const HyperDualCombination& H, long double p) {
			HyperDualCombination result = H;
			result.cells[0] += p;
			return result;
		}
		friend HyperDualCombination operator+(const rational::Rational r, const HyperDualCombination& H) {
			HyperDualCombination result = H;
			result.cells[0] += r.toLongDouble();
			return result;
		}
		friend HyperDualCombination operator+(const HyperDualCombination& H, const rational::Rational r) {
			HyperDualCombination result = H;
			result.cells[0] += r.toLongDouble();
			return result;
		}

		friend HyperDualCombination operator-(const HyperDualCombination& H1, const HyperDualCombination& H2) {
			H1.Check(H2);
			HyperDualCombination result(H1.orders);
			for (size_t i = 0; i < H1.cells.size(); ++i) {
				result.cells[i] = H1.cells[i] - H2.cells[i];
			}
			return result;
		}
		friend HyperDualCombination operator-(long double p, const HyperDualCombination& H) {
			HyperDualCombination result = -H;
			result.cells[0] += p;
			return result;
		}
		friend HyperDualCombination operator-(const HyperDualCombination& H, long double p) {
			HyperDualCombination result = H;
			result.cells[0] -= p;
			return result;
		}
		friend HyperDualCombination operator-(const rational::Rational r, const HyperDualCombination& H) {
			HyperDualCombination result = -H;
			result.cells[0] += r.toLongDouble();
			return result;
		}
		friend HyperDualCombination operator-(const HyperDualCombination& H, const rational::Rational r) {
			HyperDualCombination result = H;
			result.cells[0] -= r.toLongDouble();
			return result;
		}

		// Вспомогательная функция
		template<typename Func>
		static void all_multi_index(const std::vector<size_t>& k, size_t dimension, std::vector<size_t>& s, Func f) {
			if (dimension == k.size()) {
				f(s);
				return;
			}
			for (size_t i = 0; i <= k[dimension]; ++i) {
				s[dimension] = i;
				all_multi_index(k, dimension + 1, s, f);
			}
		}

		friend HyperDualCombination operator*(const HyperDualCombination& H1, const HyperDualCombination& H2) {
			H1.Check(H2);
			HyperDualCombination result(H1.orders);
			const auto& cells1 = H1.cells;
			const auto& cells2 = H2.cells;

			for (size_t i = 0; i < result.total_size(); ++i) {
				std::vector<size_t> k = H1.multi_index(i);
				long double sum = 0.0L;
				std::vector<size_t> s(H1.number_of_orders(), 0);
				HyperDualCombination::all_multi_index(k, 0, s, [&](const std::vector<size_t>& s) {
					std::vector<size_t> t(H1.number_of_orders());
					bool t_in_range = true;
					for (size_t dim = 0; dim < H1.number_of_orders(); ++dim) {
						t[dim] = k[dim] - s[dim];
						if (t[dim] >= H1.orders[dim]) {
							t_in_range = false;
							break;
						}
					}
					if (t_in_range) {
						size_t index_s = H1.linear_index(s);
						size_t index_t = H1.linear_index(t);
						sum += cells1[index_s] * cells2[index_t];
					}
					});
				result.cells[i] = sum;
			}
			
			return result;
		}
		friend HyperDualCombination operator*(long double p, const HyperDualCombination& H) {
			HyperDualCombination result(H.orders);
			for (size_t i = 0; i < H.cells.size(); ++i) {
				result.cells[i] = p * H.cells[i];
			}
			return result;
		}
		friend HyperDualCombination operator*(const HyperDualCombination& H, long double p) {
			HyperDualCombination result(H.orders);
			for (size_t i = 0; i < H.cells.size(); ++i) {
				result.cells[i] = p * H.cells[i];
			}
			return result;
		}
		friend HyperDualCombination operator*(const rational::Rational& r, const HyperDualCombination& H) {
			HyperDualCombination result(H.orders);
			for (size_t i = 0; i < H.cells.size(); ++i) {
				result.cells[i] = r * H.cells[i];
			}
			return result;
		}
		friend HyperDualCombination operator*(const HyperDualCombination& H, const rational::Rational& r) {
			HyperDualCombination result(H.orders);
			for (size_t i = 0; i < H.cells.size(); ++i) {
				result.cells[i] = r * H.cells[i];
			}
			return result;
		}

		friend HyperDualCombination operator/(const HyperDualCombination& H1, const HyperDualCombination& H2) {
			if (H2[0] == 0.0L) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/");
			}
			H1.Check(H2);
			HyperDualCombination result(H1.orders);
			const auto& cells1 = H1.cells;
			const auto& cells2 = H2.cells;

			for (size_t i = 0; i < result.total_size(); ++i) {
				std::vector<size_t> k = H1.multi_index(i);
				long double sum = 0.0L;
				std::vector<size_t> s(H1.number_of_orders());
				HyperDualCombination::all_multi_index(k, 0, s, [&](const std::vector<size_t>& s) {
					bool s_is_zero = true;
					for (size_t dim = 0; dim < H1.number_of_orders(); ++dim) {
						if (s[dim] != 0) {
							s_is_zero = false;
							break;
						}
					}
					if (s_is_zero) {
						return;
					}
					std::vector<size_t> t(H1.number_of_orders());
					bool t_in_range = true;
					for (size_t dim = 0; dim < H1.number_of_orders(); ++dim) {
						t[dim] = k[dim] - s[dim];
						if (t[dim] >= H1.orders[dim]) {
							t_in_range = false;
							break;
						}
					}
					if (t_in_range) {
						size_t index_s = H1.linear_index(s);
						size_t index_t = H1.linear_index(t);
						if (index_t >= i) {
							throw std::logic_error("Internal error: division order violation");
						}
						sum += cells2[index_s] * result.cells[index_t];
					}
				});
				result.cells[i] = (cells1[i] - sum) / cells2[0];
			}

			return result;
		}
		friend HyperDualCombination operator/(long double p, const HyperDualCombination& H) {
			if (H[0] == 0.0L) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/");
			}
			return HyperDualCombination(H.orders, p) / H;
		}
		friend HyperDualCombination operator/(const HyperDualCombination& H, long double p) {
			if (p == 0.0L) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/");
			}
			HyperDualCombination result(H.orders);
			for (size_t i = 0; i < H.cells.size(); ++i) {
				result.cells[i] = H.cells[i] / p;
			}
			return result;
		}
		friend HyperDualCombination operator/(const rational::Rational& r, const HyperDualCombination& H) {
			if (H[0] == 0.0L) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/");
			}
			return HyperDualCombination(H.orders, r.toLongDouble()) / H;
		}
		friend HyperDualCombination operator/(const HyperDualCombination& H, const rational::Rational& r) {
			if (r == rational::Rational(0,1)) {
				throw std::domain_error("Division by zero in HyperDualCombination::operator/");
			}
			HyperDualCombination result(H.orders);
			for (size_t i = 0; i < H.cells.size(); ++i) {
				result.cells[i] = H.cells[i] / r;
			}
			return result;
		}
		
		// Максимальная сумма показателей степеней нильпотентных элементов, при которых степени не обнуляются
		size_t max_degree_sum() const {
			size_t sum = 0;
			for (size_t order : orders) {
				sum += order - 1;
			}
			return sum;
		}
	};

	// Прибавление нильпотентного элемента с индексом nilpotent_index, приводящее к образованию гипердуальной комбинации с порядками nilpotent_orders
	static HyperDualCombination Nilpotent_Add(long double value, const std::vector<size_t>& nilpotent_orders, size_t nilpotent_index) {
		HyperDualCombination result(nilpotent_orders);
		result[0] = value;
		std::vector<size_t> multi(nilpotent_orders.size(), 0);
		multi[nilpotent_index] = 1;
		size_t index = result.linear_index(multi);
		result[index] = 1.0L;
		return result;
	}
	static HyperDualCombination Nilpotent_Add(const rational::Rational& value, const std::vector<size_t>& nilpotent_orders, size_t nilpotent_index) {
		HyperDualCombination result(nilpotent_orders);
		result[0] = value.toLongDouble();
		std::vector<size_t> multi(nilpotent_orders.size(), 0);
		multi[nilpotent_index] = 1;
		size_t index = result.linear_index(multi);
		result[index] = 1.0L;
		return result;
	}

	// Вспомогательная функция
	HyperDualCombination simple_pow(const HyperDualCombination& H, size_t k) {
		HyperDualCombination result(H.get_orders(), 1.0L);
		for (size_t i = 0; i < k; ++i) {
			result *= H;
		}
		return result;
	}

	// Экспонента
	HyperDualCombination exp(const HyperDualCombination& H) {
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += std::exp(H[0]) * simple_pow(H - H[0], k) / factorial(k);
		}
		return result;
	}

	// Логарифмы
	HyperDualCombination log(long double p, const HyperDualCombination& H) {
		// Проверка ОДЗ логарифма
		if ((p <= 0.0) || (p == 1.0) || (H[0] <= 0.0)) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); k++) {
			if (k == 0) {
				result += ::log(p, H[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial(k - 1) * std::pow(H[0], -1 * k) * simple_pow(H - H[0], k) / (::ln(p) * factorial(k));
				}
				else {
					result += -1 * factorial(k - 1) * std::pow(H[0], -1 * k) * simple_pow(H - H[0], k) / (::ln(p) * factorial(k));
				}
			}
		}
		return result;
	}
	HyperDualCombination log(const rational::Rational& r, const HyperDualCombination& H) {
		// Проверка ОДЗ логарифма
		if ((r <= rational::Rational(0, 1)) || (r == rational::Rational(1, 1)) || (H[0] <= 0.0)) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); k++) {
			if (k == 0) {
				result += ::log(r.toLongDouble(), H[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial(k - 1) * std::pow(H[0], -1 * k) * simple_pow(H - H[0], k) / (ln(r) * factorial(k));
				}
				else {
					result += -1 * factorial(k - 1) * std::pow(H[0], -1 * k) * simple_pow(H - H[0], k) / (ln(r) * factorial(k));
				}
			}
		}
		return result;
	}
	HyperDualCombination ln(const HyperDualCombination& H) {
		// Проверка ОДЗ логарифма
		if (H[0] <= 0.0) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); k++) {
			if (k == 0) {
				result += ::ln(H[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial(k - 1) * std::pow(H[0], -1 * k) * simple_pow(H - H[0], k) / factorial(k);
				}
				else {
					result += -1 * factorial(k - 1) * std::pow(H[0], -1 * k) * simple_pow(H - H[0], k) / factorial(k);
				}
			}
		}
		return result;
	}
	HyperDualCombination log(const HyperDualCombination& H) {
		// Проверка ОДЗ логарифма
		if (H[0] <= 0.0) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); k++) {
			if (k == 0) {
				result += ::ln(H[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial(k - 1) * std::pow(H[0], -1 * k) * simple_pow(H - H[0], k) / factorial(k);
				}
				else {
					result += -1 * factorial(k - 1) * std::pow(H[0], -1 * k) * simple_pow(H - H[0], k) / factorial(k);
				}
			}
		}
		return result;
	}
	HyperDualCombination lg(const HyperDualCombination& H) {
		// Проверка ОДЗ логарифма
		if (H[0] <= 0.0) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); k++) {
			if (k == 0) {
				result += ::lg(H[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial(k - 1) * std::pow(H[0], -1 * k) * simple_pow(H - H[0], k) / (::ln(10) * factorial(k));
				}
				else {
					result += -1 * factorial(k - 1) * std::pow(H[0], -1 * k) * simple_pow(H - H[0], k) / (::ln(10) * factorial(k));
				}
			}
		}
		return result;
	}
	HyperDualCombination lb(const HyperDualCombination& H) {
		// Проверка ОДЗ логарифма
		if (H[0] <= 0.0) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); k++) {
			if (k == 0) {
				result += ::lb(H[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial(k - 1) * std::pow(H[0], -1 * k) * simple_pow(H - H[0], k) / (::ln(2) * factorial(k));
				}
				else {
					result += -1 * factorial(k - 1) * std::pow(H[0], -1 * k) * simple_pow(H - H[0], k) / (::ln(2) * factorial(k));
				}
			}
		}
		return result;
	}
	HyperDualCombination log(const HyperDualCombination& H1, const HyperDualCombination& H2) {
		if ((H1[0] <= 0.0) || (H1[0] == 1) || (H2[0] <= 0.0)) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		return ln(H2) / ln(H1);
	}
	HyperDualCombination log(const HyperDualCombination& H, long double p) {
		if ((H[0] <= 0.0) || (H[0] == 1) || (p <= 0.0)) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		return ::ln(p) / ln(H);
	}
	HyperDualCombination log(const HyperDualCombination& H, const rational::Rational& r) {
		if ((H[0] <= 0.0) || (H[0] == 1) || (r <= rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the range of allowable values of the HyperDualCombination logarithm");
		}
		return ln(r) / ln(H);
	}

	// Возведение в степень
	HyperDualCombination sq(const HyperDualCombination& H) {
		return H * H;
	}
	HyperDualCombination cb(const HyperDualCombination& H) {
		return H * H * H;
	}
	HyperDualCombination pow(const HyperDualCombination& H, int p) {
		if ((H[0] == 0.0) && (p <= 0)) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		if (p >= 0) {
			HyperDualCombination result(H.get_orders());
			HyperDualCombination term(H.get_orders());
			for (size_t k = 0; k < H.max_degree_sum(); ++k) {
				if (k <= p) {
					term = std::pow(H[0], p - k) * simple_pow(H - H[0], k) / factorial(k);
					for (size_t i = p; i >= p - k + 1; --i) {
						term *= i;
					}
					result += term;
				}
			}
			return result;
		}
		else {
			return HyperDualCombination(H.get_orders(), 1.0L) / pow(H, -p);
		}
	}
	HyperDualCombination pow(const HyperDualCombination& H, int64_t p) {
		if ((H[0] == 0.0) && (p <= 0)) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		if (p >= 0) {
			HyperDualCombination result(H.get_orders());
			HyperDualCombination term(H.get_orders());
			for (size_t k = 0; k < H.max_degree_sum(); ++k) {
				if (k <= p) {
					term = std::pow(H[0], p - k) * simple_pow(H - H[0], k) / factorial(k);
					for (size_t i = p; i >= p - k + 1; --i) {
						term *= i;
					}
					result += term;
				}
			}
			return result;
		}
		else {
			return HyperDualCombination(H.get_orders(), 1.0L) / pow(H, -p);
		}
	}
	HyperDualCombination pow(const HyperDualCombination& H, long double p) {
		if (H[0] <= 0.0) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		return exp(p * ln(H));
	}
	HyperDualCombination pow(const HyperDualCombination& H, const rational::Rational& r) {
		if (H[0] == 0.0) {
			if (r <= rational::Rational(0, 1)) {
				throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
			}
			return HyperDualCombination(H.get_orders());
		}
		if (H[0] < 0.0) {
			if (r.denominatorIsEven()) {
				throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
			}
			HyperDualCombination result = pow(-H, r.toLongDouble());
			if (r.numeratorIsOdd()) {
				result *= -1;
			}
			return result;
		}
		else {
			return pow(H, r.toLongDouble());
		}
	}
	HyperDualCombination pow(const HyperDualCombination& H1, const rational::Rational& H2_Real, std::vector<long double> H2_Nilpotent, std::vector<size_t> Nilpotent_orders) {
		if ((H1[0] == 0.0) && (H2_Real <= rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		if ((H1[0] < 0.0) && H2_Real.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		std::vector<long double> coeffs(H2_Nilpotent.size() + 1);
		for (size_t i = 0; i < coeffs.size(); ++i) {
			if (i == 0) {
				coeffs[i] = H2_Real.toLongDouble();
			}
			else {
				coeffs[i] = H2_Nilpotent[i - 1];
			}
		}
		return exp(HyperDualCombination(Nilpotent_orders, coeffs) * ln(H1));
	}
	HyperDualCombination pow(long double p, const HyperDualCombination& H) {
		if (p <= 0.0) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		return exp(::ln(p) * H);
	}
	HyperDualCombination pow(const HyperDualCombination& H1, const HyperDualCombination& H2) {
		if (H1[0] <= 0.0) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		return exp(ln(H1) * H2);
	}
	HyperDualCombination pow(const rational::Rational& r, const HyperDualCombination& H) {
		if (r <= rational::Rational(0, 1)) {
			throw std::domain_error("Exceeding the domain of allowable values of a HyperDualCombination power function");
		}
		return exp(ln(r) * H);
	}

	// Извлечение корня
	HyperDualCombination sqrt(const HyperDualCombination& H) {
		if (H[0] < 0) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination root function");
		}
		if (H[0] == 0.0) {
			if (H.total_size() > 1) {
				for (size_t i = 1; i < H.total_size(); ++i) {
					if (H[i] != 0.0) {
						throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
					}
				}
			}
			return H;
		}
		HyperDualCombination result(H.get_orders());
		HyperDualCombination term(H.get_orders());
		for (size_t k = 0; k < H.total_size(); ++k) {
			term = std::pow(H[0], 0.5L - k) * simple_pow(H - H[0], k) / factorial(k);
			if (k > 0) {
				for (size_t i = 0; i <= k - 1; ++i) {
					term *= 0.5L - i;
				}
			}
			result += term;
		}
		return result;
	}
	HyperDualCombination cbrt(const HyperDualCombination& H) {
		HyperDualCombination result(H.get_orders());
		HyperDualCombination term(H.get_orders());
		if (H[0] == 0.0) {
			if (H.total_size() > 1) {
				for (size_t i = 1; i < H.total_size(); ++i) {
					if (H[i] != 0.0) {
						throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
					}
				}
			}
			return H;
		}
		for (size_t k = 0; k < H.total_size(); ++k) {
			term = std::pow(H[0], 1.0L / 3.0L - k) * simple_pow(H - H[0], k) / factorial(k);
			if (k > 0) {
				for (size_t i = 0; i <= k - 1; ++i) {
					term *= 1.0L / 3.0L - i;
				}
			}
			result += term;
		}
		return result;
	}
	HyperDualCombination rt(const HyperDualCombination& H, const rational::Rational& r) {
		if ((r == rational::Rational(0, 1)) || ((H[0] < 0) && r.numeratorIsEven())) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if (H[0] == 0.0) {
			if (r < rational::Rational(0, 1)) {
				throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
			}
			if (r == rational::Rational(1, 1)) {
				return H;
			}
			else {
				if (H.total_size() > 1) {
					for (size_t i = 1; i < H.total_size(); ++i) {
						if (H[i] != 0.0) {
							throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
						}
					}
				}
				return H;
			}
		}
		if (r >= rational::Rational(0, 1)) {
			HyperDualCombination result(H.get_orders());
			HyperDualCombination term(H.get_orders());
			for (size_t k = 0; k < H.total_size(); ++k) {
				term = std::pow(H[0], 1 / r.toLongDouble() - k) * simple_pow(H - H[0], k) / factorial(k);
				if (k > 0) {
					for (size_t i = 0; i <= k - 1; ++i) {
						term *= 1 / r.toLongDouble() - i;
					}
				}
				result += term;
			}
			return result;
		}
		else {
			return 1.0 / rt(H, -r);
		}
	}
	HyperDualCombination rt(const HyperDualCombination& H, long double p) {
		if ((p == 0.0) || (H[0] < 0)) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if (H[0] == 0.0) {
			if (p < 0.0) {
				throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
			}
			if (p == 1.0) {
				return H;
			}
			else {
				if (H.total_size() > 1) {
					for (size_t i = 1; i < H.total_size(); ++i) {
						if (H[i] != 0.0) {
							throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
						}
					}
				}
				return H;
			}
		}
		if (p >= 0) {
			HyperDualCombination result(H.get_orders());
			HyperDualCombination term(H.get_orders());
			for (size_t k = 0; k < H.total_size(); ++k) {
				term = std::pow(H[0], 1 / p - k) * simple_pow(H - H[0], k) / factorial(k);
				if (k > 0) {
					for (size_t i = 0; i <= k - 1; ++i) {
						term *= 1 / p - i;
					}
				}
				result += term;
			}
			return result;
		}
		else {
			return 1.0 / rt(H, -p);
		}
	}
	HyperDualCombination rt(const HyperDualCombination& H1, const HyperDualCombination& H2) {
		if ((H2[0] == 0.0) || (H1[0] < 0)) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		if ((H1[0] == 0) && (H2[0] < 0.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		return pow(H1, 1 / H2);
	}
	HyperDualCombination rt(long double p, const HyperDualCombination& H) {
		if ((H[0] == 0.0) || (p < 0)) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		if ((p == 0) && (H[0] < 0.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		return pow(p, 1 / H);
	}
	HyperDualCombination rt(const HyperDualCombination& H1, const rational::Rational& H2_Real, std::vector<long double> H2_Nilpotent, std::vector<size_t> Nilpotent_orders) {
		if (H2_Real == rational::Rational(0, 1)) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		if ((H1[0] == 0) && (H2_Real < rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		if ((H1[0] < 0) && H2_Real.numeratorIsEven()) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		std::vector<long double> coeffs(H2_Nilpotent.size() + 1);
		for (size_t i = 0; i < coeffs.size(); ++i) {
			if (i == 0) {
				coeffs[i] = H2_Real.toLongDouble();
			}
			else {
				coeffs[i] = H2_Nilpotent[i - 1];
			}
		}
		return exp(ln(H1) / HyperDualCombination(Nilpotent_orders, coeffs));
	}

	// Тригонометрические функции
	HyperDualCombination sin(const HyperDualCombination& H) {
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += std::sin(H[0] + ::pi() * k / 2) * simple_pow(H - H[0], k) / factorial(k);
		}
		return result;
	}
	HyperDualCombination cos(const HyperDualCombination& H) {
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += std::cos(H[0] + ::pi() * k / 2) * simple_pow(H - H[0], k) / factorial(k);
		}
		return result;
	}
	HyperDualCombination tg(const HyperDualCombination& H) {
		if (std::cos(H[0]) == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination tangent function");
		}
		return sin(H) / cos(H);
	}
	HyperDualCombination ctg(const HyperDualCombination& H) {
		if (std::sin(H[0]) == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination cotangent function");
		}
		return cos(H) / sin(H);
	}
	HyperDualCombination sec(const HyperDualCombination& H) {
		if (std::cos(H[0]) == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination secant function");
		}
		return 1 / cos(H);
	}
	HyperDualCombination cosec(const HyperDualCombination& H) {
		if (std::sin(H[0]) == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination cosecant function");
		}
		return 1 / sin(H);
	}

	// Обратные тригонометрические функции
	HyperDualCombination arcsin(const HyperDualCombination& H) {
		if ((H[0] < -1.0) || (H[0] > 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination arcsine function");
		}
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += arcsin_derivative(H[0], k) * simple_pow(H - H[0], k) / factorial(k);
		}
		return result;
	}
	HyperDualCombination arccos(const HyperDualCombination& H) {
		if ((H[0] < -1.0) || (H[0] > 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination arccosine function");
		}
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += arccos_derivative(H[0], k) * simple_pow(H - H[0], k) / factorial(k);
		}
		return result;
	}
	HyperDualCombination arctg(const HyperDualCombination& H) {
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += arctg_derivative(H[0], k) * simple_pow(H - H[0], k) / factorial(k);
		}
		return result;
	}
	HyperDualCombination arcctg(const HyperDualCombination& H) {
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			result += arcctg_derivative(H[0], k) * simple_pow(H - H[0], k) / factorial(k);
		}
		return result;
	}
	HyperDualCombination arcsec(const HyperDualCombination& H) {
		if ((H[0] > -1.0) && (H[0] < 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination arcsecant function");
		}
		return arccos(1 / H);
	}
	HyperDualCombination arccosec(const HyperDualCombination& H) {
		if ((H[0] > -1.0) && (H[0] < 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination arccosecant function");
		}
		return arcsin(1 / H);
	}

	// Гиперболические функции
	HyperDualCombination sh(const HyperDualCombination& H) {
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			if (k % 2 == 0) {
				result += ::sh(H[0]) * simple_pow(H - H[0], k) / factorial(k);
			}
			else {
				result += ::ch(H[0]) * simple_pow(H - H[0], k) / factorial(k);
			}
		}
		return result;
	}
	HyperDualCombination ch(const HyperDualCombination& H) {
		HyperDualCombination result(H.get_orders());
		for (size_t k = 0; k < H.max_degree_sum(); ++k) {
			if (k % 2 == 0) {
				result += ::ch(H[0]) * simple_pow(H - H[0], k) / factorial(k);
			}
			else {
				result += ::sh(H[0]) * simple_pow(H - H[0], k) / factorial(k);
			}
		}
		return result;
	}
	HyperDualCombination th(const HyperDualCombination& H) {
		return sh(H) / ch(H);
	}
	HyperDualCombination cth(const HyperDualCombination& H) {
		if (H[0] == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination hyperbolic cotangent function");
		}
		return ch(H) / sh(H);
	}
	HyperDualCombination sch(const HyperDualCombination& H) {
		return 1 / ch(H);
	}
	HyperDualCombination csch(const HyperDualCombination& H) {
		if (H[0] == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination hyperbolic cosecant function");
		}
		return 1 / sh(H);
	}

	// Смешанные производные
	template<typename Func>
	long double D(Func f, const std::vector<long double>& X, const std::vector<size_t>& derivatives_orders) {
		size_t n = X.size(); // Количество аргументов функции
		if (derivatives_orders.size() != n) {
			throw std::invalid_argument("D: derivatives_orders size must match X size");
		}

		std::vector<size_t> nilpotent_orders(n);
		for (size_t i = 0; i < n; ++i) {
			nilpotent_orders[i] = derivatives_orders[i] + 1;
		}

		std::vector<HyperDualCombination> X_hyperdual(n);
		for (size_t i = 0; i < n; ++i) {
			if (derivatives_orders[i] > 0) {
				X_hyperdual[i] = Nilpotent_Add(X[i], nilpotent_orders, i);
			}
			else {
				X_hyperdual[i] = HyperDualCombination(nilpotent_orders, X[i]);
			}
		}

		HyperDualCombination F = f(X_hyperdual);

		std::vector<size_t> multi_index(n);
		for (size_t i = 0; i < n; ++i) {
			multi_index[i] = derivatives_orders[i];
		}
		size_t linear_index = F.linear_index(multi_index);
		long double result = F[linear_index];
		
		for (size_t derivative_order : derivatives_orders) {
			result *= factorial(derivative_order);
		}

		return result;
	}

	// Гессиан
	template<typename Func>
	std::vector<std::vector<long double>> Hess(Func f, const std::vector<long double>& X) {
		size_t n = X.size();
		std::vector<size_t> nilpotent_orders(n, 3);
		std::vector<HyperDualCombination> X_hyperdual(n);
		for (size_t i = 0; i < n; ++i) {
			X_hyperdual[i] = Nilpotent_Add(X[i], nilpotent_orders, i);
		}
		HyperDualCombination F = f(X_hyperdual);
		std::vector<std::vector<long double>> H(n, std::vector<long double>(n));
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = i; j < n; ++j) {
				std::vector<size_t> multi_index(n, 0);
				if (i == j) {
					multi_index[i] = 2;
				}
				else {
					multi_index[i] = 1;
					multi_index[j] = 1;
				}
				size_t linear_index = F.linear_index(multi_index);
				long double value = F[linear_index];
				if (i == j) {
					value *= 2; // Учёт факториала для вычисления производных второго порядка
				}
				H[i][j] = value;
				H[j][i] = value;
			}
		}
		return H;
	}

	// Градиент порядка n, являющийся тензором этого порядка
	template<typename Func>
	std::map<std::vector<size_t>, long double> grad(Func f, const std::vector<long double>& X, size_t n) {
		size_t m = X.size();
		std::vector<size_t> nilpotent_orders(m, n + 1);
		std::vector<HyperDualCombination> X_hyperdual(m);
		for (size_t i = 0; i < m; ++i) {
			X_hyperdual[i] = Nilpotent_Add(X[i], nilpotent_orders, i);
		}
		HyperDualCombination F = f(X_hyperdual);
		std::map<std::vector<size_t>, long double> T; // Выходной тензор
		std::vector<size_t> multi_index(m);
		std::function<void(size_t, size_t, size_t)> rec = [&](size_t index, size_t start, size_t sum) {
			if (index == m) {
				if (sum == n) {
					size_t linear_index = F.linear_index(multi_index);
					T[multi_index] = F[linear_index];
					for (size_t i = 0; i < m; ++i) {
						T[multi_index] *= factorial(multi_index[i]);
					}
				}
				return;
			}
			size_t maximum = n - sum;
			for (size_t v = start; v <= maximum; ++v) {
				multi_index[index] = v;
				rec(index + 1, v, sum + v);
			}
		};
		rec(0, 0, 0);
		return T;
	}

}