#pragma once // Защита от множественного включения

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include "real_utils.h"
#include "rational.h" // Рациональные числа

namespace dual {

	// Класс дуальных комбинаций (порядок определяется автоматически конструктором)
	class DualCombination {
	public:
		std::vector<long double> cells; // Ячейки 

		// Конструкторы
		explicit DualCombination(size_t order = 1) : cells(order, 0.0L) {}
		explicit DualCombination(long double value, size_t order = 1) : cells(order, 0.0L) { cells[0] = value; }
		explicit DualCombination(const rational::Rational& value, size_t order = 1) : cells(order, 0.0L) { cells[0] = value.toLongDouble(); }
		DualCombination(const std::vector<long double>& values) : cells(values) {}

		// Возвращение порядка дуальной комбинации при заданных ячейках
		size_t order() const { return cells.size(); }

		// Обращение к ячейкам
		long double& operator[](size_t j) { return cells[j]; }
		const long double& operator[](size_t j) const { return cells[j]; }

		// Составные операторы присваивания
		DualCombination& operator+= (const DualCombination& other) {
			if (order() != other.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator+=");
			}
			for (size_t i = 0; i < order(); ++i) cells[i] += other.cells[i];
			return *this;
		}
		DualCombination& operator+= (long double other) {
			cells[0] += other;
			return *this;
		}
		DualCombination& operator+= (const rational::Rational& other) {
			cells[0] += other.toLongDouble();
			return *this;
		}

		DualCombination& operator-= (const DualCombination& other) {
			if (order() != other.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator-=");
			}
			for (size_t i = 0; i < order(); ++i) cells[i] -= other.cells[i];
			return *this;
		}
		DualCombination& operator-= (long double other) {
			cells[0] -= other;
			return *this;
		}
		DualCombination& operator-= (const rational::Rational& other) {
			cells[0] -= other.toLongDouble();
			return *this;
		}

		DualCombination& operator*= (const DualCombination& other) {
			if (order() != other.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator*=");
			}
			*this = *this * other;
			return *this;
		}
		DualCombination& operator*= (long double other) {
			for (auto& cell : cells) {
				cell *= other;
			}
			return *this;
		}
		DualCombination& operator*= (const rational::Rational& other) {
			for (auto& cell : cells) {
				cell *= other.toLongDouble();
			}
			return *this;
		}

		DualCombination& operator/= (const DualCombination& other) {
			if (order() != other.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator/=");
			}
			*this = *this / other;
			return *this;
		}
		DualCombination& operator/= (long double other) {
			if (other == 0.0L) {
				throw std::domain_error("Division by zero in DualCombination::operator/=");
			}
			for (auto& cell : cells) {
				cell /= other;
			}
			return *this;
		}
		DualCombination& operator/= (const rational::Rational& other) {
			if (other == rational::Rational(0, 1)) {
				throw std::domain_error("Division by zero in DualCombination::operator/=");
			}
			for (auto& cell : cells) {
				cell /= other.toLongDouble();
			}
			return *this;
		}

		// Дружественные операторы
		friend DualCombination operator+(const DualCombination& D1, const DualCombination& D2) {
			if (D1.order() != D2.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator+");
			}
			DualCombination Sum(D1.order());
			for (size_t i = 0; i < D1.order(); ++i) {
				Sum[i] = D1[i] + D2[i];
			}
			return Sum;
		}
		friend DualCombination operator+(long double p, const DualCombination& D) {
			DualCombination Sum = D;
			Sum[0] += p;
			return Sum;
		}
		friend DualCombination operator+(const DualCombination& D, long double p) {
			DualCombination Sum = D;
			Sum[0] += p;
			return Sum;
		}
		friend DualCombination operator+(const rational::Rational r, const DualCombination& D) {
			DualCombination Sum = D;
			Sum[0] += r.toLongDouble();
			return Sum;
		}
		friend DualCombination operator+(const DualCombination& D, const rational::Rational r) {
			DualCombination Sum = D;
			Sum[0] += r.toLongDouble();
			return Sum;
		}

		friend DualCombination operator-(const DualCombination& D1, const DualCombination& D2) {
			if (D1.order() != D2.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator-");
			}
			DualCombination Sum(D1.order());
			for (size_t i = 0; i < D1.order(); ++i) {
				Sum[i] = D1[i] - D2[i];
			}
			return Sum;
		}
		friend DualCombination operator-(long double p, const DualCombination& D) {
			DualCombination Sum = D;
			for (size_t i = 0; i < Sum.order(); ++i) {
				Sum.cells[i] *= -1;
			}
			Sum[0] += p;
			return Sum;
		}
		friend DualCombination operator-(const DualCombination& D, long double p) {
			DualCombination Sum = D;
			Sum[0] -= p;
			return Sum;
		}
		friend DualCombination operator-(const rational::Rational r, const DualCombination& D) {
			DualCombination Sum = D;
			for (size_t i = 0; i < Sum.order(); ++i) {
				Sum.cells[i] *= -1;
			}
			Sum[0] += r.toLongDouble();
			return Sum;
		}
		friend DualCombination operator-(const DualCombination& D, const rational::Rational r) {
			DualCombination Sum = D;
			Sum[0] -= r.toLongDouble();
			return Sum;
		}

		friend DualCombination operator*(const DualCombination& D1, const DualCombination& D2) {
			if (D1.order() != D2.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator*");
			}
			DualCombination Product(D1.order());
			for (size_t k = 0; k < Product.order(); ++k) {
				for (size_t s = 0; s <= k; ++s) {
					Product[k] += D1.cells[s] * D2.cells[k - s];
				}
			}
			return Product;
		}
		friend DualCombination operator*(long double p, const DualCombination& D) {
			DualCombination Product(D.order());
			for (size_t i = 0; i < D.order(); ++i) {
				Product[i] = D[i] * p;
			}
			return Product;
		}
		friend DualCombination operator*(const DualCombination& D, long double p) {
			DualCombination Product(D.order());
			for (size_t i = 0; i < D.order(); ++i) {
				Product[i] = D[i] * p;
			}
			return Product;
		}
		friend DualCombination operator*(const rational::Rational& r, const DualCombination& D) {
			DualCombination Product(D.order());
			for (size_t i = 0; i < D.order(); ++i) {
				Product[i] = D[i] * r;
			}
			return Product;
		}
		friend DualCombination operator*(const DualCombination& D, const rational::Rational& r) {
			DualCombination Product(D.order());
			for (size_t i = 0; i < D.order(); ++i) {
				Product[i] = D[i] * r;
			}
			return Product;
		}

		friend DualCombination operator/(const DualCombination& D1, const DualCombination& D2) {
			if (D1.order() != D2.order()) {
				throw std::invalid_argument("DualCombination order mismatch in DualCombination::operator/");
			}
			if (D2[0] == 0.0L) {
				throw std::domain_error("Division by zero in DualCombination::operator/");
			}
			DualCombination Quotient(D1.order());
			Quotient[0] = D1[0] / D2[0];
			for (size_t i = 1; i < Quotient.order(); ++i) {
				long double Sum = 0.0L;
				for (size_t j = 1; j <= i; ++j) {
					Sum += D2.cells[j] * Quotient[i - j];
				}
				Quotient[i] = (D1.cells[i] - Sum) / D2.cells[0];
			}
			return Quotient;
		}
		friend DualCombination operator/(long double p, const DualCombination& D) {
			if (D[0] == 0.0L) {
				throw std::domain_error("Division by zero in DualCombination::operator/");
			}
			DualCombination Quotient(D.order());
			Quotient[0] = p / D[0];
			for (size_t i = 1; i < Quotient.order(); ++i) {
				long double Sum = 0.0L;
				for (size_t j = 1; j <= i; ++j) {
					Sum += D.cells[j] * Quotient[i - j];
				}
				Quotient[i] = -Sum / D.cells[0];
			}
			return Quotient;
		}
		friend DualCombination operator/(const DualCombination& D, long double p) {
			if (p == 0.0L) {
				throw std::domain_error("Division by zero in DualCombination::operator/");
			}
			DualCombination Quotient(D.order());
			for (size_t i = 0; i < D.order(); ++i) {
				Quotient[i] = D[i] / p;
			}
			return Quotient;
		}
		friend DualCombination operator/(const rational::Rational& r, const DualCombination& D) {
			if (D[0] == 0.0L) {
				throw std::domain_error("Division by zero in DualCombination::operator/");
			}
			DualCombination Quotient(D.order());
			Quotient[0] = r / D[0];
			for (size_t i = 1; i < Quotient.order(); ++i) {
				long double Sum = 0.0L;
				for (size_t j = 1; j <= i; ++j) {
					Sum += D.cells[j] * Quotient[i - j];
				}
				Quotient[i] = -Sum / D.cells[0];
			}
			return Quotient;
		}
		friend DualCombination operator/(const DualCombination& D, const rational::Rational& r) {
			if (r == rational::Rational(0,1)) {
				throw std::domain_error("Division by zero in DualCombination::operator/");
			}
			DualCombination Quotient(D.order());
			for (size_t i = 0; i < D.order(); ++i) {
				Quotient[i] = D[i] / r;
			}
			return Quotient;
		}

		// Унарный минус
		DualCombination operator-() const {
			DualCombination NewComb(order());
			for (size_t i = 0; i < order(); ++i) {
				NewComb[i] = -cells[i];
			}
			return NewComb;
		}
	};

	// Функция для прибавления к вещественному числу нильпотентного элемента заданного порядка n в первой степени
	DualCombination Nilpotent_Add(long double value, size_t n) {
		DualCombination result(value, n);
		if (n >= 2) {
			result[1] = 1.0L;
		}
		return result;
	}
	DualCombination Nilpotent_Add(const rational::Rational& value, size_t n) {
		DualCombination result(value, n);
		if (n >= 2) {
			result[1] = 1.0L;
		}
		return result;
	}

	// Вспомогательная функция
	DualCombination simple_pow(const DualCombination& D, size_t k) {
		if (D[0] != 0.0L) {
			throw std::domain_error("Exceeding the domain of allowable values of a multidual simple power function");
		}
		DualCombination result(1.0L, D.order());
		for (size_t i = 0; i < k; ++i) {
			result *= D;
		}
		return result;
	}

	// Факториал
	long double factorial(size_t k) {
		if ((k == 0) || (k == 1)) {
			return 1.0L;
		}
		else {
			long double result = 1.0L;
			for (size_t i = 2; i <= k; ++i) {
				result *= i;
			}
			return result;
		}
	}

	// Экспонента
	DualCombination exp(const DualCombination& D) {
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += std::exp(D[0]) * simple_pow(D - D[0], k) / factorial(k);
		}
		return result;
	}

	// Логарифмы
	DualCombination log(long double p, const DualCombination& D) {
		// Проверка ОДЗ логарифма
		if ((p <= 0.0) || (p == 1.0) || (D[0] <= 0.0)) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); k++) {
			if (k == 0) {
				result += ::log(p, D[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial(k - 1) * std::pow(D[0], -1 * k) * simple_pow(D - D[0], k) / (::ln(p) * factorial(k));
				}
				else {
					result += -1 * factorial(k - 1) * std::pow(D[0], -1 * k) * simple_pow(D - D[0], k) / (::ln(p) * factorial(k));
				}
			}
		}
		return result;
	}
	DualCombination log(const rational::Rational& r, const DualCombination& D) {
		// Проверка ОДЗ логарифма
		if ((r <= rational::Rational(0, 1)) || (r == rational::Rational(1, 1)) || (D[0] <= 0.0)) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); k++) {
			if (k == 0) {
				result += ::log(r.toLongDouble(), D[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial(k - 1) * std::pow(D[0], -1 * k) * simple_pow(D - D[0], k) / (ln(r) * factorial(k));
				}
				else {
					result += -1 * factorial(k - 1) * std::pow(D[0], -1 * k) * simple_pow(D - D[0], k) / (ln(r) * factorial(k));
				}
			}
		}
		return result;
	}
	DualCombination ln(const DualCombination& D) {
		// Проверка ОДЗ логарифма
		if (D[0] <= 0.0) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); k++) {
			if (k == 0) {
				result += ::ln(D[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial(k - 1) * std::pow(D[0], -1 * k) * simple_pow(D - D[0], k) / factorial(k);
				}
				else {
					result += -1 * factorial(k - 1) * std::pow(D[0], -1 * k) * simple_pow(D - D[0], k) / factorial(k);
				}
			}
		}
		return result;
	}
	DualCombination log(const DualCombination& D) {
		// Проверка ОДЗ логарифма
		if (D[0] <= 0.0) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); k++) {
			if (k == 0) {
				result += ::ln(D[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial(k - 1) * std::pow(D[0], -1 * k) * simple_pow(D - D[0], k) / factorial(k);
				}
				else {
					result += -1 * factorial(k - 1) * std::pow(D[0], -1 * k) * simple_pow(D - D[0], k) / factorial(k);
				}
			}
		}
		return result;
	}
	DualCombination lg(const DualCombination& D) {
		// Проверка ОДЗ логарифма
		if (D[0] <= 0.0) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); k++) {
			if (k == 0) {
				result += ::lg(D[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial(k - 1) * std::pow(D[0], -1 * k) * simple_pow(D - D[0], k) / (::ln(10) * factorial(k));
				}
				else {
					result += -1 * factorial(k - 1) * std::pow(D[0], -1 * k) * simple_pow(D - D[0], k) / (::ln(10) * factorial(k));
				}
			}
		}
		return result;
	}
	DualCombination lb(const DualCombination& D) {
		// Проверка ОДЗ логарифма
		if (D[0] <= 0.0) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); k++) {
			if (k == 0) {
				result += ::lb(D[0]);
			}
			else {
				if ((k - 1) % 2 == 0) {
					result += factorial(k - 1) * std::pow(D[0], -1 * k) * simple_pow(D - D[0], k) / (::ln(2) * factorial(k));
				}
				else {
					result += -1 * factorial(k - 1) * std::pow(D[0], -1 * k) * simple_pow(D - D[0], k) / (::ln(2) * factorial(k));
				}
			}
		}
		return result;
	}
	DualCombination log(const DualCombination& D1, const DualCombination& D2) {
		if ((D1[0] <= 0.0) || (D1[0] == 1) || (D2[0] <= 0.0)) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		return ln(D2) / ln(D1);
	}
	DualCombination log(const DualCombination& D, long double p) {
		if ((D[0] <= 0.0) || (D[0] == 1) || (p <= 0.0)) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		return ::ln(p) / ln(D);
	}
	DualCombination log(const DualCombination& D, const rational::Rational& r) {
		if ((D[0] <= 0.0) || (D[0] == 1) || (r <= rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the range of allowable values of the DualCombination logarithm");
		}
		return ln(r) / ln(D);
	}

	// Возведение в степень
	DualCombination sq(const DualCombination& D) {
		return D * D;
	}
	DualCombination cb(const DualCombination& D) {
		return D * D * D;
	}
	DualCombination pow(const DualCombination& D, int p) {
		if ((D[0] == 0.0) && (p <= 0)) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		if (p >= 0) {
			DualCombination result(D.order());
			DualCombination term(D.order());
			for (size_t k = 0; k < D.order(); ++k) {
				if (k <= p) {
					term = std::pow(D[0], p - k) * simple_pow(D - D[0], k) / factorial(k);
					for (size_t i = p; i >= p - k + 1; --i) {
						term *= i;
					}
					result += term;
				}
			}
			return result;
		}
		else {
			return DualCombination(1.0L, D.order()) / pow(D, -p);
		}
	}
	DualCombination pow(const DualCombination& D, int64_t p) {
		if ((D[0] == 0.0) && (p <= 0)) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		if (p >= 0) {
			DualCombination result(D.order());
			DualCombination term(D.order());
			for (size_t k = 0; k < D.order(); ++k) {
				if (k <= p) {
					term = std::pow(D[0], p - k) * simple_pow(D - D[0], k) / factorial(k);
					for (size_t i = p; i >= p - k + 1; --i) {
						term *= i;
					}
					result += term;
				}
			}
			return result;
		}
		else {
			return DualCombination(1.0L, D.order()) / pow(D, -p);
		}
	}
	DualCombination pow(const DualCombination& D, long double p) {
		if (D[0] <= 0.0) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		return exp(p * ln(D));
	}
	DualCombination pow(const DualCombination& D, const rational::Rational& r) {
		if (D[0] == 0.0) {
			if (r <= rational::Rational(0, 1)) {
				throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
			}
			return DualCombination(D.order());
		}
		if (D[0] < 0.0) {
			if (r.denominatorIsEven()) {
				throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
			}
			DualCombination result = pow(-D, r.toLongDouble());
			if (r.numeratorIsOdd()) {
				result *= -1;
			}
			return result;
		}
		else {
			return pow(D, r.toLongDouble());
		}
	}
	DualCombination pow(const DualCombination& D1, const rational::Rational& D2_Real, const rational::Rational& D2_Nilpotent) {
		if ((D1[0] == 0.0) && (D2_Real <= rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		if ((D1[0] < 0.0) && D2_Real.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		std::vector<long double> coeffs { D2_Real.toLongDouble(), D2_Nilpotent.toLongDouble() };
		return exp(DualCombination(coeffs) * ln(D1));
	}
	DualCombination pow(const DualCombination& D1, const rational::Rational& D2_Real, long double D2_Nilpotent) {
		if ((D1[0] == 0.0) && (D2_Real <= rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		if ((D1[0] < 0.0) && D2_Real.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		std::vector<long double> coeffs{ D2_Real.toLongDouble(), D2_Nilpotent };
		return exp(DualCombination(coeffs) * ln(D1));
	}
	DualCombination pow(long double p, const DualCombination& D) {
		if (p <= 0.0) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		return exp(::ln(p) * D);
	}
	DualCombination pow(const DualCombination& D1, const DualCombination& D2) {
		if (D1[0] <= 0.0) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		return exp(ln(D1) * D2);
	}
	DualCombination pow(const rational::Rational& r, const DualCombination& D) {
		if (r <= rational::Rational(0, 1)) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination power function");
		}
		return exp(ln(r) * D);
	}

	// Извлечение корня
	DualCombination sqrt(const DualCombination& D) {
		if (D[0] < 0) {
			throw std::domain_error("Exceeding the domain of allowable values of a DualCombination root function");
		}
		if (D[0] == 0.0) {
			if (D.order() > 1) {
				for (size_t i = 1; i < D.order(); ++i) {
					if (D[i] != 0.0) {
						throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
					}
				}
			}
			return D;
		}
		DualCombination result(D.order());
		DualCombination term(D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			term = std::pow(D[0], 0.5L - k) * simple_pow(D - D[0], k) / factorial(k);
			if (k > 0) {
				for (size_t i = 0; i <= k - 1; ++i) {
					term *= 0.5L - i;
				}
			}
			result += term;
		}
		return result;
	}
	DualCombination cbrt(const DualCombination& D) {
		DualCombination result(D.order());
		DualCombination term(D.order());
		if (D[0] == 0.0) {
			if (D.order() > 1) {
				for (size_t i = 1; i < D.order(); ++i) {
					if (D[i] != 0.0) {
						throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
					}
				}
			}
			return D;
		}
		for (size_t k = 0; k < D.order(); ++k) {
			term = std::pow(D[0], 1.0L / 3.0L - k) * simple_pow(D - D[0], k) / factorial(k);
			if (k > 0) {
				for (size_t i = 0; i <= k - 1; ++i) {
					term *= 1.0L / 3.0L - i;
				}
			}
			result += term;
		}
		return result;
	}
	DualCombination rt(const DualCombination& D, const rational::Rational& r) {
		if ((r == rational::Rational(0, 1)) || ((D[0] < 0) && r.numeratorIsEven())) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if (D[0] == 0.0) {
			if (r < rational::Rational(0, 1)) {
				throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
			}
			if (r == rational::Rational(1, 1)) {
				return D;
			}
			else {
				if (D.order() > 1) {
					for (size_t i = 1; i < D.order(); ++i) {
						if (D[i] != 0.0) {
							throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
						}
					}
				}
				return D;
			}
		}
		if (r >= rational::Rational(0, 1)) {
			DualCombination result(D.order());
			DualCombination term(D.order());
			for (size_t k = 0; k < D.order(); ++k) {
				term = std::pow(D[0], 1 / r.toLongDouble() - k) * simple_pow(D - D[0], k) / factorial(k);
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
			return DualCombination(1.0L, D.order()) / rt(D, -r);
		}
	}
	DualCombination rt(const DualCombination& D, long double p) {
		if ((p == 0.0) || (D[0] < 0)) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if (D[0] == 0.0) {
			if (p < 0.0) {
				throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
			}
			if (p == 1.0) {
				return D;
			}
			else {
				if (D.order() > 1) {
					for (size_t i = 1; i < D.order(); ++i) {
						if (D[i] != 0.0) {
							throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
						}
					}
				}
				return D;
			}
		}
		if (p >= 0) {
			DualCombination result(D.order());
			DualCombination term(D.order());
			for (size_t k = 0; k < D.order(); ++k) {
				term = std::pow(D[0], 1 / p - k) * simple_pow(D - D[0], k) / factorial(k);
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
			return DualCombination(1.0L, D.order()) / rt(D, -p);
		}
	}
	DualCombination rt(const DualCombination& D1, const DualCombination& D2) {
		if ((D2[0] == 0.0) || (D1[0] < 0)) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if ((D1[0] == 0) && (D2[0] < 0.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		return pow(D1, 1 / D2);
	}
	DualCombination rt(long double p, const DualCombination& D) {
		if ((D[0] == 0.0) || (p < 0)) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if ((p == 0) && (D[0] < 0.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		return pow(p, 1 / D);
	}
	DualCombination rt(const DualCombination& D1, const rational::Rational& D2_Real, long double D2_Nilpotent) {
		if (D2_Real == rational::Rational(0, 1)) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if ((D1[0] == 0) && (D2_Real < rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		if ((D1[0] < 0) && D2_Real.numeratorIsEven()) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination root function");
		}
		return pow(D1, rational::Rational(D2_Real.denominator(), D2_Real.numerator()), -D2_Nilpotent / (D2_Real * D2_Real));
	}
	DualCombination rt(const DualCombination& D1, const rational::Rational& D2_Real, std::vector<long double> D2_Nilpotent) {
		if (D2_Real == rational::Rational(0, 1)) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		if ((D1[0] == 0) && (D2_Real < rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		if ((D1[0] < 0) && D2_Real.numeratorIsEven()) {
			throw std::domain_error("Exceeding the domain of valid values of the HyperDualCombination root function");
		}
		std::vector<long double> coeffs(D2_Nilpotent.size() + 1);
		for (size_t i = 0; i < coeffs.size(); ++i) {
			if (i == 0) {
				coeffs[i] = D2_Real.toLongDouble();
			}
			else {
				coeffs[i] = D2_Nilpotent[i - 1];
			}
		}
		return exp(ln(D1) / DualCombination(coeffs));
	}

	// Тригонометрические функции
	DualCombination sin(const DualCombination& D) {
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += std::sin(D[0] + ::pi() * k / 2) * simple_pow(D - D[0], k) / factorial(k);
		}
		return result;
	}
	DualCombination cos(const DualCombination& D) {
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += std::cos(D[0] + ::pi() * k / 2) * simple_pow(D - D[0], k) / factorial(k);
		}
		return result;
	}
	DualCombination tg(const DualCombination& D) {
		if (std::cos(D[0]) == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination tangent function");
		}
		return sin(D) / cos(D);
	}
	DualCombination ctg(const DualCombination& D) {
		if (std::sin(D[0]) == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination cotangent function");
		}
		return cos(D) / sin(D);
	}
	DualCombination sec(const DualCombination& D) {
		if (std::cos(D[0]) == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination secant function");
		}
		return 1 / cos(D);
	}
	DualCombination cosec(const DualCombination& D) {
		if (std::sin(D[0]) == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination cosecant function");
		}
		return 1 / sin(D);
	}

	// Обратные тригонометрические функции
	long double arcsin_derivative(long double p, size_t k) {
		if ((p < -1.0) || (p > 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the arcsine_derivative function");
		}
		if (k == 0) {
			return ::arcsin(p);
		}
		else if (k == 1) {
			return 1 / std::sqrt(1 - ::sq(p));
		}
		else {
			return ((2 * k - 3) * p * arcsin_derivative(p, k - 1) + ::sq(k - 2) * arcsin_derivative(p, k - 2)) / (1 - ::sq(p));
		}
	}
	DualCombination arcsin(const DualCombination& D) {
		if ((D[0] < -1.0) || (D[0] > 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination arcsine function");
		}
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += arcsin_derivative(D[0], k) * simple_pow(D - D[0], k) / factorial(k);
		}
		return result;
	}
	long double arccos_derivative(long double p, size_t k) {
		if ((p < -1.0) || (p > 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the arccosine_derivative function");
		}
		if (k == 0) {
			return ::arccos(p);
		}
		else if (k == 1) {
			return -1 / std::sqrt(1 - ::sq(p));
		}
		else {
			return ((2 * k - 3) * p * arccos_derivative(p, k - 1) + ::sq(k - 2) * arccos_derivative(p, k - 2)) / (1 - ::sq(p));
		}
	}
	DualCombination arccos(const DualCombination& D) {
		if ((D[0] < -1.0) || (D[0] > 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination arccosine function");
		}
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += arccos_derivative(D[0], k) * simple_pow(D - D[0], k) / factorial(k);
		}
		return result;
	}
	long double arctg_derivative(long double p, size_t k) {
		if (k == 0) {
			return ::arctg(p);
		}
		else if (k == 1) {
			return 1 / (1 + ::sq(p));
		}
		else {
			return -(2 * (k - 1) * p * arctg_derivative(p, k - 1) + (k - 1) * (k - 2) * arctg_derivative(p, k - 2)) / (1 + ::sq(p));
		}
	}
	DualCombination arctg(const DualCombination& D) {
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += arctg_derivative(D[0], k) * simple_pow(D - D[0], k) / factorial(k);
		}
		return result;
	}
	long double arcctg_derivative(long double p, size_t k) {
		if (k == 0) {
			return ::arcctg(p);
		}
		else if (k == 1) {
			return -1 / (1 + ::sq(p));
		}
		else {
			return -(2 * (k - 1) * p * arcctg_derivative(p, k - 1) + (k - 1) * (k - 2) * arcctg_derivative(p, k - 2)) / (1 + ::sq(p));
		}
	}
	DualCombination arcctg(const DualCombination& D) {
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			result += arcctg_derivative(D[0], k) * simple_pow(D - D[0], k) / factorial(k);
		}
		return result;
	}
	DualCombination arcsec(const DualCombination& D) {
		if ((D[0] > -1.0) && (D[0] < 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination arcsecant function");
		}
		return arccos(1 / D);
	}
	DualCombination arccosec(const DualCombination& D) {
		if ((D[0] > -1.0) && (D[0] < 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination arccosecant function");
		}
		return arcsin(1 / D);
	}

	// Гиперболические функции
	DualCombination sh(const DualCombination& D) {
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			if (k % 2 == 0) {
				result += ::sh(D[0]) * simple_pow(D - D[0], k) / factorial(k);
			}
			else {
				result += ::ch(D[0]) * simple_pow(D - D[0], k) / factorial(k);
			}
		}
		return result;
	}
	DualCombination ch(const DualCombination& D) {
		DualCombination result(D.order());
		for (size_t k = 0; k < D.order(); ++k) {
			if (k % 2 == 0) {
				result += ::ch(D[0]) * simple_pow(D - D[0], k) / factorial(k);
			}
			else {
				result += ::sh(D[0]) * simple_pow(D - D[0], k) / factorial(k);
			}
		}
		return result;
	}
	DualCombination th(const DualCombination& D) {
		return sh(D) / ch(D);
	}
	DualCombination cth(const DualCombination& D) {
		if (D[0] == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination hyperbolic cotangent function");
		}
		return ch(D) / sh(D);
	}
	DualCombination sch(const DualCombination& D) {
		return 1 / ch(D);
	}
	DualCombination csch(const DualCombination& D) {
		if (D[0] == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the DualCombination hyperbolic cosecant function");
		}
		return 1 / sh(D);
	}

	// Полные производные высших порядков
	template<typename Func>
	long double D(Func f, long double x, size_t k) {
		DualCombination x_dual = Nilpotent_Add(x, k + 1);
		DualCombination F = f(x_dual);
		return F[k] * factorial(k);
	}

	// Частные производные высших порядков
	template<typename Func>
	long double D(Func f, const std::vector<long double>& X, size_t x_index, size_t k) {
		size_t n = X.size();
		if (x_index >= n) throw std::out_of_range("Invalid variable index");
		std::vector<DualCombination> X_dual(n);

		for (size_t j = 0; j < n; ++j) {
			if (j == x_index) {
				X_dual[j] = Nilpotent_Add(X[j], k + 1);
			}
			else {
				X_dual[j] = DualCombination(X[j], k + 1);
			}
		}

		DualCombination F = f(X_dual);
		return F[k] * factorial(k);
	}

}