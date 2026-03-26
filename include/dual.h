#pragma once // Защита от множественного включения

#include <iostream> // Ввод-вывод
#include <stdexcept>
#include <cmath>
#include <array>
#include <vector>
#include "rational.h" // Рациональные числа

namespace dual {

	class Dual {
	public:
		long double a; // Действительная часть
		long double b; // Нильпотентная часть

		// Конструктор с параметрами по умолчанию
		Dual(long double a = 0.0, long double b = 0.0);

		// Операторы сравнения
		bool operator==(const Dual& other) const;
		bool operator!=(const Dual& other) const;

		// Декларации дружественных операторов
		friend Dual operator+(const Dual& w1, const Dual& w2);
		friend Dual operator+(long double p, const Dual& w);
		friend Dual operator+(const Dual& w, long double p);
		friend Dual operator+(const rational::Rational& r, const Dual& w);
		friend Dual operator+(const Dual& w, const rational::Rational& r);

		friend Dual operator-(const Dual& w1, const Dual& w2);
		friend Dual operator-(long double p, const Dual& w);
		friend Dual operator-(const Dual& w, long double p);
		friend Dual operator-(const rational::Rational& r, const Dual& w);
		friend Dual operator-(const Dual& w, const rational::Rational& r);

		friend Dual operator*(const Dual& w1, const Dual& w2);
		friend Dual operator*(long double p, const Dual& w);
		friend Dual operator*(const Dual& w, long double p);
		friend Dual operator*(const rational::Rational& r, const Dual& w);
		friend Dual operator*(const Dual& w, const rational::Rational& r);

		friend Dual operator/(const Dual& w1, const Dual& w2);
		friend Dual operator/(long double p, const Dual& w);
		friend Dual operator/(const Dual& w, long double p);
		friend Dual operator/(const rational::Rational& r, const Dual& w);
		friend Dual operator/(const Dual& w, const rational::Rational& r);

		// Составные операторы присваивания
		Dual& operator+=(const Dual& other);
		Dual& operator+=(long double other);

		Dual& operator-=(const Dual& other);
		Dual& operator-=(long double other);

		Dual& operator*=(const Dual& other);
		Dual& operator*=(long double other);

		Dual& operator/=(const Dual& other);
		Dual& operator/=(long double other);

		// Унарный минус
		Dual operator-() const;
	};

	// Экспонента
	Dual exp(const Dual& w);

	// Логарифмы
	Dual log(long double p, const Dual& w);
	Dual log(const rational::Rational& r, const Dual& w);
	Dual ln(const Dual& w);
	Dual log(const Dual& w); // Интерпретация логарифма без указанного основания как натурального
	Dual lg(const Dual& w);
	Dual lb(const Dual& w);
	Dual log(const Dual& w1, const Dual& w2);
	Dual log(const Dual& w, long double p);
	Dual log(const Dual& w, const rational::Rational& r);

	// Возведение в степень
	Dual sq(const Dual& w);
	Dual cb(const Dual& w);
	Dual pow(const Dual& w, int p);
	Dual pow(const Dual& w, int64_t p);
	Dual pow(const Dual& w, const rational::Rational& r);
	Dual pow(const Dual& w1, const rational::Rational& w2_Real, const rational::Rational& w2_Nilpotent); // На случай, если вещественная часть показателя рациональна
	Dual pow(const Dual& w1, const rational::Rational& w2_Real, long double w2_Nilpotent);
	Dual pow(const Dual& w, long double p);
	Dual pow(long double p, const Dual& w);
	Dual pow(const Dual& w1, const Dual& w2);
	Dual pow(const rational::Rational& r, const Dual& w);

	// Извлечение корня
	Dual sqrt(const Dual& w);
	Dual cbrt(const Dual& w);
	Dual rt(const Dual& w, const rational::Rational& r);
	Dual rt(const Dual& w, long double p);
	Dual rt(const Dual& w1, const Dual& w2);
	Dual rt(long double p, const Dual& w);
	Dual rt(const Dual& w1, const rational::Rational& w2_Real, const rational::Rational& w2_Nilpotent);
	Dual rt(const Dual& w1, const rational::Rational& w2_Real, long double w2_Nilpotent);

	// Тригонометрические функции
	Dual sin(const Dual& w);
	Dual cos(const Dual& w);
	Dual tg(const Dual& w);
	Dual ctg(const Dual& w);
	Dual sec(const Dual& w);
	Dual cosec(const Dual& w);

	// Обратные тригонометрические функции
	Dual arcsin(const Dual& w);
	Dual arccos(const Dual& w);
	Dual arctg(const Dual& w);
	Dual arcctg(const Dual& w);
	Dual arcsec(const Dual& w);
	Dual arccosec(const Dual& w);

	// Гиперболические функции
	Dual sh(const Dual& w);
	Dual ch(const Dual& w);
	Dual th(const Dual& w);
	Dual cth(const Dual& w);
	Dual sch(const Dual& w);
	Dual csch(const Dual& w);

	// Произвольная функция и оператор дифференцирования (реализованный в заголовке, так как является шаблонным)
	template<typename Func>
	long double D(Func f, long double x) {
		Dual x_dual(x, 1.0);
		Dual F = f(x_dual);
		return F.b;
	}

	// Частные производные
	template<size_t N, typename Func>
	long double D(Func f, const std::array<long double, N>& X, size_t x_index) {
		if (x_index >= N) throw std::out_of_range("Invalid variable index");
		std::array<Dual, N> X_dual;

		for (size_t j = 0; j < N; ++j) {
			if (j == x_index) {
				X_dual[j] = Dual(X[j], 1.0);
			}
			else {
				X_dual[j] = Dual(X[j], 0.0);
			}
		}

		Dual F = f(X_dual);
		return F.b;
	}

	// Градиент (вектор, а не массив)
	template<typename Func>
	std::vector<long double> grad(Func f, const std::vector<long double>& X) {
		size_t n = X.size();
		std::vector<long double> gradient(n);
		std::vector<Dual> X_dual(n);
		Dual F;

		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				X_dual[j] = Dual(X[j], (j == i) ? 1.0 : 0.0);
			}
			F = f(X_dual);
			gradient[i] = F.b;
		}

		return gradient;
	}

}
