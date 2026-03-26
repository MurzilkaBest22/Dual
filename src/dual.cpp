#include "real_utils.h"
#include "dual.h" // Подключение заголовка с объявлениями
#include <stdexcept> // Для обработки ошибок

namespace dual {

	// Реализация конструктора
	Dual::Dual(long double a, long double b) : a(a), b(b) {}

	// Операторы сравнения
	bool Dual::operator== (const Dual& other) const {
		return a == other.a && b == other.b; // Так как при инициализации оба числа приводятся к виду несократимой дроби
	}
	bool Dual::operator!= (const Dual& other) const {
		return a != other.a || b != other.b; // Так как при инициализации оба числа приводятся к виду несократимой дроби
	}

	// Арифметические операторы для Dual (ключевое слово friend не нужно, так как используется только при объявлении)
	Dual operator+(const Dual& w1, const Dual& w2) {
		return Dual(w1.a + w2.a, w1.b + w2.b);
	}
	Dual operator+(long double p, const Dual& w) {
		return Dual(w.a + p, w.b);
	}
	Dual operator+(const Dual& w, long double p) {
		return Dual(w.a + p, w.b);
	}
	Dual operator+(const rational::Rational& r, const Dual& w) {
		return Dual(w.a + r, w.b);
	}
	Dual operator+(const Dual& w, const rational::Rational& r) {
		return Dual(w.a + r, w.b);
	}

	Dual operator-(const Dual& w1, const Dual& w2) {
		return Dual(w1.a - w2.a, w1.b - w2.b);
	}
	Dual operator-(long double p, const Dual& w) {
		return Dual(p - w.a, -w.b);
	}
	Dual operator-(const Dual& w, long double p) {
		return Dual(w.a - p, w.b);
	}
	Dual operator-(const rational::Rational& r, const Dual& w) {
		return Dual(r - w.a, -w.b);
	}
	Dual operator-(const Dual& w, const rational::Rational& r) {
		return Dual(w.a - r, w.b);
	}

	Dual operator*(const Dual& w1, const Dual& w2) {
		return Dual(w1.a * w2.a, w1.a * w2.b + w2.a * w1.b);
	}
	Dual operator*(long double p, const Dual& w) {
		return Dual(p * w.a, p * w.b);
	}
	Dual operator*(const Dual& w, long double p) {
		return Dual(p * w.a, p * w.b);
	}
	Dual operator*(const rational::Rational& r, const Dual& w) {
		return Dual(r * w.a, r * w.b);
	}
	Dual operator*(const Dual& w, const rational::Rational& r) {
		return Dual(r * w.a, r * w.b);
	}

	Dual operator/(const Dual& w1, const Dual& w2) {
		// Проверка возможности провести деление
		if (w2.a == 0.0) {
			throw std::runtime_error("Division by zero in Dual::operator/");
		}
		return Dual(w1.a / w2.a, (w2.a * w1.b - w1.a * w2.b) / (w2.a * w2.a));
	}
	Dual operator/(long double p, const Dual& w) {
		// Проверка возможности провести деление
		if (w.a == 0.0) {
			throw std::runtime_error("Division by zero in Dual::operator/");
		}
		return Dual(p / w.a, (-p * w.b) / (w.a * w.a));
	}
	Dual operator/(const Dual& w, long double p) {
		// Проверка возможности провести деление
		if (p == 0.0) {
			throw std::runtime_error("Division by zero in Dual::operator/");
		}
		return Dual(w.a / p, w.b / p);
	}
	Dual operator/(const rational::Rational& r, const Dual& w) {
		// Проверка возможности провести деление
		if (w.a == 0.0) {
			throw std::runtime_error("Division by zero in Dual::operator/");
		}
		return Dual(r / w.a, (-r.toLongDouble() * w.b) / (w.a * w.a));
	}
	Dual operator/(const Dual& w, const rational::Rational& r) {
		// Проверка возможности провести деление
		if (r == rational::Rational(0, 1)) {
			throw std::runtime_error("Division by zero in Dual::operator/");
		}
		return Dual(w.a / r, w.b / r);
	}

	// Составные операторы присваивания
	Dual& Dual::operator+=(const Dual& other) {
		a += other.a;
		b += other.b;
		return *this;
	}
	Dual& Dual::operator+=(long double other) {
		a += other;
		return *this;
	}

	Dual& Dual::operator-=(const Dual& other) {
		a -= other.a;
		b -= other.b;
		return *this;
	}
	Dual& Dual::operator-=(long double other) {
		a -= other;
		return *this;
	}

	Dual& Dual::operator*=(const Dual& other) {
		b = a * other.b + b * other.a;
		a *= other.a;
		return *this;
	}
	Dual& Dual::operator*=(long double other) {
		a *= other;
		b *= other;
		return *this;
	}

	Dual& Dual::operator/=(const Dual& other) {
		if (other.a == 0) {
			throw std::runtime_error("Division by zero in Dual::operator/=");
		}
		b = (b * other.a - a * other.b) / (other.a * other.a);
		a = a / other.a;
		return *this;
	}
	Dual& Dual::operator/=(long double other) {
		if (other == 0) {
			throw std::runtime_error("Division by zero in Dual::operator/=");
		}
		a /= other;
		b /= other;
		return *this;
	}

	// Унарный минус
	Dual Dual::operator-() const {
		return Dual(-a, -b);
	}

	// Экспонента
	Dual exp(const Dual& w) {
		return Dual(std::exp(w.a), w.b * std::exp(w.a));
	}

	// Логарифмы
	Dual log(long double p, const Dual& w) {
		// Проверка ОДЗ логарифма
		if ((p <= 0.0) || (p == 1.0) || (w.a <= 0.0)) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return Dual(::log(p, w.a), w.b / (w.a * ::ln(p)));
	}
	Dual log(const rational::Rational& r, const Dual& w) {
		// Проверка ОДЗ логарифма
		if ((r <= rational::Rational(0, 1)) || (r == rational::Rational(1, 1)) || (w.a <= 0.0)) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return Dual(log(r, w.a), w.b / (w.a * ln(r)));
	}
	Dual ln(const Dual& w) {
		// Проверка ОДЗ логарифма
		if (w.a <= 0.0) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return Dual(::ln(w.a), w.b / w.a);
	}
	Dual log(const Dual& w) {
		// Проверка ОДЗ логарифма
		if (w.a <= 0.0) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return Dual(::ln(w.a), w.b / w.a);
	}
	Dual lg(const Dual& w) {
		// Проверка ОДЗ логарифма
		if (w.a <= 0.0) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return Dual(::lg(w.a), w.b / (w.a * ::ln(10)));
	}
	Dual lb(const Dual& w) {
		// Проверка ОДЗ логарифма
		if (w.a <= 0.0) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return Dual(::lb(w.a), w.b / (w.a * ::ln(2)));
	}
	Dual log(const Dual& w1, const Dual& w2) {
		if ((w1.a <= 0.0) || (w1.a == 1) || (w2.a <= 0.0)) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return ln(w2) / ln(w1);
	}
	Dual log(const Dual& w, long double p) {
		if ((w.a <= 0.0) || (w.a == 1) || (p <= 0.0)) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return ::ln(p) / ln(w);
	}
	Dual log(const Dual& w, const rational::Rational& r) {
		if ((w.a <= 0.0) || (w.a == 1) || (r <= rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the range of allowable values of the dual logarithm");
		}
		return ln(r) / ln(w);
	}

	// Возведение в степень
	Dual sq(const Dual& w) {
		return Dual(::sq(w.a), 2 * w.a * w.b);
	}
	Dual cb(const Dual& w) {
		return Dual(::cb(w.a), 3 * ::sq(w.a) * w.b);
	}
	Dual pow(const Dual& w, int p) { // Всегда принимает ещё и отрицательные основания степени
		if ((w.a == 0.0) && (p <= 0)) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return Dual(std::pow(w.a, p), p * std::pow(w.a, p - 1) * w.b);
	}
	Dual pow(const Dual& w, int64_t p) { // Всегда принимает ещё и отрицательные основания степени
		if ((w.a == 0.0) && (p <= 0)) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return Dual(std::pow(w.a, p), p * std::pow(w.a, p - 1) * w.b);
	}
	Dual pow(const Dual& w, const rational::Rational& r) {
		if ((w.a == 0.0) && (r <= rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		if ((w.a < 0.0) && r.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return Dual(rational::pow(w.a, r), r * std::pow(w.a, r - 1.0L) * w.b);
	}
	Dual pow(const Dual& w1, const rational::Rational& w2_Real, const rational::Rational& w2_Nilpotent) { // На случай, если вещественная часть показателя рациональна
		if ((w1.a == 0.0) && (w2_Real <= rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		if ((w1.a < 0.0) && w2_Real.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return exp(Dual(w2_Real.toLongDouble(), w2_Nilpotent.toLongDouble()) * ln(w1));
	}
	Dual pow(const Dual& w1, const rational::Rational& w2_Real, long double w2_Nilpotent) {
		if ((w1.a == 0.0) && (w2_Real <= rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		if ((w1.a < 0.0) && w2_Real.denominatorIsEven()) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return exp(Dual(w2_Real.toLongDouble(), w2_Nilpotent) * ln(w1));
	}
	Dual pow(const Dual& w, long double p) { // Показатель степени воспринимается как иррациональный, в силу чего ОДЗ является положительными числами
		if (w.a <= 0.0) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return Dual(std::pow(w.a, p), p * std::pow(w.a, p - 1) * w.b);
	}
	Dual pow(long double p, const Dual& w) {
		if (p <= 0.0) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return Dual(std::pow(p, w.a), w.b * ::ln(p) * std::pow(p, w.a));
	}
	Dual pow(const Dual& w1, const Dual& w2) {
		if (w1.a <= 0.0) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return exp(w2 * ln(w1));
	}
	Dual pow(const rational::Rational& r, const Dual& w) {
		if (r <= rational::Rational(0, 1)) {
			throw std::domain_error("Exceeding the domain of allowable values of a dual power function");
		}
		return Dual(pow(r, w.a), w.b * ln(r) * pow(r, w.a));
	}

	// Извлечение корня
	Dual sqrt(const Dual& w) {
		if (w.a < 0) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return Dual(std::sqrt(w.a), 0.5 * w.b / std::sqrt(w.a));
	}
	Dual cbrt(const Dual& w) {
		return Dual(std::cbrt(w.a), w.b / (3 * std::pow(w.a, 2.0L / 3.0L)));
	}
	Dual rt(const Dual& w, const rational::Rational& r) {
		if ((r == rational::Rational(0, 1)) || ((w.a < 0) && r.numeratorIsEven())) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if ((w.a == 0) && (r < rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return pow(w, 1 / r);
	}
	Dual rt(const Dual& w, long double p) {
		if ((p == 0.0) || (w.a < 0)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if ((w.a == 0) && (p < 0.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return pow(w, 1 / p);
	}
	Dual rt(const Dual& w1, const Dual& w2) {
		if ((w2.a == 0.0) || (w1.a < 0)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if ((w1.a == 0) && (w2.a < 0.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return pow(w1, 1 / w2);
	}
	Dual rt(long double p, const Dual& w) {
		if ((w.a == 0.0) || (p < 0)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if ((p == 0) && (w.a < 0.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return pow(p, 1 / w);
	}
	Dual rt(const Dual& w1, const rational::Rational& w2_Real, const rational::Rational& w2_Nilpotent) {
		if (w2_Real == rational::Rational(0, 1)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if ((w1.a == 0) && (w2_Real < rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if ((w1.a < 0) && w2_Real.numeratorIsEven()) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return pow(w1, 1 / Dual(w2_Real.toLongDouble(), w2_Nilpotent.toLongDouble()));
	}
	Dual rt(const Dual& w1, const rational::Rational& w2_Real, long double w2_Nilpotent) {
		if (w2_Real == rational::Rational(0, 1)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if ((w1.a == 0) && (w2_Real < rational::Rational(0, 1))) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		if ((w1.a < 0) && w2_Real.numeratorIsEven()) {
			throw std::domain_error("Exceeding the domain of valid values of the dual root function");
		}
		return pow(w1, 1 / Dual(w2_Real.toLongDouble(), w2_Nilpotent));
	}

	// Тригонометрические функции
	Dual sin(const Dual& w) {
		return Dual(std::sin(w.a), w.b * std::cos(w.a));
	}
	Dual cos(const Dual& w) {
		return Dual(std::cos(w.a), -w.b * std::sin(w.a));
	}
	Dual tg(const Dual& w) {
		if (std::cos(w.a) == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the dual tangent function");
		}
		return Dual(::tg(w.a), w.b / std::pow(std::cos(w.a), 2));
	}
	Dual ctg(const Dual& w) {
		if (std::sin(w.a) == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the dual cotangent function");
		}
		return Dual(::ctg(w.a), -w.b / std::pow(std::sin(w.a), 2));
	}
	Dual sec(const Dual& w) {
		if (std::cos(w.a) == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the dual secant function");
		}
		return Dual(::sec(w.a), w.b * ::tg(w.a) / std::cos(w.a));
	}
	Dual cosec(const Dual& w) {
		if (std::sin(w.a) == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the dual cosecant function");
		}
		return Dual(::cosec(w.a), -w.b * ::ctg(w.a) / std::sin(w.a));
	}

	// Обратные тригонометрические функции
	Dual arcsin(const Dual& w) {
		if ((w.a < -1.0) || (w.a > 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual arcsine function");
		}
		return Dual(::arcsin(w.a), w.b / std::sqrt(1 - std::pow(w.a, 2)));
	}
	Dual arccos(const Dual& w) {
		if ((w.a < -1.0) || (w.a > 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual arccosine function");
		}
		return Dual(::arccos(w.a), -w.b / std::sqrt(1 - std::pow(w.a, 2)));
	}
	Dual arctg(const Dual& w) {
		return Dual(::arctg(w.a), w.b / (1 + std::pow(w.a, 2)));
	}
	Dual arcctg(const Dual& w) {
		return Dual(::arcctg(w.a), -w.b / (1 + std::pow(w.a, 2)));
	}
	Dual arcsec(const Dual& w) {
		if ((w.a > -1.0) && (w.a < 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual arcsecant function");
		}
		return Dual(::arcsec(w.a), w.b / std::sqrt(w.a * w.a * (w.a * w.a - 1)));
	}
	Dual arccosec(const Dual& w) {
		if ((w.a > -1.0) && (w.a < 1.0)) {
			throw std::domain_error("Exceeding the domain of valid values of the dual arccosecant function");
		}
		return Dual(::arccosec(w.a), -w.b / std::sqrt(w.a * w.a * (w.a * w.a - 1)));
	}

	// Гиперболические функции
	Dual sh(const Dual& w) {
		return Dual(::sh(w.a), w.b * ::ch(w.a));
	}
	Dual ch(const Dual& w) {
		return Dual(::ch(w.a), w.b * ::sh(w.a));
	}
	Dual th(const Dual& w) {
		return Dual(::th(w.a), w.b * std::pow(::sch(w.a), 2));
	}
	Dual cth(const Dual& w) {
		if (w.a == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the dual hyperbolic cotangent function");
		}
		return Dual(::cth(w.a), -w.b * std::pow(::csch(w.a), 2));
	}
	Dual sch(const Dual& w) {
		return Dual(::sch(w.a), -w.b * ::th(w.a) / ::ch(w.a));
	}
	Dual csch(const Dual& w) {
		if (w.a == 0.0) {
			throw std::domain_error("Exceeding the domain of valid values of the dual hyperbolic cosecant function");
		}
		return Dual(::csch(w.a), -w.b * ::cth(w.a) / ::sh(w.a));
	}

}
