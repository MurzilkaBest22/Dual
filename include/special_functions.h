#pragma once // Защита от множественного включения

#include <cmath>
#include <vector>
#include <stdexcept>
#include <limits>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include "real_utils.h"
#include "rational.h"
#include "dual.h"
#include "multidual.h"
#include "hyperdual.h"

namespace binomial {

    static long long binomial_coefficient(size_t n, size_t k) {
        if (k > n) {
            return 0;
        }
        if ((k == 0) || (k == n)) {
            return 1;
        }
        long long result = 1;
        for (size_t i = 1; i <= k; ++i) {
            result = result * (n - k + i) / i;
        }
        return result;
    }

}

namespace rational {

    // Гамма-функция
    long double tgamma(const Rational& r) {
        if ((r <= Rational(0, 1)) && (r.denominator() == 1)) {
            throw std::domain_error("rational::tgamma: pole at non-positive integer");
        }
        return boost::math::tgamma(r.toLongDouble());
    }

    // Бета-функция
    long double beta(const Rational& r1, const Rational& r2) {
        if ((r1 <= Rational(0, 1)) || (r2 <= Rational(0, 1))) {
            throw std::domain_error("rational::beta: arguments must be positive");
        }
        return boost::math::beta(r1.toLongDouble(), r2.toLongDouble());
    }

    // Дзета-функция Римана
    long double zeta(const Rational& r) {
        if (r == Rational(1, 1)) {
            throw std::domain_error("rational::zeta: arguments must be > 1");
        }
        return boost::math::zeta(r.toLongDouble());
    }

}

namespace dual {

    // Гамма-функция
    Dual tgamma(const Dual& w) {
        if ((w.a <= 0.0L) && (w.a == std::floor(w.a))) {
            throw std::domain_error("dual::tgamma: pole at non-positive integer");
        }
        return Dual(boost::math::tgamma(w.a), w.b * std::tgamma(w.a) * boost::math::digamma(w.a));
    }

    std::vector<long double> tgamma_derivatives(long double x, size_t n) {
        if (n == 0) {
            return {};
        }

        std::vector<long double> result(n);
        result[0] = boost::math::tgamma(x);

        for (size_t m = 1; m < n; ++m) {
            long double sum = 0.0L;
            for (size_t k = 0; k < m; ++k) {
                sum += binomial::binomial_coefficient(m - 1, k) * result[k] * boost::math::polygamma(static_cast<int>(m - 1 - k), x);
            }
            result[m] = sum;
        }
        return result;
    }

    DualCombination tgamma(const DualCombination& D) {
        if ((D[0] <= 0.0L) && (D[0] == std::floor(D[0]))) {
            throw std::domain_error("dual::tgamma: pole at non-positive integer");
        }
        DualCombination result(D.order());
        std::vector<long double> derivatives = tgamma_derivatives(D[0], D.order());
        for (size_t k = 0; k < D.order(); ++k) {
            result += derivatives[k] * simple_pow(D - D[0], k) / factorial(k);
        }
        return result;
    }

    HyperDualCombination tgamma(const HyperDualCombination& H) {
        if ((H[0] <= 0.0L) && (H[0] == std::floor(H[0]))) {
            throw std::domain_error("dual::tgamma: pole at non-positive integer");
        }
        HyperDualCombination result(H.get_orders());
        std::vector<long double> derivatives = tgamma_derivatives(H[0], H.max_degree_sum());
        for (size_t k = 0; k < H.max_degree_sum(); ++k) {
            result += derivatives[k] * simple_pow(H - H[0], k) / factorial(k);
        }
        return result;
    }

    // Бета-функция
    Dual beta(const Dual& w1, const Dual& w2) {
        if ((w1.a <= 0) || (w2.a <= 0)) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(w1) * tgamma(w2) / tgamma(w1 + w2);
    }
    Dual beta(long double p, const Dual& w) {
        if ((w.a <= 0) || (p <= 0)) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(w) * boost::math::tgamma(p) / tgamma(w + p);
    }
    Dual beta(const Dual& w, long double p) {
        if ((w.a <= 0) || (p <= 0)) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(w) * boost::math::tgamma(p) / tgamma(w + p);
    }
    Dual beta(const rational::Rational& r, const Dual& w) {
        if ((w.a <= 0) || (r <= rational::Rational(0,1))) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(w) * rational::tgamma(r) / tgamma(w + r);
    }
    Dual beta(const Dual& w, const rational::Rational& r) {
        if ((w.a <= 0) || (r <= rational::Rational(0,1))) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(w) * rational::tgamma(r) / tgamma(w + r);
    }

    DualCombination beta(const DualCombination& D1, const DualCombination& D2) {
        if ((D1[0] <= 0) || (D2[0] <= 0)) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(D1) * tgamma(D2) / tgamma(D1 + D2);
    }
    DualCombination beta(long double p, const DualCombination& D) {
        if ((D[0] <= 0) || (p <= 0)) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(D) * boost::math::tgamma(p) / tgamma(D + p);
    }
    DualCombination beta(const DualCombination& D, long double p) {
        if ((D[0] <= 0) || (p <= 0)) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(D) * boost::math::tgamma(p) / tgamma(D + p);
    }
    DualCombination beta(const rational::Rational& r, const DualCombination& D) {
        if ((D[0] <= 0) || (r <= rational::Rational(0, 1))) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(D) * rational::tgamma(r) / tgamma(D + r);
    }
    DualCombination beta(const DualCombination& D, const rational::Rational& r) {
        if ((D[0] <= 0) || (r <= rational::Rational(0, 1))) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(D) * rational::tgamma(r) / tgamma(D + r);
    }

    HyperDualCombination beta(const HyperDualCombination& H1, const HyperDualCombination& H2) {
        if ((H1[0] <= 0) || (H2[0] <= 0)) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(H1) * tgamma(H2) / tgamma(H1 + H2);
    }
    HyperDualCombination beta(long double p, const HyperDualCombination& H) {
        if ((H[0] <= 0) || (p <= 0)) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(H) * boost::math::tgamma(p) / tgamma(H + p);
    }
    HyperDualCombination beta(const HyperDualCombination& H, long double p) {
        if ((H[0] <= 0) || (p <= 0)) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(H) * boost::math::tgamma(p) / tgamma(H + p);
    }
    HyperDualCombination beta(const rational::Rational& r, const HyperDualCombination& H) {
        if ((H[0] <= 0) || (r <= rational::Rational(0, 1))) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(H) * rational::tgamma(r) / tgamma(H + r);
    }
    HyperDualCombination beta(const HyperDualCombination& H, const rational::Rational& r) {
        if ((H[0] <= 0) || (r <= rational::Rational(0, 1))) {
            throw std::domain_error("dual::beta: arguments must be positive");
        }
        return tgamma(H) * rational::tgamma(r) / tgamma(H + r);
    }

    // Дзета-функция Римана
    std::vector<long double> zeta_derivatives(long double x, size_t n, size_t h) {
        if (n == 0) {
            return {};
        }

        if (h < 2) {
            h = 2;
        }

        std::vector<long double> result(n);
        result[0] = boost::math::zeta(x);
        if (n == 1) {
            return result;
        }
        
        std::vector<long double> sums(n - 1, 0.0L);
        std::vector<long double> terms(n - 1);

        size_t k = 2;
        bool enough = false;

        size_t current_h = h;

        while (!enough) {
            long double ln_k = ::ln(k);
            long double pow_k = std::pow(k, -x);

            std::vector<long double> pow_ln(n - 1);
            pow_ln[0] = ln_k;

            for (size_t i = 1; i < n - 1; ++i) {
                pow_ln[i] = pow_ln[i - 1] * ln_k;
            }

            bool the_error_is_permissible = true;
            for (size_t i = 0; i < n - 1; ++i) {
                terms[i] = pow_ln[i] * pow_k;
                sums[i] += terms[i];
                if ((std::abs(terms[i]) > 1e-12L * std::abs(sums[i])) && (std::abs(terms[i]) > 1e-14L)) {
                    the_error_is_permissible = false;
                }
            }

            ++k;

            if (the_error_is_permissible) {
                enough = true;
                break;
            }

            if (!enough) {
                current_h *= 10;
            }
        }


        for (size_t i = 0; i < n - 1; ++i) {
            result[i + 1] = ((i + 1) % 2 == 0) ? sums[i] : -sums[i];
        }

        return result;
    }
   
    Dual zeta(const Dual& w, size_t h = 100000) {
        if (w.a <= 1) {
            throw std::domain_error("dual::zeta: argument must be > 1");
        }
        auto derivatives = zeta_derivatives(w.a, 2, h);
        return Dual(derivatives[0], w.b * derivatives[1]);
    }

    DualCombination zeta(const DualCombination& D, size_t h = 100000) {
        if (D[0] <= 1) {
            throw std::domain_error("dual::zeta: argument must be > 1");
        }
        DualCombination result(D.order());
        auto derivatives = zeta_derivatives(D[0], D.order(), h);
        for (size_t k = 0; k < D.order(); ++k) {
            result += derivatives[k] * simple_pow(D - D[0], k) / factorial(k);
        }
        return result;
    }

    HyperDualCombination zeta(const HyperDualCombination& H, size_t h = 100000) {
        if (H[0] <= 1) {
            throw std::domain_error("dual::zeta: argument must be > 1");
        }
        HyperDualCombination result(H.get_orders());
        auto derivatives = zeta_derivatives(H[0], H.max_degree_sum(), h);
        for (size_t k = 0; k < H.max_degree_sum(); ++k) {
            result += derivatives[k] * simple_pow(H - H[0], k) / factorial(k);
        }
        return result;
    }

}