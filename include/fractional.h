#pragma once // Защита от множественного включения

#include <cmath>
#include <vector>
#include <stdexcept>
#include <limits>
#include <functional>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include "real_utils.h"
#include "rational.h"
#include "dual.h"
#include "multidual.h"
#include "hyperdual.h"

namespace dual {

    // Для полных производных функций одного вещественного переменного

    template<typename Func>
    long double I(const Func& f, long double x, long double order, long double a) {
        if (order <= 0.0L || order >= 1.0L) {
            throw std::domain_error("Riemann-Liouville_integral:: wrong order");
        }
        if (x <= a) {
            throw std::domain_error("Riemann-Liouville_integral:: must be x > a");
        }

        long double upper = std::pow(x - a, order);
        if (upper <= 0.0L) {
            upper = 1e-20L;
        }
        upper *= (1.0L - 10 * std::numeric_limits<long double>::epsilon());

        auto integral_expression = [&](long double u) -> long double {
            long double t = x - std::pow(u, 1.0L / order);

            if (t < a) {
                t = a;
            }
            if (t > x) {
                t = x;
            }

            return f(t) / order;
        };

        boost::math::quadrature::tanh_sinh<long double> integrator;
        long double error;
        return integrator.integrate(integral_expression, 0.0L, upper, 1e-10L, &error) / boost::math::tgamma(order);
    }

    template<typename Func>
    DualCombination I(const Func& f, DualCombination x, long double order, long double a) {
        if (order <= 0.0L || order >= 1.0L) {
            throw std::domain_error("Riemann-Liouville_integral:: wrong order");
        }
        if (x[0] <= a) {
            throw std::domain_error("Riemann-Liouville_integral:: must be x > a");
        }

        long double upper = std::pow(x[0] - a, order);
        if (upper <= 0.0L) {
            upper = 1e-20L;
        }
        upper *= (1.0L - 10 * std::numeric_limits<long double>::epsilon());

        size_t n = x.order();
        DualCombination result(n);

        boost::math::quadrature::tanh_sinh<long double> integrator;
        long double error;
        for (size_t j = 0; j < n; ++j) {
            auto integral_expression = [&](long double u) -> long double {
                long double t = x[0] - std::pow(u, 1.0L / order);

                if (t < a) {
                    t = a;
                }
                if (t > x[0]) {
                    t = x[0];
                }

                DualCombination t_dual(t, n);
                DualCombination f_t = f(t_dual);
                return f_t[j] / order;
                };
            result[j] = integrator.integrate(integral_expression, 0.0L, upper, 1e-10L, &error) / boost::math::tgamma(order);
        }

        return result;
    }

    template<typename Func>
    long double RL_D(const Func& f, long double x, long double order, long double a = 0.0L) {
        if (x <= a) {
            throw std::domain_error("Riemann-Liouville_derivative:: must be x > a");
        }
        if (order < 0.0L) {
            throw std::domain_error("Riemann-Liouville_derivative:: order must be non-negative");
        }
        if (order == 0.0L) {
            DualCombination x_dual = Nilpotent_Add(x, 1);
            DualCombination F = f(x_dual);
            return F[0];
        }

        int n = static_cast<int>(std::ceil(order));
        long double beta = static_cast<long double>(n) - order;

        auto integral = [&](DualCombination t) -> DualCombination {
            return I(f, t, beta, a);
        };

        return D(integral, x, n);
    }

    template<typename Func>
    long double C_D(const Func& f, long double x, long double order, long double a = 0.0L) {
        if (x <= a) {
            throw std::domain_error("Caputo_derivative:: must be x > a");
        }
        if (order < 0.0L) {
            throw std::domain_error("Caputo_derivative:: order must be non-negative");
        }
        if (order == 0.0L) {
            DualCombination x_dual = Nilpotent_Add(x, 1);
            DualCombination F = f(x_dual);
            return F[0];
        }

        int n = static_cast<int>(std::ceil(order));
        long double beta = static_cast<long double>(n) - order;

        long double upper = std::pow(x - a, beta);
        if (upper <= 0.0L) {
            upper = 1e-20L;
        }

        auto derivative = [&](long double t) -> long double {
            return D(f, t, n);
        };

        return I(derivative, x, beta, a);
    }

    // Для частных производных функций нескольких вещественных переменных

    template<typename Func>
    long double I(const Func& f, const std::vector<long double>& X, size_t x_index, long double order, long double a) {
        if (order <= 0.0L || order >= 1.0L) {
            throw std::domain_error("Riemann-Liouville_integral:: wrong order");
        }
        if (X[x_index] <= a) {
            throw std::domain_error("Riemann-Liouville_integral:: must be x > a");
        }

        long double upper = std::pow(X[x_index] - a, order);
        if (upper <= 0.0L) {
            upper = 1e-20L;
        }
        upper *= (1.0L - 10 * std::numeric_limits<long double>::epsilon());

        auto integral_expression = [&](long double u) -> long double {
            long double t = X[x_index] - std::pow(u, 1.0L / order);

            if (t < a) {
                t = a;
            }
            if (t > X[x_index]) {
                t = X[x_index];
            }

            std::vector<long double> X_t = X;
            X_t[x_index] = t;
            return f(X_t) / order;
        };

        boost::math::quadrature::tanh_sinh<long double> integrator;
        long double error;
        return integrator.integrate(integral_expression, 0.0L, upper, 1e-10L, &error) / boost::math::tgamma(order);
    }

    template<typename Func>
    DualCombination I(const Func& f, const std::vector<DualCombination>& X, size_t x_index, long double order, long double a) {
        DualCombination x = X[x_index];

        if (order <= 0.0L || order >= 1.0L) {
            throw std::domain_error("Riemann-Liouville_integral:: wrong order");
        }
        if (x[0] <= a) {
            throw std::domain_error("Riemann-Liouville_integral:: must be x > a");
        }

        long double upper = std::pow(x[0] - a, order);
        if (upper <= 0.0L) {
            upper = 1e-20L;
        }
        upper *= (1.0L - 10 * std::numeric_limits<long double>::epsilon());

        size_t n = x.order();
        DualCombination result(n);

        boost::math::quadrature::tanh_sinh<long double> integrator;
        long double error;
        for (size_t j = 0; j < n; ++j) {
            auto integral_expression = [&](long double u) -> long double {
                long double t = x[0] - std::pow(u, 1.0L / order);

                if (t < a) {
                    t = a;
                }
                if (t > x[0]) {
                    t = x[0];
                }

                std::vector<DualCombination> X_t = X;
                X_t[x_index] = DualCombination(t, n);
                DualCombination f_t = f(X_t);
                return f_t[j] / order;
            };
            result[j] = integrator.integrate(integral_expression, 0.0L, upper, 1e-10L, &error) / boost::math::tgamma(order);
        }

        return result;
    }

    template<typename Func>
    long double RL_D(const Func& f, const std::vector<long double>& X, size_t x_index, long double order, long double a = 0.0L) {
        if (X[x_index] <= a) {
            throw std::domain_error("Riemann-Liouville_derivative:: must be x > a");
        }
        if (order < 0.0L) {
            throw std::domain_error("Riemann-Liouville_derivative:: order must be non-negative");
        }
        if (order == 0.0L) {
            std::vector<DualCombination> X_dual(X.size());
            for (int i = 0; i < X_dual.size(); ++i) {
                X_dual[i] = Nilpotent_Add(X[i], 1);
            }
            DualCombination F = f(X_dual);
            return F[0];
        }

        int n = static_cast<int>(std::ceil(order));
        long double beta = static_cast<long double>(n) - order;

        auto integral = [&](const std::vector<DualCombination>& X_t) -> DualCombination {
            return I(f, X_t, x_index, beta, a);
        };

        return D(integral, X, x_index, n);
    }

    template<typename Func>
    long double C_D(const Func& f, const std::vector<long double>& X, size_t x_index, long double order, long double a = 0.0L) {
        if (X[x_index] <= a) {
            throw std::domain_error("Caputo_derivative:: must be x > a");
        }
        if (order < 0.0L) {
            throw std::domain_error("Caputo_derivative:: order must be non-negative");
        }
        if (order == 0.0L) {
            std::vector<DualCombination> X_dual(X.size());
            for (int i = 0; i < X_dual.size(); ++i) {
                X_dual[i] = Nilpotent_Add(X[i], 1);
            }
            DualCombination F = f(X_dual);
            return F[0];
        }

        int n = static_cast<int>(std::ceil(order));
        long double beta = static_cast<long double>(n) - order;

        long double upper = std::pow(X[x_index] - a, beta);
        if (upper <= 0.0L) {
            upper = 1e-20L;
        }
        upper *= (1.0L - 10 * std::numeric_limits<long double>::epsilon());

        auto derivative = [&](const std::vector<long double>& X_t) -> long double {
            return D(f, X_t, x_index, n);
        };

        return I(derivative, X, x_index, beta, a);
    }

    // Для смешанных производных функций нескольких вещественных переменных

    template<typename Func>
    HyperDualCombination I(const Func& f, const std::vector<HyperDualCombination>& X, const std::vector<long double>& orders, const std::vector<long double>& a = std::vector<long double>(X.size(), 0.0L)) {
        size_t s = X.size(); // Количество аргументов функции

        if (orders.size() != s) {
            throw std::invalid_argument("Riemann-Liouville_integral:: derivatives_orders size must match X size");
        }
        if (a.size() != s) {
            throw std::invalid_argument("Riemann-Liouville_integral:: a size must match X size");
        }

        std::function<HyperDualCombination(const std::vector<HyperDualCombination>&)> F = f;

        for (size_t j = 0; j < s; ++j) {
            if (orders[j] <= 0.0L || orders[j] >= 1.0L) {
                throw std::domain_error("Riemann-Liouville_integral:: wrong order");
            }
            if (orders[j] <= 1e-12L) {
                continue;
            }
            F = [F, j, order_j = orders[j], a_j = a[j]](const std::vector<HyperDualCombination>& X_arguments) -> HyperDualCombination {
                long double x = X_arguments[j][0];

                long double upper = std::pow(x - a_j, order_j);
                if (upper <= 0.0L) {
                    upper = 1e-20L;
                }
                upper *= (1.0L - 10 * std::numeric_limits<long double>::epsilon());

                if (x <= a_j) {
                    throw std::domain_error("Riemann-Liouville_integral:: must be x > a");
                }

                HyperDualCombination result(X_arguments[0].get_orders());

                boost::math::quadrature::tanh_sinh<long double> integrator;
                long double error;

                for (size_t i = 0; i < result.total_size(); ++i) {
                    auto integral_expression = [&](long double u) -> long double {
                        long double t = x - std::pow(u, 1.0L / order_j);

                        if (t < a_j) {
                            t = a_j;
                        }
                        if (t > x) {
                            t = x;
                        }

                        std::vector<HyperDualCombination> X_t = X_arguments;
                        X_t[j] = HyperDualCombination(X_arguments[0].get_orders(), t);
                        HyperDualCombination F_t = F(X_t);
                        return F_t[i] / order_j;
                    };
                    result[i] = integrator.integrate(integral_expression, 0.0L, upper, 1e-10L, &error) / boost::math::tgamma(order_j);
                }

                return result;
            };
        }

        return F(X);
    }

    template<typename Func>
    long double RL_D(const Func& f, const std::vector<long double>& X, const std::vector<long double>& derivatives_orders, const std::vector<long double>& a = std::vector<long double>(X.size(), 0.0L)) {
        size_t s = X.size(); // Количество аргументов функции
        if (derivatives_orders.size() != s) {
            throw std::invalid_argument("Riemann-Liouville_derivative:: derivatives_orders size must match X size");
        }
        if (a.size() != s) {
            throw std::invalid_argument("Riemann-Liouville_derivative:: a size must match X size");
        }

        std::vector<int> n(s);
        std::vector<long double> beta(s);

        for (size_t i = 0; i < s; ++i) {
            if (X[i] <= a[i]) {
                throw std::domain_error("Riemann-Liouville_derivative:: must be x > a");
            }
            if (derivatives_orders[i] < 0.0L) {
                throw std::domain_error("Riemann-Liouville_derivative:: order must be non-negative");
            }

            n[i] = static_cast<int>(std::ceil(derivatives_orders[i]));
            beta[i] = static_cast<long double>(n[i]) - derivatives_orders[i];
        }

        std::function<HyperDualCombination(const std::vector<HyperDualCombination>&)> integral = [&](const std::vector<HyperDualCombination>& X_t) -> HyperDualCombination {
            return I(f, X_t, beta, a);
        };

        std::vector<size_t> n_size_t(n.begin(), n.end());

        return D(integral, X, n_size_t);
    }

    template<typename Func>
    long double C_D(const Func& f, const std::vector<long double>& X, const std::vector<long double>& derivatives_orders, const std::vector<long double>& a = std::vector<long double>(X.size(), 0.0L)) {
        size_t s = X.size(); // Количество аргументов функции
        if (derivatives_orders.size() != s) {
            throw std::invalid_argument("Caputo_derivative:: derivatives_orders size must match X size");
        }
        if (a.size() != s) {
            throw std::invalid_argument("Caputo_derivative:: a size must match X size");
        }

        std::vector<int> n(s);
        std::vector<long double> beta(s);
        for (size_t i = 0; i < s; ++i) {
            if (X[i] <= a[i]) {
                throw std::domain_error("Caputo_derivative:: must be x > a");
            }
            if (derivatives_orders[i] < 0.0L) {
                throw std::domain_error("Caputo_derivative:: order must be non-negative");
            }

            n[i] = static_cast<int>(std::ceil(derivatives_orders[i]));
            beta[i] = static_cast<long double>(n[i]) - derivatives_orders[i];
        }

        std::vector<size_t> n_size_t(n.begin(), n.end());

        auto derivative = [&](const std::vector<long double>& X_t) -> long double {
            return D(f, X_t, n_size_t);
        };

        // Последовательное применение интегральных операторов
        std::function<long double(const std::vector<long double>&)> integral = derivative;
        for (size_t j = 0; j < s; ++j) {
            if (beta[j] > 1e-12L) {
                integral = [integral, j, order_j = beta[j], a_j = a[j]](const std::vector<long double>& X_t) -> long double {
                    return I(integral, X_t, j, order_j, a_j);
                };
            }
        }

        return integral(X);
    }

    // Для рациональных порядков производных

    template<typename Func>
    long double I(const Func& f, long double x, const rational::Rational& order, long double a) {
        return I(f, x, order.toLongDouble(), a);
    }

    template<typename Func>
    long double RL_D(const Func& f, long double x, const rational::Rational& order, long double a = 0.0L) {
        return RL_D(f, x, order.toLongDouble(), a);
    }

    template<typename Func>
    long double C_D(const Func& f, long double x, const rational::Rational& order, long double a = 0.0L) {
        return C_D(f, x, order.toLongDouble(), a);
    }

    template<typename Func>
    long double I(const Func& f, const std::vector<long double>& X, size_t x_index, const rational::Rational& order, long double a) {
        return I(f, X, x_index, order.toLongDouble(), a);
    }

    template<typename Func>
    long double RL_D(const Func& f, const std::vector<long double>& X, size_t x_index, const rational::Rational& order, long double a = 0.0L) {
        return RL_D(f, X, x_index, order.toLongDouble(), a);
    }

    template<typename Func>
    long double C_D(const Func& f, const std::vector<long double>& X, size_t x_index, const rational::Rational& order, long double a = 0.0L) {
        return C_D(f, X, x_index, order.toLongDouble(), a);
    }

    template<typename Func>
    long double RL_D(const Func& f, const std::vector<long double>& X, const std::vector<rational::Rational>& derivatives_orders, const std::vector<long double>& a = std::vector<long double>(X.size(), 0.0L)) {
        std::vector<long double> long_double_derivatives_orders(derivatives_orders.size());
        for (size_t i = 0; i < derivatives_orders.size(); ++i) {
            long_double_derivatives_orders[i] = derivatives_orders[i].toLongDouble();
        }
        return RL_D(f, X, long_double_derivatives_orders, a);
    }

    template<typename Func>
    long double C_D(const Func& f, const std::vector<long double>& X, const std::vector<rational::Rational>& derivatives_orders, const std::vector<long double>& a = std::vector<long double>(X.size(), 0.0L)) {
        std::vector<long double> long_double_derivatives_orders(derivatives_orders.size());
        for (size_t i = 0; i < derivatives_orders.size(); ++i) {
            long_double_derivatives_orders[i] = derivatives_orders[i].toLongDouble();
        }
        return C_D(f, X, long_double_derivatives_orders, a);
    }

}