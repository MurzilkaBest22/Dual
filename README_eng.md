DESCRIPTION
Dual is a full‑featured automatic differentiation library based on dual and hyper‑dual numbers. It provides:
Complete, partial, and mixed derivatives of arbitrary natural order.
Higher‑order gradients (including Hessians).
Fractional derivatives (Riemann‑Liouville and Caputo).
Support for rational numbers (when using the Rational library).
Special functions (gamma, beta, Riemann zeta).
Classes for dual combinations and hyper‑dual combinations.

FEATURES
Classes: Dual, DualCombination, HyperDualCombination.
Differentiation operators for various argument types.
Hybrid approach for fractional derivatives (AD + numerical integration).
Full domain checking for all functions.

INSTALLATION AND BUILD
The library consists of:
dual.h, dual.cpp - dual numbers
multidual.h, multidual.cpp - dual combinations
hyperdual.h, hyperdual.cpp - hyperdual combinations
special_functions.h, special_functions.cpp - special functions
fractional.h, fractional.cpp - fractional orders
alldiff.h, alldiff.cpp - the main file

TO USE:
Copy the Dual folder into your project.
Add all .cpp files (dual.cpp, special_functions.cpp, multidual.cpp, hyperdual.cpp, fractional.cpp) to your build.
Include the main header:

#include "Dual/alldiff.h"

(this includes all other headers).
Use the namespace dual.
If you don't need all features, you can include only the required headers.

DEPENDENCIES
RealUtils (required)
Rational (for rational numbers and fractional orders)
Boost.Math (for special functions and numerical integration). Make sure Boost is installed and accessible.

EXAMPLE 

#include <iostream>
#include "Dual/alldiff.h"

auto f = [](dual::Dual x) {
    return dual::sin(x) * dual::exp(x);
};

int main() {
    double x = 1.0;
    double df = dual::D(f, x);
    std::cout << "f'(1) = " << df << std::endl;
    return 0;
}

LICENSE
MIT License. See the LICENSE file.

AUTHOR
Mikhail D. Sychev
Email: murzilkabest@icloud.com
Telegram: @Murz1k22
