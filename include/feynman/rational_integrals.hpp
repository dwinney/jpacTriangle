// Analytic forms of anti-derivatives of ratios of rational functions.
// These show up in integrals over feynman parameters
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include <complex>
#include <vector>
#include <cmath>
#include <iomanip>

#include "constants.hpp"

// redefine arctangent in terms of logarithms for complex argument
// to better control singularity structure
std::complex<double> c_atan(std::complex<double> z);

// Antiderivative of 1/(a y^2 + b y + c)
std::complex<double> ri_poly1(double y,
    std::complex<double> a,
    std::complex<double> b,
    std::complex<double> c);

// Antiderivative of (e y^2 + f y + g) / (a y ^2 + b y + c)
std::complex<double> ri_poly2(double y,
    std::complex<double> a,
    std::complex<double> b,
    std::complex<double> c,
    std::complex<double> e,
    std::complex<double> f,
    std::complex<double> g);

// Antiderivative of Log(1/(a y^2 + b y + c))
std::complex<double> ri_log1(double y,
    std::complex<double> a,
    std::complex<double> b,
    std::complex<double> c);
