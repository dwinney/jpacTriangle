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

#include <iostream>

// redefine arctangent in terms of logarithms for complex argument
// to better control singularity structure
std::complex<double> c_atan(std::complex<double> z);

// Antiderivative of 1/(a y^2 + b y + c)
std::complex<double> ri_poly0(double y,
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

// Antiderivative of (e y^4 + f y^3 + g y^2 + h y + i) / (a y ^2 + b y + c)
std::complex<double> ri_poly4(double y,
    std::complex<double> a,
    std::complex<double> b,
    std::complex<double> c,
    std::complex<double> e,
    std::complex<double> f,
    std::complex<double> g,
    std::complex<double> h,
    std::complex<double> i);

// Antiderivative of Log(1/(a y^2 + b y + c))
std::complex<double> ri_log0(double y,
    std::complex<double> a,
    std::complex<double> b,
    std::complex<double> c);

// Antiderivative of (e y^2 + f y + g) * Log(1/(a y^2 + b y + c))
std::complex<double> ri_log2(double y,
    std::complex<double> a,
    std::complex<double> b,
    std::complex<double> c,
    std::complex<double> e,
    std::complex<double> f,
    std::complex<double> g);
