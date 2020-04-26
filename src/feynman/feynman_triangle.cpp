// This object defines a Triangle amplitude, i.e. the rescattering
// diagram associated with an intermediate t-channel exchange.
//
// Evaluated specifically with the Feynman method
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman/feynman_triangle.hpp"

// ---------------------------------------------------------------------------
// Evaluate the convolution of the LHC function with triangle function
std::complex<double> feynman_triangle::eval(int j, int jp, double s)
{
    auto F = [&](double tp)
    {
      return lhc_func->disc(tp) *kernel.eval(1, j, jp, s, tp);
    };

    double error;
    std::complex<double> result;
    result = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(F, r_thresh, std::numeric_limits<double>::infinity(), 5, 1.E-9, &error);

    result /= M_PI;

    return 1. + result;
};
