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
    integ.check_weights();

    // Calculate the dispersion piece coming from the triangle
    std::complex<double> sum = 0.;
    for (int i = 0; i < integ.xN; i++)
    {
      double tp = r_thresh + tan(M_PI * integ.abscissas[i] / 2.);

      std::complex<double> temp;
      temp = lhc_func->disc(tp);
      temp *= triangle_kernel(1, j, jp, s, tp);

      temp *= (M_PI / 2.);
      temp /= pow(cos(M_PI * integ.abscissas[i] / 2.), 2.); // jacobian

      sum += integ.weights[i] * temp;
    }
    sum /= M_PI;

    return 1. + sum;
};
