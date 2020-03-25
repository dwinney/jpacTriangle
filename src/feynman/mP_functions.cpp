// Functions to evaluate the polynomial contribution to the inhomogeneity
// associated with the two-point "forward" rescattering diagram
//
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman/feynman_triangle.hpp"

// Once subtracted polynomial contribution (P^0(s))
std::complex<double> feynman_triangle::mP1(double s)
{
  integ.check_weights();

  std::complex<double> result, sum = 0.;
  for (int i = 0; i < integ.xN; i++)
  {
    double x_i = integ.abscissas[i];
    std::complex<double> temp = mPi2 - x_i * (1. - x_i) * (s + ieps) ;

    sum += integ.weights[i] * log(mPi2 / temp);
  }

  result = sum / M_PI;

  return result;
};
