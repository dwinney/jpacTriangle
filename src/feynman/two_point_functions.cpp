// Functions to evaluate the polynomial contribution to the inhomogeneity
// associated with the two-point "forward" rescattering diagram
//
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman_triangle.hpp"

std::complex<double> feynman_triangle::two_point_0(double s)
{
  check_weights();

  std::complex<double> result, sum = 0.;
  for (int i = 0; i < xN; i++)
  {
    double x_i = abscissas[i];

    std::complex<double> temp = x_i * (1. - x_i) * (s + ieps) - m1sq;
    sum += weights[i] * log(4. * M_PI / temp);
  }

  result = M_EULER - sum;
  result /= 16. * M_PI * M_PI;

  return result;
};
