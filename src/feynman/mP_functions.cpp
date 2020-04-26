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
  auto F = [&](double x)
  {
    std::complex<double> temp = mPi2 - x * (1. - x) * (s + ieps) ;
    return log(mPi2 / temp);
  };

  double error;
  std::complex<double> result;
  result = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(F, 0, 1, 5, 1.E-20, &error);

  result /= M_PI;

  return result;
};
