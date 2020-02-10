// Kernel functions for the triangle convolution integral. Higher ordered
// kernels correspond to higher number of subtractions.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman_triangle.hpp"

// ---------------------------------------------------------------------------
// Once-subtracted
// ---------------------------------------------------------------------------

// Evaluate the kernel function by integrating over x
std::complex<double> feynman_triangle::mT1(double s, double t)
{
  check_weights();

  //integrate over x
  std::complex<double> sum1 = 0., sum2 = 0.;
  for (int i = 0; i < xN; i++)
  {
    double x_i = abscissas[i];

    // the integrand is the analytic form of the y_integral evaluated at:
    // y = 1 - x and y = 0
    sum1 += weights[i] * ( mT1_yintegral1(s, t, x_i) - mT1_yintegral1(0., t, x_i) );
    sum2 += weights[i] * ( mT1_yintegral2(s, t, x_i) - mT1_yintegral2(0., t, x_i) );
  }

  return sum2;
};

// ---------------------------------------------------------------------------
// Intermediate functions for the mT1 y_integral

// Analytic form of the first integral over y
std::complex<double> feynman_triangle::mT1_yintegral1(double s, double t, double x)
{
  std::complex<double> a, b, c, d;
  std::complex<double> e, f, g;

  // coeffs of denominator polynomial
  a = mPi2;
  b = mPi2 + (x - 1.) * mPi2 + x * (mDec2 + ieps) - x * s - (t - xi*.15);
  c = (1. - x) * t + x * mPi2 + x*(x-1.)* (mDec2 + ieps);

  // coeffs of numerator polynomial
  e = mPi2;
  f = x * (mDec2 + ieps + mPi2 - s);
  g = x*x * (mDec2 + ieps);

  // Evaluate definite integral with bounds y = [0, 1-x]
  return ri_poly2(1. - x, a, b, c, e, f, g) - ri_poly2(0., a, b, c, e, f, g);
};

// Analytic form of the second integral over y
std::complex<double> feynman_triangle::mT1_yintegral2(double s, double t, double x)
{
  std::complex<double> a, b, c, d;

  // coeffs of denominator polynomial
  a = mPi2;
  b = mPi2 + (x - 1.) * mPi2 + x * (mDec2 + ieps) - x * s - (t - xi*.15);
  c = (1. - x) * (t - xi*.15) + x * mPi2 + x * (x-1.) * (mDec2 + ieps);

  // Evaluate definite integral with bounds y = [0, 1-x]
  return ri_log1(1. - x, a, b, c) - ri_log1(0., a, b, c);
}
