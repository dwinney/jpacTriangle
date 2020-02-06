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
    sum1 += weights[i] * (mT1_yintegral1(s,t,x_i,1.-x_i) - mT1_yintegral1(0.,t,x_i,1.-x_i));
    sum2 += weights[i] * (mT1_yintegral2(s,t,x_i,1.-x_i) - mT1_yintegral2(0.,t,x_i,1.-x_i));
  }

  std::complex<double> result;
  result = sum1 + sum2 / (4. * M_PI);

  return result / (16. * M_PI);
};

// ---------------------------------------------------------------------------
// Intermediate functions for the mT1 y_integral

// Analytic form of the first integral over y
std::complex<double> feynman_triangle::mT1_yintegral1(double s, double t, double x, double y)
{
  std::complex<double> a, b, c, d;
  std::complex<double> e, f, g;

  // coeffs of denominator polynomial
  a = mPi2;
  b = mPi2 + (x - 1.) * mPi2 + x * (mDec2 + ieps) - x * s - t;
  c = (1. - x) * t + x * mPi2 + x*(x-1.)* (mDec2 + ieps);
  d = b * b - 4. * a * c; // discriminant

  // coeffs of numerator polynomial
  e = mPi2;
  f = x * (mDec2 + ieps + mPi2 - s);
  g = x*x * (mDec2 + ieps);

  // Antiderivative done by mathematica
  std::complex<double> term1;
  term1 = atan((2. * a * y + b) / sqrt(-d));
  term1 *= b*b*e - a*f*b + 2.*a*(a*g-c*e);
  term1 /= a*a * sqrt(-d);

  std::complex<double> term2;
  term2 = log(a*y*y + b*y + c);
  term2 *= a*f - b*e;
  term2 /= 2.*a*a;

  std::complex<double> term3;
  term3 = e * y / a;

  return term1 + term2 + term3;
};

// Analytic form of the second integral over y
std::complex<double> feynman_triangle::mT1_yintegral2(double s, double t, double x, double y)
{
  std::complex<double> a, b, c, d;

  // coeffs of denominator polynomial
  a = mPi2;
  b = mPi2 + (x - 1.) * mPi2 + x * (mDec2 + ieps) - x * s - t;
  c = (1. - x) * t + x * mPi2 + x * (x-1.) * (mDec2 + ieps);
  d = b * b - 4. * a * c; // discriminant

  std::complex<double> term1;
  term1 = - atan((2.*a*y + b) / sqrt(-d));
  term1 *= sqrt(-d) / a;

  std::complex<double> term2;
  term2 = y * log(4.*M_PI / (a*y*y + b*y + c));

  std::complex<double> term3;
  term3 = log(a*y*y + b*y + c);
  term3 *= -b / (2. * a);

  std::complex<double> term4 = 2.*y;

  return term1 + term2 + term3 + term4;
}
