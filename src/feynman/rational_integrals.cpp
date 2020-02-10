// Analytic forms of integrals over rational functions which show up in Feynman integrals
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "rational_integrals.hpp"

// // ---------------------------------------------------------------------------
// // Integral over 1/(a y^2 + b y + c)
// std::complex<double> ri_poly1(double y,
//   std::complex<double> a,
//   std::complex<double> b,
//   std::complex<double> c)
// {
//   std::complex<double> d = b * b - 4. * a * c; // discriminant
//
//   // Roots of the polynomial
//   std::complex<double> y_plus = (-b + sqrt(xr * d)) / (2. * a);
//   std::complex<double> y_minus = (-b - sqrt(xr * d)) / (2. * a);
//
//   std::complex<double> result;
//   result = log(y_plus + x - xr) - log(y_minus + x - xr);
//   result -= log(y_plus) - log(y_minus);
//   result /= sqrt(xr * d);
//
//   return result;
// };

// ---------------------------------------------------------------------------
// Integral over (e y^2 + f y + g) / (a y ^2 + b y + c)
std::complex<double> ri_poly2(double y,
  std::complex<double> a,
  std::complex<double> b,
  std::complex<double> c,
  std::complex<double> e,
  std::complex<double> f,
  std::complex<double> g)
{
  std::complex<double> d = b * b - 4. * a * c; // discriminant

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

// ---------------------------------------------------------------------------
// Integral over Log(1/(a y^2 + b y + c))
std::complex<double> ri_log1(double y,
    std::complex<double> a,
    std::complex<double> b,
    std::complex<double> c)
{
  std::complex<double> d = b * b - 4. * a * c; // discriminant

  std::complex<double> term1;
  term1 = - atan((2.*a*y + b) / sqrt(-d * xr));
  term1 *= sqrt(-d * xr) / a;

  std::complex<double> term2;
  term2 = y * log(xr / (a*y*y + b*y + c - ieps));

  std::complex<double> term3;
  term3 = log(a*y*y + b*y + c + ieps);
  term3 *= -b / (2. * a);

  std::complex<double> term4 = 2.*y;

  return term1 + term2 + term3 + term4;
};
