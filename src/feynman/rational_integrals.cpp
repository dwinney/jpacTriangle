// Analytic forms of integrals over rational functions which show up in Feynman integrals
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman/rational_integrals.hpp"

// redefine arctangent in terms of logarithms for complex argument to better control
// singularity structure
std::complex<double> c_atan(std::complex<double> z)
{
  std::complex<double> result = log((z + xi) / (z - xi));
  result /= -2.*xi;

  return result;
};

// ---------------------------------------------------------------------------
// Integral over 1/(a y^2 + b y + c)
std::complex<double> ri_poly0(double y,
  std::complex<double> a,
  std::complex<double> b,
  std::complex<double> c)
{
  std::complex<double> d = b*b - 4. * a * c; // discriminant

  std::complex<double> result = c_atan((2.*a*y+b) / sqrt(-d*xr));;
  result *= 2. / sqrt(-d*xr);

  return result;
};

// ---------------------------------------------------------------------------
// Integral over (e + f * y + g * y^2) / (a y ^2 + b y + c)
std::complex<double> ri_poly2(double y,
  std::complex<double> a,
  std::complex<double> b,
  std::complex<double> c,
  std::complex<double> e,
  std::complex<double> f,
  std::complex<double> g)
{
  std::complex<double> atan_term;
  atan_term  = 2. * c_atan((2.*a*y+b) / sqrt(4.*a*c - b*b));
  atan_term /= sqrt(4.*a*c - b*b);
  atan_term *= -a*(b*f+2.*c*g) + b*b*g;

  std::complex<double> log_term;
  log_term  = log(a*y*y + b*y + c);
  log_term *= a*f - b*g;

  std::complex<double> poly_term;
  poly_term  = 2.*a*g*y;

  std::complex<double> result = atan_term + log_term + poly_term;
  result /= 2.*a*a;
  result += e * ri_poly0(y, a, b, c);

  return result;
};

// ---------------------------------------------------------------------------
// Antiderivative of (e + f * y + g * y^2 + h * y^3 + i * y^4) / (a y ^2 + b y + c)
std::complex<double> ri_poly4(double y,
    std::complex<double> a,
    std::complex<double> b,
    std::complex<double> c,
    std::complex<double> e,
    std::complex<double> f,
    std::complex<double> g,
    std::complex<double> h,
    std::complex<double> i)
{

  std::complex<double> atan_term;
  atan_term  = c_atan((2.*a*y+b) / sqrt(4.*a*c*xr - b*b));
  atan_term *= 6.*(a*a*c*(3.*b*h+2.*c*i) - a*b*b*(b*h+4.*c*i) + i*b*b*b*b);
  atan_term /= sqrt(4.*a*c - b*b);

  std::complex<double> log_term;
  log_term  = log(a*y*y + b*y + c);
  log_term *= -3. * (a*a*c*h - a*b*(b*h+2.*c*i) + b*b*b*i);

  std::complex<double> poly_term;
  poly_term  = 2.*a*a*a*i * y*y*y;
  poly_term += 3.*a*a*(a*h - b*i) * y*y;
  poly_term += -6. * a*(a*b*h+a*c*i-i*b*b) * y;

  std::complex<double> result = poly_term;
  result += log_term;
  result += atan_term;
  result /= 6.*a*a*a*a;

  result += ri_poly2(y, a, b, c, e, f, g);

  return result;
};

// ---------------------------------------------------------------------------
// Integral over Log(1/(a y^2 + b y + c))
std::complex<double> ri_log0(double y,
    std::complex<double> a,
    std::complex<double> b,
    std::complex<double> c)
{
  std::complex<double> d = b * b - 4. * a * c; // discriminant

  std::complex<double> term1;
  term1  = - c_atan((2.*a*y + b) / sqrt(-d * xr));
  term1 *= sqrt(-d * xr) / a;

  std::complex<double> term2;
  term2 = y * log(xr / (a*y*y + b*y + c * xr));

  std::complex<double> term3;
  term3  = log(a*y*y + b*y + c * xr);
  term3 *= -b / (2. * a);

  std::complex<double> term4 = 2.*y;

  return term1 + term2 + term3 + term4;
};

// ---------------------------------------------------------------------------
// Antiderivative of (e + f * y + g * y^2) * Log(1/(a y^2 + b y + c))
std::complex<double> ri_log2(double y,
    std::complex<double> a,
    std::complex<double> b,
    std::complex<double> c,
    std::complex<double> e,
    std::complex<double> f,
    std::complex<double> g)
{
    std::complex<double> log_term;
    log_term  = -6.*a*a*a*y*y * (3.*f + 2.*g*y);
    log_term += -3.*(6.*a*a*c*f - 3.*a*b*(b*f + 2.*c*g) + 2.*b*b*b*g);
    log_term *= log(a*y*y + b*y + c);

    std::complex<double> atan_term;
    atan_term  = -6. * c_atan((2.*a*y+b) / sqrt(4.*a*c*xr - b*b));
    atan_term *= sqrt(4.*a*c - b*b);
    atan_term *= -3.*a*b*f - 2.*a*c*g + 2.*b*b*g;

    std::complex<double> poly_term;
    poly_term  = 2.*a*y;
    poly_term *= a*a*y*(9.*f+ 4.*g*y) - 3.*a*(3.*b*f+b*g*y + 4.*c*g) + 6.*b*b*g;

    std::complex<double> result = log_term + poly_term + atan_term;
    result /= 36.*a*a*a;

    result += e * ri_log0(y, a, b, c);

    return result;
};
