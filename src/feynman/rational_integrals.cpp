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
std::complex<double> ri_poly1(double y,
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
  term1 = c_atan((2. * a * y + b) / sqrt(-d));
  term1 *= b*b*e - a*f*b + 2.*a*(a*g-c*e);
  term1 /= a*a * sqrt(-d);

  std::complex<double> term2;
  term2 = log(a*y*y + b*y + c * xr);
  term2 *= a*f - b*e;
  term2 /= 2.*a*a;

  std::complex<double> term3;
  term3 = e * y / a;

  return term1 + term2 + term3;
};

// ---------------------------------------------------------------------------
// Antiderivative of (e y^4 + f y^3 + g y^2 + h y + i) / (a y ^2 + b y + c)
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
  std::complex<double> d = 4. * a * c - b * b; // discriminant

  std::complex<double> term1;
  term1 = -0.5 * log(a*y*y + b*y + c);
  term1 *= a*a*c*f - a*b*b*f - 2.*a*b*c*e + b*b*b*e;
  term1 /= a*a*a*a;

  std::complex<double> term2;
  term2 = c_atan((2.*a*y+b) / sqrt(d));
  term2 *= 3.*a*a*b*c*f + 2.*a*a*c*c*e - a*b*b*b*f - 4.*a*b*b*c*e + b*b*b*b*e;
  term2 /= a*a*a*a * sqrt(d);

  std::complex<double> term3;
  term3 = a*y;
  term3 *= -3.*a*b*(e*y + 2.*f) + a*a*y*(2.*e*y+3.*f) - 6.*a*c*e + 6.*b*b*e;
  term3 /= 6.* a*a*a*a;

  return term1 + term2 + term3 + ri_poly2(y, a, b, c, g, h, i);
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
  term1 = - c_atan((2.*a*y + b) / sqrt(-d * xr));
  term1 *= sqrt(-d * xr) / a;

  std::complex<double> term2;
  term2 = y * log(xr / (a*y*y + b*y + c * xr));

  std::complex<double> term3;
  term3 = log(a*y*y + b*y + c);
  term3 *= -b / (2. * a);

  std::complex<double> term4 = 2.*y;

  return term1 + term2 + term3 + term4;
};

// ---------------------------------------------------------------------------
// Antiderivative of (e y^2 + f y + g) * Log(1/(a y^2 + b y + c))
std::complex<double> ri_log2(double y,
    std::complex<double> a,
    std::complex<double> b,
    std::complex<double> c,
    std::complex<double> e,
    std::complex<double> f,
    std::complex<double> g)
{
    std::complex<double> d = 4. * a * c - b * b; // discriminant

    std::complex<double> term1;
    term1 = 6. * a*a*a*y*y*(2.*e*y + 3.*f);
    term1 *= log(xr / (a*y*y + b*y + c * xr));
    term1 /= 36.*a*a*a;

    std::complex<double> term2;
    term2 = -3. * (6.*a*a*c*f - 3.*a*b*b*f - 6.*a*b*c*e + 2.*b*b*b*e);
    term2 *= log(a*y*y + b*y + c);
    term2 /= 36.*a*a*a;

    std::complex<double> term3;
    term3 = 2.*a*y;
    term3 *= -3.*a*b*(e*y+3.*f) + a*a*y*(4.*e*y+9.*f) - 12.*c*e*a + 6. * b*b*e;
    term3 /= 36.*a*a*a;

    std::complex<double> term4;
    term4 = -.6 * sqrt(d) * (-3.*a*b*f - 2.*a*c*e + 2.*b*b*e);
    term4 *= c_atan((2.*a*y + b) / sqrt(d));
    term4 /= 36.*a*a*a;

    return term1 + term2 + term3 + term4 + g * ri_log0(y,a,b,c);
};
