// Kernel functions for the triangle convolution integral. Higher ordered
// kernels correspond to higher number of subtractions.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman/triangle_function.hpp"

// ---------------------------------------------------------------------------
// This is the final triangle function.
// Combines appropriate combinations of basis functions to correspond to exchanges with spin
std::complex<double> triangle_function::eval(int n, int j, int jp, double s, double tp)
{
  std::complex<double> result;

  switch (jp)
  {
    // scalar exchange
    case 0: return mT(n, 1, j, s, tp) / tp;

    // vector exchange
    case 1:
    {
      result = mT(n, 1, 0, s, tp);
      result += 2. * s * mT(n - 1, 0, j, s, tp);
      result -= (mDec2 + 3. * mPi2) * mT(n, 0, j, s, tp);
      break;
    }
    default:
    {
      std::cout << "\n j and jp combination not available! Quitting... \n";
      exit(0);
    }
  }

  return result;
};

// ---------------------------------------------------------------------------
// Takes selected basis function (specified by orders of k^2 and projection order j)
// and subtracts it n times at s = 0
// Also actually does the numerical integration over the final parameter x
std::complex<double> triangle_function::mT(int n, int k, int j, double s, double t)
{
  auto F = [&](double x)
  {
    return mT_integrand(n, k, j, s, t, x);
  };

  double error;
  std::complex<double> result;
  result = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(F, 0, 1, 5, 1.E-9, &error);

  return result / M_PI;
};

// ---------------------------------------------------------------------------
// Basis functions based on order of divergence.
// Filters number of powers of k^2 in the numerator of the triangle function
std::complex<double> triangle_function::mT_integrand(int n, int k, int j, double s, double t, double x)
{
  switch (k)
  {
    case 0:
    {
      if (n == 0)
      {
        return int_mT0(j, s, t, x);
      }
      else if (n == 1)
      {
        return int_mT0(j, s, t, x) - int_mT0(j, 0., t, x);
      }
    }
    case 1:
    {
      if (n == 0)
      {
        return int_mT1(j, s, t, x);
      }
      else if (n == 1)
      {
        return int_mT1(j, s, t, x) - int_mT1(j, 0., t, x);
      }
    }
    default:
    {
      std::cout << "\n Error finding mT integrand. Quitting... \n";
      exit(0);
    }
  }
};


// ---------------------------------------------------------------------------
// Vanilla Triangle function: mT_j{1}(s,t)
// ---------------------------------------------------------------------------

std::complex<double> triangle_function::int_mT0(int j, double s, double t, double x)
{
  std::complex<double> a, b, c, d;
  std::complex<double> e, f, g;

  // coeffs of denominator polynomial
  a = mPi2;
  b = mPi2 + (x - 1.) * mPi2 + x * mDec2 - x * s - t;
  c = (1. - x) * t + x * mPi2 + x*(x-1.)* mDec2 - ieps;

  switch (j)
  {
    case 0:
    {
      // 1
      g = 0.;
      f = 0.;
      e = 1.;
      break;
    }
    case 1:
    {
      // z
      g = 0.;
      f = -1.;
      e = 1. - x;
      break;
    }
    default:
    {
      std::cout << " \n mT0: only j <= 1 implemented. Quiting... \n";
      exit(0);
    }
  }

  // Evaluate definite integral with bounds y = [0, 1-x]
  std::complex<double> result;
  result = ri_poly2(1. - x, a, b, c, e, f, g) - ri_poly2(0., a, b, c, e, f, g);

  return result;
};

// ---------------------------------------------------------------------------
// mT_j{k^2}(s,t)
// ---------------------------------------------------------------------------

std::complex<double> triangle_function::int_mT1(int j, double s, double t, double x)
{
  std::complex<double> a, b, c, d; // Denominator polynomial
  std::complex<double> e, f, g, h, i; // Numerator polynomial
  std::complex<double> l, n, m; // Log polynomial

  // coeffs of denominator polynomial (always the same)
  a = mPi2;
  b = mPi2 + (x - 1.) * mPi2 + x * mDec2 - x * s - t;
  c = (1. - x) * t + x * mPi2 + x*(x-1.)* mDec2 - ieps;

  // the relative momenta (pX - p3)^2 / 4
  double p = (mDec2 + mPi2) / 2. - s / 4.;

  std::complex<double> result = 0.;
  switch (j)
  {
    case 0:
    {
      // (x+y)^2
      g = 1.;
      f = 2. * x;
      e = x * x;

      result = p * (ri_poly2(1. - x, a, b, c, e, f, g) - ri_poly2(0., a, b, c, e, f, g));
      result -= 2. * (ri_log0(1. - x, a, b, c) - ri_log0(0., a, b, c));
      break;
    }
    case 1:
    {
      // z (x+y)^2
      i = 0.;
      h = -1.; // problem child
      g = 1. - 3.*x;
      f = 2.*x - 3.*x*x;
      e = x*x* (1. - x);

      // z
      m = 0.;
      n = -1.; // problem child
      l = 1. - x;

      result = p * (ri_poly4(1. - x, a, b, c, e, f, g, h, i) - ri_poly4(0., a, b, c, e, f, g, h, i));
      result -= 2. * (ri_log2(1. - x, a, b, c, l, n, m) - ri_log2(0., a, b, c, l, n, m));
      break;
    }
    default:
    {
      std::cout << " \n mT1: only j <= 1 implemented. Quiting... \n";
      exit(0);
    }
  }

  return result;
};
