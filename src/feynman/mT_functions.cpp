// Kernel functions for the triangle convolution integral. Higher ordered
// kernels correspond to higher number of subtractions.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman_triangle.hpp"

// ---------------------------------------------------------------------------
// Triangle kernels
// built from basis functions in appropriate combinations for spin polynomials
std::complex<double> feynman_triangle::triangle_kernel(int n, int j, int jp, double s, double tp)
{
  std::complex<double> result;

  switch (j)
  {
    // s - wave
    case 0:
    {
        switch (jp)
        {
          // scalar exchange
          case 0: return mT(n, 0, s, tp);

          // vector exchange
          case 1:
          {
            result = mT(n, 1, s, tp);
            result += 2. * s * mT(n - 1, 0, s, tp) - (mDec2 + 3. * mPi2) * mT(n, 0, s, tp);
            break;
          }
        }
        break;
    }

    // // p - wave
    // case 1:
    // {
    //   std::complex<double> result;
    //   result = 2. * mT1(s,t) + (s - mDec2 - 3.*mPi2) * mT0(s,t);
    //
    //   return result;
    // }
  default: std::cout << "\n j and jp combination not available! Quitting... \n"; exit(0);
  }

  return result;
};

// Takes selected basis function and subtracts it n times
std::complex<double> feynman_triangle::mT(int n, int m, double s, double t)
{
  check_weights();

  //integrate over x
  std::complex<double> sum = 0.;
  for (int i = 0; i < xN; i++)
  {
    double x_i = abscissas[i];

    switch (n)
    {
      // One subtraction
      case 1:
      {
        sum += weights[i] * (mT_integrand(m, s, t, x_i) - mT_integrand(m, 0., t, x_i));
        break;
      }

      // No subtractions
      case 0:
      {
        sum += weights[i] * (mT_integrand(m, s, t, x_i));
        break;
      }
      default:
      {
        std::cout << "\n Error! Weird number of subtractions: " << n << " recieved. Quitting... \n";
        exit(0);
      }
    }
  }

  return sum / M_PI;
};

// Basis functions based on order of divergence.
// Filters number of powers of k^2 in the numerator of the triangle function
std::complex<double> feynman_triangle::mT_integrand(int m, double s, double t, double x)
{
  switch (m)
  {
    case 0:
    {
      return mT0_yintegrand(s, t, x);
    }
    case 1:
    {
      return mT1_yintegrand(s, t, x);
    }
    default:
    {
      std::cout << "\n Error finding mT integrand. Quitting... \n";
      exit(0);
    }
  }
};


// ---------------------------------------------------------------------------
// Vanilla Triangle function (no powers of k^2 in numerator)
// ---------------------------------------------------------------------------

std::complex<double> feynman_triangle::mT0_yintegrand(double s, double t, double x)
{
  std::complex<double> a, b, c, d;

  // coeffs of denominator polynomial
  a = mPi2;
  b = mPi2 + (x - 1.) * mPi2 + x * mDec2 - x * s - t;
  c = (1. - x) * t + x * mPi2 + x*(x-1.)* mDec2 - ieps;

  // Evaluate definite integral with bounds y = [0, 1-x]
  std::complex<double> result;
  result = ri_poly1(1. - x, a, b, c) - ri_poly1(0., a, b, c);

  return result;
};


// ---------------------------------------------------------------------------
// k^2 in the numerator
// ---------------------------------------------------------------------------

std::complex<double> feynman_triangle::mT1_yintegrand(double s, double t, double x)
{
  std::complex<double> a, b, c, d;
  std::complex<double> e, f, g;

  // coeffs of denominator polynomial
  a = mPi2;
  b = mPi2 + (x - 1.) * mPi2 + x * mDec2 - x * s - t;
  c = (1. - x) * t + x * mPi2 + x*(x-1.)* mDec2 - ieps;

  // coeffs of numerator polynomial
  e = mPi2;
  f = x * (mDec2 + mPi2 - s);
  g = x*x * mDec2;

  // Evaluate definite integral with bounds y = [0, 1-x]
  std::complex<double> result;
  result = ri_poly2(1. - x, a, b, c, e, f, g) - ri_poly2(0., a, b, c, e, f, g);
  result -= 2. * (ri_log1(1. - x, a, b, c) - ri_log1(0., a, b, c));

  return result;
};
