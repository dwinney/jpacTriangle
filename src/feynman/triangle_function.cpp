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
// Combines appropriate combinations of basis functions to correspond to spin exchanges
std::complex<double> triangle_function::eval(int n, int j, int jp, double s, double tp)
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
          case 0: return mT(n, 0, 0, s, tp);

          // vector exchange
          case 1:
          {
            result = mT(n, 1, 0, s, tp);
            result += 2. * s * mT(n - 1, 0, j, s, tp) - (mDec2 + 3. * mPi2) * mT(n, 0, j, s, tp);
            break;
          }
        }
        break;
    }

    // p - wave
    case 1:
    {
      switch (jp)
      {
        // scalar exchange
        case 0: return mT(n, 0, j, s, tp);
      }
      break;
    }

    // d - wave projection
    case 2:
    {
      if (jp == 0)
      {
        return mT(n, 0, j, s, tp);
      }
    }
  default: std::cout << "\n j and jp combination not available! Quitting... \n"; exit(0);
  }

  return result;
};

// ---------------------------------------------------------------------------
// Takes selected basis function (specified by orders of k^2 and z) and subtracts it n times
std::complex<double> triangle_function::mT(int n, int k, int z, double s, double t)
{
  integ.check_weights();

  //integrate over x
  std::complex<double> sum = 0.;
  for (int i = 0; i < integ.xN; i++)
  {
    double x_i = integ.abscissas[i];

    switch (n)
    {
      // One subtraction
      case 1:
      {
        sum += integ.weights[i] * (mT_integrand(k, z, s, t, x_i) - mT_integrand(k, z, 0., t, x_i));
        break;
      }

      // No subtractions
      case 0:
      {
        sum += integ.weights[i] * (mT_integrand(k, z, s, t, x_i));
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
std::complex<double> triangle_function::mT_integrand(int k, int z, double s, double t, double x)
{
  switch (k)
  {
    case 0:
    {
      return int_mT0(z, s, t, x);
    }
    case 1:
    {
      if (z == 0)
      {
      return int_mT1(z, s, t, x);
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

std::complex<double> triangle_function::int_mT0(int z, double s, double t, double x)
{
  std::complex<double> a, b, c, d;
  std::complex<double> e, f, g;

  // coeffs of denominator polynomial
  a = mPi2;
  b = mPi2 + (x - 1.) * mPi2 + x * mDec2 - x * s - t;
  c = (1. - x) * t + x * mPi2 + x*(x-1.)* mDec2 - ieps;

  switch (z)
  {
    case 0:
    {
      e = 0.;
      f = 0.;
      g = 1.;
      break;
    }
    case 1:
    {
      e = 0.;
      f = -1.;
      g = 1. - x;
      break;
    }
    case 2:
    {
      e = 1.;
      f = 2. *(x-1.);
      g = x*x - 2.*x + 1.;
      break;
    }
    default:
    {
      std::cout << " \n mT0: only j <= 2 implemented. Quiting... \n";
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

std::complex<double> triangle_function::int_mT1(int z, double s, double t, double x)
{
  std::complex<double> a, b, c, d; // Denominator polynomial
  std::complex<double> e, f, g, h, i; // Numerator polynomial
  std::complex<double> j, k, l; // Log polynomial

  // coeffs of denominator polynomial
  a = mPi2;
  b = mPi2 + (x - 1.) * mPi2 + x * mDec2 - x * s - t;
  c = (1. - x) * t + x * mPi2 + x*(x-1.)* mDec2 - ieps;

  switch (z)
  {
    case 0:
    {
      e = 0.;
      f = 0.;
      g = mPi2;
      h = x * (mDec2 + mPi2 - s);
      i = x*x * mDec2;

      j = 0.;
      k = 0.;
      l = 1.;
      break;
    }
    default:
    {
      std::cout << " \n mT1: only j <= 1 implemented. Quiting... \n";
      exit(0);
    }
  }

  // Evaluate definite integral with bounds y = [0, 1-x]
  std::complex<double> result;
  result = ri_poly4(1. - x, a, b, c, e, f, g, h, i) - ri_poly4(0., a, b, c, e, f, g, h, i);
  result -= 2. * (ri_log2(1. - x, a, b, c, j, k, l) - ri_log2(0., a, b, c, j, k, l));

  return result;
};
