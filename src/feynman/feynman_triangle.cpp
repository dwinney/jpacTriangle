// This object defines a Triangle amplitude, i.e. the rescattering
// diagram associated with an intermediate t-channel exchange.
//
// Evaluated specifically with the Feynman method
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman_triangle.hpp"

// Evaluate the convolution of the LHC function with triangle function
std::complex<double> feynman_triangle::eval(double s)
{
    check_weights();

    return mT1(s, .770);

    // std::complex<double> sum = 0.;
    // for (int i = 0; i < xN; i++)
    // {
    //   double tp = r_thresh + tan(M_PI * abscissas[i] / 2.);
    //
    //   std::complex<double> temp;
    //   temp = lhc_func->disc(tp) / tp;
    //   temp *= mT1(s, tp);
    //   temp *= (M_PI / 2.);
    //   temp /= pow(cos(M_PI * abscissas[i] / 2.), 2.); // jacobian
    //
    //   sum += weights[i] * temp;
    // }
    //
    // sum /= M_PI;
    //
    // return sum;
};

// Triangle function from the perturbation theory result
std::complex<double> feynman_triangle::kernel(double s, double t)
{
  check_weights();

  // integrate over x
  std::complex<double> sum = 0.;
  for (int i = 0; i < xN; i++)
  {
    double x_i = abscissas[i];
    sum += weights[i] * kernel_integrand(s, t, x_i);
  }

  return sum;
};

// Logarithm from integrating over y and z
std::complex<double> feynman_triangle::kernel_integrand(double s, double t, double x)
{
  std::complex<double> a, b, c, d;

  a = mPi2;
  b = mPi2 + (x - 1.) * mPi2 + x * (mDec2 + ieps) - x * s - t;
  c = (1. - x) * t + x * mPi2 + x*(x-1.)* (mDec2 + ieps);
  d = b * b - 4. * a * c; // discriminant

  // Roots of the polynomial
  std::complex<double> y_plus = (-b + sqrt(xr * d)) / (2. * a);
  std::complex<double> y_minus = (-b - sqrt(xr * d)) / (2. * a);

  std::complex<double> result;
  result = log(y_plus + x - xr) - log(y_minus + x - xr);
  result -= log(y_plus) - log(y_minus);
  result /= sqrt(xr * d);

  return result / M_PI;
};

// ---------------------------------------------------------------------------
// UTILITY FUNCTIONS

// -----------------------------------------------------------------------------
// Check whether or not the integration weights are already saved.
void feynman_triangle::check_weights()
{
  if (WG_GENERATED == false)
  {

    weights.clear(); abscissas.clear();

    double w[xN + 1], a[xN + 1];
    gauleg(0., 1., a, w, xN + 1);

    for (int i = 1; i < xN + 1; i++)
    {
      weights.push_back(w[i]);
      abscissas.push_back(a[i]);
    }

    if (weights.size() != xN || abscissas.size() != xN)
    {
      cout << "ERROR: wrong number of weights generated for some reason. Quitting... \n";
      exit(0);
    }
    else
    {
      WG_GENERATED = true;
    }
  }
};
