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

// ---------------------------------------------------------------------------
// Evaluate the convolution of the LHC function with triangle function
std::complex<double> feynman_triangle::eval(int n, int j, double s)
{
    check_weights();

    // Calculate the dispersion piece coming from the triangle
    std::complex<double> sum = 0.;
    for (int i = 0; i < xN; i++)
    {
      double tp = r_thresh + tan(M_PI * abscissas[i] / 2.);

      std::complex<double> temp;
      temp = lhc_func->disc(tp);
      temp *= triangle_kernel(n, j, s, tp);

      temp *= (M_PI / 2.);
      temp /= pow(cos(M_PI * abscissas[i] / 2.), 2.); // jacobian

      sum += weights[i] * temp;
    }
    sum /= M_PI;

    // if needed add the polynomial contribution
    if (n == 1)
    {
      sum += 1. + mP1(s);
    }

    return sum;
};

// ---------------------------------------------------------------------------
// Triangle kernel which encodes all the spin stuff
std::complex<double> feynman_triangle::triangle_kernel(int n, int j, double s, double t)
{
  switch (n)
  {
    // unsubtracted
    case 0:
    {
      switch (j)
      {
        // s - wave
        case 0: return mT0(s, t);

        // p - wave
        case 1:
        {
          std::complex<double> result;
          result = 2. * mT1(s,t) + (s - mDec2 - 3.*mPi2) * mT0(s,t);

          return result;
        }
      }
    }
    // once-subtracted
    case 1:
    {
      return mT1(s, t) / t;
    }
    default: std::cout << "Number of subtractions not implemented yet! Quitting... \n"; exit(0);
  }
};

// ---------------------------------------------------------------------------
// UTILITY FUNCTIONS

// ---------------------------------------------------------------------------
// Kallen triangle function
std::complex<double> feynman_triangle::Kallen(double x, double y, double z)
{
  return x * x + y * y + z * z - 2. * (x * z + y * z + x * y);
};

// ---------------------------------------------------------------------------
// Kacser function which includes the correct analytic structure of
// product of breakup momenta, 4 * p(s) * q(s)
std::complex<double> feynman_triangle::Kacser(double s)
{
  std::complex<double> result;

  result = sqrt(pow(sqrt(s) + mPi, 2.) - mDec2 - ieps);
  result *= sqrt(pow(sqrt(s) - mPi, 2.) - mDec2 - ieps);
  result *= sqrt(Kallen(s, mPi2, mPi2)) / s;

  return result;
};

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
