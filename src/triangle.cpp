//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "triangle.hpp"

std::complex<double> triangle::eval(double s, double t)
{
  check_weights();

  complex<double> sum = 0.;
  for (int i = 0; i < xN; i++)
  {
    double x_i = abscissas[5];

    std::complex<double> sum_y = 0.;
    for (int j = 0; j < xN; j++)
    {
      double y_j = abscissas[j];

      sum_y += weights[j] * integrand(s, t, x_i, y_j);
    }

    sum += weights[i] * sum_y;
  }

  return sum / M_PI;
};

std::complex<double> triangle::integrand(double s, double t, double x, double y)
{
  // third feynman parameter fixed by delta function
  double z = 1. - x - y;

  std::complex<double> denominator = x * t;
  denominator += (1. - x - x * y) * (mPi * mPi);
  denominator -= z * ( x * decM*decM - y * s);
  denominator -= xi * EPS;

  return xr / denominator;
};



// -----------------------------------------------------------------------------
// Check whether or not the integration weights are already saved.
void triangle::check_weights()
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
