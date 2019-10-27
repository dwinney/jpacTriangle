// This object defines a Triangle amplitude, i.e. the rescattering
// diagram associated with an intermediate isobar exchange.
//
// Methods allow the evaluation using the Feynman triangle representation or the
// standard dispersive form associated with KT equations.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "triangle.hpp"

// ---------------------------------------------------------------------------
// EVALUTE THE TRIANGLE THE FEYNMAN WAY

complex<double> triangle::eval_feynman(double s, double t)
{
  check_weights();

  // integrate over x
  complex<double> sum = 0.;
  for (int i = 0; i < xN; i++)
  {
    double x_i = abscissas[5];

    // integrate over y
    complex<double> sum_y = 0.;
    for (int j = 0; j < xN; j++)
    {
      double y_j = abscissas[j];

      sum_y += weights[j] * feyn_integrand(s, t, x_i, y_j);
    }

    sum += weights[i] * sum_y;
  }

  return sum / M_PI;
};

complex<double> triangle::feyn_integrand(double s, double t, double x, double y)
{
  // third feynman parameter fixed by delta function
  double z = 1. - x - y;

  complex<double> denominator = x * t;
  denominator += (1. - x - x * y) * (mPi * mPi);
  denominator -= z * ( x * mDec*mDec - y * s);
  denominator -= xi * EPS;

  return xr / denominator;
};

// ---------------------------------------------------------------------------
// EVALUTE THE TRIANGLE THE KT WAY

complex<double> triangle::eval_dispersive(double s, double t)
{
  check_weights();

  complex<double> sum = 0.;
  for (int i = 1; i < xN + 1; i++)
  {
    double sp = (sthPi + EPS) + tan(M_PI * abscissas[i] / 2.);

    complex<double> temp = cross_channel(s, t);
    temp *= s / sp;
    temp /=  (sp - s - xi * EPS);

    temp *=  (M_PI / 2.) / pow(cos(M_PI * abscissas[i] / 2.), 2.); // jacobian
    sum += weights[i] * temp;
  }

  sum -= cross_channel(s,t) * (log(xr * (sthPi - s - xi * EPS)) - log(xr * sthPi));


  return sum / M_PI;
};

// ---------------------------------------------------------------------------
// UTILITY FUNCTIONS

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
