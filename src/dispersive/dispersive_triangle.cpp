// This object defines a Triangle amplitude, i.e. the rescattering
// diagram associated with an intermediate t-channel exchange.
//
// Evaluated specifically with the dispersion relation method
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dispersive_triangle.hpp"

std::complex<double> dispersive_triangle::eval(int j, int jp, double s)
{
  // if pseudo threshold is in the bounds of integration, exclude a small interval around it
  // this is alwys the case for decays, but to keep the code general I consider the other case
  if (p_thresh < s_thresh)
  {
    return 1. + s_dispersion_inf(j, jp, s, s_thresh);
  }
  else
  {
    std::complex<double> temp;
    temp = 1. + s_dispersion(j, jp, s, s_thresh + EPS, p_thresh - exc);
    temp += s_dispersion_inf(j, jp, s, p_thresh + exc);

    return temp;
  }
};

// ---------------------------------------------------------------------------
// calculate the dispersion integral over s with finite bounds of integration
std::complex<double> dispersive_triangle::s_dispersion(int j, int jp, double s, double low, double high)
{
  double w[xN + 1], x[xN + 1];
  gauleg(low, high, x, w, xN);

  // Integrate
  std::complex<double> sum = 0.;
  for (int i = 1; i <= xN; i++)
  {
    double sp = x[i];
    std::complex<double> temp;

    temp = rho(sp) * b(j, jp, sp) - rho(s) * b(j, jp, s);
    temp *= s / sp;
    temp /= (sp - s - ieps);

    sum += w[i] * temp;
  }

  // Log term from subtracted singularity
  std::complex<double> log_term;
  log_term = log(high - s * xr) - log(high);
  log_term -= log(low - s * xr) - log(low);
  log_term *= rho(s) * b(j, jp, s);

  return (sum + log_term) / M_PI;
};

// calculate the dispersion integral over s up to infinity
std::complex<double> dispersive_triangle::s_dispersion_inf(int j, int jp, double s, double low)
{
  check_weights();

  // Integrate
  std::complex<double> sum = 0.;
  for (int i = 0; i < xN; i++)
  {
    double sp = low + tan(M_PI * abscissas[i] / 2.);

    std::complex<double> temp;
    temp = rho(sp) * b(j, jp, sp) - rho(s) * b(j, jp, s);
    temp *= s / sp;
    temp /= (sp - s - ieps);
    temp *= (M_PI / 2.) / pow(cos(M_PI * abscissas[i] / 2.), 2.); // jacobian

    sum += weights[i] * temp;
  }

  // subtracted point
  std::complex<double> log_term = - rho(s) * b(j, jp, s);
  log_term *= log(low - s * xr) - log(low);

  return (sum + log_term) / M_PI;
};

// ---------------------------------------------------------------------------
// calculate the dispersion intgral over t
std::complex<double> dispersive_triangle::b(int j, int jp, double s)
{
  check_weights();

  std::complex<double> sum = 0.;
  for (int i = 0; i < xN; i++)
  {
    double tp = (r_thresh + EPS) + tan(M_PI * abscissas[i] / 2.);

    std::complex<double> temp;
    temp = lhc_func->disc(tp);
    temp *= projector(0, j, jp, s, tp);
    temp *= (M_PI / 2.) / pow(cos(M_PI * abscissas[i] / 2.), 2.);

    sum += weights[i] * temp;
  }

  sum /= M_PI;

  return sum;
};

// ---------------------------------------------------------------------------
// Kacser function which includes the correct analytic structure of
// product of breakup momenta, p(s) * q(s)
std::complex<double> dispersive_triangle::Kacser(double s)
{
  std::complex<double> result;

  result = sqrt(pow(sqrt(s) + mPi, 2.) - mDec2 - ieps);
  result *= sqrt(pow(sqrt(s) - mPi, 2.) - mDec2 - ieps);
  result *= rho(s);

  return result;
};

// Two particle phase-space function
std::complex<double> dispersive_triangle::rho(double s)
{
  return sqrt(Kallen(s, mPi2, mPi2)) / s;
};

// ---------------------------------------------------------------------------
// complex Bounds of integration
std::complex<double> dispersive_triangle::t_minus(double s)
{
  return (mDec2 + ieps) + mPi2 - (s + mDec2 + ieps - mPi2) / 2. - Kacser(s) / 2.;
};

std::complex<double> dispersive_triangle::t_plus(double s)
{
  return (mDec2 + ieps) + mPi2 - (s + mDec2 + ieps - mPi2) / 2. + Kacser(s) / 2.;
};

// ---------------------------------------------------------------------------
// UTILITY FUNCTIONS

// -----------------------------------------------------------------------------
// Check whether or not the integration weights are already generated
void dispersive_triangle::check_weights()
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
