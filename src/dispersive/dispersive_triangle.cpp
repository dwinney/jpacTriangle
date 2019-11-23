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

std::complex<double> dispersive_triangle::eval(int jp, double s)
{
  // if pseudo threshold is in the bounds of integration, exclude a small interval around it
  // this is alwys the case for decays, but to keep the code general I consider the other case
  if (p_thresh < s_thresh)
  {
    return 1. + s_dispersion_inf(jp, s, s_thresh);
  }
  else
  {
    std::complex<double> temp;
    temp = s_dispersion(jp, s, s_thresh + EPS, p_thresh - exc);
    temp += s_dispersion_inf(jp, s, p_thresh + exc);

    return 1. + temp;
  }
};

// ---------------------------------------------------------------------------
// calculate the dispersion integral over s with finite bounds of integration
std::complex<double> dispersive_triangle::s_dispersion(int jp, double s, double low, double high)
{
  double w[xN + 1 ], x[xN + 1];
  gauleg(low, high, x, w, xN);

  // Subtract off the pole at s = sp
  std::complex<double> sub_point = phase_space(s) * t_dispersion(jp, s);

  std::complex<double> sum = 0.;
  for (int i = 1; i <= xN; i++)
  {
    std::complex<double> temp;
    temp = phase_space(x[i]) * t_dispersion(jp, x[i]) - sub_point;
    temp *= s / x[i];
    temp /= (x[i] - s - ieps);

    sum += w[i] * temp;
  }

  std::complex<double> log_term;
  log_term = log(high - s - ieps) - log(high);
  log_term -= log(low - s - ieps) - log(low);
  log_term *= sub_point;

  return (sum + log_term) / M_PI;
};

// calculate the dispersion integral over s up to infinity
std::complex<double> dispersive_triangle::s_dispersion_inf(int jp, double s, double low)
{
  check_weights();

  std::complex<double> sub_point = phase_space(s) * t_dispersion(jp, s);

  std::complex<double> sum = 0.;
  for (int i = 0; i < xN; i++)
  {
    double sp = low + tan(M_PI * abscissas[i] / 2.);

    std::complex<double> temp;
    temp = phase_space(s) * t_dispersion(jp, sp) - sub_point;
    temp *= s / sp;
    temp /= (sp - s - ieps);
    temp *= (M_PI / 2.) / pow(cos(M_PI * abscissas[i] / 2.), 2.); // jacobian

    sum += weights[i] * temp;
  }

  std::complex<double> log_term = - sub_point;
  log_term *= log(low - s - ieps) - log(low);

  return (sum + log_term) / M_PI;
};

// ---------------------------------------------------------------------------
// calculate the dispersion intgral over t
std::complex<double> dispersive_triangle::t_dispersion(int jp, double s)
{
  check_weights();

  std::complex<double> sum = 0.;
  for (int i = 0; i < xN; i++)
  {
    double tp = (t_thresh + EPS) + tan(M_PI * abscissas[i] / 2.);

    std::complex<double> temp;
    temp = lhc_func->disc(tp) * projection(jp, s, tp);
    temp *= (M_PI / 2.);
    temp /= pow(cos(M_PI * abscissas[i] / 2.), 2.); // jacobian

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

  result = sqrt(r_thresh - s - ieps);
  result *= sqrt(p_thresh - s + ieps);
  result *= sqrt(Kallen(s, m1sq, m2sq));
  result /= s;

  return result;
};

std::complex<double> dispersive_triangle::phase_space(double s)
{
  return sqrt(Kallen(s, m1sq, m2sq)) / (2. * s);
};

// ---------------------------------------------------------------------------
// complex Bounds of integration
std::complex<double> dispersive_triangle::t_minus(double s)
{
  return p2*p2 + m1*m1 - (s + p2*p2 - p1*p1) * (s + m1*m1 - m2*m2) / (2. * s) - Kacser(s) / 2.;
};

std::complex<double> dispersive_triangle::t_plus(double s)
{
  return p2*p2 + m1*m1 - (s + p2*p2 - p1*p1) * (s + m1*m1 - m2*m2) / (2. * s) + Kacser(s) / 2.;
};

// ---------------------------------------------------------------------------
// UTILITY FUNCTIONS

// -----------------------------------------------------------------------------
// Check whether or not the integration weights are already saved.
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
