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

// Evaluate the convolution of the LHC function with triangle function
std::complex<double> triangle::eval_feynman(double s)
{
    check_weights();

    std::complex<double> sum = 0.;
    for (int i = 0; i < xN; i++)
    {
      double tp = (t_thresh + EPS) + tan(M_PI * abscissas[i] / 2.);

      std::complex<double> temp;
      temp = lhc_func->disc(tp) * triangle_kernel(s, tp) / tp;
      temp *= (M_PI / 2.);
      temp /= pow(cos(M_PI * abscissas[i] / 2.), 2.); // jacobian

      sum += weights[i] * temp;
    }

    sum /= M_PI;

    return sum;
};

// Triangle function from the perturbation theory result
std::complex<double> triangle::triangle_kernel(double s, double t)
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
std::complex<double> triangle::kernel_integrand(double s, double t, double x)
{
  std::complex<double> a, b, c, d;
  
  a = p2 * p2;
  b =  (t - ieps) + x * s - m1*m1 - (x - 1.) * p2*p2 - x * p1*p1 ;
  c = (1. - x) * (t - ieps) + x * m1*m1 + x*(x-1.)* p1*p1;
  d = b * b - 4. * a * c; // discriminant

  // Roots of the polynomial
  std::complex<double> y_plus = (-b + sqrt(xr * d + ieps)) / (2. * a);
  std::complex<double> y_minus = (-b - sqrt(xr * d + ieps)) / (2. * a);

  std::complex<double> result;
  result = log(xr - x + y_minus) - log(xr - x + y_plus);
  result -= log(y_minus) - log(y_plus);
  result /= sqrt(xr * d + ieps);

  return result;
};

// ---------------------------------------------------------------------------
// EVALUTE THE TRIANGLE THE KT WAY

std::complex<double> triangle::eval_dispersive(double s)
{
  // if pseudo threshold is in the bounds of integration, exclude a small interval around it
  // this is alwys the case for decays, but to keep the code general I consider the other case
  if (p_thresh < s_thresh)
  {
    return s_dispersion_inf(s, s_thresh);
  }
  else
  {
    std::complex<double> temp;
    temp = s_dispersion(s, s_thresh + EPS, p_thresh - exc);
    temp += s_dispersion_inf(s, p_thresh + exc);

    return temp;
  }
};

// ---------------------------------------------------------------------------
// calculate the dispersion integral over s with finite bounds of integration
std::complex<double> triangle::s_dispersion(double s, double low, double high)
{
  double w[xN + 1 ], x[xN + 1];
  gauleg(low, high, x, w, xN);

  // Subtract off the pole at s = sp
  std::complex<double> sub_point = t_dispersion(s);

  std::complex<double> sum = 0.;
  for (int i = 1; i <= xN; i++)
  {
    std::complex<double> temp;
    temp = t_dispersion(x[i]) - sub_point;
    temp /= x[i] * (x[i] - s - ieps);

    sum += w[i] * temp;
  }

  std::complex<double> log_term;
  log_term = log(high - s - ieps) - log(high);
  log_term -= log(low - s - ieps) - log(low);
  log_term *= sub_point / s;

  return (sum + log_term) / M_PI;
};

// calculate the dispersion integral over s up to infinity
std::complex<double> triangle::s_dispersion_inf(double s, double low)
{
  check_weights();

  std::complex<double> sub_point = t_dispersion(s);

  std::complex<double> sum = 0.;
  for (int i = 0; i < xN; i++)
  {
    double sp = low + tan(M_PI * abscissas[i] / 2.);

    std::complex<double> temp;
    temp = t_dispersion(sp) - sub_point;
    temp /= sp * (sp - s - ieps);
    temp *= (M_PI / 2.) / pow(cos(M_PI * abscissas[i] / 2.), 2.); // jacobian

    sum += weights[i] * temp;
  }

  std::complex<double> log_term = - sub_point / s;
  log_term *= log(low - s - ieps) - log(low);

  return (sum + log_term) / M_PI;
};

// ---------------------------------------------------------------------------
// calculate the dispersion intgral over t
std::complex<double> triangle::t_dispersion(double s)
{
  check_weights();

  std::complex<double> sum = 0.;
  for (int i = 0; i < xN; i++)
  {
    double tp = (t_thresh + EPS) + tan(M_PI * abscissas[i] / 2.);

    std::complex<double> temp;
    temp = lhc_func->disc(tp) * projection(s, tp);
    temp *= (M_PI / 2.);
    temp /= pow(cos(M_PI * abscissas[i] / 2.), 2.); // jacobian

    sum += weights[i] * temp;
  }

  sum /= M_PI;

  return sqrt(Kallen(s, m1*m1, m2*m2)) * sum;
};

// ---------------------------------------------------------------------------
// Kacser function which includes the correct analytic structure of
// product of breakup momenta, p(s) * q(s)
std::complex<double> triangle::Kacser(double s)
{
  std::complex<double> result;

  result = sqrt(r_thresh - s - ieps);
  result *= sqrt(p_thresh - s + ieps);
  result *= sqrt(Kallen(s, m1*m1, m2*m2));
  result /= s;

  return result;
};

// ---------------------------------------------------------------------------
// complex Bounds of integration
std::complex<double> triangle::t_minus(double s)
{
  return p2*p2 + m1*m1 - (s + p2*p2 - p1*p1) * (s + m1*m1 - m2*m2) / (2. * s) - Kacser(s) / 2.;
};

std::complex<double> triangle::t_plus(double s)
{
  return p2*p2 + m1*m1 - (s + p2*p2 - p1*p1) * (s + m1*m1 - m2*m2) / (2. * s) + Kacser(s) / 2.;
};

// ---------------------------------------------------------------------------
// Angular projection kernel
std::complex<double> triangle::projection(double s, double tp)
{
  std::complex<double> result;
  result = log((tp - ieps) - t_minus(s));
  result -= log((tp - ieps) - t_plus(s));
  result /= Kacser(s);

  return result;
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
