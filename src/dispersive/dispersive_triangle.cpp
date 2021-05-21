// Class to output the evaluation of the triangle diagram using
// dispersive / spectral function representation.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dispersive/dispersive_triangle.hpp"

std::complex<double> jpacTriangle::dispersive_triangle::eval(double s, double t)
{
  // Store s and t so i dont have to keep passing them around
  fix_energies(s, t);

  // Pseudo threshold
  double p_thresh = (qns->mDec - mPi) * (qns->mDec - mPi);

  std::complex<double> result;
  result  = s_dispersion(4.*mPi2, p_thresh - exc);
  result += s_dispersion(p_thresh + exc, std::numeric_limits<double>::infinity());

  return result + (sum_rule() + qns->extra_subtract) * s;
};

// ---------------------------------------------------------------------------
// Two particle phase-space function
std::complex<double> jpacTriangle::dispersive_triangle::rho(double s)
{
  return sqrt(Kallen(s, mPi2, mPi2)) / s;
};

// ---------------------------------------------------------------------------
// calculate the dispersion integral over s with finite bounds of integration
std::complex<double> jpacTriangle::dispersive_triangle::s_dispersion(double low, double high)
{
  auto dsprime = [&](double sp)
  {
    std::complex<double> temp;
    temp = rho(sp) * projector.eval(sp, t) * (s / sp);
    temp -= rho(s) * projector.eval(s, t);
    temp *= s / sp;
    temp /= (sp - s - ieps);
    return temp;
  };

  std::complex<double> result;
  result = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(dsprime, low, high, 5, 1.E-9, NULL);

  std::complex<double> log_term;
  if (high == std::numeric_limits<double>::infinity())
  {
    // subtracted point
    log_term  = - rho(s) * projector.eval(s, t);
    log_term *= log(low - s * xr) - log(low);
  }
  else
  {
    // Log term from subtracted singularity
    log_term = log(high - s * xr) - log(high);
    log_term -= log(low - s * xr) - log(low);
    log_term *= rho(s) * projector.eval(s, t);
  }

  return (result + log_term) / M_PI;
};

std::complex<double> jpacTriangle::dispersive_triangle::sum_rule()
{
  auto dsprime = [&](double sp)
  {
    std::complex<double> temp;
    temp = rho(sp) * projector.eval(sp, t);
    temp /= sp * sp;
    return temp;
  };

  std::complex<double> result;
  result = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(dsprime, 4.*mPi2, std::numeric_limits<double>::infinity(), 5, 1.E-9, NULL);
  
  return result / M_PI;
};
