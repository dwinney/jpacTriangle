// Class to define a simple breit-wigner propgator with fixed masses
// taken from my kt_3pi repository for KT
//
// Methods allow the evaluation using the Feynman triangle representation
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "breit_wigner.hpp"

// ---------------------------------------------------------------------------
// Evaluate the Breit-Wigner
std::complex<double> breit_wigner::eval(double s)
{
  double real_part, imag_part;
  real_part = s - res_mass * res_mass;
  imag_part = res_width * res_mass;

  return 1. / (real_part * xr - imag_part * xi);
};

std::complex<double> breit_wigner::disc(double s)
{
  double im = std::imag(eval(s));
  return im * xr;
};
