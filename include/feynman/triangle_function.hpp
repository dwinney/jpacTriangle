// This object defines a Triangle amplitude, i.e. the rescattering
// diagram associated with an intermediate t-channel exchange.
//
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TRIFUNC_
#define _TRIFUNC_

#include <complex>
#include <vector>
#include <cmath>
#include <iomanip>

#include "constants.hpp"
#include "rational_integrals.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>

class triangle_function
{
public:
  // Constructor
  // All we need is the decaying particle mass (so far)
  triangle_function(double m)
  : mDec(m), mDec2(m*m)
  {};

  // Feynman triangle kernel for n subtractions spin j projection and spin jp exchange
  std::complex<double> eval(int n, int j, int jp, double s, double t);

private:
  double mDec, mDec2;

  // We build the above from a basis of functions which we evaluate with these
  std::complex<double> mT(int n, int k, int j, double s, double t);
  std::complex<double> mT_integrand(int n, int k, int j, double s, double t, double x);

  // Analytically calculated 4-dimensional (dimensionally-regularized) integrals
  // int_mTell = mT_j{(k^2)^ell}(s,t)
  std::complex<double> int_mT0(int j, double s, double t, double x);
  std::complex<double> int_mT1(int j, double s, double t, double x);
};

#endif
