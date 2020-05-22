// This object defines a Triangle amplitude, i.e. the rescattering
// diagram associated with an intermediate t-channel exchange.
//
// Evaluated specifically with the Feynman method. That is, the perturbative
// triangle convoluted with an isobar for the exchange.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _FEYN_TRI_
#define _FEYN_TRI_

#include <complex>
#include <vector>
#include <cmath>
#include <iomanip>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "cubature.h"

#include "constants.hpp"
#include "quantum_numbers.hpp"
#include "lefthand_cut.hpp"
#include "feynman/dF3_integrand.hpp"

class feynman_triangle
{
public:
  feynman_triangle(quantum_numbers * xqn, lefthand_cut * a_x)
  : lhc_func(a_x), qns(xqn), integrand(qns)
  {};

  // Evalate the diagram
  std::complex<double> eval(double s);

  // Evaluate using fixed-mass exchange with an infinitesimal width
  std::complex<double> kernel(double s, double t);

// ---------------------------------------------------------------------------
private:
  // Isobar lineshape and triangle ampltiude to convolute
  lefthand_cut * lhc_func;

  // All the associated quantum numbers and parameters for the amplitude
  quantum_numbers * qns;

  // Feynman parameter integrand
  dF3_integrand integrand;

  // Wrapper for interfacing the integrand with hcubature routine
  static int wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);
};

#endif
