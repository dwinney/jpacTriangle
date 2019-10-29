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

#ifndef _TRIANGLE_
#define _TRIANGLE_

#include <complex>
#include <vector>
#include <cmath>

#include "t_integral.hpp"
#include "constants.hpp"
#include "utilities.hpp"

using namespace std;

class triangle
{
public:
  // Constructor
  // The diagram is entirely determined by the decaying particle mass
  triangle(double mass)
  : mDec(mass), cross_channel(mass)
  {};

  // Evaluate the triangle amplitude
  complex<double> eval_feynman(complex<double> s, complex<double> t);
  complex<double> eval_dispersive(double s, double t);

private:
  double mDec = 0.; // Decay mass squared
  double EPS = 1.e-6; // small epsilon for i epsilon :)
  complex<double> ieps = xi * EPS;

  // Integration quantities
  int xN = 30;
  bool WG_GENERATED = false;
  vector<double> weights, abscissas;
  void check_weights();

  // Feynman triangle kernel
  // function of energies s and t as well as two Feynman parameters
  complex<double> feyn_integrand(complex<double> s, complex<double> t, double x);

  // Dispersive KT integral
  t_integral cross_channel;
};

#endif
