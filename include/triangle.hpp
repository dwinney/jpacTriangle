//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TRIANGLE_
#define _TRIANGLE_

#include <complex>
#include <vector>

#include "constants.hpp"
#include "utilities.hpp"

class triangle
{
public:
  // Constructor
  triangle(double mass)
  : decM(mass)
  {};

  // Evaluate the triangle amplitude
  std::complex<double> eval(double s, double t);

private:
  double decM = 0.; // Decay mass squared
  double EPS = 1.e-6; // small epsilon for i epsilon :)

  // Integration quantities
  int xN = 60;
  bool WG_GENERATED = false;
  std::vector<double> weights, abscissas;
  void check_weights();

  // Triangle kernel, function of s and t as well as two Feynman parameters
  std::complex<double> integrand(double s, double t, double x, double y);

};

#endif
