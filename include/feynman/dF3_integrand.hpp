// This class formulates the kernels in terms of feynman parameters.
// Subtractions, and spin combinations are applied BEFORE integrating to
// save on integration calls
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _INTEGRAND_
#define _INTEGRAND_

#include "constants.hpp"
#include "quantum_numbers.hpp"

class dF3_integrand
{
public:
  dF3_integrand(quantum_numbers* xqns)
  : qns(xqns), mDec2(xqns->mDec * xqns->mDec)
  {};

  // Evaluate the feynman parameters
  std::complex<double> eval(double x, double y, double z);

  // Fix the energies s and t
  void set_energies(double s, double t);

private:
  // All the associated quantum numbers and parameters for the amplitude
  quantum_numbers* qns;

  double mDec2, s, t, psqr; // center of mass energies, t is the exchange particle mass
};

#endif
