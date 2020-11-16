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
  inline void set_energies(double xs, double xt)
  {
    s = xs; t = xt;
  };

private:
  // All the associated quantum numbers and parameters for the amplitude
  quantum_numbers* qns;

  double denom, delta;
  double denom0, delta0;
  double mDec2;
  double s, t; // center of mass energies, t is the exchange particle mass

  // Currently stored feynman parameters
  double x, y ,z;
  inline void update_fparams(double xx, double yy, double zz)
  {
    // Save the feynman parameters so we dont have to keep passing them
    x = xx; y = yy; z = zz;

    // Update the denominators
    denom0 = z*t + (1.-z)*mPi2 - x*z*mDec2 - y*z* mPi2;

    // update P's
    delta0 = x*(1.-z)*mDec2 + y*(1.-z)*mPi2;
  };

  // Dimensionally regularized integrals
  // denom is the combined denominators of all the propagators
  // ell is the degree of divergence
  std::complex<double> T(int ell);

  // The dimensionally regularized integrals but reparameterized
  // in terms of the shifted loop momentum relevant for the triangle
  // optional bool if True evaluates mT at s = 0
  std::complex<double> mT(int id, bool SUB = false);

};

#endif
