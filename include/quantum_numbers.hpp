// Small struct to carry along all the parameters the amplitude depends on
// In a KT context this may be replaced with a decay_kinematic object.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _QNS_
#define _QNS_

struct quantum_numbers
{
  int n = 1; // number of Subtractions
  int j = 0, lam = 0; // spin and helicity projection in the s-channel
  int jp = 0, lamp = 0; // spin and helicity projections in the t-channel

  double mDec = 0.;
};

#endif
