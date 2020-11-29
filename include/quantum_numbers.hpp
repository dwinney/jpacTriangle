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
  int  n = 1; // number of Subtractions in s integral
  int  l = 0; // number of subtractions in t integral
  int  J = 0; // spin of the decaying particle
  int  P = 1; // Parity of decaying particle
  int  j = 0,  lam = 0;   // spin and helicity projection in the s-channel
  int jp = 0, lamp = 0; // spin and helicity projections in the t-channel

  double mDec = 0.;

  inline int id()
  {
    int id = 10000 * J + 1000 * lam + 100 * lamp + 10 * j + jp;
    return id * P;
  };

  inline void set_id(int x)
  {
    (x > 0) ? (P = 1) : (P = -1);
    x *= P;
    J    = x % 10000; x -= 10000 * J;
    lam  = x % 1000;  x -= 1000  * lam;
    lamp = x % 100;   x -= 100   * lamp;
    j    = x % 10;    x -= 10    * j;
    jp   = x;
  };
};

#endif
