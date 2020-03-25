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

#include "rational_integrals.hpp"
#include "constants.hpp"
#include "utilities.hpp"
#include "lefthand_cut.hpp"
#include "integration.hpp"

class feynman_triangle
{
public:
  feynman_triangle(lefthand_cut * a_x)
  : lhc_func(a_x), integ()
  {};

  feynman_triangle(lefthand_cut * a_x, int n)
  : lhc_func(a_x), integ(n)
  {};


  // Evalate the diagram
  std::complex<double> eval(int j, int jp, double s);

  // ---------------------------------------------------------------------------
  // Utilities
  void set_decayMass(double m)
  {
    mDec = m; mDec2 = m * m;
    update_thresholds();
  };

// ---------------------------------------------------------------------------
private:
  lefthand_cut * lhc_func;

  // Masses
  double mDec, mDec2;

  // Physical thresholds
  double s_thresh, p_thresh, r_thresh;
  void update_thresholds()
  {
    // s final-state thresholds
    s_thresh = sthPi;

    // regular and psueodo threshold
    p_thresh = (mDec - mPi) * (mDec - mPi);
    r_thresh = (mDec + mPi) * (mDec + mPi);
  };

  // Integration stuff
  gauleg integ;

  // Two-point functions, from polynomial contribution
  std::complex<double> mP1(double s);

  // Feynman triangle kernel for n subtractions spin j projection and spin jp exchange
  std::complex<double> triangle_kernel(int n, int j, int jp, double s, double t);

  // We build the above from a basis of functions which we evaluate with these
  std::complex<double> mT(int n, int k, int z, double s, double t);
  std::complex<double> mT_integrand(int k, int z, double s, double t, double x);
  std::complex<double> int_mT00(double s, double t, double x);
  std::complex<double> int_mT10(double s, double t, double x);
  std::complex<double> int_mT20(double s, double t, double x);
  std::complex<double> int_mT01(double s, double t, double x);

};


#endif
