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

#include "triangle_function.hpp"
#include "constants.hpp"
#include "lefthand_cut.hpp"
#include "integration.hpp"

class feynman_triangle
{
public:
  feynman_triangle(lefthand_cut * a_x, double m)
  : lhc_func(a_x), integ(), kernel(m),
    mDec(m), mDec2(m*m), r_thresh((mDec + mPi) * (mDec + mPi))
  {};

  feynman_triangle(lefthand_cut * a_x, double m, int n)
  : lhc_func(a_x), integ(n), kernel(m),
    mDec(m), mDec2(m*m), r_thresh((mDec + mPi) * (mDec + mPi))
  {};

  // Evalate the diagram
  std::complex<double> eval(int j, int jp, double s);

// ---------------------------------------------------------------------------
private:
  // Isobar lineshape and triangle ampltiude to convolute
  lefthand_cut * lhc_func;
  triangle_function kernel;

  // Mass of the decaying particle and the threshold X Pi threshold
  double mDec, mDec2;
  double r_thresh;

  // Integration stuff
  gauleg integ;

  // Two-point functions, from polynomial contribution
  std::complex<double> mP1(double s);
};


#endif
