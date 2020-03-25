// This object defines a Triangle amplitude, i.e. the rescattering
// diagram associated with an intermediate t-channel exchange.
//
// Evaluated specifically with the dispersion relation method
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DISP_TRI_
#define _DISP_TRI_

#include <complex>
#include <vector>
#include <cmath>
#include <iomanip>

#include "constants.hpp"
#include "utilities.hpp"
#include "lefthand_cut.hpp"
#include "integration.hpp"

class dispersive_triangle
{
public:
  dispersive_triangle(lefthand_cut * a_x)
  : lhc_func(a_x), integ()
  {};

  dispersive_triangle(lefthand_cut * a_x, int n)
  : lhc_func(a_x), integ(n)
  {};

  // Evalate the diagram
  std::complex<double> eval(int j, int jp, double s);

  // Setting utility
  void set_decayMass(double m1)
  {
    mDec = m1;
    mDec2 = m1 * m1;
    update_thresholds();
  };

// ---------------------------------------------------------------------------
private:
  lefthand_cut * lhc_func;

  double mDec; // external decay mass
  double mDec2, mPi2 = mPi * mPi; // masses squared

  // Physical thresholds
  double s_thresh, p_thresh, r_thresh;
  void update_thresholds()
  {
    // s final-state thresholds
    s_thresh = 4. * mPi2;
    // regular and psueodo threshold
    p_thresh = (mDec - mPi) * (mDec - mPi);
    r_thresh = (mDec + mPi) * (mDec + mPi);
  };

  // Integration stuff
  gauleg integ;

  // ---------------------------------------------------------------------------
  // Kacser function analytically continues momenta between s and t channels
  std::complex<double> Kacser(double s);
  std::complex<double> Kallen(std::complex<double> x, std::complex<double> y, std::complex<double> z)
  {
    return x * x + y * y + z * z - 2. * (x * z + y * z + x * y);
  };

  // Ratio of barrier factors that need to be taken out
  std::complex<double> barrier_ratio(int l, double s);

  // Two-body phase space
  std::complex<double> rho(double s);

  // Complex bounds of integtion
  std::complex<double> t_minus(double s);
  std::complex<double> t_plus(double s);

  // Calculation of dispersion integrals
  double exc = EPS; // small interval around pseudo-threshold to exclude
  std::complex<double> s_dispersion(int j, int jp, double s, double low, double high);
  std::complex<double> s_dispersion_inf(int j, int jp, double s, double low);

  // Cross-channel projected amplitude, b(s).
  std::complex<double> b(int j, int jp, double s);

  // Angular kernel functions
  std::complex<double> Q_0(double s, double tp);
  std::complex<double> Q(int n, double s, double tp);

  std::complex<double> projector(int n, int j, int jp, double s, double tp);

};


#endif
