// This object defines a Triangle amplitude, i.e. the rescattering
// diagram associated with an intermediate t-channel exchange.
//
// Evaluated specifically with the Feynman method
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

#include "constants.hpp"
#include "utilities.hpp"
#include "lefthand_cut.hpp"

class feynman_triangle
{
public:
  feynman_triangle(lefthand_cut * b_x)
  : lhc_func(b_x)
  {};

  // Evalate the diagram
  std::complex<double> eval(double s);

  // ---------------------------------------------------------------------------
  // Utilities
  void set_Nint(int n)
  {
    xN = n;
    WG_GENERATED = false;
  }

  void set_externalMasses(double m1, double m2)
  {
    p1 = m1; p2 = m2;
    update_thresholds();
  };

  void set_internalMasses(double q1, double q2)
  {
    m1 = q1; m2 = q2;
    update_thresholds();
  };

// ---------------------------------------------------------------------------
private:
  lefthand_cut * lhc_func;

  double p1 = 0., p2 = 0.; // external masses
  double m1 = 0., m2 = 0.; // internal loop mass

  // Physical thresholds
  double s_thresh, t_thresh, p_thresh, r_thresh;
  void update_thresholds()
  {
    // s & t final-state thresholds
    s_thresh =  (m1 + m2) * (m1 + m2);
    t_thresh = (p2 + m1) * (p2 + m1);

    // regular and psueodo threshold
    p_thresh = (p1 - p2) * (p1 - p2);
    r_thresh = (p1 + p2) * (p1 + p2);
  };

  // Integration quantities
  int xN = 800;
  bool WG_GENERATED = false;
  std::vector<double> weights, abscissas;
  void check_weights();

  // Feynman triangle kernel
  // function of energies s and one feynman parameter x
  std::complex<double> kernel_integrand(double s, double t, double x);
  std::complex<double> kernel(double s, double t);
};


#endif
