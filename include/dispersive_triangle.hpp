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

class dispersive_triangle
{
public:
  dispersive_triangle(lefthand_cut * b_x)
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
    p1sq = m1 * m1; p2sq = m2*m2;

    update_thresholds();
  };

  void set_internalMasses(double q1, double q2)
  {
    m1 = q1; m2 = q2;
    m1sq = q1 * q1; m2sq = q2*q2;
    update_thresholds();
  };

// ---------------------------------------------------------------------------
private:
  lefthand_cut * lhc_func;

  double p1 = 0., p2 = 0.; // external masses
  double m1 = 0., m2 = 0.; // internal loop mass

  double p1sq, p2sq, m1sq, m2sq; // masses squared

  // Physical thresholds
  double s_thresh, t_thresh, p_thresh, r_thresh;
  void update_thresholds()
  {
    // s & t final-state thresholds
    s_thresh = (m1 + m2) * (m1 + m2);
    t_thresh = (p1 + m1) * (p1 + m1);

    // regular and psueodo threshold
    p_thresh = (p1 - p2) * (p1 - p2);
    r_thresh = (p1 + p2) * (p1 + p2);
  };

  // Integration quantities
  int xN = 600;
  bool WG_GENERATED = false;
  std::vector<double> weights, abscissas;
  void check_weights();

  // Kacser function analytically continues momenta between s and t channels
  std::complex<double> Kacser(double s);
  std::complex<double> Kallen(std::complex<double> x, std::complex<double> y, std::complex<double> z)
  {
    return x * x + y * y + z * z - 2. * (x * z + y * z + x * y);
  };

  // Two-body phase space
  std::complex<double> rho(double s);

  // Complex bounds of integtion
  std::complex<double> t_minus(double s);
  std::complex<double> t_plus(double s);

  // Calculation of dispersion integrals
  double exc = EPS; // small interval around pseudo-threshold to exclude
  std::complex<double> s_dispersion(double s, double low, double high);
  std::complex<double> s_dispersion_inf(double s, double low);

  // Cross-channel projected amplitude, b(s).
  std::complex<double> b(double s);

  // Angular kernel functions
  std::complex<double> projection(int jp, double s, double tp);
  std::complex<double> Q_0(double s, double tp);
  std::complex<double> Q_1(double s, double tp);


};


#endif
