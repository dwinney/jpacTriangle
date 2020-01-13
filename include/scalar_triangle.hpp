// This object defines a Triangle amplitude, i.e. the rescattering
// diagram associated with an intermediate t-channel exchange.
//
// Methods allow the evaluation using the Feynman triangle representation
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TRIANGLE_
#define _TRIANGLE_

#include <complex>
#include <vector>
#include <cmath>
#include <iomanip>

#include "constants.hpp"
#include "utilities.hpp"
#include "lefthand_cut.hpp"

class scalar_triangle
{
public:
  // Empty Constructor
  scalar_triangle(lefthand_cut * b_x)
  : lhc_func(b_x)
  {};

  // Parameterized Constructor
  scalar_triangle(lefthand_cut * b_x, double x1, double x2, double q1, double q2)
  : lhc_func(b_x), p1(x1), p2(x2), m1(q1), m2(q2)
  {
    update_thresholds();
  }

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
  // Evaluate the triangle amplitude
  std::complex<double> eval_feynman(double s);
  std::complex<double> eval_dispersive(double s);

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
  int xN = 200;
  bool WG_GENERATED = false;
  std::vector<double> weights, abscissas;
  void check_weights();

  // ---------------------------------------------------------------------------
  // Feynman FUNCTIONS

  // Feynman triangle kernel
  // function of energies s and one feynman parameter x
  std::complex<double> kernel_integrand(double s, double t, double x);
  std::complex<double> triangle_kernel(double s, double t);

  // ---------------------------------------------------------------------------
  // Dispersive FUNCTIONS

  // Kacser function analytically continues momenta between s and t channels
  std::complex<double> Kacser(double s);
  std::complex<double> Kallen(std::complex<double> x, std::complex<double> y, std::complex<double> z)
  {
    return x * x + y * y + z * z - 2. * (x * z + y * z + x * y);
  };

  // Complex bounds of integtion
  std::complex<double> t_minus(double s);
  std::complex<double> t_plus(double s);

  // Calculation of dispersion integrals
  double exc = EPS; // small interval around pseudo-threshold to exclude
  std::complex<double> s_dispersion(double s, double low, double high);
  std::complex<double> s_dispersion_inf(double s, double low);


  // Dispersion kernel functions
  std::complex<double> projection(double s, double tp);
  std::complex<double> t_dispersion(double s);
};

#endif
