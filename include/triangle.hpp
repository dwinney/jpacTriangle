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

#include "constants.hpp"
#include "utilities.hpp"

class triangle
{
public:
  // Empty Constructor
  triangle()
  {};

  // Parameterized Constructor
  // triangle(double m1, double m1, double q1, double q2)

  // ---------------------------------------------------------------------------
  // Utilities
  void set_Nint(int n)
  {
    xN = n;
  }

  void set_externalMasses(double m1, double m2)
  {
    p1 = m1; p2 = m2;

  };

  void set_internalMass(double q1, double q2)
  {
    m1 = q1; m2 = q2;
  };

  void set_exchangeMass(double mEx, double gamEx = 0.0001)
  {
    t = mEx * mEx - xi * gamEx;
  };

  // ---------------------------------------------------------------------------
  // Evaluate the triangle amplitude
  std::complex<double> eval_feynman(double s);
  std::complex<double> eval_dispersive(double s);

private:
  double p1 = 0., p2 = 0.; // external masses
  double m1 = 0., m2 = 0.; // internal loop mass

  // t = complex invariant mass of exchange particle
  std::complex<double> t = 0.;

  // Integration quantities
  int xN = 800;
  bool WG_GENERATED = false;
  std::vector<double> weights, abscissas;
  void check_weights();

  // ---------------------------------------------------------------------------
  // Feynman FUNCTIONS

  // Feynman triangle kernel
  // function of energies s and one feynman parameter x
  std::complex<double> feyn_integrand(double s, double x);

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
  double exc = 0.005; // small interval around pseudo-threshold to exclude
  std::complex<double> s_dispersion(double s, double low, double high);
  std::complex<double> s_dispersion_inf(double s, double low);


  // Dispersion kernel functions
  std::complex<double> projection(double s, double tp);
  std::complex<double> propagator(double tp);
  std::complex<double> t_dispersion(double s);
};

#endif
