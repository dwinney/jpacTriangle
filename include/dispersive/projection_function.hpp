// Spin projection functions, Q_j(s,t)
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PROJECTORS_
#define _PROJECTORS_

#include "constants.hpp"
#include "quantum_numbers.hpp"

std::complex<double> Kallen(std::complex<double> x, std::complex<double> y, std::complex<double> z);

class projection_function
{
public:
  projection_function(quantum_numbers * xqn)
  : qns(xqn), mDec2(xqn->mDec*xqn->mDec)
  {};

  // Evalate the diagram at fixed CoM energy^2, s, and exchange mass^2, t
  std::complex<double> eval(double s, double t);

private:
  quantum_numbers * qns;
  double mDec2;

  // Save s and t so dont need to pass them around
  double s, t;
  void set_energies(double xs, double xt)
  {
      s = xs; t = xt;
  };

  // Kacser function analytically continues momenta between s and t channels
  std::complex<double> Kacser();

  // Complex bounds of integtion
  std::complex<double> t_minus();
  std::complex<double> t_plus();

  // Ratio of momenta q(s) / p(s)
  std::complex<double> barrier_ratio(int ell);

  // Angular kernel functions
  std::complex<double> Q_0();
  std::complex<double> Q(int n);

};

#endif
