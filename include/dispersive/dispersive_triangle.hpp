// Class to output the evaluation of the triangle diagram using
// dispersive / spectral function representation.
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _DISP_TRI_
#define _DISP_TRI_

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "constants.hpp"
#include "quantum_numbers.hpp"
#include "projection_function.hpp"

namespace jpacTriangle
{
  class dispersive_triangle
  {
  public:
    dispersive_triangle(quantum_numbers * xqn)
    : qns(xqn), projector(qns), mDec2(xqn->mDec*xqn->mDec)
    {};

    // Evalate the diagram at fixed CoM energy^2, s, and exchange mass^2, t
    std::complex<double> eval(double s, double t);

  // ---------------------------------------------------------------------------
  private:
    // All the associated quantum numbers and parameters for the amplitude
    quantum_numbers * qns;
    double mDec2;

    double s, t;
    void fix_energies(double xs, double xt)
    {
      s = xs; t = xt;
    };

    // Two-body phase space
    std::complex<double> rho(double s);

    // Q-function
    projection_function projector;

    // Calculation of dispersion integrals
    double exc = 0.; // small interval around pseudo-threshold to exclude
    std::complex<double> s_dispersion(double low, double high);
    
    std::complex<double> sum_rule();
  };
};

#endif
