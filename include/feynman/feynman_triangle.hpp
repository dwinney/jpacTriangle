// Class to output the evaluation of the feynman triangle diagram using
// integration over feynman parameters.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _FEYN_TRI_
#define _FEYN_TRI_

#include "cubature.h"

#include "constants.hpp"
#include "quantum_numbers.hpp"
#include "feynman/dF3_integrand.hpp"

namespace jpacTriangle
{
  class feynman_triangle
  {
  public:
    feynman_triangle(quantum_numbers * xqn)
    : qns(xqn), integrand(xqn)
    {};

    // Evalate the diagram at fixed CoM energy^2, s, and exchange mass^2, t
    std::complex<double> eval(double s, double t);

  // ---------------------------------------------------------------------------
  private:
    // All the associated quantum numbers and parameters for the amplitude
    quantum_numbers * qns;

    // Feynman parameter integrand
    dF3_integrand integrand;

    // Wrapper for interfacing the integrand with hcubature routine
    static int wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);
  };
}

#endif
