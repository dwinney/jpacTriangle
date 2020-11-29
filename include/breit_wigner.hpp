// Class to define a simple breit-wigner propgator with fixed masses
// taken from my kt_3pi repository for KT
//
// Methods allow the evaluation using the Feynman triangle representation
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _BW_
#define _BW_

#include "lefthand_cut.hpp"
#include "constants.hpp"

namespace jpacTriangle
{
  class breit_wigner : public lefthand_cut
  {
  public:
    breit_wigner(double mass, double width = 0.001)
    : res_mass(mass), res_width(width)
    {};

    // Evaluate the BW
    std::complex<double> eval(double s);

    // Evaluate the discontinuity which in this case is just the imaginary part
    std::complex<double> disc(double s);

  private:
    double res_mass, res_width;
};
}

#endif
