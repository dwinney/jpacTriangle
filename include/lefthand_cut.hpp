// Abstract class to define a dispersive model for the left-hand cut
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _LHC_
#define _LHC_

#include <complex>

namespace jpacTriangle
{
  class lefthand_cut
  {
  public:
    // Empty constructor
    lefthand_cut(){};

    // ---------------------------------------------------------------------------
    // THESE FUNCTIONS MUST BE OVERWRITTEN IN ANY SPECIFIC IMPLEMENTATION

    // Method to evalute the whole model for the lefthand cut
    virtual std::complex<double> eval(double s) = 0;

    // Method to evaluate only the discontinuity
    virtual std::complex<double> disc(double s) = 0;
  };
};

#endif
