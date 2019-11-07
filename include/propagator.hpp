// Abstract class to define a dispersive popogator.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PROPAGATOR_
#define _PROPAGATOR_

#include <complex>


class propagator
{
public:
  // Empty constructor
  propagator(){};

  // ---------------------------------------------------------------------------
  // THESE FUNCTIONS MUST BE OVERWRITTEN IN ANY SPECIFIC IMPLEMENTATION

  // Method to evalute the whole propagator
  virtual std::complex<double> eval(double s) = 0;

  // Method to evaluate only the discontinuity
  virtual std::complex<double> disc(double s) = 0;
};

#endif
