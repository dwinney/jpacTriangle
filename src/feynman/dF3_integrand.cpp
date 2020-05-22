// This class formulates the kernels in terms of feynman parameters.
// Subtractions, and spin combinations are applied BEFORE integrating to
// save on integration calls
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman/dF3_integrand.hpp"

std::complex<double> dF3_integrand::eval(double x, double y, double z)
{
  // Combined denominator for the triangle
 std::complex<double> D, result;

 D = z*t + (1.-z)*mPi2 - x*z* (mDec2 - ieps) - y*z* mPi2 - x*y* s;

 result = xr / (D - ieps);

 return result / (2. * M_PI);
};

void dF3_integrand::set_energies(double xs, double xt)
{
  s = xs; t = xt;
  psqr = (2. * (mDec2 + mPi2) - s) / 4.;
};
