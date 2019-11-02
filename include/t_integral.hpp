// This object calculates the integration along the complex plane
// required to analytically continue the projection of cross-channel quantities
// to the direct channel
//
// Seperated as its own object for easier debugging
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _COMP_INT_
#define _COMP_INT_

#include "constants.hpp"
#include "utilities.hpp"

#include <complex>

using namespace std;

class t_integral
{
public:
  // Constructor
  t_integral(double mass)
  : mDec(mass)
  {
    a = (mDec - mPi) * (mDec - mPi); // Pseudo-threshold
    b = (mDec + mPi) * (mDec + mPi); // threshold
    c = 0.5 *  (mDec * mDec - mPi * mPi); // turnaround point above cut
  };

  complex<double> operator() (complex<double> s, complex<double> t);
private:
  double mDec;
  double a, b, c;

  //
  int xN = 60;
  double EPS = 1.e-6;
  complex<double> ieps = xi * EPS;
  complex<double> integrand(complex<double> t, complex<double> tp);


  // Dispersive functions
  complex<double> Kacser(complex<double> s);
  complex<double> t_minus(complex<double> s);
  complex<double> t_plus(complex<double> s);

  // complex integration
  complex<double> integ_sthPi_c(complex<double> s, complex<double> t);
  complex<double> integ_c_a(complex<double> s, complex<double> t);
  complex<double> integ_a_b(complex<double> s, complex<double> t);
  complex<double> integ_b(complex<double> s, complex<double> t);
};

#endif
