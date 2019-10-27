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
    c = 0.5 *  (mDec * mDec - mPi * mPi);
  };

  complex<double> operator() (double s, double t);
private:
  double mDec;
  double a, b, c;

  //
  int xN = 30;
  double EPS = 1.e-6;
  complex<double> integrand(double t, complex<double> tp);


  // Dispersive functions
  complex<double> Kacser(complex<double> s);
  complex<double> t_minus(double s);
  complex<double> t_plus(double s);

  // complex integration
  complex<double> integ_sthPi_c(double s, double t);
  complex<double> integ_c_a(double s, double t);
  complex<double> integ_a_b(double s, double t);
  complex<double> integ_b(double s, double t);
};

#endif
