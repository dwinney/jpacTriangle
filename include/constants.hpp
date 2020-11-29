// Header file with global phyiscal constants.
// Everything is in GeV unless explicitly stated otherwise.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------


#ifndef _DEBUG_
#define _DEBUG_

#include <iostream>
#include <iomanip>

namespace jpacTriangle
{
  // ---------------------------------------------------------------------------
  // Generic functions
  template<typename T>
  void debug(T x)
  {
    std::cout << x << std::endl;
  };

  template<typename T, typename F>
  void debug(T x, F y)
  {
    std::cout << std::left;
    std::cout << std::setw(30) << x;
    std::cout << std::setw(30) << y << std::endl;
  };

  template<typename T, typename F, typename G>
  void debug(T x, F y, G z)
  {
    std::cout << std::left;
    std::cout << std::setw(30) << x;
    std::cout << std::setw(30) << y;
    std::cout << std::setw(30) << z << std::endl;
  };
}

#endif

#ifndef _CONSTANT_
#define _CONSTANT_

#include <cmath>
#include <complex>
#include <iostream>

namespace jpacTriangle
{
  //-----------------------------------------------------------------------------
  const double conv = (M_PI / 180.);
  const double EPS = 1.e-6;
  const double M_EULER = 0.577215;

  //Masses
  const double mPi = 0.13957061;
  const double mPi2 = mPi*mPi;
  const double mK = 0.496;
  const double mEta = 0.54753;

  const double mRho = .77545;
  const double mRho2 = mRho * mRho;
  const double mF2 = 1.2754;

  //Thresholds for pi, eta, and K
  const double sthPi = 4.*mPi2;
  const double sthK = 4.*mK*mK;
  const double sthEta = 4.*mEta*mEta;

  //Unit imaginary and real
  const std::complex<double> xr(1., 0.);
  const std::complex<double> xi(0., 1.);
  const std::complex<double> ieps(0., EPS);
};


#endif
