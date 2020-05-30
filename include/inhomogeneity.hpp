// This object defines a Triangle amplitude, i.e. the rescattering
// diagram associated with an intermediate t-channel exchange given by a
// lefhand_cut object
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _INHOM_
#define _INHOM_

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "cubature.h"

#include "constants.hpp"
#include "quantum_numbers.hpp"
#include "lefthand_cut.hpp"
#include "feynman/dF3_integrand.hpp"

template <class T>
class inhomogeneity
{
public:
  inhomogeneity(quantum_numbers * xqn, lefthand_cut * a_x)
  : lhc_func(a_x), qns(xqn), triangle_kernel(xqn)
  {};

  // Evalate the diagram
  std::complex<double> eval(double s)
  {
    double r_thresh = (qns->mDec + mPi) * (qns->mDec + mPi);

    auto dTp = [&](double tp)
    {
      return lhc_func->disc(tp) * triangle_kernel.eval(s, tp);
    };

    std::complex<double> result;
    result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(dTp, r_thresh, std::numeric_limits<double>::infinity(), 0, 1.E-4, NULL);

    result /= M_PI;

    return 1. + result;
  };

// ---------------------------------------------------------------------------
private:
  // Class which defines the triangle function to convolute.
  // Can be either <dispersive_triangle> or <feynman_triangle>
  T triangle_kernel;

  // Isobar lineshape and triangle ampltiude to convolute
  lefthand_cut * lhc_func;

  // All the associated quantum numbers and parameters for the amplitude
  quantum_numbers * qns;
};

#endif
