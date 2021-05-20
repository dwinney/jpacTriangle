// This class formulates the kernels in terms of feynman parameters.
// Subtractions, and spin combinations are applied BEFORE integrating to
// save on integration calls
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman/dF3_integrand.hpp"

// ---------------------------------------------------------------------------
// Return the value of the integrand in terms of the feynman parameters
std::complex<double> jpacTriangle::dF3_integrand::eval(double x, double y, double z)
{
    // Store the feynman parameters so to not have to keep passing them
    update_fparams(x, y, z);

    std::complex<double> result = 0.;
    
    // check if theres sufficiently many subtractions applied
    if (qns->n < 0)
    {
      std::cout << "\nError! Insufficient subtractions!\n";
      std::cout << "j = " << qns->j << ", j' = " << qns->jp;
      std::cout << " integral does not converge with n = " << qns->n << " subtractions.";
      std::cout << " Quitting... \n";
      exit(0);
    }

    // apply the necessary subtractions in s
    int id = qns->id();

    switch (qns->n)
    {
      // No subtractions
      case 0:
      {
          return mT(id, s);
      }
      // One subtraction
      case 1:
      {
          return mT(id, s) - mT(id, 0.);
      }
      default:
      {
        std::cout << "\nError! n = " << qns->n << " times subtracted integrands";
        std::cout << " not yet implimented. Quitting... \n";
        exit(0);
      }
    }
};

// ---------------------------------------------------------------------------
// Triangle kernels
// optional bool if True evaluates mT at s = 0
std::complex<double> jpacTriangle::dF3_integrand::mT(int id, double _s)
{
  // Whether or not to evaluate at s or at s = 0
  denom = denom0 - x*y* _s,  delta = delta0 - x*y* _s;

  auto error = [&] ()
  {
      std::cout << "\nError! projection_function:";
      std::cout << " j = " << std::to_string(qns->j);
      std::cout << " and j' = " << std::to_string(qns->jp);
      std::cout << " (code " << std::to_string(qns->id()) << ")";
      std::cout << " combination not available. Quitting... \n";
      exit(0);
  };

  std::complex<double> result = 0.;
  switch (id)
  {
      // S-wave, scalar exchange
      case 0:
      {
          result = T(0);
          break;
      }

      // S-wave, vector exchange
      case 1:
      {
          result  = T(1);
          result += (delta + 2.*_s - mDec2 - 3.*mPi2) * T(0);
          break;
      }

      // P-wave, scalar exchange
      case 10:
      {
          result = z * T(0);
          break;
      };

      // P-wave, vector exchange
      case 11:
      {
          result  = (3.*z - 1.) * T(1) / 2.;
          result += z * (delta + 2.*_s - mDec2 - 3.*mPi2) * T(0);
          break;
      };

      // D-wave, scalar exchange
      case 20: 
      {
          result  = z*z * T(0);
          break;
      };

      case 10000:
      { 
          result  = (_s + mDec2 - mPi2) * (T(1) + delta * T(0));
          result += (_s - mDec2 - mPi2)*(mDec2 - mPi2) * T(0);
          break;
      };

      // Omega case
      case -11111:
      {
          result  = - 2. * T(1);
          break;
      };

      case 11010:
      {
          result  = sqrt(mDec2) * 2. * T(1);
          break;
      };

      default: error();
  };
  return result;
};

// ---------------------------------------------------------------------------
// Dimensionally regularized integral of divergence order k
// D is the combined denominators of all the propagators
std::complex<double> jpacTriangle::dF3_integrand::T(int ell)
{
  std::complex<double> result;
  switch (ell)
  {
    // Convergent
    case 0:
    {
      result = xr / (denom - ieps);
      break;
    }
    // linearly divergent in ell^2
    case 1:
    {
      result = 2. * log(denom - ieps);
      break;
    }
    default:
    {
      std::cout << "\nError! Feynman integrand T of divergence order l =" << ell;
      std::cout << " not yet implimented. Quitting... \n";
      exit(0);
    }
  }
  return result / (2. * M_PI);
}
