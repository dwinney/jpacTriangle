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
std::complex<double> dF3_integrand::eval(double x, double y, double z)
{
  // Store the feynman parameters so to not have to keep passing them
  update_fparams(x, y, z);

  std::complex<double> result = 0.;
  switch (qns->jp)
  {
    case 0:
    {
      result = sub_mT(qns->n, 0);
      break;
    }
    case 1:
    {
      result  = sub_mT(qns->n, 1);
      result -= (mDec2 + 3. * mPi2) * sub_mT(qns->n, 0);
      result += 2.*s * sub_mT(qns->n-1, 0);
      break;
    }
    default:
    {
      std::cout << "\nError! j' = " << qns->jp << " integrands";
      std::cout << " not yet implimented. Quitting... \n";
      exit(1);
    }
  }

  return result;
};

// ---------------------------------------------------------------------------
// mT with n subtractions applied
std::complex<double> dF3_integrand::sub_mT(int n, int k)
{
  // check if theres sufficiently many subtractions applied
  if (n < 0)
  {
    std::cout << "\nError! Insufficient subtractions!\n";
    std::cout << "j = " << qns->j << ", j' = " << qns->jp;
    std::cout << " integral does not converge with n = " << qns->n << " subtractions.";
    std::cout << " Quitting... \n";
    exit(1);
  }

  switch (n)
  {
    // No subtractions
    case 0:
    {
      return mT(k);
    }
    // One subtraction
    case 1:
    {
      return mT(k) - mT(k, true);
    }
    default:
    {
      std::cout << "\nError! n = " << n << " times subtracted integrands";
      std::cout << " not yet implimented. Quitting... \n";
      exit(1);
    }
  }
};


// ---------------------------------------------------------------------------
// The dimensionally regularized integrals but reparameterized
// in terms of the shifted loop momentum relevant for the triangle
// optional bool if True evaluates mT at s = 0
std::complex<double> dF3_integrand::mT(int k, bool SUB)
{
  // Whether or not to evaluate at s or at s = 0
  double denom, psqr, xs;
  if (SUB == true)
  {
    denom = D0; psqr = P0_sqr; xs = 0.;
  }
  else
  {
    denom = D, psqr = P_sqr; xs =  s;
  }

  // Error message
  auto error = [&](int k)
  {
    std::cout << "\nError!\n";
    std::cout << "Feynman integrand mT of divergence order k =" << k << "\n";
    std::cout << "and s-channel projection j = " << qns->j << "\n";
    std::cout << "not yet implimented. Quitting... \n";
    exit(1);
  };

  switch (qns->j)
  {
    // S-wave projection
    case 0:
    {
      switch (k)
      {
        // triangle transform of 1
        case 0:
        {
          return T(0, denom);
        }
        // Triangle transform of k^2
        case 1:
        {
          return T(1, denom) + (x*(1.-z)*mDec2 + y*(1.-z)*mPi2 - x*y*xs)*T(0, denom);
        }
        default: error(k);
      }
      break;
    }

    // P-wave projection
    case 1:
    {
      switch (k)
      {
        // triangle transform of 1
        case 0:
        {
          return z * T(0, denom);
        }
        // Triangle transform of k^2
        case 1:
        {
          return z * T(1, denom) + z*(x+y)*(x+y) * psqr * T(0, denom);
        }
        default: error(k);
      }
      break;
    }
    default:
    {
      std::cout << "\nError! Feynman integrands for j > 1 not yet implimented. Quitting... \n";
      exit(1);
    }
  }
};

// ---------------------------------------------------------------------------
// Dimensionally regularized integral of divergence order k
// D is the combined denominators of all the propagators
std::complex<double> dF3_integrand::T(int ell, double denom)
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
      exit(1);
    }
  }
  return result / (2. * M_PI);
}

// ---------------------------------------------------------------------------
// Update the stored values of the energies s and t
void dF3_integrand::set_energies(double xs, double xt)
{
  s = xs; t = xt;
  P_sqr = P0_sqr - xs / 4.;
};

void dF3_integrand::update_fparams(double xx, double yy, double zz)
{
  // Save the feynman parameters so we dont have to keep passing them
  x = xx; y = yy; z = zz;

  // Update the denominators
  D0 = z*t + (1.-z)*mPi2 - x*z*mDec2 - y*z* mPi2;
  D  = D0 - x*y* s;
};
