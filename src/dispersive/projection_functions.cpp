// Spin projection functions
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dispersive_triangle.hpp"

std::complex<double> dispersive_triangle::projector(int n, int j, int jp, double s, double tp)
{
  std::complex<double> result = 0.;

  switch (j)
  {
    // s-wave projection
    case 0:
    {
      switch (jp)
      {
        // s wave, scalar exchange
        case 0: result = Q(n, s, tp); break;

        // s wave, vector exchange
        case 1: //
        {
          result = Q(n+1, s, tp);
          result += (2.*s - mDec2 - 3.*mPi2) * Q(0, s, tp);
          break;
        }
        default: std::cout << "j and j' combination not available. Quitting... \n"; exit(0);
      }
      break;
    }

    // p-wave projection
    case 1:
    {
      switch (jp)
      {
        // p - wave, scalar exchange
        case 0:
        {
          result = 2. * Q(n+1, s, tp);
          result += (s - mDec2 - 3.*mPi2) * Q(n, s, tp);
          break;
        }

        // p - wave, vector exchange
        case 1:
        {
          result = 2. * Q(n+2, s, tp);
          result += (5.*s + - 3. * mDec2 - 9. * mPi2) * Q(n+1, s, tp);
          result += (2.*s*s - 3.*mDec2*s - 9.*mPi2*s + mDec2*mDec2 + 6.*mDec2*mPi2 + 9.*mPi2*mPi2) * Q(n,s,tp);
          break;
        }
        default: std::cout << "j and j' combination not available. Quitting... \n"; exit(0);
      }
      result /= Kacser(s);
      break;
    }

    // d-wave projection
    case 2:
    {
      if (jp == 0)
      {
        result = 12.* Q(n+2, s, tp);
        result += (12. * s - 12. * mDec2 - 36. * mPi2) * Q(n+1, s, tp);
        result += (3.*s*s - 6.*mDec2*s + 3.*mDec2*mDec2 - 18.*mPi2*s + 18.*mPi2*mDec2 + 27.*mPi2*mPi2 - 1.) * Q(n,s,tp);
        result /= 2.;

        result /= pow(Kacser(s), 2.);
      }
      break;
    }

    default: std::cout << "j and j' combination not available. Quitting... \n"; exit(0);
  }
  result *= barrier_ratio(j, s);
  result /= pow(tp, double(n));
  return result;
};

// ---------------------------------------------------------------------------
// Angular projection Q kernel functions
// These are of the form:
// 1/Kacser(s) * \int_{t_minus}^{t_plus} x^n / (tp - tp - ieps)
std::complex<double> dispersive_triangle::Q_0(double s, double tp)
{
  std::complex<double> result;
  result = log(tp - ieps - t_minus(s));
  result -= log(tp - ieps - t_plus(s));

  result /= Kacser(s);
  return result;
};

std::complex<double> dispersive_triangle::Q(int n, double s, double tp)
{
  switch (n)
  {
    case 0: return Q_0(s,tp);
    case 1: return tp * Q_0(s,tp) - 1.;
    case 2: return tp * tp * Q_0(s,tp) - tp - 0.5* (pow(t_plus(s), 2.) - pow(t_minus(s), 2.));
    default: std::cout << "Not enough Q's!!! \n"; exit(0);
  }
};
