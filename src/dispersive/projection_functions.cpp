// Spin projection functions
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dispersive_triangle.hpp"

// ---------------------------------------------------------------------------
// Angular projection Q kernel functions
// These are of the form
//
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

// ---------------------------------------------------------------------------
// Angular projection P_{nj} kernel functions
// combinations of the above with n subtractions
//  and spin j projection onto s channel

std::complex<double> dispersive_triangle::P_0(int n, double s, double tp)
{
  std::complex<double> result;
  result = Q(n, s, tp);
  result /= pow(tp, n);

  return result;
};

// P-wave projector with n subtractions
std::complex<double> dispersive_triangle::P_1(int n, double s, double tp)
{
  std::complex<double> result;
  result = 2. * Q(n+1, s, tp);
  result += (s - mDec2 - 3.*mPi2) * Q(n, s,tp);

  result /= Kacser(s);

  result /= pow(tp, n);
  return -result;
};
