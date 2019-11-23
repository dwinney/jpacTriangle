// Spin projection functions
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "dispersive_triangle.hpp"

std::complex<double> dispersive_triangle::projection(int jp, double s, double tp)
{
  switch(jp)
  {
    case 0: return Q_0(s,tp);
    case 1:
    {
      std::complex<double> result = (2.*s - p1sq - p2sq - m1sq - m2sq) * Q_0(s,tp);
      result += Q_1(s,tp);
      return result;
    }
    default:  std::cout << "invalid jp in projection! Quitting... \n";
              exit(0);
  }

  return 0.;
};

// ---------------------------------------------------------------------------
// Angular projection Q kernel functions
std::complex<double> dispersive_triangle::Q_0(double s, double tp)
{
  std::complex<double> result;
  result = log((tp - ieps) - t_minus(s));
  result -= log((tp - ieps) - t_plus(s));

  result /= Kacser(s);

  return result;
};

std::complex<double> dispersive_triangle::Q_1(double s, double tp)
{
  std::complex<double> result;
  result = tp * Q_0(s, tp);
  result -= (t_plus(s) - t_minus(s)) / Kacser(s);

  return result;
};
