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
    case 1: return Q_1(s,tp);
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
  result = log(tp - t_minus(s));
  result -= log(tp - t_plus(s));

  result /= Kacser(s);

  return result;
};

std::complex<double> dispersive_triangle::Q_1(double s, double tp)
{
  std::complex<double> result;
  result = tp * Q_0(s, tp);
  result -= 1.;

  return result;
};
