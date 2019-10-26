#include "triangle.hpp"
#include <iostream>

int main(int argc, char** argv)
{
  if (argc != 3)
  {
    std::cout << "Usage: \t";
    std::cout << "./plot_triangle <decay_mass> <exchange_mass>" << std::endl;
    exit(0);
  }
  double dec = std::stod(argv[1]);
  double ex = std::stod(argv[2]);

  // Plotting bounds
  double low = 0;
  double high = 1.;

  if (ex * ex + 3. * mPi * mPi > dec * dec)
  {
    std::cout << "ERROR! Input exchange is not kinematically allowed. Quitting... \n";
    exit(0);
  }

  triangle tri(dec * dec);

  std::vector<double> s;
  vector< std::complex<double> > fx;

  for (int i = 0; i < 160; i++)
  {
    double si = low + double(i) * (high - low) / 160.;
    complex<double> fx_i = tri.eval(si, ex * ex);

    s.push_back(si);
    fx.push_back(fx_i);
  }

  quick_plot(s, fx, "triangle");

  return 1.;
};
