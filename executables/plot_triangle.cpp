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
  double low = sthPi + 1.e-6;
  double high = 1.;

  if (ex * ex + 3. * mPi * mPi > dec * dec)
  {
    std::cout << "ERROR! Input exchange is not kinematically allowed. Quitting... \n";
    exit(0);
  }

  triangle tri(dec);

  std::vector<double> s;
  vector< std::complex<double> > feyn, disp;

  for (int i = 0; i < 160; i++)
  {
    double si = low + double(i) * (high - low) / 160.;
    complex<double> fx_f = tri.eval_feynman(si, ex * ex);
    complex<double> fx_d = tri.eval_dispersive(si, ex * ex);

    s.push_back(si);
    feyn.push_back(fx_f);
    disp.push_back(fx_d);
  }

  quick_plot(s, feyn, "feynman");
  quick_plot(s, disp, "dispersive");

  return 1.;
};
