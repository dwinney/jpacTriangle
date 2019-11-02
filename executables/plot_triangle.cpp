#include "triangle.hpp"
#include <iostream>

int main(int argc, char** argv)
{
  // if (argc != 4)
  // {
  //   std::cout << "Usage: \t";
  //   std::cout << "./plot_triangle <decay_mass> <exchange_mass> <exchange_width>" << std::endl;
  //   exit(0);
  // }
  // double dec = std::stod(argv[1]);
  // double ex = std::stod(argv[2]);
  // double Gamma = std::stod(argv[3]);

  // Plotting bounds
  double low = 1.e-3;
  double high = 10. * 0.14;
  int Np = 360;

  triangle tri;
  tri.set_exchange(0.7755, .002);
  tri.set_internal(.14, .14);
  tri.set_external(1.02, .14);

  std::vector<double> s;
  vector< std::complex<double> > feyn, disp;

  for (int i = 0; i < Np; i++)
  {
    double si = low + double(i) * (high - low) / double(Np);
    complex<double> fx_f = tri.eval_feynman(si);
    // complex<double> fx_d = tri.eval_dispersive(si);

    s.push_back(sqrt(si) / 0.14);
    feyn.push_back(fx_f);
    // disp.push_back(fx_d);
  }

  quick_plot(s, feyn, "feynman");
  // quick_plot(s, disp, "dispersive");

  return 1.;
};
