#include "triangle.hpp"
#include <iostream>
#include <string>

int main()
{

  // Initialize a triangle object
  triangle tri;

  // Set the exchange rho particle mass and width,
  tri.set_exchangeMass(0.770, .145);

  // Set the two other intermediate particles, in this case both pions
  tri.set_internalMass(mPi, mPi);

  // Set the two external particles: omega and pion
  tri.set_externalMasses(.780, mPi);

  // Choose the name for the output files to have (sans and extentions)
  std::string filename = "test";

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

  // Plotting bounds
  double low = sthPi+ EPS;
  double high = 1.;

  int Np = 100; // Number of points to plot

  // Normalization
  // complex<double> fxf_0 = tri.eval_feynman(low);
  complex<double> fxd_0 = tri.eval_dispersive(low);

  std::vector<double> s;
  vector< std::complex<double> > feyn, disp;
  for (int i = 0; i < Np; i++)
  {
    double si = low + double(i) * (high - low) / double(Np);

    // complex<double> fx_f = tri.eval_feynman(si) / fxf_0;
    complex<double> fx_d = tri.eval_dispersive(si) / fxd_0;

    s.push_back(sqrt(si));
    // feyn.push_back(fx_f);
    disp.push_back(fx_d);

    if ((i % 25) == 0)
    {
      cout << i << "/" << Np << endl;
    }
  }

  // // Plot the vectors and save them in .pdf and .dat
  // quick_plot(s, feyn, filename + "_feyn");
  // quick_print(s, feyn, filename + "_feyn");

  quick_plot(s, disp, filename + "_disp");
  quick_print(s, disp, filename + "_disp");

  return 1.;
};
