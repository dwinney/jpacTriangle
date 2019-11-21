// Test code for dispersive triangle with spin
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "breit_wigner.hpp"
#include "dispersive_triangle.hpp"

#include <iostream>
#include <string>
#include <ctime>

int main()
{
  // Start a breit_wigner object for rho exchange.
  // This will serve as our LHC being convoluted in the triangle diagram
  breit_wigner left_hand_cut(.770, .145);

  // Initialize a triangle object passing the above propogator
  dispersive_triangle tri(&left_hand_cut);

  // Set the two external particles: omega and pion
  // here argument 1 >= argument 2
  tri.set_externalMasses(0.780, mPi);

  // Set the two other intermediate particles, in this case both pions
  tri.set_internalMasses(mPi, mPi);

  // Choose the name for the output files to have (sans and extentions)
  std::string filename = "omega";

  // Plotting bounds
  double low = 1.e-3;
  double high = 81. * mPi * mPi;

  int Np = 100; // Number of points to plot

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Calculate the triangle function first with the feynman representation

  std::cout << "\n";
  std::cout << "Calculating Dispersive triangle... \n";

  std::vector<double> s;
  std::vector< std::complex<double> > disp;

  clock_t begin = clock();

  std::complex<double> fxd_0 = tri.eval(low);   // Normalization

  for (int i = 0; i < Np; i++)
  {
    double si = low + double(i) * (high - low) / double(Np);

    std::complex<double> fx_d = tri.eval(si) / fxd_0;

    s.push_back(sqrt(si) / mPi);
    disp.push_back(fx_d);
  }

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  std::cout << "Done in " << elapsed_secs << " seconds. \n";
  std::cout << "\n";

  // Plot the vectors and save them in .pdf and .dat
  quick_plot(s, disp, filename + "_disp");
  quick_print(s, disp, filename + "_disp");

  std::cout << "\n";

  return 1.;
};
