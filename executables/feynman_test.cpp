// Test code for feynman triangle with spin
//
// Produces .dat and .pdf plotting the triangle for both the feynman evaluation
// or dispersive evaluation for comparison.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "breit_wigner.hpp"
#include "feynman_triangle.hpp"

#include <iostream>
#include <string>
#include <ctime>

int main()
{
  // Start a breit_wigner object for rho exchange.
  // This will serve as our LHC being convoluted in the triangle diagram
  breit_wigner left_hand_cut(.770, .145);

  // Initialize a triangle object passing the above propogator
  feynman_triangle tri(&left_hand_cut);

  // Set the two external particles: omega and pion
  // here argument 1 >= argument 2
  tri.set_decayMass(0.780);

  // Choose the name for the output files to have (sans and extentions)
  std::string filename = "omega";

  // Plotting bounds
  double low = 1.e-3;
  double high = 16 * mPi * mPi;

  int Np = 100; // Number of points to plot

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Calculate the triangle function first with the feynman representation

  std::cout << "\n";
  std::cout << "Calculating Feynman triangle... \n";

  std::vector<double> s;
  std::vector< std::complex<double> > feyn;

  clock_t begin = clock();

  for (int i = 0; i < Np; i++)
  {
    double si = low + double(i) * (high - low) / double(Np);

    std::complex<double> fx_f = tri.eval(si);

    s.push_back(sqrt(si) / mPi);
    feyn.push_back(fx_f);
  }

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  std::cout << "Done in " << elapsed_secs << " seconds. \n";
  std::cout << "\n";

  // Plot the vectors and save them in .pdf and .dat
  quick_plot(s, feyn, filename + "_feyn");
  quick_print(s, feyn, filename + "_feyn");

  std::cout << "\n";

  return 1.;
};
