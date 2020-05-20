// Test code for dispersive triangle with spin
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "breit_wigner.hpp"
#include "dispersive_triangle.hpp"

#include "jpacGraph1Dc.hpp"

#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <ctime>

int main(int argc, char** argv)
{
  // Desired quantum numbers
  int j = 0, jp = 0;

  // Parse inputs
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-j")==0) j = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-jp")==0) jp = atof(argv[i+1]);
  }

  // Start a breit_wigner object for rho exchange.
  // This will serve as our LHC being convoluted in the triangle diagram
  breit_wigner left_hand_cut(.770, .145);

  // Initialize a triangle object passing the above propogator
  dispersive_triangle tri(&left_hand_cut);

  // Set the decaying particle
  tri.set_decayMass(0.780);

  // Choose the name for the output files to have (sans and extentions)
  std::string filename = "omega_disp.pdf";

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

  for (int i = 0; i < Np; i++)
  {
    double si = low + double(i) * (high - low) / double(Np);

    std::complex<double> fx_d = tri.eval(j, jp, si);

    s.push_back(sqrt(si) / mPi);
    disp.push_back(fx_d);
  }

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  std::cout << "Done in " << elapsed_secs << " seconds. \n";
  std::cout << "\n";

  jpacGraph1Dc* plotter = new jpacGraph1Dc();
  plotter->AddEntry(s, disp, "");
  plotter->SetLegend(false);

  plotter->SetXaxis("#sqrt{s} / m_{#pi}", 0, 10);
  plotter->Plot(filename);

  delete plotter;
  return 1.;
};
