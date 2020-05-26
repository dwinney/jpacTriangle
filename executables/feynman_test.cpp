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
#include "feynman/feynman_triangle.hpp"
#include "constants.hpp"

#include "jpacGraph1Dc.hpp"

#include <cstring>
#include <string>
#include <ctime>

int main( int argc, char** argv )
{
  // Desired quantum numbers
  int j = 0, jp = 0;

  // Parse inputs
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-j")==0) j = atof(argv[i+1]);
    if (std::strcmp(argv[i],"-jp")==0) jp = atof(argv[i+1]);
  }

  // All the associated quantum numbers for the amplitude
  quantum_numbers qns;
  qns.n = 1;
  qns.j = j;
  qns.jp = jp;
  qns.mDec = .780; // The decaying particle mass

  // Start a breit_wigner object for rho exchange.
  // This will serve as our LHC being convoluted in the triangle diagram
  breit_wigner left_hand_cut(.770, .145);

  // Initialize a triangle object passing the above propogator
  // and the mass of the decaying particle
  feynman_triangle tri(&qns, &left_hand_cut);


  // Choose the name for the output files to have (sans and extentions)
  std::string filename = "omega_feyn.pdf";

  // Plotting bounds
  double low = 0.;
  double high = 81. * mPi2;

  int Np = 50; // Number of points to plot

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Calculate the triangle function first with the feynman representation

  std::cout << "\n";
  std::cout << "Calculating Feynman triangle... \n\n";

  std::vector<double> s;
  std::vector< std::complex<double> > feyn;

  clock_t begin = clock();

  for (int i = 0; i <= Np; i++)
  {
    double si = low + EPS + double(i) * (high - low) / double(Np);
    // std::complex<double> fx_f = tri.eval(si);
    std::complex<double> fx_f = tri.kernel(si, mRho2);

    s.push_back(sqrt(si) / mPi);
    feyn.push_back(fx_f);

    std::cout << std::left;
    std::cout << std::setw(7)  << i;
    std::cout << std::setw(15) << si;
    std::cout << std::setw(30) << fx_f << std::endl;
  }

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  std::cout << "\nDone in " << elapsed_secs << " seconds. \n";
  std::cout << "\n";

  jpacGraph1Dc* plotter = new jpacGraph1Dc();
  plotter->AddEntry(s, feyn, "");
  plotter->SetLegend(false);

  plotter->SetXaxis("#sqrt{s} / m_{#pi}", 0, 10);
  plotter->Plot(filename);

  delete plotter;

  return 1.;
};
