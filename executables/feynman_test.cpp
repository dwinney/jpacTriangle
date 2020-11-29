// Test code for the feynman evaluation
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman/feynman_triangle.hpp"
#include "dispersive/dispersive_triangle.hpp"
#include "quantum_numbers.hpp"

#include "jpacGraph1Dc.hpp"

#include <cstring>
#include <string>
#include <ctime>

using namespace jpacTriangle;

int main( int argc, char** argv )
{
  // Desired quantum numbers
  int id = 0;
  int Np = 20;
  
  // Parse inputs
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-id")==0) id = atoi(argv[i+1]);
  }

  // All the associated quantum numbers for the amplitude
  quantum_numbers qns;
  qns.n = 1;
  qns.set_id(id);
  qns.mDec = .780; // The decaying particle mass

  // Initialize a triangle object
  feynman_triangle tri_feyn(&qns);

  // Choose the name for the output files to have (sans and extentions)
  std::string filename = "omega_feyn.pdf";

  // Plotting bounds
  double low = EPS;
  double high = 81. * mPi2;

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Calculate the triangle function first with the feynman representation

  std::cout << "\n";
  std::cout << "Calculating Feynman triangle... \n\n";

  std::vector<double> s;
  std::vector< std::complex<double> > feyn, disp;

  clock_t begin = clock();

  for (int i = 0; i <= Np; i++)
  {
    double si = low + EPS + double(i) * (high - low) / double(Np);
    s.push_back(sqrt(si) / mPi);

    std::complex<double> fx_f = tri_feyn.eval(si, mRho2);
    feyn.push_back(fx_f);

    std::cout << std::left;
    std::cout << std::setw(7)  << i;
    std::cout << std::setw(15) << sqrt(si) / mPi;
    std::cout << std::setw(30) << fx_f;
    std::cout << std::endl;
  }

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  std::cout << "\nDone in " << elapsed_secs << " seconds. \n";
  std::cout << "\n";

  jpacGraph1Dc* plotter = new jpacGraph1Dc();
  plotter->AddEntry(s, feyn, "feynman");

  plotter->SetLegend(false);

  plotter->SetXaxis("#sqrt{s} / m_{#pi}", 0, 10);
  plotter->Plot(filename);

  delete plotter;

  return 1.;
};
