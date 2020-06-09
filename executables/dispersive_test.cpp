// Test code for dispersive triangle with spin
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "quantum_numbers.hpp"
#include "constants.hpp"
#include "dispersive/dispersive_triangle.hpp"
#include "dispersive/projection_function.hpp"

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

  // All the associated quantum numbers for the amplitude
  quantum_numbers qns;
  qns.n = 1;
  qns.j = j;
  qns.jp = jp;
  qns.mDec = .780; // The decaying particle mass

  // Initialize a triangle object passing the above propogator
  dispersive_triangle tri(&qns);

  projection_function test(&qns);

  // Choose the name for the output files to have (sans and extentions)
  std::string filename = "omega_disp.pdf";

  // Plotting bounds
  double low = EPS;
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

    std::complex<double> fx_d = tri.eval(si, mRho2);

    s.push_back(sqrt(si) / mPi);
    disp.push_back(fx_d);

    std::cout << std::left;
    std::cout << std::setw(7)  << i;
    std::cout << std::setw(15) << sqrt(si) / mPi;
    std::cout << std::setw(30) << fx_d << std::endl;
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
