// Test code for plotting both at the same time
//
// Produces .dat and .pdf plotting the triangle for both the feynman evaluation
// or dispersive evaluation for comparison.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman/feynman_triangle.hpp"
#include "dispersive/dispersive_triangle.hpp"
#include "quantum_numbers.hpp"

#include "triangle_plotter.hpp"
#include "jpacGraph1Dc.hpp"

#include <cstring>
#include <string>
#include <ctime>

using namespace jpacTriangle;

int main( int argc, char** argv )
{
  // Desired quantum numbers
  int id = 0;
  int Np = 100;
  
  // Parse inputs
  for (int i = 0; i < argc; i++)
  {
    if (std::strcmp(argv[i],"-id")==0) id = atoi(argv[i+1]);
    if (std::strcmp(argv[i],"-n")==0)  Np = atoi(argv[i+1]);
  }

  // All the associated quantum numbers for the amplitude
  quantum_numbers qns;
  qns.n = 1;
  qns.set_id(id);
  qns.mDec = .780; // The decaying particle mass
//   qns.mDec = .5478; // eta decay

  // Initialize a triangle object
  feynman_triangle tri_feyn(&qns);
  dispersive_triangle tri_disp(&qns);

  // Choose the name for the output files to have (sans and extentions)
  std::string filename = "compare.pdf";

  // Plotting bounds
  double low = EPS;
  double high = 81. * mPi2;

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Calculate the triangle function first with the feynman representation

  std::cout << "\n";
  std::cout << "Calculating triangles... id = " << id << "\n\n";

  std::cout << std::left;
  std::cout << std::setw(7)  << "i";
  std::cout << std::setw(15) << "sqrt(s)/mPi";
  std::cout << std::setw(30) << "disp";
  std::cout << std::setw(30) << "feynman";
  std::cout << std::setw(15) << "abs(disp - feyn)" << std::endl;

  std::vector<double> s;
  std::vector< std::complex<double> > tri[2];

  clock_t begin = clock();

  for (int i = 0; i <= Np; i++)
  {
    double si = low + EPS + double(i) * (high - low) / double(Np);
    s.push_back(sqrt(si) / mPi);

    std::complex<double> fx_f = tri_feyn.eval(si, mRho2);
    tri[0].push_back(fx_f);

    std::complex<double> fx_d = tri_disp.eval(si, mRho2);
    tri[1].push_back(fx_d);

    std::cout << std::left;
    std::cout << std::setw(7)  << i;
    std::cout << std::setw(15) << sqrt(si) / mPi;
    std::cout << std::setw(30) << fx_d;
    std::cout << std::setw(30) << fx_f;
    std::cout << std::setw(15) << std::abs(fx_d - fx_f) << std::endl;
  }

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  std::cout << "\nDone in " << elapsed_secs << " seconds. \n";
  std::cout << "\n";

  triangle_plotter* plotter = new triangle_plotter();
  plotter->SetVerticals(2., qns.mDec / mPi - 1., qns.mDec / mPi + 1.);
  plotter->AddEntry(s, tri);
  plotter->SetXaxis("#sqrt{s} / m_{#pi}", 0, 9);

  switch (id)
  {
    case 0: 
    {
      plotter->SetLegend(0.7, 0.5);
      plotter->SetYaxis("", -1.62, 2.);
      break;
    };
    case 1: 
    {
      plotter->SetLegend(0.7, 0.2);
      plotter->SetYaxis("", -0.1, 2.3);
      break;
    };
    case 10:
    {
      plotter->SetLegend(0.7, 0.15);
      plotter->SetYaxis("", -0.15, 0.3);
      break;
    };
    case 11:
    {
      plotter->SetLegend(0.64, 0.7);
      plotter->SetYaxis("", -0.02, 0.8);
      break;
    };
    case -11111:
    {
      plotter->SetLegend(0.34, 0.65);
      plotter->SetYaxis("", -0.1, 1.2);
      break;
    };
  }
  plotter->Plot(filename);

  delete plotter;

  return 1.;
};
