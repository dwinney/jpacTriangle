// Example executable to print out the real and imaginary parts of the triangle diagram
// for the omega to 3 pi kinematics with a fixed mass rho exchange LHC.
//
// Produces .dat and .pdf plotting the triangle for both the feynman evaluation
// or dispersive evaluation for comparison.
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "triangle.hpp"
#include "breit_wigner.hpp"

#include <iostream>
#include <string>
#include <ctime>

int main()
{
  // Start a breit_wigner object for rho exchange.
  // This will serve as our LHC being convoluted in the triangle diagram
  breit_wigner left_hand_cut(.770, .145);

  // Initialize a triangle object passing the above propogator
  triangle tri(&left_hand_cut);

  // Set the two external particles: omega and pion
  // here argument 1 >= argument 2
  tri.set_externalMasses(0.780, mPi);

  // Set the two other intermediate particles, in this case both pions
  tri.set_internalMasses(mPi, mPi);

  // Choose the name for the output files to have (sans and extentions)
  std::string filename = "omega";

  // Plotting bounds
  double low = 4. * mPi * mPi;
  double high = 81. * mPi * mPi;

  int Np = 100; // Number of points to plot

// ---------------------------------------------------------------------------
// You shouldnt need to change anything below this line
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Calculate the triangle function first with the feynman representation

  std::cout << "\n";
  std::cout << "Calulating Feynman triangle... \n";

  std::vector<double> s;
  std::vector< std::complex<double> > feyn;

  clock_t begin = clock();

  std::complex<double> fxf_0 = tri.eval_feynman(low);   // Normalization

  for (int i = 0; i < Np; i++)
  {
    double si = low + double(i) * (high - low) / double(Np);

    std::complex<double> fx_f = tri.eval_feynman(si) / fxf_0;

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
//
// // ---------------------------------------------------------------------------
// // Calculate the same thing now with the dispersive representation
//
//   std::cout << "Calulating Dispersive triangle... \n";
//
//   begin = clock();
//
//   std::vector<std::complex<double>> disp;
//   std::complex<double> fxd_0 = tri.eval_dispersive(low);
//
//   s.clear();
//   for (int i = 0; i < Np; i++)
//   {
//     double si = low + double(i) * (high - low) / double(Np);
//
//     std::complex<double> fx_d = tri.eval_dispersive(si) / fxd_0;
//
//     s.push_back(sqrt(si) / mPi);
//     disp.push_back(fx_d);
//   }
//
//   end = clock();
//   elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//
//   std::cout << "Done in " << elapsed_secs << " seconds. \n";
//   std::cout << "\n";
//
//   quick_plot(s, disp, filename + "_disp");
//   quick_print(s, disp, filename + "_disp");
//
//   std::cout << "\n";

// ---------------------------------------------------------------------------
  return 1.;
};
