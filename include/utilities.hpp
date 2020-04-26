// Misc utility functions and classes
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// -----------------------------------------------------------------------------

#ifndef _UTIL_
#define _UTIL_

#include <complex>
#include <vector>
#include <iterator>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include <TH2.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TError.h>

using std::string;
using std::vector;
using std::cout;
using std::complex;

//-----------------------------------------------------------------------------
// given int n, outputs right english string
// "1st" vs "2nd" vs "3rd" etc...
// :)
std::string st_nd_rd(int n);

//-----------------------------------------------------------------------------
// Utility functions that take in a vector of complex<doubles> and return vectors of the same size
// containing only the real or imaginary parts.
std::vector<double> vec_real( std::vector<std::complex<double>> fx);
std::vector<double> vec_imag( std::vector<std::complex<double>> fx);

//-----------------------------------------------------------------------------
// Simple function to call ROOT to print a plot of real and imag parts
void quick_plot(vector<double> s, vector<complex<double>> fx, string filename);

// make a symmetric dalitz plot from a file
// optional parameter is whether to plot using white color scheme to see % difference
void quick_dplot(string file, bool DEV = false);

// Print a complex function to a .dat file
void quick_print(vector<double> s, vector<complex<double>> fx, string file);

#endif
