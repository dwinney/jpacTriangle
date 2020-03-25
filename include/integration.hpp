// Methods and utilities for integrating functions
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _INTEG_
#define _INTEG_

#include <vector>
#include <cmath>
#include <iostream>

//-----------------------------------------------------------------------------
// Numerical integration routine. Generates Gauss-Legendre weights and abscissas
void NR_gauleg(double x1, double x2, double x[], double w[], int n);

// Self contained class that generates weights and abscissas using NR_gauleg,
// stores them in a vector from 0...xN-1 instead of from 1...N
class gauleg
{
public:
  // Constructor
  gauleg()
  {
    check_weights();
  };

  gauleg(int n)
  : xN(n)
  {
    check_weights();
  };

  // Number of Gaussian points
  int xN = 200;
  
  // Stored weights and abscissas
  std::vector<double> weights, abscissas;

  // Check if weights are saved
  void check_weights();
private:
  // Check
  bool WG_GENERATED = false;
};

#endif
