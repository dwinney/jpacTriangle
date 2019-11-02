// This object calculates the integration along the complex plane
// required to analytically continue the projection of cross-channel quantities
// to the direct channel
//
// Seperated as its own object for easier debugging
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "t_integral.hpp"

// ---------------------------------------------------------------------------
// Evaluate the cross channel projection integral.
// This function alls the approriate integration routine depending
// on the value of s which determines which part of the complex path and therefore
// whether the bounds of integration are complex or real
complex<double> t_integral::operator() (complex<double> s, complex<double> t)
{
  complex<double> result = 0.;
  double sr = real(s);
  if (sr >= sthPi + EPS && sr < c)
  {
    result = integ_sthPi_c(s, t);
  }
  else if (sr >= c && sr <= a)
  {
    result = integ_c_a(s, t);
  }
  else if (sr > a && sr < b)
  {
    result = integ_a_b(s, t);
  }
  else if (sr >= b)
  {
    result = integ_b(s, t);
  }
  else
  {
    cout << "s = " << sr << "\t aaaaaaaaaaaa! " << endl;
    exit(0);
  }

  result /= Kacser(s);
  return result;
};

// ---------------------------------------------------------------------------
// Dynamical left-hand cut function
complex<double> t_integral::integrand(complex<double>  t, complex<double> tp)
{
  return xr / (t - xi * .15 - tp);
};

// Kacser function which includes the correct analytic structure of
// product of breakup momenta and phase space factor.
complex<double> t_integral::Kacser( complex<double> s)
{
  complex<double> result;
  result = sqrt(xr * (s - sthPi - EPS + ieps)/ s); //includes phase-space factor
  result *= sqrt(xr * (b - s)) * sqrt(xr * (a - s));

  return result;
};


// ---------------------------------------------------------------------------
// Complex Bounds of integration
complex<double> t_integral::t_minus(complex<double> s)
{
  return (mDec * mDec + 3. * mPi * mPi - s) / 2. - Kacser(s) / 2.;
};

complex<double> t_integral::t_plus(complex<double> s)
{
  return (mDec * mDec + 3. * mPi * mPi - s) / 2. + Kacser(s) / 2.;
};

// ---------------------------------------------------------------------------
// Integration path in the complex plane

// ---------------------------------------------------------------------------
// Both limits are above the unitarity cut
complex<double> t_integral::integ_sthPi_c(complex<double> s, complex<double> t)
{
  if (real(s) <= sthPi || real(s) >= c)
  {
    cout << "integ_s0_a0: Integration out of range! Quitting... \n";
    exit(1);
  }

  if (real(t) >= real(t_minus(s)) && real(t) <= real(t_plus(s)))
  {
    complex<double> temp;
    temp = log(t - t_plus(s));
    temp -= log(t - t_minus(s));
    return temp;
  }
  else
  {
    double w[xN + 1], x[xN + 1];
    gauleg(real(t_minus(s)), real(t_plus(s)), x, w, xN);

    complex<double> sum = 0.;
    for (int i = 1; i <= xN; i++)
    {
      sum += w[i] * integrand(t, x[i]);
    }
    return sum;
  }
};

// ---------------------------------------------------------------------------
// both limits are purely real but one is above the other below
complex<double> t_integral::integ_c_a(complex<double> s, complex<double> t)
{
  if (real(s) < c || real(s) > a )
  {
    cout << "integ_c_a: Integration out of range! Quitting... \n";
    exit(1);
  }

  complex<double> sumM, sumP;

  if (std::abs(t_minus(s) - (sthPi + EPS)) > 0.001)
  {
    if (real(t) >= real(t_minus(s)) && real(t) <= sthPi + EPS)
    {
      sumM = log(t - (sthPi + EPS));
      sumM -= log(t - t_minus(s));
    }
    else
    {
      double wM[xN + 1], xM[xN + 1];
      gauleg(real(t_minus(s)), sthPi + EPS, xM, wM, xN);

      for(int i = 1; i <= xN; i++)
      {
        sumM += wM[i] * integrand(t, xM[i]);
      }
    }
  }

  if (std::abs(t_plus(s) - (sthPi + EPS)) > 0.001)
  {
    if (real(t) >= sthPi + EPS && real(t) <= real(t_plus(s)))
    {
      sumP = log(t - t_plus(s));
      sumP -= log(t - (sthPi + EPS));
    }
    else
    {
      double wP[xN + 1], xP[xN + 1];
      gauleg(sthPi + EPS, real(t_plus(s)), xP, wP, xN);

      for(int i = 1; i <= xN; i++)
      {
        sumP += wP[i] * integrand(t, xP[i]);
      }
    }
  }

  return sumM + sumP;
};

// ---------------------------------------------------------------------------
// this is the unphysical region, the bounds of integration are complex to avoid singularities
complex<double> t_integral::integ_a_b(complex<double> s, complex<double> t)
{
  if (real(s) < a || real(s) > b)
  {
    cout << "integ_a_b: Integration out of range! Quitting... \n";
    exit(1);
  }

  double w[xN + 1], x[xN + 1];
  gauleg(0., 1., x, w, xN);

  complex<double> sumM = 0., sumP = 0.;
  for(int i = 1; i <= xN; i++)
  {
    complex<double> z1_i = (1. - x[i]) * t_minus(s) + x[i] * (2. * t_minus(b));
    sumM +=  w[i] * (2. * t_minus(b) - t_minus(s)) * integrand(t, z1_i);

    complex<double> z2_i = (1. - x[i]) *  (2. * t_plus(b)) + x[i] * t_plus(s);
    sumP +=  w[i] * (t_plus(s) - 2. * t_plus(b)) * integrand(t, z2_i);
  }

  return (sumP + sumM);
};


// t-channel scattering region, limits are real again
complex<double> t_integral::integ_b(complex<double> s, complex<double> t)
{
  if (real(s) < b)
  {
    cout << "integ_b: Integration out of range! Quitting... \n";
    exit(1);
  }

  double w[xN + 1], x[xN + 1];
  gauleg(real(t_minus(s)), real(t_plus(s)), x, w, xN);

  complex<double> sum = 0.;
  for (int i = 1; i <= xN; i++)
  {
    sum += w[i] * integrand(t, x[i]);
  }

  return sum;
};
