#include "TMath.h"
#include <iostream>

#include "DPWignerFunctionJ1over2.hh"

// parameter theta is actually cos(theta)
double DPWignerFunctionJ1over2::function(double cosTheta, double mm, double nn)
{

  int m=(int)2*mm;
  int n=(int)2*nn;
  double theta=TMath::ACos(cosTheta);
  if ( m == 1 )
  {
    switch (n)
    {
      case 1: return dp(theta);
              break;
      case -1: return dm(theta);
              break;
    }
  }
  else if ( m == -1 )
  {
    switch (n)
    {
      case 1: return -dm(theta);
              break;
      case -1: return dp(theta);
              break;
    }
  }

  return -1000; // Give crazy number, alternatively we can exit or throw exception
}

double DPWignerFunctionJ1over2::dp(double theta)
{
  return TMath::Cos(theta/2);
}

double DPWignerFunctionJ1over2::dm(double theta)
{
  return -TMath::Sin(theta/2);
}


