#include "TMath.h"
#include <iostream>

#include "DPWignerFunctionJ1.hh"

// parameter theta is actually cos(theta)
double DPWignerFunctionJ1::function(double theta, double mm, double nn)
{

  int m=(int)mm;
  int n=(int)nn;
  if ( m == 0 )
  {
    switch (n)
    {
      case 0: return d00(theta);
              break;
      case 1: return -dp10(theta);
              break;
      case -1: return dp10(theta);
              break;
    }    
  }
  else if ( m == 1 )
  {
    switch (n)
    {
      case 0: return dp10(theta);
              break;
      case 1: return dp1p1(theta);
              break;
      case -1: return dp1m1(theta);
              break;
    }    
  }
  else if ( m == -1 )
  {
    switch (n)
    {
      case 0: return -dp10(theta);
              break;
      case 1: return dp1m1(theta);
              break;
      case -1: return dp1p1(theta);
              break;
    }    
  }

  std::cerr<<"What is going on? For spin 1, m and n has to be 0, 1 or -1\n";

  return -1000; // Give crazy number, alternatively we can exit or throw exception
}

double DPWignerFunctionJ1::d00(double cosTheta)
{
  return cosTheta;
}
 
double DPWignerFunctionJ1::dp10(double cosTheta)
{
  double sinTheta=TMath::Sqrt(1-cosTheta*cosTheta);
  return -sinTheta/TMath::Sqrt(2);
}

double DPWignerFunctionJ1::dp1p1(double cosTheta)
{
  return 0.5*(1+cosTheta);
}

double DPWignerFunctionJ1::dp1m1(double cosTheta)
{
  return 0.5*(1-cosTheta);
}

