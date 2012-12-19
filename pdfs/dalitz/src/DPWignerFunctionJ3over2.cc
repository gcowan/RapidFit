#include "TMath.h"
#include <iostream>

#include "DPWignerFunctionJ3over2.hh"

// parameter theta is actually cos(theta)
double DPWignerFunctionJ3over2::function(double cosTheta, double mm, double nn)
{

  int m=(int)2*mm;
  int n=(int)2*nn;
  double theta=TMath::ACos(cosTheta);
  if ( m == -3 )
  {
    switch (n)
    {
      case 3: return -dp3m3(theta);
              break;
      case 1: return -dp3m1(theta);
              break;
      case -1: return -dp3p1(theta);
              break;
      case -3: return dp3p3(theta);
              break;
    }    
  }
  else if ( m == -1 )
  {
    switch (n)
    {
      case 3: return -dp3m1(theta);
              break;
      case 1: return -dp1m1(theta);
              break;
      case -1: return dp1p1(theta);
              break;
      case -3: return dp3p1(theta);
              break;
    }    
  }
  else if ( m == 1 )
  {
    switch (n)
    {
      case 3: return -dp3p1(theta);
              break;
      case 1: return dp1p1(theta);
              break;
      case -1: return dp1m1(theta);
              break;
      case -3: return dp3m1(theta);
              break;
    }    
  }
  else if ( m == 3 )
  {
    switch (n)
    {
      case 3: return dp3p3(theta);
              break;
      case 1: return dp3p1(theta);
              break;
      case -1: return dp3m1(theta);
              break;
      case -3: return dp3m3(theta);
              break;
    }    
  }

  return -1000; // Give crazy number, alternatively we can exit or throw exception
}

double DPWignerFunctionJ3over2::dp3p3(double theta)
{
  return (1.+TMath::Cos(theta))/2.*TMath::Cos(theta/2);
}

double DPWignerFunctionJ3over2::dp3m3(double theta)
{
  return -(1.-TMath::Cos(theta))/2.*TMath::Sin(theta/2);
}

double DPWignerFunctionJ3over2::dp3p1(double theta)
{
  return -TMath::Sqrt(3)*(1.+TMath::Cos(theta))/2.*TMath::Sin(theta/2);
}

double DPWignerFunctionJ3over2::dp3m1(double theta)
{
  return TMath::Sqrt(3)*(1.-TMath::Cos(theta))/2.*TMath::Cos(theta/2);
}

double DPWignerFunctionJ3over2::dp1p1(double theta)
{
  return (3*TMath::Cos(theta)-1)/2.*TMath::Cos(theta/2);
}

double DPWignerFunctionJ3over2::dp1m1(double theta)
{
  return -(3*TMath::Cos(theta)+1)/2.*TMath::Sin(theta/2);
}


