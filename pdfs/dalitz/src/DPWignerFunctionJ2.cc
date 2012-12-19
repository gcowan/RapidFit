#include "TMath.h"
#include <iostream>

#include "DPWignerFunctionJ2.hh"

double DPWignerFunctionJ2::function(double cosTheta, double mm, double nn)
{

  int m=(int)mm;
  int n=(int)nn;
  if ( m == 0 )
  {
    switch (n)
    {
      case 0: return d00(cosTheta);
              break;
      case 1: return -dp10(cosTheta);
              break;
      case -1: return dp10(cosTheta);
              break;
      case 2: return dp20(cosTheta);
              break;
      case -2: return dp20(cosTheta);
              break;
    }    
  }
  else if ( m == 1 )
  {
    switch (n)
    {
      case 0: return dp10(cosTheta);
              break;
      case 1: return dp1p1(cosTheta);
              break;
      case -1: return dp1m1(cosTheta);
              break;
      case 2: return -dp2p1(cosTheta);
              break;
      case -2: return dp2m1(cosTheta);
              break;
    }    
  }
  else if ( m == -1 )
  {
    switch (n)
    {
      case 0: return -dp10(cosTheta);
              break;
      case 1: return dp1m1(cosTheta);
              break;
      case -1: return dp1p1(cosTheta);
              break;
      case 2: return -dp2m1(cosTheta);
              break;
      case -2: return dp2p1(cosTheta);
              break;
    }    
  }
  else if ( m == 2 )
  {
    switch (n)
    {
      case 0: return dp20(cosTheta);
              break;
      case 1: return dp2p1(cosTheta);
              break;
      case -1: return dp2m1(cosTheta);
              break;
      case 2: return dp2p2(cosTheta);
              break;
      case -2: return dp2m2(cosTheta);
              break;
    }    
  }
  else if ( m == -2 )
  {
    switch (n)
    {
      case 0: return dp20(cosTheta);
              break;
      case 1: return -dp2m1(cosTheta);
              break;
      case -1: return -dp2p1(cosTheta);
              break;
      case 2: return dp2m2(cosTheta);
              break;
      case -2: return dp2p2(cosTheta);
              break;
    }    
  }

  std::cerr<<"What is going on? For spin 2, m and n has to be 0, 1, -1, 2 or -2\n";

  return -1000; // Give crazy number, alternatively we can exit or throw exception
}

double DPWignerFunctionJ2::d00(double cosTheta)
{
  return (1.5*cosTheta*cosTheta-0.5);
}
 
double DPWignerFunctionJ2::dp10(double cosTheta)
{
  double sinTheta=TMath::Sqrt(1-cosTheta*cosTheta);
  return -TMath::Sqrt(1.5)*sinTheta*cosTheta;
}

double DPWignerFunctionJ2::dp1p1(double cosTheta)
{
  return 0.5*(1+cosTheta)*(2*cosTheta-1);
}

double DPWignerFunctionJ2::dp1m1(double cosTheta)
{
  return 0.5*(1-cosTheta)*(2*cosTheta+1);
}

double DPWignerFunctionJ2::dp2p2(double cosTheta)
{
  double result = 0.5*(1+cosTheta);
  
  return result*result;
}

double DPWignerFunctionJ2::dp2p1(double cosTheta)
{
  double sinTheta=TMath::Sqrt(1-cosTheta*cosTheta);

  return -0.5*(1+cosTheta)*sinTheta;
}

double DPWignerFunctionJ2::dp20(double cosTheta)
{
  double sinTheta=TMath::Sqrt(1-cosTheta*cosTheta);

  return TMath::Sqrt(6)/4.*sinTheta*sinTheta;
}

double DPWignerFunctionJ2::dp2m1(double cosTheta)
{
  double sinTheta=TMath::Sqrt(1-cosTheta*cosTheta);

  return -0.5*(1-cosTheta)*sinTheta;
}

double DPWignerFunctionJ2::dp2m2(double cosTheta)
{
  double result = 0.5*(1-cosTheta);

  return result*result;
}
