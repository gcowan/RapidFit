#include "DPBarrierL3.hh"
#include "TMath.h"

DPBarrierL3::DPBarrierL3(double radius) : DPBarrierFactor(radius)
{
}

double DPBarrierL3::barrier(double p0, double p)
{
  double z  = p*p*R()*R();
  double z0 = p0*p0*R()*R();

  double factor = function(z0)/function(z); 
  return factor; 
}

double DPBarrierL3::function(double z)
{
  return TMath::Sqrt(z*z*z+6*z*z+45*z+225);
}
