#include "DPBarrierL2.hh"
#include "TMath.h"

DPBarrierL2::DPBarrierL2(double radius) : DPBarrierFactor(radius)
{
}

double DPBarrierL2::barrier(double p0, double p)
{
  double z  = p*p*R()*R();
  double z0 = p0*p0*R()*R();

  double factor = function(z0)/function(z); 
  return factor; 
}

double DPBarrierL2::function(double z)
{
  return TMath::Sqrt(z*z+3*z+9);
}
