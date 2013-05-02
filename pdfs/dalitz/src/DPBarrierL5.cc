#include "DPBarrierL5.hh"
#include "TMath.h"

DPBarrierL5::DPBarrierL5(double radius) : DPBarrierFactor(radius)
{
}

double DPBarrierL5::barrier(double p0, double p)
{
  double z  = p*p*R()*R();
  double z0 = p0*p0*R()*R();

  double factor = function(z0)/function(z);
  return factor;
}

double DPBarrierL5::function(double z)
{
  return TMath::Sqrt(z*z*z*z*z+15*z*z*z*z+315*z*z*z+6300*z*z+99225*z+893025);
}
