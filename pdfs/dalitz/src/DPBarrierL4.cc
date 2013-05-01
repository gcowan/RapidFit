#include "DPBarrierL4.hh"
#include "TMath.h"

DPBarrierL4::DPBarrierL4(double radius) : DPBarrierFactor(radius)
{
}

double DPBarrierL4::barrier(double p0, double p)
{
  double z  = p*p*R()*R();
  double z0 = p0*p0*R()*R();

  double factor = function(z0)/function(z);
  return factor;
}

double DPBarrierL4::function(double z)
{
  return TMath::Sqrt(z*z*z*z+10*z*z*z+135*z*z+1575*z+11025);
}
