#include "DPBarrierL1.hh"
#include "TMath.h"

DPBarrierL1::DPBarrierL1(double radius) : DPBarrierFactor(radius)
{
 //std::cout << "DPBarrierL1 const " << std::endl;
}

double DPBarrierL1::barrier(double p0, double p)
{
  double z  = p*p*R()*R();
  double z0 = p0*p0*R()*R();

  double factor = function(z0)/function(z);
  return factor;
}

double DPBarrierL1::function(double z)
{
  return TMath::Sqrt(1+z);
}
