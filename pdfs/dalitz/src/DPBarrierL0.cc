#include "DPBarrierL0.hh"

DPBarrierL0::DPBarrierL0(double radius) : DPBarrierFactor(radius)
{
}

double DPBarrierL0::barrier(double p0, double p)
{
  return 1; 
}
