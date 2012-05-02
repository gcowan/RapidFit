#ifndef DP_BARRIER_L2
#define DP_BARRIER_L2

#include "DPBarrierFactor.hh"

class DPBarrierL2: public virtual DPBarrierFactor
{ 
  public:
    DPBarrierL2(double radius);
    ~DPBarrierL2() {};
 
    double barrier(double p0, double p);

  private:

    static double function(double z);
};

#endif 
