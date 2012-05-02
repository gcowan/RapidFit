#ifndef DP_BARRIER_L3
#define DP_BARRIER_L3

#include "DPBarrierFactor.hh"

class DPBarrierL3: public virtual DPBarrierFactor
{ 
  public:
    DPBarrierL3(double radius);
    ~DPBarrierL3() {};
 
    double barrier(double p0, double p);

  private:

    static double function(double z);
};

#endif 
