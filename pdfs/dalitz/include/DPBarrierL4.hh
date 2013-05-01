#ifndef DP_BARRIER_L4
#define DP_BARRIER_L4

#include "DPBarrierFactor.hh"

class DPBarrierL4: public virtual DPBarrierFactor
{
  public:
    DPBarrierL4(double radius);
    ~DPBarrierL4() {};

    double barrier(double p0, double p);

  private:

    static double function(double z);
};

#endif
