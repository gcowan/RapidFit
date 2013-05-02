#ifndef DP_BARRIER_L5
#define DP_BARRIER_L5

#include "DPBarrierFactor.hh"

class DPBarrierL5: public virtual DPBarrierFactor
{
  public:
    DPBarrierL5(double radius);
    ~DPBarrierL5() {};

    double barrier(double p0, double p);

  private:

    static double function(double z);
};

#endif
