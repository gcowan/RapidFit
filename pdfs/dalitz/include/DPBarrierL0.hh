#ifndef DP_BARRIER_L0
#define DP_BARRIER_L0

#include "DPBarrierFactor.hh"

class DPBarrierL0: public virtual DPBarrierFactor
{
  public:
    DPBarrierL0(double radius);
    ~DPBarrierL0() {};//std::cout << "DPBarrierL0 dest" << std::endl;};

    double barrier(double p0, double p);
};

#endif
