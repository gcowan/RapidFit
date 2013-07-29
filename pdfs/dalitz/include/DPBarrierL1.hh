#ifndef DP_BARRIER_L1
#define DP_BARRIER_L1

#include "DPBarrierFactor.hh"

class DPBarrierL1: public virtual DPBarrierFactor
{
  public:
    DPBarrierL1(double radius);
    ~DPBarrierL1() {};//std::cout << "DPBarrierL1 dest" << std::endl;};

    double barrier(double p0, double p);

  private:

    static double function(double z);
};

#endif
