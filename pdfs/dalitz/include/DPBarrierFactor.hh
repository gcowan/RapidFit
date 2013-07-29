#ifndef DP_BARRIER_FUNCTION
#define DP_BARRIER_FUNCTION
#include <iostream>

class DPBarrierFactor
{
  public:

    DPBarrierFactor(double radius);
    virtual ~DPBarrierFactor() {};//std::cout << "DPBarrierFactor dest" << std::endl;};

    virtual double barrier(double p0, double p) = 0;

  protected:

    double R() {return radius;};

  private:

    double radius;  // Blatt-Weiskopf radius

};

#endif
