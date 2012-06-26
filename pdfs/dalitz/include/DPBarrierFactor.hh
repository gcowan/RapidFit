#ifndef DP_BARRIER_FUNCTION
#define DP_BARRIER_FUNCTION


class DPBarrierFactor
{
  public:

    DPBarrierFactor(double radius);
    ~DPBarrierFactor() {};

    virtual double barrier(double p0, double p) = 0;

  protected:

    double R() {return radius;};

  private:

    double radius;  // Blatt-Weiskopf radius

};

#endif
