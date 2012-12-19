#ifndef DP_WIGNER_FUNCTION_J1
#define DP_WIGNER_FUNCTION_J1


#include "DPWignerFunction.hh"

class DPWignerFunctionJ1 : public virtual DPWignerFunction
{
  public:

    DPWignerFunctionJ1() {};
    ~DPWignerFunctionJ1() {};

    double function(double cosTheta, double m, double n);

  protected:

  private:

    static double d00(double cosTheta);
    static double dp10(double cosTheta);
    static double dp1p1(double cosTheta);
    static double dp1m1(double cosTheta);
   
};

#endif
