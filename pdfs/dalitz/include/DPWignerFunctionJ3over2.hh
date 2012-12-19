#ifndef DP_WIGNER_FUNCTION_J3over2
#define DP_WIGNER_FUNCTION_J3over2


#include "DPWignerFunction.hh"

class DPWignerFunctionJ3over2 : public virtual DPWignerFunction
{
  public:

    DPWignerFunctionJ3over2() {};
    ~DPWignerFunctionJ3over2() {};

    double function(double cosTheta, double m, double n);

  protected:

  private:

    static double dp3p3(double cosTheta);
    static double dp3m3(double cosTheta);
    static double dp3p1(double cosTheta);
    static double dp3m1(double cosTheta);
    static double dp1p1(double cosTheta);
    static double dp1m1(double cosTheta);
   
};

#endif
