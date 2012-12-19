#ifndef DP_WIGNER_FUNCTION_GENERAL
#define DP_WIGNER_FUNCTION_GENERAL


#include "DPWignerFunction.hh"

class DPWignerFunctionGeneral : public virtual DPWignerFunction
{
  public:

    DPWignerFunctionGeneral() {};
    DPWignerFunctionGeneral(double J) {j=J;};
    ~DPWignerFunctionGeneral() {};

    double function(double cosTheta, double m, double n);

  protected:

  private:

    static long factorials[20];
    static long fact(int n);
   
    double j;
};

#endif
