#ifndef DP_WIGNER_FUNCTION_J1over2
#define DP_WIGNER_FUNCTION_J1over2


#include "DPWignerFunction.hh"

class DPWignerFunctionJ1over2 : public virtual DPWignerFunction
{
  public:

    DPWignerFunctionJ1over2() {};
    ~DPWignerFunctionJ1over2() {};

// This will take helicity*2 rather than helicity used by integer spins
    double function(double cosTheta, double m, double n);

  protected:

  private:

    static double dp(double cosTheta);
    static double dm(double cosTheta);
   
};

#endif
