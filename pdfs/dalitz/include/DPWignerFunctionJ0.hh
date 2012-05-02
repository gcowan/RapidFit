#ifndef DP_WIGNER_FUNCTION_J0
#define DP_WIGNER_FUNCTION_J0

#include "DPWignerFunction.hh"

class DPWignerFunctionJ0 : public virtual DPWignerFunction
{
  public:

    DPWignerFunctionJ0() {};
    ~DPWignerFunctionJ0() {};

    double function(double cosTheta, int m, int n);

  protected:

  private:

};

#endif
