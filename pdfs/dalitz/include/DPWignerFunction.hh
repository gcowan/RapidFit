#ifndef DP_WIGNER_FUNCTION
#define DP_WIGNER_FUNCTION


class DPWignerFunction
{
  public:

    DPWignerFunction() {};
    virtual ~DPWignerFunction() {};

    virtual double function(double cosTheta, double m, double n) = 0;

  protected:

  private:

};

#endif
