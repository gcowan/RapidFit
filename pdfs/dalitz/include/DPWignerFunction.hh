#ifndef DP_WIGNER_FUNCTION
#define DP_WIGNER_FUNCTION


class DPWignerFunction
{
  public:

    DPWignerFunction() {};
    ~DPWignerFunction() {};

    virtual double function(double cosTheta, int m, int n) = 0;

  protected:

  private:

};

#endif
