#ifndef DP_NONRESONANT_SHAPE
#define DP_NONRESONANT_SHAPE

#include "DPMassShape.hh"

#include "TComplex.h"

class DPNonresonant: public virtual DPMassShape
{

  public:

    DPNonresonant(double mR, double gammaR, int L, double m1, 
                       double m2, double R);
    DPNonresonant( const DPNonresonant& );
    ~DPNonresonant();

    TComplex massShape(double m);

    void setParameters(double* pars){};

  private:

    double mR;
    double gammaR;
    int LR;
    double m1;
    double m2; 

    void setResonanceParameters( double mass, double sigma ){};
};

#endif
