#ifndef DP_NONRESONANT_SHAPE
#define DP_NONRESONANT_SHAPE

#include "DPMassShape.hh"

#include "TComplex.h"

class DPNonresonant: public virtual DPMassShape
{

  public:

    DPNonresonant();
    DPNonresonant(const DPNonresonant& other);
    ~DPNonresonant();

    TComplex massShape(double m);

    void setParameters(double* pars) { (void)pars; }

  private:

    void setResonanceParameters( double mass, double sigma ) { (void)mass; (void)sigma; }
};

#endif
