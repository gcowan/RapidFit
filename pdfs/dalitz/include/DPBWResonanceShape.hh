#ifndef DP_BW_RESONANCE_SHAPE
#define DP_BW_RESONANCE_SHAPE

#include "DPMassShape.hh"
#include "DPBarrierFactor.hh"

#include "TComplex.h"

class DPBWResonanceShape: public virtual DPMassShape
{

  public:

    DPBWResonanceShape(double mR, double gammaR, int L, double m1, 
                       double m2, double R);
    DPBWResonanceShape( const DPBWResonanceShape& );
    ~DPBWResonanceShape();

    TComplex massShape(double m);

    void setParameters(double* pars){};

  private:

    double mR;
    double gammaR;
    int LR;
    double m1;
    double m2; 
    double R; 
    DPBarrierFactor* barrier;
    double pR0;  // Momentum of daughters at mR

    double gamma(double m);
    double daughterMomentum(double m);
    void setResonanceParameters( double mass, double sigma );
};

#endif
