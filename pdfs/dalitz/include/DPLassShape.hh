#ifndef DP_LASS_SHAPE
#define DP_LASS_SHAPE

#include "DPMassShape.hh"
#include "DPBarrierFactor.hh"

#include "TComplex.h"

class DPLassShape: public virtual DPMassShape
{

  public:

    DPLassShape(double mR, double gammaR, int L, double m1, 
                       double m2, double R, double a, double r);
    ~DPLassShape();

    TComplex massShape(double m);

    void setResonanceParameters(double a, double r);

  private:

    double mR;
    double gammaR;
    double LR;
    double m1;
    double m2; 
    double a;
    double r;
    DPBarrierFactor* barrier;
    double pR0;  // Momentum of daughters at mR

    double gamma(double m);
    double daughterMomentum(double m);
};

#endif
