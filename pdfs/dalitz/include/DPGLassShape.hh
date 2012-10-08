#ifndef DP_GLASS_SHAPE
#define DP_GLASS_SHAPE

#include "DPMassShape.hh"
#include "DPBarrierFactor.hh"

#include "TComplex.h"

class DPGLassShape: public virtual DPMassShape
{

  public:

    DPGLassShape(double mR, double gammaR, int L, double m1, 
                       double m2, double R, double a, double r);
    DPGLassShape( const DPGLassShape& );
    ~DPGLassShape();

    TComplex massShape(double m);

    void setResonanceParameters(double a, double r);
    void setParameters(double* pars);

  private:

    double mR;
    double gammaR;
    int LR;
    double m1;
    double m2; 
    double a;
    double r;
    double R;
    DPBarrierFactor* barrier;
    double pR0;  // Momentum of daughters at mR

    double fraction;
    double phaseR;
    double phaseB;

    double gamma(double m);
    double daughterMomentum(double m);
};

#endif
