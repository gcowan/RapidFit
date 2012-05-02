#ifndef DP_JPSIKSPINONE
#define DP_JPSIKSPINONE

#include "TComplex.h"
#include "DPComponent.hh"
#include "DPWignerFunctionJ1.hh"

class DPJpsiKSpinOne: public DPComponent
{

  public:

    DPJpsiKSpinOne(int LB, int LR, double mB, double mR, double gammaR,
                   double m1, double m2, double RB, double RR, double mJpsi); 
    ~DPJpsiKSpinOne();

    TComplex amplitude(double m23, double cosTheta1, double cosTheta2, 
                       double phi, int twoLambda, int twoLambdaPsi);
  
    void setHelicityAmplitudes(double magA0, double magAplus, 
                     double magAminus, double phaseA0, double phaseAplus,
                     double phaseAminus);

    void setResonanceParameters( double mass, double sigma );

  private:
   
    double mJpsi;
    double m1;
    double m2;

    double pR0;

  // Three helicity amplitudes
  TComplex A0;
  TComplex Aplus;
  TComplex Aminus;

  // Wigner function for our calculations
  DPWignerFunctionJ1 wigner;
};

#endif
