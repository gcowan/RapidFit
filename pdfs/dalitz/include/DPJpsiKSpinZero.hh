#ifndef DP_JPSIKSPINZERO
#define DP_JPSIKSPINZERO

#include "TComplex.h"
#include "DPComponent.hh"
#include "DPWignerFunctionJ1.hh"
#include "DPWignerFunctionJ0.hh"

class DPJpsiKSpinZero: public DPComponent
{

  public:

    DPJpsiKSpinZero(int LB, int LR, double mB, double mR, double gammaR,
                   double m1, double m2, double RB, double RR, double mJpsi); 
    ~DPJpsiKSpinZero();

    TComplex amplitude(double m23, double cosTheta1, double cosTheta2, 
                       double phi, int twoLambda, int twoLambdaPsi);
  
// Only A0 with its phase have meaning here
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

  // Wigner function for our calculations
  DPWignerFunctionJ1 wigner1;
  DPWignerFunctionJ0 wigner0;
};

#endif
