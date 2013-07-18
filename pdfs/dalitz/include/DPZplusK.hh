#ifndef DP_ZPLUSK
#define DP_ZPLUSK

#include "TComplex.h"
#include "DPComponent.hh"
#include "DPWignerFunctionJ0.hh"
#include "DPWignerFunctionJ1.hh"
#include "DPWignerFunctionJ2.hh"

class DPZplusK: public DPComponent
{

  public:

    DPZplusK(int LB, int LR, double mB, double mR, double gammaR,
	     double m1, double m2, double RB, double RR,
	     double mJpsi, int spin, int resonanceIn=23);

    DPZplusK( const DPZplusK& input );
    virtual ~DPZplusK();

    TComplex amplitude(double m23, double cosTheta1, double cosTheta2,
                       double phi, int twoLambda, int twoLambdaPsi, int pionID);

    void setHelicityAmplitudes(double magA0, double magAplus,
                     double magAminus, double phaseA0, double phaseAplus,
                     double phaseAminus);

    TComplex amplitudeProperVars(double m13, double cosTheta1, double cosTheta2,
				 double phi, int pionID, int twoLambda, int twoLambdaPsi);

    void setResonanceParameters(double mass, double sigma);

  private:
   DPBarrierFactor* barrierB;
   DPBarrierFactor* barrierR;
   int LB;
   int LR;
   double RB;
   double RR;
   double mB;
   double mR;
   double gammaR;
   DPMassShape* massShape;

    double mJpsi;
    double m1;  // Should be kaon mass
    double m2;  // Should be pion mass

  // Three helicity amplitudes
  TComplex A0;
  TComplex Aplus;
  TComplex Aminus;

  //
    int spinZplus;
    int resonanceIn;

  // Wigner function for our calculations
  DPWignerFunction* wigner;
  DPWignerFunctionJ1 wignerPsi;
};

#endif
