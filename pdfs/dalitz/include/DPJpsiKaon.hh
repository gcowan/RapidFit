#ifndef DP_JPSIKAON
#define DP_JPSIKAON

#include <string>

#include "TComplex.h"
#include "DPComponent.hh"
#include "DPWignerFunctionJ0.hh"
#include "DPWignerFunctionJ1.hh"
#include "DPWignerFunctionJ2.hh"

class DPJpsiKaon: public DPComponent
{

  public:

    DPJpsiKaon(int LB, int LR, double mB, double mR, double gammaR,
                   double m1, double m2, double RB, double RR,
                   double mJpsi, int spin, std::string massShape="BW",
                   double a=1, double r=0);
    DPJpsiKaon( const DPJpsiKaon& input );
    ~DPJpsiKaon();

    TComplex amplitude(double m23, double cosTheta1, double cosTheta2,
                       double phi, int twoLambda, int twoLambdaPsi, int pionID);

    TComplex amplitudeProperVars(double m23, double cosTheta1, double cosTheta2,
				 double phi, int pionID, int twoLambda, int twoLambdaPsi);

    void setHelicityAmplitudes(double magA0, double magAplus,
                     double magAminus, double phaseA0, double phaseAplus,
                     double phaseAminus);

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
   double a;
   double r;
   DPMassShape* massShape;

    double mJpsi;
    double m1;  // Should be kaon mass
    double m2;  // Should be pion mass

    // Three helicity amplitudes
    TComplex A0;
    TComplex Aplus;
    TComplex Aminus;

    //
    int spinKaon;
    std::string mShape;

    // Wigner function for our calculations
    DPWignerFunction* wigner;
    DPWignerFunctionJ1 wignerPsi;
};

#endif
