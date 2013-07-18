#ifndef DP_COMPONENT
#define DP_COMPONENT

#include "TComplex.h"
#include "DPBarrierFactor.hh"
#include "DPMassShape.hh"

class DPComponent
{
  public:

    DPComponent(){};
    virtual ~DPComponent(){};

    virtual TComplex amplitude(double m23, double cosTheta1, double cosTheta2,
                               double phi, int twoLambda, int twoLambdaPsi=0, int pionID = -211) = 0;

    virtual TComplex amplitudeProperVars(double m13, double cosTheta1, double cosTheta2,
					 double phi, int pionID, int twoLambda, int twoLambdaPsi=0) = 0;

    virtual void setHelicityAmplitudes(double magA0, double magAplus,
                     double magAminus, double phaseA0, double phaseAplus,
                     double phaseAminus) = 0;

    virtual void setResonanceParameters(double mass, double sigma) = 0;
/*
  protected:

   DPBarrierFactor* barrierB;
   DPBarrierFactor* barrierR;
   int LB;
   int LR;
   double mB;
   double mR;
   DPMassShape* massShape;
*/
};

#endif

