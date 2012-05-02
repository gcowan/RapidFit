#ifndef DPHELPERS_HH
#define DPHELPERS_HH

#include "TLorentzVector.h"

namespace DPHelpers
{
  double daughterMomentum(double mR, double m1, double m2);
  void calculateFinalStateMomenta(double mB0, double m23, double mMuMu, double cosTheta1, 
        double cosTheta2, double phi, double mMuPlus, double mMuMinus,
        double mPi, double mK, TLorentzVector& pMuPlus,
        TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK);
  void calculateZplusAngles(TLorentzVector& pB, TLorentzVector& pMuPlus,
        TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK,
        double* cosThetaZ, double* cosThetaPsi, double* dphi); 
  double referenceAxisCosAngle(TLorentzVector& pB, TLorentzVector& pMuPlus,
        TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK);
};

#endif
