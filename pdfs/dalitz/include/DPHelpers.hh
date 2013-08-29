#ifndef DPHELPERS_HH
#define DPHELPERS_HH

#include "TLorentzVector.h"

namespace DPHelpers
{
 void calculateFinalStateMomentaBelle(double m_b, double ms, double m_psi,
                            double cos2s, double cospi, double phi2s,
                            double m_mu,
                            double m_pi, double m_k,
                            TLorentzVector& p_mu_b, TLorentzVector& p_mu1_b,
                            TLorentzVector& p_pi_b, TLorentzVector& p_k_b);

 void Belle(const TLorentzVector & _pMuPlus,
            const TLorentzVector & _pMuMinus,
            const TLorentzVector & _pPi,
            const TLorentzVector & _pK,
            double & m23,
            double & cosKPi,
            double & cosPsi,
            double & phiKPiPsi,
            double & m13,
            double & cosZ,
            double & cosPsi_Z,
            double & phiPsiZ,
            double & phiZPsiPsi);

  double daughterMomentum(double mR, double m1, double m2);
  void calculateFinalStateMomenta(double mB0, double m23, double mMuMu, double cosTheta1,
				  double cosTheta2, double phi, int pion_ID, double mMuPlus, double mMuMinus,
				  double mPi, double mK, TLorentzVector& pMuPlus,
				  TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK);
  void calculateZplusAngles(TLorentzVector& pB, TLorentzVector& pMuPlus,
			    TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK,
			    double* cosThetaZ, double* cosThetaPsi, double* dphi, int pion_ID);
  double referenceAxisCosAngle(TLorentzVector& pB, TLorentzVector& pMuPlus,
			       TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK);
  double referenceAxisCosAngleMu(TLorentzVector& pB, TLorentzVector& pMuPlus,
				 TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK);
  double planeAngle2_BIN(const TLorentzVector & _muPlus,
			 const TLorentzVector & _muMinus,
			 const TLorentzVector & _psi,
			 const TLorentzVector & _z,
			 const TLorentzVector & _b);

};

#endif
