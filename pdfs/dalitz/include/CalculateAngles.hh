//////////////////////////////////////////////////////////
//
//
//         Compute angles for Z analysis
//
//
//////////////////////////////////////////////////////////

#ifndef CalculateAngles_h
#define CalculateAngles_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>



class CalculateAngles {
public :
  double decayAngleChi ( const TLorentzVector& d1 ,
			 const TLorentzVector& d2 ,
			 const TLorentzVector& h1 ,
			 const TLorentzVector& h2 );

  double decayAngle( const TLorentzVector& P ,
		     const TLorentzVector& Q ,
		     const TLorentzVector& D );

  static void calculateZplusAngles_GOLD(TLorentzVector& pB,
					TLorentzVector& pMuPlus,
					TLorentzVector& pMuMinus,
					TLorentzVector& pPi,
					TLorentzVector& pK,
					double* cosThetaZ,
					double* cosThetaPsi,
					double* dphi,
					int pion_ID);

   double HelCos_Tomasz(const double & ms,const double &mz,const int style);

   double coshel0_Tomasz(TLorentzVector particle, TLorentzVector parent, TLorentzVector grandparent);

   double coshel1_Tomasz(TLorentzVector particle, TLorentzVector parent,
			 TLorentzVector grandparent,TLorentzVector grandgrandparent);

   double planeAngle0_Tomasz(TLorentzVector particleFrame,
			     TLorentzVector particleA, TLorentzVector particleB,
			     TLorentzVector particleC, TLorentzVector particleD);

   Float_t planeAngle1_Tomasz(TLorentzVector particleFrame,
			      TLorentzVector particleA, TLorentzVector particleB,
			      TLorentzVector particleC, TLorentzVector particleD);

   static void calculateFinalStateMomenta_GOLD(double mB0, double m23, double mMuMu, double cosTheta1,
				   double cosTheta2, double phi, double mMuPlus, double mMuMinus,
				   double mPi, double mK,
				   int pion_ID,
				   TLorentzVector& pMuPlus,
				   TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK);

   static double daughterMomentum(double m, double m1, double m2);

   static double planeAngle2_BIN_GOLD(const TLorentzVector & _muPlus,
				      const TLorentzVector & _muMinus,
				      const TLorentzVector & _psi,
				      const TLorentzVector & _z,
				      const TLorentzVector & _b);

   void Calculates_Fit_Var_EPFL(const TLorentzVector & pMplus,
				const TLorentzVector & pMminus,
				const TLorentzVector & pPiminus,
				const TLorentzVector & pKplus,
				int pion_ID,
				double& m23,
				double& cosTheta1,
				double& cosTheta2,
				double& phi);

};

#endif
