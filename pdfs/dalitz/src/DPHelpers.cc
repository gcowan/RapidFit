#include "DPHelpers.hh"
#include "TMath.h"
#include "TLorentzVector.h"
#include "CalculateAngles.hh"
#include <iostream>

double DPHelpers::daughterMomentum(double m, double m1, double m2)
{
  double momentum;

  momentum=(m*m-(m1+m2)*(m1+m2))*(m*m-(m1-m2)*(m1-m2));
  momentum=TMath::Sqrt(momentum);
  momentum/=2*m;

  return momentum;
}


/*
 * Function which takes three angular variables, invariant mass of the Kpi
 * pair, masses of all four final state particles and returns 4-momenta of
 * all final state particles in B0 rest frame. Due to the rotational
 * symmetry, J/psi momentum is along z-axis and muons are in x-z plane.
 */
void DPHelpers::calculateFinalStateMomenta(double mB0, double m23, double mMuMu, double cosTheta1,
					   double cosTheta2, double phi, 
					   int pion_ID, 
					   double mMuPlus, double mMuMinus,
					   double mPi, double mK, TLorentzVector& pMuPlus,
					   TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK)
{

  

 
// // void CalculateAngles::calculateFinalStateMomenta(double mB0, double m23, double mMuMu, 
// // 						 double cosTheta1, double cosTheta2, double phi, 
// // 						 double mMuPlus, double mMuMinus,
// // 						 double mPi, double mK, 
// // 						 int pion_ID,
// // 						 TLorentzVector& pMuPlus,
// // 						 TLorentzVector& pMuMinus, 
// // 						 TLorentzVector& pPi, 
// // 						 TLorentzVector& pK)
//    if (pion_ID > 0){ // anti-B0
//     if(phi >0) phi = TMath::Pi() - phi;
//     else phi  = -TMath::Pi() - phi;
//     cosTheta1 = -cosTheta1;
//   }


//   //std::cout << " " << pMuPlus.X() << " " << pMuPlus.Y() << " " << pMuPlus.Z() << std::endl;

//   // 4-momenta of dimuon and Kpi system in B0 rest frame aligned along
//   // z-axis.
//   double pJpsi=DPHelpers::daughterMomentum(mB0, mMuMu, m23);
//   TLorentzVector p4Jpsi(0,0,+pJpsi,TMath::Sqrt(mMuMu*mMuMu+pJpsi*pJpsi));
//   TLorentzVector p4Kpi(0,0,-pJpsi,TMath::Sqrt(m23*m23+pJpsi*pJpsi));

//   // 4-momenta of muons, first in dimuon rest frame which are then boosted
//   // using p4Jpsi to B0 rest frame
//   double pMu=DPHelpers::daughterMomentum(mMuMu, mMuPlus, mMuMinus);
//   pMuPlus.SetPxPyPzE(pMu*TMath::Sqrt(1-cosTheta1*cosTheta1),
// 		     0,
// 		     pMu*cosTheta1,
//                      TMath::Sqrt(mMuPlus*mMuPlus+pMu*pMu));

//   pMuMinus.SetPxPyPzE(-pMu*TMath::Sqrt(1-cosTheta1*cosTheta1),
// 		      0,
// 		      -pMu*cosTheta1,
// 		      TMath::Sqrt(mMuMinus*mMuMinus+pMu*pMu));
//   pMuPlus.Boost(p4Jpsi.BoostVector());
//   pMuMinus.Boost(p4Jpsi.BoostVector());

//   // Now kaon and pion to finish
//   double ppK=DPHelpers::daughterMomentum(m23, mK, mPi);
//   double pz=ppK*cosTheta2;
//   double pT=ppK*TMath::Sqrt(1-cosTheta2*cosTheta2);
//   double py=-pT*TMath::Sin(phi);
//   double px=-pT*TMath::Cos(phi);

//   pK. SetPxPyPzE(  px, -py, -pz, TMath::Sqrt(mK*mK+ppK*ppK));
//   pPi.SetPxPyPzE( -px,  py,  pz, TMath::Sqrt(mPi*mPi+ppK*ppK));
//   pK.Boost(p4Kpi.BoostVector());
//   pPi.Boost(p4Kpi.BoostVector());



  CalculateAngles::calculateFinalStateMomenta_GOLD( mB0, m23, mMuMu, cosTheta1,
					       cosTheta2, phi, 
					       mMuPlus,  mMuMinus,
					       mPi,  mK,  
					       pion_ID, 
					       pMuPlus,
					       pMuMinus,  pPi,  pK);

}

void DPHelpers::calculateZplusAngles(TLorentzVector& pB, TLorentzVector& pMuPlus,
      TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK,
				     double* cosThetaZ, double* cosThetaPsi, double* dphi, int pion_ID)
{

  CalculateAngles::calculateZplusAngles_GOLD(pB, 
					     pMuPlus,
					     pMuMinus, 
					     pPi, 
					     pK,
					     cosThetaZ, 
					     cosThetaPsi, 
					     dphi,
					     pion_ID);
// void CalculateAngles::calculateZplusAngles_GOLD(TLorentzVector& pB, 
// 						TLorentzVector& pMuPlus,
// 						TLorentzVector& pMuMinus, 
// 						TLorentzVector& pPi, 
// 						TLorentzVector& pK,
// 						double* cosThetaZ, 
// 						double* cosThetaPsi, 
// 						double* dphi,
// 						int pion_ID)

  /*
  TLorentzVector p4Jpsi=pMuPlus+pMuMinus;
  TLorentzVector p4Zplus=p4Jpsi+pPi;

  // If we are not in B0 rest frame, boost there
  if ( pB.BoostVector().Mag() > 0. )
  {
    p4Jpsi.Boost(-1.0*pB.BoostVector());
    p4Zplus.Boost(-1.0*pB.BoostVector());
    pMuPlus.Boost(-1.0*pB.BoostVector());
    pMuMinus.Boost(-1.0*pB.BoostVector());
    pPi.Boost(-1.0*pB.BoostVector());
    pK.Boost(-1.0*pB.BoostVector());
  }

  // For reminder, work with local copy of final state momenta
  TLorentzVector p4MuPlus(pMuPlus);
  TLorentzVector p4MuMinus(pMuMinus);
  TLorentzVector p4K(pK);
  TLorentzVector p4Pi(pPi);

  // Go to Z+ rest frame
  pB.Boost(-1.0*p4Zplus.BoostVector());
  p4Jpsi.Boost(-1.0*p4Zplus.BoostVector());
  p4MuMinus.Boost(-1.0*p4Zplus.BoostVector());
  p4MuPlus.Boost(-1.0*p4Zplus.BoostVector());
  p4K.Boost(-1.0*p4Zplus.BoostVector());
  p4Pi.Boost(-1.0*p4Zplus.BoostVector());
  p4Zplus.Boost(-1.0*p4Zplus.BoostVector());

  // Get cos of the angle between pi+ and B0
//   *cosThetaZ=p4Pi.Vect().Dot(pB.Vect())/p4Pi.Vect().Mag()/pB.Vect().Mag();
  *cosThetaZ=p4Pi.Vect().Dot(pK.Vect())/p4Pi.Vect().Mag()/pK.Vect().Mag();

  // Now make another copy of vectors and move to J/psi rest frame
  TLorentzVector p4MuPlusJpsiRest(p4MuPlus);
  TLorentzVector p4ZPlusJpsiRest(p4Zplus);
  p4MuPlusJpsiRest.Boost(-1.0*p4Jpsi.BoostVector());
  p4ZPlusJpsiRest.Boost(-1.0*p4Jpsi.BoostVector());
  *cosThetaPsi=p4MuPlusJpsiRest.Vect().Dot(p4ZPlusJpsiRest.Vect())/p4MuPlusJpsiRest.Vect().Mag()/p4ZPlusJpsiRest.Vect().Mag();


  TVector3 normal1=p4Jpsi.Vect().Cross(pPi.Vect());
  TVector3 normal2=p4MuMinus.Vect().Cross(p4MuPlus.Vect());


  TLorentzVector _muPlus(pMuPlus);
  TLorentzVector _muMinus(pMuMinus);
  TLorentzVector _psi(pMuPlus + pMuMinus);
  TLorentzVector _z(pMuPlus + pMuMinus + pPi);
  TLorentzVector _b(pMuPlus + pMuMinus + pPi + pK);
  
  *dphi = planeAngle2_BIN( _muPlus, _muMinus, _psi, _z, _b);

  if(pion_ID > 0)
    {
      (*cosThetaPsi) = -(*cosThetaPsi);
      if(*dphi >0) {(*dphi) = TMath::Pi() - (*dphi);}
      else{ (*dphi) = -TMath::Pi() - (*dphi);}
    }
  */
}

/*
 * Calculate angle between (Kpi)-momentum in psi rest frame (B0 --> psi K*) and pi momentum
 * in psi rest (B0 --> ZK)
 */
double DPHelpers::referenceAxisCosAngle(TLorentzVector& pB, TLorentzVector& pMuPlus,
      TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK)
{
  TLorentzVector p4Jpsi=pMuPlus+pMuMinus;
  TLorentzVector p4Zplus=p4Jpsi+pPi;

  if ( pB.BoostVector().Mag() > 0. )
  {
    p4Jpsi.Boost(-1.0*pB.BoostVector());
    p4Zplus.Boost(-1.0*pB.BoostVector());
    pMuPlus.Boost(-1.0*pB.BoostVector());
    pMuMinus.Boost(-1.0*pB.BoostVector());
    pPi.Boost(-1.0*pB.BoostVector());
    pK.Boost(-1.0*pB.BoostVector());
  }

  // First reference momentum
  TLorentzVector p4Kpi=pK+pPi;
  p4Kpi.Boost(-1.0*p4Jpsi.BoostVector());

  // Second reference momentum
  TLorentzVector p4Pi(pPi);
  p4Pi.Boost(-1.0*p4Zplus.BoostVector());
  p4Jpsi.Boost(-1.0*p4Zplus.BoostVector());
  p4Pi.Boost(-1.0*p4Jpsi.BoostVector());

  // Cos of the angle
  double cosTheta=p4Kpi.Vect().Dot(p4Pi.Vect())/p4Kpi.Vect().Mag()/p4Pi.Vect().Mag();

  return cosTheta;
}

double DPHelpers::referenceAxisCosAngleMu(TLorentzVector& pB, TLorentzVector& pMuPlus,
      TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK)
{
  TLorentzVector p4Jpsi=pMuPlus+pMuMinus;
  TLorentzVector p4Zplus=p4Jpsi+pPi;

  if ( pB.BoostVector().Mag() > 0. )
  {
    p4Jpsi.Boost(-1.0*pB.BoostVector());
    p4Zplus.Boost(-1.0*pB.BoostVector());
    pMuPlus.Boost(-1.0*pB.BoostVector());
    pMuMinus.Boost(-1.0*pB.BoostVector());
    pPi.Boost(-1.0*pB.BoostVector());
    pK.Boost(-1.0*pB.BoostVector());
  }

  // First reference momentum
  TLorentzVector ref1(pB);
  ref1.Boost(-1.0*p4Jpsi.BoostVector());

  // Second reference momentum
  TLorentzVector p4Pi(pPi);
  p4Pi.Boost(-1.0*p4Zplus.BoostVector());
  p4Jpsi.Boost(-1.0*p4Zplus.BoostVector());
  TLorentzVector ref2=p4Pi+p4Jpsi;

  // Cos of the angle
  double cosTheta=ref1.Vect().Dot(ref2.Vect())/ref1.Vect().Mag()/ref2.Vect().Mag();

  return cosTheta;
}

double DPHelpers::planeAngle2_BIN(const TLorentzVector & _muPlus,
				  const TLorentzVector & _muMinus,
				  const TLorentzVector & _psi,
				  const TLorentzVector & _z,
				  const TLorentzVector & _b)


{
  TLorentzVector * particle1     = new TLorentzVector(_muPlus);
  TLorentzVector * particle2     = new TLorentzVector(_muMinus);
  TLorentzVector * mother1       = new TLorentzVector(_psi);
  TLorentzVector * grandmother1  = new TLorentzVector(_z);
  TLorentzVector * Bcand         = new TLorentzVector(_b);


  // Boost everything to the B frame
  TVector3 boostToB = -Bcand -> BoostVector();

  particle1    -> Boost(boostToB);
  particle2    -> Boost(boostToB);
  mother1      -> Boost(boostToB);
  grandmother1 -> Boost(boostToB);

  // Boost X's children to X frame
  TVector3 boostToX = -grandmother1 -> BoostVector();

  particle1    -> Boost(boostToX);
  particle2    -> Boost(boostToX);
  mother1      -> Boost(boostToX);

  // Boost mu to Jpsi frame
  TVector3 boostToJpsi = -mother1 -> BoostVector();

  particle1    -> Boost(boostToJpsi);
  particle2    -> Boost(boostToJpsi);

  // Define some vectors that we need
  TVector3 vecA   = particle1      -> Vect().Unit();
  TVector3 vecB   = particle2      -> Vect().Unit();
  TVector3 vecM1  = mother1        -> Vect().Unit();
  TVector3 vecGM1 = grandmother1   -> Vect().Unit();

  // Get the normals to the decay planes
  TVector3 eJpsi = ( vecA.Cross( vecM1 ) ).Unit();
  TVector3 eX    = ( vecM1.Cross( vecGM1 ) ).Unit();

  // The direction of mother1 in the X/Z frame
  TVector3 ejz = vecM1.Unit();

  // Angle between X and Jpsi planes
  double cosPhi_X_Jpsi = ( eJpsi.Dot( eX ) );
  double sinPhi_X_Jpsi = ( eX.Cross( eJpsi ) ).Dot( ejz );
  double phi_X_Jpsi    = acos( cosPhi_X_Jpsi ) ;

  delete particle1;
  delete particle2;
  delete mother1;
  delete grandmother1;
  delete Bcand;


  // Resolve ambiguity
  return ( sinPhi_X_Jpsi > 0.0 ? phi_X_Jpsi : -phi_X_Jpsi );
//   return 1.;

}
