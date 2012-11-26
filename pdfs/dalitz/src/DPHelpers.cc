#include "DPHelpers.hh"
#include "TMath.h"
#include "TLorentzVector.h"
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
        double cosTheta2, double phi, double mMuPlus, double mMuMinus, 
        double mPi, double mK, TLorentzVector& pMuPlus,
        TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK)
{
  //std::cout << " " << pMuPlus.X() << " " << pMuPlus.Y() << " " << pMuPlus.Z() << std::endl;

  // 4-momenta of dimuon and Kpi system in B0 rest frame aligned along
  // z-axis.
  double pJpsi=DPHelpers::daughterMomentum(mB0, mMuMu, m23);
  TLorentzVector p4Jpsi(0,0,+pJpsi,TMath::Sqrt(mMuMu*mMuMu+pJpsi*pJpsi));
  TLorentzVector p4Kpi(0,0,-pJpsi,TMath::Sqrt(m23*m23+pJpsi*pJpsi));

  // 4-momenta of muons, first in dimuon rest frame which are then boosted
  // using p4Jpsi to B0 rest frame
  double pMu=DPHelpers::daughterMomentum(mMuMu, mMuPlus, mMuMinus);
  pMuPlus.SetPxPyPzE(pMu*TMath::Sqrt(1-cosTheta1*cosTheta1),0,pMu*cosTheta1,
                     TMath::Sqrt(mMuPlus*mMuPlus+pMu*pMu));
  pMuMinus.SetPxPyPzE(-pMu*TMath::Sqrt(1-cosTheta1*cosTheta1),0,-pMu*cosTheta1,
                     TMath::Sqrt(mMuMinus*mMuMinus+pMu*pMu));
  pMuPlus.Boost(p4Jpsi.BoostVector());
  pMuMinus.Boost(p4Jpsi.BoostVector());

  // Now kaon and pion to finish
  double ppK=DPHelpers::daughterMomentum(m23, mK, mPi);
  double pz=ppK*cosTheta2;
  double pT=ppK*TMath::Sqrt(1-cosTheta2*cosTheta2);
  double py=-pT*TMath::Sin(phi);
  double px=-pT*TMath::Cos(phi);
  pK.SetPxPyPzE(-px,-py,-pz,TMath::Sqrt(mK*mK+ppK*ppK));
  pPi.SetPxPyPzE(px,py,pz,TMath::Sqrt(mPi*mPi+ppK*ppK));
  pK.Boost(p4Kpi.BoostVector());
  pPi.Boost(p4Kpi.BoostVector());
}

void DPHelpers::calculateZplusAngles(TLorentzVector& pB, TLorentzVector& pMuPlus,
      TLorentzVector& pMuMinus, TLorentzVector& pPi, TLorentzVector& pK,
      double* cosThetaZ, double* cosThetaPsi, double* dphi)
{
  TLorentzVector p4Jpsi=pMuPlus+pMuMinus;
  TLorentzVector p4Zplus=p4Jpsi+pPi;
  //std::cout << "Z " << p4Jpsi.X() << " " << p4Jpsi.Y() << " " << p4Jpsi.Z() << std::endl;
  //std::cout << "Z " << p4Zplus.X() << " " << pPi.X() << std::endl;

  // If we are not in B0 rest frame, boost there
  if ( pB.BoostVector().Mag() != 0 )
  {
    p4Jpsi.Boost(-1.0*pB.BoostVector());
    p4Zplus.Boost(-1.0*pB.BoostVector());
    pMuPlus.Boost(-1.0*pB.BoostVector());
    pMuMinus.Boost(-1.0*pB.BoostVector());
    pPi.Boost(-1.0*pB.BoostVector());
    pK.Boost(-1.0*pB.BoostVector());
  }
  //std::cout << " " << p4Jpsi.X() << " " << pMuPlus.X() << " " << pMuMinus.X() << std::endl;
  //std::cout << " " << p4Zplus.X() << " " << pPi.X() << std::endl;

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
  *cosThetaZ=p4Pi.Vect().Dot(pB.Vect())/p4Pi.Vect().Mag()/pB.Vect().Mag();

  // Now make another copy of vectors and move to J/psi rest frame
  TLorentzVector p4MuPlusJpsiRest(p4MuPlus);
  TLorentzVector p4ZPlusJpsiRest(p4Zplus);
  p4MuPlusJpsiRest.Boost(-1.0*p4Jpsi.BoostVector());
  p4ZPlusJpsiRest.Boost(-1.0*p4Jpsi.BoostVector());
  *cosThetaPsi=p4MuPlusJpsiRest.Vect().Dot(p4ZPlusJpsiRest.Vect())/p4MuPlusJpsiRest.Vect().Mag()/p4ZPlusJpsiRest.Vect().Mag();

  //std::cout << " " << *cosThetaZ << " " << *cosThetaPsi << " " << *dphi << std::endl;


  // Finally angle between J/psi-pi plane and Mu+ Mu- plane
  // Get normals to two planes and then calculate angle between them  (0;2pi)
//  TVector3 normal1=p4Jpsi.Vect().Cross(pB.Vect());
//  TVector3 normal2=p4Jpsi.Vect().Cross(p4MuPlus.Vect());
  TVector3 normal1=pJpsi.Vect().Cross(pPi.Vect());
  TVector3 normal2=pMuMinus.Vect().Cross(pMuPlus.Vect());
  *dphi=TMath::ACos(normal1.Dot(normal2)/normal1.Mag()/normal2.Mag());
  // Now we still need to check in which quadrant we should be
  TVector3 dir(normal1.Cross(normal2));
//  float a1=dir.Dot(p4Jpsi.Vect())/dir.Mag()/p4Jpsi.Vect().Mag();
  float a1=dir.Dot(pPi.Vect())/dir.Mag()/pPi.Vect().Mag();
  a1=TMath::ACos(a1);

  if ( a1 > TMath::Pi() -0.0001 )
  {
    (*dphi)=2*TMath::Pi()-(*dphi);
  }
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

  if ( pB.BoostVector().Mag() != 0 )
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

  if ( pB.BoostVector().Mag() != 0 )
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

