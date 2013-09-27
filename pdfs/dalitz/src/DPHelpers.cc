#include "TMath.h"
#include "DPHelpers.hh"
#include "TLorentzVector.h"
#include "CalculateAngles.hh"
#include <iostream>
#include <cmath>

void DPHelpers::calculateFinalStateMomentaBelle(double m_b, double ms, double m_psi,
                            double cos2s, double cospi, double phi2s,
                            double m_mu,
                            double m_pi, double m_k,
                            TLorentzVector& p_mu_b, TLorentzVector& p_mu1_b,
                            TLorentzVector& p_pi_b, TLorentzVector& p_k_b)
{
      TLorentzVector p_b_b(0,0,0,m_b);
      const double p_mu_psi_p = sqrt( 0.25*m_psi*m_psi - m_mu*m_mu );

      const double & m = ms;
      const double q=0.5*sqrt( (m+m_k+m_pi)*(m+m_k-m_pi)*(m+m_pi-m_k)*(m-m_pi-m_k) )/m;
      const double p=0.5*sqrt( (m_b+m_psi+m)*(m_b+m_psi-m)*(m_b+m-m_psi)*(m_b-m-m_psi) )/m_b;

      const double p2 = p*p;
      const double q2 = q*q;

      const double ep = ( sqrt((m_psi*m_psi+p2)*(m*m+p2))+p2 )/m;
      const double pp = sqrt(ep*ep-m_psi*m_psi);

      const double spsipi = m_psi*m_psi + m_pi*m_pi + 2.0*( ep*sqrt(m_pi*m_pi+q2) - pp*q*cospi );

      TLorentzVector p_psi_b(0,0,p,sqrt(m_psi*m_psi+p2));
      TLorentzVector p_ks_b(0,0,-p,sqrt(m*m+p2));
      TLorentzVector p_b_psi= p_b_b;
      TLorentzVector p_b_ks= p_b_b;
      TVector3 boostfromks_b;
      TVector3 boostfrompsi_b;

      TVector3 boosttopsi_b = -(p_psi_b.BoostVector());
      p_b_psi.Boost( boosttopsi_b );
      TVector3 boosttoks_b = -(p_ks_b.BoostVector());
      p_b_ks.Boost( boosttoks_b );
      boostfromks_b = -(p_b_ks.BoostVector());
      boostfrompsi_b = -(p_b_psi.BoostVector());

      double sinpi=sqrt(1.0-cospi*cospi);

      const double & cospsi = cos2s;
      const double sinpsi = sqrt( 1.0 - cospsi*cospsi );
      TLorentzVector p_mu_psi( p_mu_psi_p * sinpsi, 0, p_mu_psi_p * cospsi, m_psi*0.5);
      p_mu_b =p_mu_psi;
      p_mu_b.Boost( boostfrompsi_b );
      p_mu1_b = p_psi_b - p_mu_b;

      TLorentzVector p_pi_ks(
                           -q*sinpi*cos(phi2s),
                           q*sinpi*sin(phi2s),
                           q*cospi,
                           sqrt(m_pi*m_pi+q2)
                           );
     p_pi_b = p_pi_ks;
     p_pi_b.Boost( boostfromks_b );

     TLorentzVector p_z_b(p_psi_b + p_pi_b);
     p_k_b = p_ks_b - p_pi_b;
}


void DPHelpers::Belle(const TLorentzVector & _pMuPlus,
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
           double & phiZPsiPsi)
{
  // inputs can be in any reference frame

  TLorentzVector pB=_pMuPlus+_pMuMinus+_pPi+_pK;

  TLorentzVector pMuPlus(_pMuPlus);
  TLorentzVector pMuMinus(_pMuMinus);
  TLorentzVector pPi(_pPi);
  TLorentzVector pK(_pK);

  // ============== B0 rest frame ========================

  if ( pB.BoostVector().Mag() != 0 )
  {
    pMuPlus.Boost(-pB.BoostVector());
    pMuMinus.Boost(-pB.BoostVector());
    pPi.Boost(-pB.BoostVector());
    pK.Boost(-pB.BoostVector());
  }
  TLorentzVector pKPi=pK+pPi;
  TLorentzVector pPsi=pMuPlus+pMuMinus;
  TLorentzVector pZ = pPsi + pPi;
  m23 = pKPi.M();
  m13 = pZ.M();

  //
  TVector3 p3KPi = pKPi.Vect();
  TVector3 p3K = pK.Vect();
  TVector3 p3Psi = pPsi.Vect();
  TVector3 p3MuPlus = pMuPlus.Vect();

  TVector3 aK = p3K - p3KPi * (p3K.Dot(p3KPi)/p3KPi.Mag2());
  TVector3 aMuPlus = p3MuPlus - p3Psi * (p3MuPlus.Dot(p3Psi)/p3Psi.Mag2());

  // angle between K* and Psi decay planes in B0 rest frame
  phiKPiPsi = atan2( 
			(p3Psi.Cross(aK)).Dot(aMuPlus)/(p3Psi.Mag()*aK.Mag()*aMuPlus.Mag()),
                        aK.Dot(aMuPlus)/(aK.Mag()*aMuPlus.Mag())
			);

  if( std::isnan( phiKPiPsi ) )
  {
    //std::cout << "phi is nan" << std::endl;
    phiKPiPsi = 0.;
  }


  // ============= K* rest frame ============================

  TLorentzVector pK_KPi(pK);
  TLorentzVector pPsi_KPi(pPsi);

  pK_KPi.Boost(-pKPi.BoostVector());
  pPsi_KPi.Boost(-pKPi.BoostVector());

  TVector3 p3K_KPi = pK_KPi.Vect();
  TVector3 p3Psi_KPi = pPsi_KPi.Vect();

  // K* helicity angle
  cosKPi = - p3Psi_KPi.Dot(p3K_KPi)/(p3Psi_KPi.Mag()*p3K_KPi.Mag());

  // ================== Z rest frame ================================

  TLorentzVector pMuPlus_Z(pMuPlus);
  TLorentzVector pPi_Z(pPi);
  TLorentzVector pK_Z(pK);
  TLorentzVector pPsi_Z(pPsi);
  //

  pMuPlus_Z.Boost(-pZ.BoostVector());
  pPi_Z.Boost(-pZ.BoostVector());
  pK_Z.Boost(-pZ.BoostVector());
  pPsi_Z.Boost(-pZ.BoostVector());

  TVector3 p3K_Z = pK_Z.Vect();
  TVector3 p3Psi_Z = pPsi_Z.Vect();

  // Z helicity angle
  cosZ = - p3K_Z.Dot(p3Psi_Z)/(p3K_Z.Mag()*p3Psi_Z.Mag());

  // ================== psi rest frame from Z ========================

  TLorentzVector pMuPlus_Z_Psi(pMuPlus_Z);
  TLorentzVector pPi_Z_Psi(pPi_Z);
  TLorentzVector pK_Z_Psi(pK_Z);

  pMuPlus_Z_Psi.Boost(-pPsi_Z.BoostVector());
  pPi_Z_Psi.Boost(-pPsi_Z.BoostVector());
  pK_Z_Psi.Boost(-pPsi_Z.BoostVector());

  TVector3 p3MuPlus_Z_Psi = pMuPlus_Z_Psi.Vect();
  TVector3 p3Pi_Z_Psi = pPi_Z_Psi.Vect();
  TVector3 p3K_Z_Psi = pK_Z_Psi.Vect();

  cosPsi_Z = - p3Pi_Z_Psi.Dot(p3MuPlus_Z_Psi)/(p3Pi_Z_Psi.Mag()*p3MuPlus_Z_Psi.Mag());


  TVector3 aK_Z_Psi = p3K_Z_Psi - p3Pi_Z_Psi * (p3K_Z_Psi.Dot(p3Pi_Z_Psi)/p3Pi_Z_Psi.Mag2());
  TVector3 aMuPlus_Z_Psi = p3MuPlus_Z_Psi - p3Pi_Z_Psi * (p3MuPlus_Z_Psi.Dot(p3Pi_Z_Psi)/p3Pi_Z_Psi.Mag2());
  phiPsiZ = atan2(
                          -(p3Pi_Z_Psi.Cross(aK_Z_Psi)).Dot(aMuPlus_Z_Psi)/(p3Pi_Z_Psi.Mag()*aK_Z_Psi.Mag()*aMuPlus_Z_Psi.Mag()),
                          aK_Z_Psi.Dot(aMuPlus_Z_Psi)/(aK_Z_Psi.Mag()*aMuPlus_Z_Psi.Mag())
                          );
  if ( std::isnan(phiPsiZ) )
  {
    //std::cout << "phiZ is nan" << std::endl;
    phiPsiZ = 0.;
  }



  // ================ psi rest frame from B (i.e. K* decay chain ) ===================

  TLorentzVector pMuPlus_Psi(pMuPlus);
  TLorentzVector pPi_Psi(pPi);
  TLorentzVector pK_Psi(pK);
  TLorentzVector pKPi_Psi(pKPi);
  //

  pMuPlus_Psi.Boost(-pPsi.BoostVector());
  pPi_Psi.Boost(-pPsi.BoostVector());
  pK_Psi.Boost(-pPsi.BoostVector());
  pKPi_Psi.Boost(-pPsi.BoostVector());

  TVector3 p3MuPlus_Psi = pMuPlus_Psi.Vect();
  TVector3 p3Pi_Psi = pPi_Psi.Vect();
  TVector3 p3K_Psi = pK_Psi.Vect();
  TVector3 p3KPi_Psi = pKPi_Psi.Vect();

  cosPsi = - p3KPi_Psi.Dot(p3MuPlus_Psi)/(p3KPi_Psi.Mag()*p3MuPlus_Psi.Mag());

  TVector3 aKPi_Psi = p3KPi_Psi - p3MuPlus_Psi * (p3KPi_Psi.Dot(p3MuPlus_Psi)/p3MuPlus_Psi.Mag2());
  TVector3 aPi_Psi = p3Pi_Psi - p3MuPlus_Psi * (p3Pi_Psi.Dot(p3MuPlus_Psi)/p3MuPlus_Psi.Mag2());
  phiZPsiPsi =  atan2(
                              (p3MuPlus_Psi.Cross(aPi_Psi)).Dot(aKPi_Psi)/(p3MuPlus_Psi.Mag()*aPi_Psi.Mag()*aKPi_Psi.Mag()),
                              aPi_Psi.Dot(aKPi_Psi)/(aPi_Psi.Mag()*aKPi_Psi.Mag())
                              );
  if ( std::isnan(phiZPsiPsi) )
  {
    //std::cout << "phiZPsiPsi is nan" << std::endl; // this only happens once per integration, due to first point being exactly on the boundary
    phiZPsiPsi = 0.;
  }

}

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
