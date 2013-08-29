#include "CalculateAngles.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TNtuple.h>
#include <iostream>
#include <fstream>

double CalculateAngles::HelCos_Tomasz(const double & ms,const double &mz,const int style)
{
  double m_b(5.2794),m_psi(3.686093),m_k(0.493677),m_pi(0.139570),m_mu(0.1056583715);
  //style = 0, K*, style =1, Z, style = 2, between axis

  if(style == 0){
    double Epi = 0.5*(pow(ms,2)+pow(m_pi,2)-pow(m_k,2))/ms;
    double Epsi = 0.5*(pow(m_b,2)-pow(m_psi,2)-pow(ms,2))/ms;
    double ppi = sqrt(Epi*Epi-m_pi*m_pi);
    double ppsi = sqrt(Epsi*Epsi - m_psi*m_psi);
    return -(pow(m_pi,2)+pow(m_psi,2)+2.*Epi*Epsi-pow(mz,2))/(2.*ppi*ppsi);
  }
  else {
    if(style == 1){
      double Epi = 0.5*(pow(mz,2)+pow(m_pi,2)-pow(m_psi,2))/mz;
      double Ek = 0.5*(pow(m_b,2)-pow(m_k,2)-pow(mz,2))/mz;
      double ppi = sqrt(Epi*Epi-m_pi*m_pi);
      double pk = sqrt(Ek*Ek - m_k*m_k);
      return -(pow(m_pi,2)+pow(m_k,2)+2.*Epi*Ek-pow(ms,2))/(2.*ppi*pk);
    }else{
      if(style == 2){
        double E_b = 0.5*(m_b*m_b+m_psi*m_psi-ms*ms)/m_psi;
        double E_pi = 0.5*(mz*mz-m_pi*m_pi-m_psi*m_psi)/m_psi;
        double p_b = sqrt(E_b*E_b-m_b*m_b);
        double p_pi = sqrt(E_pi*E_pi-m_pi*m_pi);
        return (2.*E_b*E_pi+m_psi*m_psi+m_k*m_k-ms*ms-mz*mz)/(2.*p_b*p_pi);
      } else if (style == 3 ) {
        double E_s = 0.5*(m_b*m_b-m_psi*m_psi-ms*ms)/m_psi;
        double E_pi = 0.5*(mz*mz-m_pi*m_pi-m_psi*m_psi)/m_psi;
        double p_s = sqrt(E_s*E_s-ms*ms);
        double p_pi = sqrt(E_pi*E_pi-m_pi*m_pi);
        return (2.*E_s*E_pi+m_k*m_k-ms*ms-m_pi*m_pi)/(2.*p_s*p_pi);
      } else return 0.0;
    }
  }

}


double CalculateAngles::coshel0_Tomasz(TLorentzVector particle, TLorentzVector parent, TLorentzVector grandparent) {

TVector3 boosttoparent = -(parent.BoostVector());

particle.Boost(boosttoparent);
grandparent.Boost(boosttoparent);

TVector3 particle3 = particle.Vect();
TVector3 grandparent3 = -grandparent.Vect();

Double_t numerator = particle3.Dot(grandparent3);
Double_t denominator = (particle3.Mag())*(grandparent3.Mag());
Double_t temp = numerator/denominator;

return temp;

}




double CalculateAngles::coshel1_Tomasz(TLorentzVector particle, TLorentzVector parent,
					 TLorentzVector grandparent,TLorentzVector grandgrandparent)
{

  TVector3 boosttoB = -(grandgrandparent.BoostVector());

  particle.Boost(boosttoB);
  parent.Boost(boosttoB);
  grandparent.Boost(boosttoB);

  TVector3 boosttoZ = -(grandparent.BoostVector());
  particle.Boost(boosttoZ);
  parent.Boost(boosttoZ);

  TVector3 boosttoP = -(parent.BoostVector());
  particle.Boost(boosttoP);


  TVector3 particle3 = particle.Vect();
  TVector3 parent3 = parent.Vect();

  Double_t numerator = particle3.Dot(parent3);
  Double_t denominator = (particle3.Mag())*(parent3.Mag());
  Double_t temp = numerator/denominator;

  return temp;

}



double CalculateAngles::planeAngle0_Tomasz(TLorentzVector particleFrame,
					     TLorentzVector particleA, TLorentzVector particleB,
					     TLorentzVector particleC, TLorentzVector particleD)
{
  TVector3 boostToFrame = -(particleFrame.BoostVector());

  TLorentzVector boostedAxis = particleC + particleD;
  boostedAxis.Boost(boostToFrame);
  TLorentzVector boostedA =particleA;
  boostedA.Boost(boostToFrame);
  TLorentzVector boostedB =particleB;
  boostedB.Boost(boostToFrame);
  TLorentzVector boostedC =particleC;
  boostedC.Boost(boostToFrame);
  TLorentzVector boostedD =particleD;
  boostedD.Boost(boostToFrame);

  TVector3 vecA = (boostedA.Vect()).Unit();
  TVector3 vecB = (boostedB.Vect()).Unit();
  TVector3 vecC = (boostedC.Vect()).Unit();
  TVector3 vecD = (boostedD.Vect()).Unit();

  TVector3 el = ( vecA.Cross( vecB ) ).Unit() ;
  //TVector3 ek = ( vecC.Cross( vecD ) ).Unit() ;
  TVector3 ek = ( vecD.Cross( vecC ) ).Unit() ;

  TVector3 ez = ( boostedAxis.Vect() ).Unit();

  double cosPhi = ( ek.Dot(el) );
  double sinPhi = ( el.Cross(ek) ).Dot( ez ) ;
  //cout <<"cosPhi  " << cosPhi <<endl;
  double phi    = acos( cosPhi ) ;

  return ( sinPhi > 0.0 ? phi : -phi ) ;
}

Float_t CalculateAngles::planeAngle1_Tomasz(TLorentzVector particleFrame,
					      TLorentzVector particleA, TLorentzVector particleB,
					      TLorentzVector particleC, TLorentzVector particleD)

{
  TVector3 boostToFrame = -(particleFrame.BoostVector());

  TLorentzVector boostedAxis = particleC + particleD;
  boostedAxis.Boost(boostToFrame);
  TLorentzVector boostedA =particleA;
  boostedA.Boost(boostToFrame);
  TLorentzVector boostedB =particleB;
  boostedB.Boost(boostToFrame);
  TLorentzVector boostedC =particleC;
  boostedC.Boost(boostToFrame);
  TLorentzVector boostedD =particleD;
  boostedD.Boost(boostToFrame);

  TVector3 boostToZ = -(boostedA+boostedB+boostedD).BoostVector();
  boostedAxis.Boost(boostToZ);
  boostedA.Boost(boostToZ);
  boostedB.Boost(boostToZ);
  boostedC.Boost(boostToZ);
  boostedD.Boost(boostToZ);

  TVector3 vecA = (boostedA.Vect()).Unit();
  TVector3 vecB = (boostedB.Vect()).Unit();
  TVector3 vecC = (boostedC.Vect()).Unit();
  TVector3 vecD = (boostedD.Vect()).Unit();

  TVector3 el = ( vecA.Cross( vecB ) ).Unit() ;
  //TVector3 ek = ( vecC.Cross( vecD ) ).Unit() ;
  TVector3 ek = ( vecD.Cross( vecC ) ).Unit() ;

  TVector3 ez = ( boostedAxis.Vect() ).Unit();

  double cosPhi = ( ek.Dot(el) );
  double sinPhi = ( el.Cross(ek) ).Dot( ez ) ;
  //cout <<"cosPhi  " << cosPhi <<endl;
  double phi    = acos( cosPhi ) ;

  return ( sinPhi > 0.0 ? phi : -phi ) ;
}


double CalculateAngles::decayAngle
( const TLorentzVector& P ,
  const TLorentzVector& Q ,
  const TLorentzVector& D )

{
  const double pd  = P.Dot  ( D ) ;       // P * D
  const double pq  = P.Dot  ( Q ) ;       // P * Q
  const double qd  = Q.Dot  ( D ) ;       // D * Q
  const double mq2 = Q.M2   (   ) ;       // Q^2
  const double mp2 = P.M2   (   ) ;       // P^2
  const double md2 = D.M2   (   ) ;       // D^2

  const double value =  ( pq * pq - mq2 * mp2 ) * ( qd * qd - mq2 * md2 ) ;
  if ( 0 > value )
    {
      std::cout << "LoKi::Kinematics::decayAngle(P,Q,D):: invalid 4-momenta, return InvalidAngle"  << std::endl;
    }
  //
  return ( pd * mq2 - pq * qd ) /(sqrt( value )) ;
//   return 2;
}




double CalculateAngles::decayAngleChi
( const TLorentzVector& d1 ,
  const TLorentzVector& d2 ,
  const TLorentzVector& h1 ,
  const TLorentzVector& h2 )
{

  TLorentzVector B = h1 + h2 + d1 + d2;


  TVector3 boostToB = B.BoostVector();

  TLorentzVector pKstar_KpiMuMu   = h1 + h2;
  TLorentzVector pKplus_KpiMuMu   = h1;
  TLorentzVector pPiminus_KpiMuMu = h2;
  TLorentzVector pPsi_KpiMuMu	  = d1 + d2;
  TLorentzVector pMplus_KpiMuMu   = d1;
  TLorentzVector pMminus_KpiMuMu  = d2;

  pKstar_KpiMuMu  .Boost( -1.*boostToB );
  pKplus_KpiMuMu  .Boost( -1.*boostToB );
  pPiminus_KpiMuMu.Boost( -1.*boostToB );
  pPsi_KpiMuMu    .Boost( -1.*boostToB );
  pMplus_KpiMuMu  .Boost( -1.*boostToB );
  pMminus_KpiMuMu .Boost( -1.*boostToB );

  TVector3 boostToKstar = pKstar_KpiMuMu.BoostVector();
  TVector3 boostToPsi   = pPsi_KpiMuMu.BoostVector();
  TLorentzVector pKplus_MuMu	(pKplus_KpiMuMu);
  TLorentzVector pPiminus_MuMu	(pPiminus_KpiMuMu);
  TLorentzVector pKplus_Kpi	(pKplus_KpiMuMu);
  TLorentzVector pPiminus_Kpi	(pPiminus_KpiMuMu);

  TLorentzVector pMplus_MuMu	(pMplus_KpiMuMu);
  TLorentzVector pMminus_MuMu	(pMminus_KpiMuMu);
  TLorentzVector pMplus_Kpi	(pMplus_KpiMuMu);
  TLorentzVector pMminus_Kpi	(pMminus_KpiMuMu);


  pKplus_Kpi   .Boost( -1.*boostToKstar );
  pPiminus_Kpi .Boost( -1.*boostToKstar );
  pKplus_MuMu  .Boost( -1.*boostToPsi );
  pPiminus_MuMu.Boost( -1.*boostToPsi );

  pMplus_Kpi  .Boost( -1.*boostToKstar );
  pMminus_Kpi .Boost( -1.*boostToKstar );
  pMplus_MuMu .Boost( -1.*boostToPsi );
  pMminus_MuMu.Boost( -1.*boostToPsi );


  // Now calculate the unit vectors etc
  TVector3 e_z_KpiMuMu =     ( pMplus_KpiMuMu.Vect() + pMminus_KpiMuMu.Vect() ).Unit();
  TVector3 e_z_Kpi     = -1.*( pMplus_Kpi.Vect()     + pMminus_Kpi.Vect()     ).Unit();
  TVector3 e_z_MuMu    = -1.*( pKplus_MuMu.Vect()    + pPiminus_MuMu.Vect()   ).Unit();

  TVector3 n_KPi = ((pKplus_KpiMuMu.Vect()).Cross(pPiminus_KpiMuMu.Vect())).Unit();
  TVector3 n_MuMu = ((pMplus_KpiMuMu.Vect()).Cross(pMminus_KpiMuMu.Vect())).Unit();

  // Calculate polar angles
  double cos_thetaK = ((pKplus_Kpi.Vect()).Unit()).Dot(e_z_Kpi );
  double cos_thetaL = ((pMplus_MuMu.Vect()).Unit()).Dot(e_z_MuMu);

  // Calculate phi
  double cos_phi = ( n_KPi.Dot( n_MuMu ) );
  double sin_phi = ( n_KPi.Cross(n_MuMu)).Dot(e_z_KpiMuMu);
  double _atan2 = atan2(sin_phi, cos_phi);
  double _phi = _atan2;

//   if (pion_ID > 0){ // anti-B0
//     if(_phi >0) _phi = TMath::Pi() - _phi;
//     else _phi = -TMath::Pi() - _phi;
//     // 		  if(_phi >0) _phi =  - _phi;
//     // 		  else _phi = +TMath::Pi() - _phi;
//     // 		  _phi = TMath::Pi() - _phi;
//     cos_thetaL = -cos_thetaL;
//   }


  return _phi;

}


void CalculateAngles::calculateZplusAngles_GOLD(TLorentzVector& pB,
						TLorentzVector& pMuPlus,
						TLorentzVector& pMuMinus,
						TLorentzVector& pPi,
						TLorentzVector& pK,
						double* cosThetaZ,
						double* cosThetaPsi,
						double* dphi,
						int pion_ID)
{
  TLorentzVector _muPlus(pMuPlus);
  TLorentzVector _muMinus(pMuMinus);
  TLorentzVector _KPlus(pK);
  TLorentzVector _piMinus(pPi);
  TLorentzVector _psi(pMuPlus + pMuMinus);
  TLorentzVector _z(pMuPlus + pMuMinus + pPi);
  TLorentzVector _b(pMuPlus + pMuMinus + pPi + pK);

  // If we are not in B0 rest frame, boost there
  if ( _b.BoostVector().Mag() != 0 )
  {
    _psi.    Boost(-1.0*_b.BoostVector());
    _z.      Boost(-1.0*_b.BoostVector());
    _muPlus. Boost(-1.0*_b.BoostVector());
    _muMinus.Boost(-1.0*_b.BoostVector());
    _piMinus.Boost(-1.0*_b.BoostVector());
    _KPlus.  Boost(-1.0*_b.BoostVector());
  }

  // Go to Z+ rest frame
  _b.      Boost(-1.0*_z.BoostVector());
  _psi.    Boost(-1.0*_z.BoostVector());
  _muMinus.Boost(-1.0*_z.BoostVector());
  _muPlus. Boost(-1.0*_z.BoostVector());
  _KPlus.  Boost(-1.0*_z.BoostVector());
  _piMinus.Boost(-1.0*_z.BoostVector());
  _z.      Boost(-1.0*_z.BoostVector());

  // Get cos of the angle between pi+ and B0
  *cosThetaZ=_piMinus.Vect().Dot(_KPlus.Vect())/_piMinus.Vect().Mag()/_KPlus.Vect().Mag();

  TLorentzVector p4MuPlusJpsiRest(_muPlus);
  TLorentzVector p4ZPlusJpsiRest(_z);
  p4MuPlusJpsiRest.Boost(-1.0*_psi.BoostVector());
  p4ZPlusJpsiRest.Boost(-1.0*_psi.BoostVector());
  *cosThetaPsi=p4MuPlusJpsiRest.Vect().Dot(-p4ZPlusJpsiRest.Vect())/p4MuPlusJpsiRest.Vect().Mag()/p4ZPlusJpsiRest.Vect().Mag();

  //Bin's version
  *dphi = planeAngle2_BIN_GOLD( _muPlus, _muMinus, _psi, _z, _b);

/*
  if(pion_ID > 0)
    {
//       (*cosThetaZ) = -(*cosThetaZ);
      (*cosThetaPsi) = -(*cosThetaPsi);
      if(*dphi >0) {(*dphi) = TMath::Pi() - (*dphi);}
      else{ (*dphi) = -TMath::Pi() - (*dphi);}
    }
*/
}


/*
 * Function which takes three angular variables, invariant mass of the Kpi
 * pair, masses of all four final state particles and returns 4-momenta of
 * all final state particles in B0 rest frame. Due to the rotational
 * symmetry, J/psi momentum is along z-axis and muons are in x-z plane.
 */
void CalculateAngles::calculateFinalStateMomenta_GOLD(double mB0, double m23, double mMuMu,
						 double cosTheta1, double cosTheta2, double phi,
						 double mMuPlus, double mMuMinus,
						 double mPi, double mK,
						 int pion_ID,
						 TLorentzVector& pMuPlus,
						 TLorentzVector& pMuMinus,
						 TLorentzVector& pPi,
						 TLorentzVector& pK)
{
/*
  if (pion_ID > 0){ // anti-B0
    if(phi >0) phi = TMath::Pi() - phi;
    else phi  = -TMath::Pi() - phi;
    cosTheta1 = -cosTheta1;
  }
*/
  double pJpsi=CalculateAngles::daughterMomentum(mB0, mMuMu, m23);
  TLorentzVector p4Jpsi(0, 0, +pJpsi, TMath::Sqrt(mMuMu*mMuMu+pJpsi*pJpsi));
  TLorentzVector p4Kpi (0, 0, -pJpsi, TMath::Sqrt(m23*m23+pJpsi*pJpsi));
  TLorentzVector pB(0,0,0,mB0);
  TVector3 BBoost = -1.*pB.BoostVector();

  // 4-momenta of muons, first in dimuon rest frame which are then boosted
  // using p4Jpsi to B0 rest frame
  double pMu=CalculateAngles::daughterMomentum(mMuMu, mMuPlus, mMuMinus);

  pMuPlus.SetPxPyPzE(pMu*TMath::Sqrt(1-cosTheta1*cosTheta1),
		     0,
		     pMu*cosTheta1,
                     TMath::Sqrt(mMuPlus*mMuPlus+pMu*pMu));

  pMuMinus.SetPxPyPzE(-pMu*TMath::Sqrt(1-cosTheta1*cosTheta1),
		      0,
		      -pMu*cosTheta1,
		      TMath::Sqrt(mMuMinus*mMuMinus+pMu*pMu));
/*

  pMuPlus.SetPxPyPzE(-pMu*TMath::Sqrt(1-cosTheta1*cosTheta1),
		     0,
		     -pMu*cosTheta1,
                     TMath::Sqrt(mMuPlus*mMuPlus+pMu*pMu));

  pMuMinus.SetPxPyPzE(pMu*TMath::Sqrt(1-cosTheta1*cosTheta1),
		      0,
		      pMu*cosTheta1,
		      TMath::Sqrt(mMuMinus*mMuMinus+pMu*pMu));
*/
  pMuPlus.Boost(+1.*p4Jpsi.BoostVector());
  pMuMinus.Boost(+1.*p4Jpsi.BoostVector());

  // Now kaon and pion to finish
  double ppK=CalculateAngles::daughterMomentum(m23, mK, mPi);
  double pz=ppK*cosTheta2;
  double pT=ppK*TMath::Sqrt(1-cosTheta2*cosTheta2);
  //double py=-pT*TMath::Sin(phi);
  //double px=-pT*TMath::Cos(phi);
  double py=pT*TMath::Sin(phi);
  double px=pT*TMath::Cos(phi);

  //pK.SetPxPyPzE(   px, -py, -pz, TMath::Sqrt(mK*mK+ppK*ppK));
  //pPi.SetPxPyPzE( -px,  py,  pz, TMath::Sqrt(mPi*mPi+ppK*ppK));
  pK.SetPxPyPzE(  -px, -py, -pz, TMath::Sqrt(mK*mK+ppK*ppK));
  pPi.SetPxPyPzE( px, py, pz, TMath::Sqrt(mPi*mPi+ppK*ppK));
  pK.Boost(+1.*p4Kpi.BoostVector());
  pPi.Boost(+1.*p4Kpi.BoostVector());


}


double CalculateAngles::daughterMomentum(double m, double m1, double m2)

{
  double momentum;

  momentum=(m*m-(m1+m2)*(m1+m2))*(m*m-(m1-m2)*(m1-m2));
  momentum=TMath::Sqrt(momentum);
  momentum/=(2*m);

  return momentum;
}


double CalculateAngles::planeAngle2_BIN_GOLD(const TLorentzVector & _muPlus,
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

}

void CalculateAngles::Calculates_Fit_Var_EPFL(const TLorentzVector & pMplus,
					      const TLorentzVector & pMminus,
					      const TLorentzVector & pPiminus,
					      const TLorentzVector & pKplus,
					      int pion_ID,
					      double& m23,
					      double& cosTheta1,
					      double& cosTheta2,
					      double& phi)

{



  TLorentzVector pKstar(pKplus.Px()+pPiminus.Px(), pKplus.Py()+pPiminus.Py(), pKplus.Pz()+pPiminus.Pz(), pKplus.E()+pPiminus.E());
  TLorentzVector pPsi  (pMplus.Px()+pMminus.Px(),pMplus.Py()+pMminus.Py(),pMplus.Pz()+pMminus.Pz(),pMplus.E()+pMminus.E());
  TLorentzVector pZ    (pPsi.Px()+pPiminus.Px(), pPsi.Py()+pPiminus.Py(), pPsi.Pz()+pPiminus.Pz(), pPsi.E()+pPiminus.E());
  TLorentzVector B     (pKstar.Px()+pPsi.Px(), pKstar.Py()+pPsi.Py(), pKstar.Pz()+pPsi.Pz(), pKstar.E()+pPsi.E());

  // First boost everything to the B frame
  TVector3 boostToB = B.BoostVector();
  TLorentzVector pKstar_KpiMuMu   (pKstar);
  TLorentzVector pKplus_KpiMuMu   (pKplus);
  TLorentzVector pPiminus_KpiMuMu (pPiminus);
  TLorentzVector pPsi_KpiMuMu	  (pPsi);
  TLorentzVector pMplus_KpiMuMu   (pMplus);
  TLorentzVector pMminus_KpiMuMu  (pMminus);

  pKstar_KpiMuMu  .Boost( -1.*boostToB );
  pKplus_KpiMuMu  .Boost( -1.*boostToB );
  pPiminus_KpiMuMu.Boost( -1.*boostToB );
  pPsi_KpiMuMu    .Boost( -1.*boostToB );
  pMplus_KpiMuMu  .Boost( -1.*boostToB );
  pMminus_KpiMuMu .Boost( -1.*boostToB );

  // Now boost the final state particles to their corresponding mother frames
  TVector3 boostToKstar = pKstar_KpiMuMu.BoostVector();
  TVector3 boostToPsi   = pPsi_KpiMuMu.BoostVector();
  TLorentzVector pKplus_MuMu	(pKplus_KpiMuMu);
  TLorentzVector pPiminus_MuMu	(pPiminus_KpiMuMu);
  TLorentzVector pKplus_Kpi	(pKplus_KpiMuMu);
  TLorentzVector pPiminus_Kpi	(pPiminus_KpiMuMu);

  TLorentzVector pMplus_MuMu	(pMplus_KpiMuMu);
  TLorentzVector pMminus_MuMu	(pMminus_KpiMuMu);
  TLorentzVector pMplus_Kpi	(pMplus_KpiMuMu);
  TLorentzVector pMminus_Kpi	(pMminus_KpiMuMu);

  pKplus_Kpi   .Boost( -1.*boostToKstar );
  pPiminus_Kpi .Boost( -1.*boostToKstar );
  pKplus_MuMu  .Boost( -1.*boostToPsi );
  pPiminus_MuMu.Boost( -1.*boostToPsi );

  pMplus_Kpi  .Boost( -1.*boostToKstar );
  pMminus_Kpi .Boost( -1.*boostToKstar );
  pMplus_MuMu .Boost( -1.*boostToPsi );
  pMminus_MuMu.Boost( -1.*boostToPsi );


  // Now calculate the unit vectors etc
  TVector3 e_z_KpiMuMu =     ( pMplus_KpiMuMu.Vect() + pMminus_KpiMuMu.Vect() ).Unit();
  TVector3 e_z_Kpi     = -1.*( pMplus_Kpi.Vect()     + pMminus_Kpi.Vect()     ).Unit();
  TVector3 e_z_MuMu    = -1.*( pKplus_MuMu.Vect()    + pPiminus_MuMu.Vect()   ).Unit();

  TVector3 n_KPi = ((pKplus_KpiMuMu.Vect()).Cross(pPiminus_KpiMuMu.Vect())).Unit();
  TVector3 n_MuMu = ((pMplus_KpiMuMu.Vect()).Cross(pMminus_KpiMuMu.Vect())).Unit();

  // Calculate polar angles
  double cos_thetaK = ((pKplus_Kpi.Vect()).Unit()).Dot(e_z_Kpi );
  double cos_thetaL = ((pMplus_MuMu.Vect()).Unit()).Dot(e_z_MuMu);

  // Calculate phi
  double cos_phi = ( n_KPi.Dot( n_MuMu ) );
  double sin_phi = ( n_KPi.Cross(n_MuMu)).Dot(e_z_KpiMuMu);
  double _atan2 = atan2(sin_phi, cos_phi);
  // 		double _phi = _atan2 > 0 ? _atan2 : _atan2 + 2.*TMath::Pi();
  // 		double _phi = _atan2 > 0 ? _atan2 : _atan2 + TMath::Pi();
  double _phi = _atan2;


  if (pion_ID > 0){ // anti-B0
    if(_phi >0) _phi = TMath::Pi() - _phi;
    else _phi  = -TMath::Pi() - _phi;
    cos_thetaL = -cos_thetaL;
//     cos_thetaK = -cos_thetaK;
  }


  m23 = pKstar.M();
  cosTheta1 = cos_thetaL;
  cosTheta2 = cos_thetaK;
  phi = _phi;

}



// void CalculateAngles::Loop()
// {
// 	if (fChain == 0) return;

// 	Long64_t nentries = fChain->GetEntriesFast();

// 	TFile * file = TFile::Open("selected_candidates_with_correct_angles.root","RECREATE");
// // 	TNtuple * tuple = new TNtuple("tuple", "tuple",  "phi:cosTheta2:cosTheta1:pion_ID:phi_Tomasz:cosTheta2_Tomasz:cosTheta1_Tomasz");
// 	TNtuple * tuple = new TNtuple("tuple", "tuple",  "m23:B0_mass:cosThetaPsi_michal:cosThPsi_Z:cosThetaZ_michal:cosThZ:dphi_michal:phiZPsi_bin");
// // 		tuple->Fill(cosThetaPsi_michal, cosThPsi_Z, cosThetaZ_michal, cosThZ, dphi_michal,  phiZPsi_bin );
// // 	tuple->Fill( _phi, cos_thetaK, cos_thetaL, pion_ID,
// // 			    , phiKsPsi, cosThKs, cosThPsi);

// 	ofstream myfile;
// 	myfile.open ("output.txt");

// 	Long64_t nbytes = 0, nb = 0;
// 	for (Long64_t jentry=0; jentry<nentries;jentry++)
// 	{
// 		Long64_t ientry = LoadTree(jentry);
// 		if (ientry < 0) break;
// 		nb = fChain->GetEntry(jentry);   nbytes += nb;



// 		if ( !(
// // 		       MuPlus_pt > 700 && MuMinus_pt > 700 &&
// // 		       sqrt(MuPlus_px * MuPlus_px + MuPlus_py * MuPlus_py + MuPlus_pz * MuPlus_pz) > 10000 &&
// // 		       sqrt(MuMinus_px * MuMinus_px + MuMinus_py * MuMinus_py + MuMinus_pz * MuMinus_pz) > 10000 &&
// 		       (mKK > 1030 || mKK < 1010) &&
// 		        B0_mass > 5200 && 5350 > B0_mass &&                            //NOMINAL SEL
//                         B0_BPVLTIME > 0.25e-12 &&                                      //NOMINAL SEL
//                         Psi_mass > 3625  && Psi_mass < 3750 &&                         //NOMINAL SEL
//                         Pi_pid < 5 && K_pid > 0  &&  K_pid-Pi_pid > -10 &&             //NOMINAL SEL
//                         5 > B0_BPVIPCHI2 &&                                            //NOMINAL SEL
//                         Pi_pt >200 && K_pt >200 && B0_pt > 2000 &&                     //NOMINAL SEL
//                         vchi2dof_DTF < 4 &&                                            //NOMINAL SEL
// 		        Psi_BPVIPCHI2 > 5 &&                                           //NOMINAL SEL
//                         mincos_jpsi_h > 0.98)                                          //NOMINAL SEL
//                         ) continue;                                                   //NOMINAL SEL


// 		TLorentzVector B	(B0_px, B0_py, B0_pz, B0_E);
// 		TLorentzVector pKstar 	(KPi_px, KPi_py, KPi_pz, KPi_E);
// 		TLorentzVector pKplus 	(K_px, K_py, K_pz, K_E);
// 		TLorentzVector pPiminus	(Pi_px, Pi_py, Pi_pz, Pi_E);
// 		TLorentzVector pPsi   	(Psi_px, Psi_py, Psi_pz, Psi_E);
// 		TLorentzVector pZ   	(Psi_px+Pi_px, Psi_py+Pi_py, Psi_pz+Pi_pz, Psi_E+Pi_E);
// 		TLorentzVector pMplus 	(MuPlus_px, MuPlus_py, MuPlus_pz, MuPlus_E);
// 		TLorentzVector pMminus 	(MuMinus_px, MuMinus_py, MuMinus_pz, MuMinus_E);


// 		double m23_EPFL;
// 		double cosTheta1_EPFL;
// 		double cosTheta2_EPFL;
// 		double phi_EPFL;

// 		Calculates_Fit_Var_EPFL( pMplus,
// 					 pMminus,
// 					 pPiminus,
// 					 pKplus,
// 					 m23_EPFL,
// 					 cosTheta1_EPFL,
// 					 cosTheta2_EPFL,
// 					 phi_EPFL);

// // 		cout << "m23_EPFL      : " << m23_EPFL<< endl;
// // 		cout << "cosTheta1_EPFL: " << cosTheta1_EPFL<< endl;
// // 		cout << "cosTheta2_EPFL: " << cosTheta2_EPFL<< endl;
// // 		cout << "phi_EPFL      : " << phi_EPFL<< endl;
// // 		cout << "" << endl;

// 		double phiZPsi_bin_EPFL = planeAngle2_BIN(pMplus, pMminus, pPsi, pZ, B);


// 		double cosThetaZ_michal;
// 		double cosThetaPsi_michal;
// 		double dphi_michal;

// 		calculateZplusAngles(B,
// 				     pMplus,
// 				     pMminus,
// 				     pPiminus,
// 				     pKplus,
// 				     &cosThetaZ_michal,
// 				     &cosThetaPsi_michal,
// 				     &dphi_michal,
// 				     pion_ID);

// 		double cos_Z_joel = decayAngle( pPiminus, pMplus+pMminus+pPiminus,  B);

// 		// 		Checks suggested by Tomasz: compute Z angles from Ks angles as done in the fitter
// 		// First get 4vet from angles with function in the fit




// 		TLorentzVector  pMuPlus_fromAngles,  pMuMinus_fromAngles,  pPi_fromAngles,  pK_fromAngles;



// 		calculateFinalStateMomenta(5279.58, // double mB0
// 					   m23_EPFL,
// 					   3686,
// 					   cosTheta1_EPFL,
// 					   cosTheta2_EPFL,
// 					   phi_EPFL,
// 					   105.658,
// 					   105.658,
// 					   139.570,
// 					   493.677,
// 					   pion_ID,
// 					   pMuPlus_fromAngles,
// 					   pMuMinus_fromAngles,
// 					   pPi_fromAngles,
// 					   pK_fromAngles);


// // 		TLorentzVector B_fromAngles	= pMuPlus_fromAngles + pMuMinus_fromAngles + pPi_fromAngles + pK_fromAngles;
// // 		TLorentzVector pKstar_fromAngles= pPi_fromAngles + pK_fromAngles;
// // 		TLorentzVector pPsi_fromAngles  = pMuPlus_fromAngles + pMuMinus_fromAngles;
// // 		TLorentzVector pZ_fromAngles  	= pMuPlus_fromAngles + pMuMinus_fromAngles + pPi_fromAngles;


// 		double m23_EPFL_fromAngles, cosTheta1_EPFL_fromAngles, cosTheta2_EPFL_fromAngles, phi_EPFL_fromAngles;

// 		TLorentzVector  pMuPlus_fromAngles2,  pMuMinus_fromAngles2,  pPi_fromAngles2,  pK_fromAngles2;

// 		Calculates_Fit_Var_EPFL( pMuPlus_fromAngles,
// 					 pMuMinus_fromAngles,
// 					 pPi_fromAngles,
// 					 pK_fromAngles,
// 					 m23_EPFL_fromAngles,
// 					 cosTheta1_EPFL_fromAngles,
// 					 cosTheta2_EPFL_fromAngles,
// 					 phi_EPFL_fromAngles);

// 		calculateFinalStateMomenta(5279.58, // double mB0
// 					   m23_EPFL_fromAngles,
// 					   3686,
// 					   cosTheta1_EPFL_fromAngles,
// 					   cosTheta2_EPFL_fromAngles,
// 					   phi_EPFL_fromAngles,
// 					   105.658,
// 					   105.658,
// 					   139.570,
// 					   493.677,
// 					   pion_ID,
// 					   pMuPlus_fromAngles2,
// 					   pMuMinus_fromAngles2,
// 					   pPi_fromAngles2,
// 					   pK_fromAngles2);

// // 		double m_kpi_fromAngles   = pKstar_fromAngles.M();
// // 		double m_psipi_fromAngles = pZ_fromAngles.M();

// // 		cout << "m_kpi_fromAngles  = " <<  m_kpi_fromAngles << endl;
// // 		cout << "m_psipi_fromAngles= " <<  m_psipi_fromAngles << endl;
// // 		cout << " " <<  endl;

// // 		double cosThKs = HelCos_Tomasz(m_kpi,m_psipi,0);

// // 		double cos_K_joel_fromAngles  =  HelCos_Tomasz(m_kpi_fromAngles, m_psipi_fromAngles, 0);
// // 		double cos_Mu_joel_fromAngles =  coshel0_Tomasz( pMuPlus_fromAngles, pPsi_fromAngles, B_fromAngles);
// // 		double phi_joel_fromAngles    =  planeAngle0_Tomasz( B_fromAngles, pMuPlus_fromAngles, pMuMinus_fromAngles,
// // 								     pK_fromAngles, pPi_fromAngles );



// // 		double cosThetaZ_fromAngles = HelCos_Tomasz(m23, m_psipi_fromAngles, 1);

// // 		double cosThetaPsi_Z_fromAngles = coshel1_Tomasz( pMuPlus_fromAngles, pPsi_fromAngles,
// // 								  pZ_fromAngles, B_fromAngles );

// // 		double dphi_fromAngles = planeAngle2_BIN( pMuPlus_fromAngles, pMuMinus_fromAngles,
// // 							  pPsi_fromAngles, pZ_fromAngles, B_fromAngles);



// // 		if(pion_ID > 0){
// // 		  cos_Mu_joel_fromAngles = -cos_Mu_joel_fromAngles;
// // 		  cosThetaPsi_Z_fromAngles = -cosThetaPsi_Z_fromAngles;

// // 		  if(phi_joel_fromAngles >0) phi_joel_fromAngles = TMath::Pi() - phi_joel_fromAngles;
// // 		  else phi_joel_fromAngles = -TMath::Pi() - phi_joel_fromAngles;

// // 		  if(dphi_fromAngles >0) dphi_fromAngles = TMath::Pi() - dphi_fromAngles;
// // 		  else dphi_fromAngles = -TMath::Pi() - dphi_fromAngles;
// // 		}

// // 		TLorentzVector  pMuPlus_fromAngles2,  pMuMinus_fromAngles2,  pPi_fromAngles2,  pK_fromAngles2;


// // 		calculateFinalStateMomenta(5.27958, // double mB0
// // 					   m_kpi_fromAngles,
// // 					   3.686,
// // 					   cos_Mu_joel_fromAngles,
// // 					   cos_K_joel_fromAngles,
// // 					   phi_joel_fromAngles,
// // 					   0.105658,
// // 					   0.105658,
// // 					   0.139570,
// // 					   0.493677,
// // 					   pMuPlus_fromAngles2,
// // 					   pMuMinus_fromAngles2,
// // 					   pPi_fromAngles2,
// // 					   pK_fromAngles2);






// 		//---------------------------------------------------------------------------------------------------------------


// 		//                                       T O M A S Z ' S   C O D E


// 		//---------------------------------------------------------------------------------------------------------------


// 		// ---------------- input variables ----------------------

// 		TLorentzVector B_Tomasz	(B0_px/1000, B0_py/1000, B0_pz/1000, B0_E/1000);
// 		TLorentzVector pKstar_Tomasz 	(KPi_px/1000, KPi_py/1000, KPi_pz/1000, KPi_E/1000);
// 		TLorentzVector pKplus_Tomasz 	(K_px/1000, K_py/1000, K_pz/1000, K_E/1000);
// 		TLorentzVector pPiminus_Tomasz	(Pi_px/1000, Pi_py/1000, Pi_pz/1000, Pi_E/1000);
// 		TLorentzVector pPsi_Tomasz   	(Psi_px/1000, Psi_py/1000, Psi_pz/1000, Psi_E/1000);
// 		TLorentzVector pZ_Tomasz   	((Psi_px+Pi_px)/1000, (Psi_py+Pi_py)/1000, (Psi_pz+Pi_pz)/1000, (Psi_E+Pi_E)/1000);
// 		TLorentzVector pMplus_Tomasz 	(MuPlus_px/1000, MuPlus_py/1000, MuPlus_pz/1000, MuPlus_E/1000);
// 		TLorentzVector pMminus_Tomasz 	(MuMinus_px/1000, MuMinus_py/1000, MuMinus_pz/1000, MuMinus_E/1000);

// 		double phiZPsi_bin = planeAngle2_BIN(pMplus_Tomasz, pMminus_Tomasz, pPsi_Tomasz, pZ_Tomasz, B_Tomasz);

// 		double m_kpi = pKstar.M()/1000;
// 		double m_psipi = pZ.M()/1000;



// // 		cout << "Tomasz m_kpi  = " <<  m_kpi << endl;
// // 		cout << "Tomasz m_psipi= " <<  m_psipi << endl;
// // 		cout << " "  << endl;
// // 		cout << " "  << endl;


// 		//  K* decay chain
// 		double cosThKs = HelCos_Tomasz(m_kpi,m_psipi,0);

// 		double deg2rad = TMath::Pi()/180;

// 		double cosThPsi = coshel0_Tomasz( pMplus_Tomasz, pPsi_Tomasz, B_Tomasz);


// // 		double phiKsPsi = planeAngle0_Tomasz( B, pMplus, pMminus, pKplus, pPiminus )*deg2rad;
// 		double phiKsPsi = planeAngle0_Tomasz( B_Tomasz, pMplus_Tomasz, pMminus_Tomasz, pKplus_Tomasz, pPiminus_Tomasz );

// 		//  Z decay chain
// 		double cosThZ = HelCos_Tomasz(m_kpi,m_psipi,1);

// 		double cosThPsi_Z = coshel1_Tomasz( pMplus_Tomasz, pPsi_Tomasz, pZ_Tomasz, B_Tomasz );
// // 		double phiZPsi = planeAngle1_Tomasz( B, pMplus, pMminus, pKplus, pPiminus )*deg2rad;
// 		double phiZPsi = planeAngle1_Tomasz( B_Tomasz, pMplus_Tomasz, pMminus_Tomasz, pKplus_Tomasz, pPiminus_Tomasz );


// // // // 		if(B_ID == -511){
// 		if(pion_ID > 0){
// 		  cosThPsi = -cosThPsi;
// 		  cosThPsi_Z = -cosThPsi_Z;
// 		  if(phiKsPsi >0) phiKsPsi = TMath::Pi() - phiKsPsi;
// 		  else phiKsPsi = -TMath::Pi() - phiKsPsi;
// 		  if(phiZPsi >0) phiZPsi = TMath::Pi() - phiZPsi;
// 		  else phiZPsi = -TMath::Pi() - phiZPsi;
// 		}





// 		//Bin's Z angle



// // 		double CalculateAngles::planeAngle2(const TLorentzVector & _muPlus,
// // 				      const TLorentzVector & _muMinus,
// // 				      const TLorentzVector & _psi,
// // 				      const TLorentzVector & _z,
// // 				      const TLorentzVector & _b) const


// 		//---------------------------------------------------------------------------------------------------------------


// 		//                                    end of    T O M A S Z ' S   C O D E


// 		//---------------------------------------------------------------------------------------------------------------



// // 		tuple->Fill( _phi, cos_thetaK, cos_thetaL, pion_ID, 	     phiKsPsi, cosThKs, cosThPsi);
// 		tuple->Fill(m23, B0_mass, cosThetaPsi_michal, cosThPsi_Z, cosThetaZ_michal, cosThZ, dphi_michal,  phiZPsi_bin );






// // 		pMplus .Boost( -1.*boostToB );
// // 		pMminus.Boost( -1.*boostToB );
// // 		pKplus .Boost( -1.*boostToB );
// // 		pMminus.Boost( -1.*boostToB );

// // 		pPi_fromAngles .Boost( -1.*boostToB );
// // 		pK_fromAngles.Boost( -1.*boostToB );
// // 		pMuPlus_fromAngles .Boost( -1.*boostToB );
// // 		pMuMinus_fromAngles.Boost( -1.*boostToB );


// 		myfile << "*********EVENT************" << endl;
// 		myfile << " " << endl;
// 		myfile << "pion ID: " << pion_ID << endl;
// // // 		myfile << "(B0_px, B0_py, B0_pz, B0_E) = (" << B0_px  << ", " << B0_py << ", " << B0_pz << ", " << B0_E <<  ") " << endl;
// // 		myfile << "(Pi_px, Pi_py, Pi_pz, Pi_E)= (" <<  Pi_px  << ", " << Pi_py << ", " << Pi_pz << ", " << Pi_E <<  ") "<< endl;
// // 		myfile << "(K_px, K_py, K_pz, K_E)= (" <<  K_px  << ", " << K_py << ", " << K_pz << ", " << K_E <<  ") "<< endl;
// // 		myfile << "(Psi_px, Psi_py, Psi_pz, Psi_E)= (" << Psi_px  << ", " << Psi_py << ", " << Psi_pz << ", " << Psi_E <<  ") " << endl;
// // 		myfile << "(MuPlus_px, MuPlus_py, MuPlus_pz, MuPlus_E)= (" << MuPlus_px  << ", " << MuPlus_py << ", " << MuPlus_pz << ", " << MuPlus_E <<  ") " << endl;
// // 		myfile << "(MuMinus_px, MuMinus_py, MuMinus_pz, MuMinus_E)= (" << MuMinus_px  << ", " << MuMinus_py << ", " << MuMinus_pz << ", " << MuMinus_E <<  ") " << endl;


// // 		myfile << "            (Pi_px, Pi_py, Pi_pz, Pi_E)= (" <<  Pi_px  << ", " << Pi_py << ", " << Pi_pz << ", " << Pi_E <<  ") "<< endl;
// 		myfile << "FROM ANGLES :(Pi_px, Pi_py, Pi_pz, Pi_E)= (" <<  pPi_fromAngles.Px()  << ", " << pPi_fromAngles.Py() << ", " << pPi_fromAngles.Pz() << ", " << pPi_fromAngles.E() <<  ") "<< endl;
// 		myfile << "FROM ANGLES2:(Pi_px, Pi_py, Pi_pz, Pi_E)= (" <<  pPi_fromAngles2.Px()  << ", " << pPi_fromAngles2.Py() << ", " << pPi_fromAngles2.Pz() << ", " << pPi_fromAngles2.E() <<  ") "<< endl;
// 		myfile << " " << endl;

// // 		myfile << "            (K_px, K_py, K_pz, K_E)= (" <<  K_px  << ", " << K_py << ", " << K_pz << ", " << K_E <<  ") "<< endl;
// 		myfile << "FROM ANGLES :(K_px, K_py, K_pz, K_E)= (" <<  pK_fromAngles.Px()  << ", " << pK_fromAngles.Py() << ", " << pK_fromAngles.Pz() << ", " << pK_fromAngles.E() <<  ") "<< endl;
// 		myfile << "FROM ANGLES2:(K_px, K_py, K_pz, K_E)= (" <<  pK_fromAngles2.Px()  << ", " << pK_fromAngles2.Py() << ", " << pK_fromAngles2.Pz() << ", " << pK_fromAngles2.E() <<  ") "<< endl;
// 		myfile << " " << endl;

// // // 		myfile << "            (Psi_px, Psi_py, Psi_pz, Psi_E)= (" << Psi_px  << ", " << Psi_py << ", " << Psi_pz << ", " << Psi_E <<  ") " << endl;
// // 		myfile << "FROM ANGLES:(Psi_px, Psi_py, Psi_pz, Psi_E)= (" << Psi_px  << ", " << Psi_py << ", " << Psi_pz << ", " << Psi_E <<  ") " << endl;
// // // 		myfile << " " << endl;

// // 		myfile << "            (MuPlus_px, MuPlus_py, MuPlus_pz, MuPlus_E)= (" << MuPlus_px  << ", " << MuPlus_py << ", " << MuPlus_pz << ", " << MuPlus_E <<  ") " << endl;
// 		myfile << "FROM ANGLES :(MuPlus_px, MuPlus_py, MuPlus_pz, MuPlus_E)= (" << pMuPlus_fromAngles.Px()  << ", " << pMuPlus_fromAngles.Py() << ", " << pMuPlus_fromAngles.Pz() << ", " << pMuPlus_fromAngles.E() <<  ") " << endl;
// 		myfile << "FROM ANGLES2:(MuPlus_px, MuPlus_py, MuPlus_pz, MuPlus_E)= (" << pMuPlus_fromAngles2.Px()  << ", " << pMuPlus_fromAngles2.Py() << ", " << pMuPlus_fromAngles2.Pz() << ", " << pMuPlus_fromAngles2.E() <<  ") " << endl;
// 		myfile << " " << endl;

// // 		myfile << "            (MuMinus_px, MuMinus_py, MuMinus_pz, MuMinus_E)= (" << MuMinus_px  << ", " << MuMinus_py << ", " << MuMinus_pz << ", " << MuMinus_E <<  ") " << endl;
// 		myfile << "FROM ANGLES :(MuMinus_px, MuMinus_py, MuMinus_pz, MuMinus_E)= (" << pMuMinus_fromAngles.Px()  << ", " << pMuMinus_fromAngles.Py() << ", " << pMuMinus_fromAngles.Pz()<< ", " << pMuMinus_fromAngles.E() <<  ") " << endl;
// 		myfile << "FROM ANGLES2:(MuMinus_px, MuMinus_py, MuMinus_pz, MuMinus_E)= (" << pMuMinus_fromAngles2.Px()  << ", " << pMuMinus_fromAngles2.Py() << ", " << pMuMinus_fromAngles2.Pz()<< ", " << pMuMinus_fromAngles2.E() <<  ") " << endl;


// //              double cos_K_joel_fromAngles =decayAngle( pK_fromAngles, pKstar_fromAngles , B_fromAngles );

// // 		double phi_joel_fromAngles = decayAngleChi(  pMuPlus_fromAngles, pMuMinus_fromAngles, pK_fromAngles, pPi_fromAngles);

// // 		m23_EPFL,
// // 		  cosTheta1_EPFL,
// // 		  cosTheta2_EPFL,
// // 		  phi_EPFL);
// 		myfile << " " << endl;
// 		myfile << " " << endl;
// // 		myfile << "cos_theta_EPFL            = " << cos_thetaK << endl;
// 		myfile << "m23                 = " << m23 << endl;
// 		myfile << "m23_EPFL            = " << m23_EPFL << endl;
// 		myfile << "m23_CHECK_fromAngles= " << m23_EPFL_fromAngles << endl;


// 		myfile << " " << endl;
// 		myfile << " " << endl;
// // 		myfile << "cos_thetaK_EPFL            = " << cos_thetaK << endl;
// 		myfile << "cos_thetaK_EPFL            = " << cosTheta2_EPFL << endl;
// 		myfile << "cos_thetaK_CHECK_fromAngles= " << cosTheta2_EPFL_fromAngles << endl;
// 		myfile << "cos_thetaK_SYRA            = " << cosThKs <<  endl;
// 		myfile << " " << endl;
// // 		myfile << "cos_thetaPsi_EPFL            = " << cos_thetaL  << endl;
// 		myfile << "cos_thetaPsi_EPFL            = " << cosTheta1_EPFL  << endl;
// // 		myfile << "cos_thetaPsi_EPFL_fromAngles = " << cos_thetaL_fromAngles  << endl;
// 		myfile << "cos_thetaPsi_CHECK_fromAngles= " << cosTheta1_EPFL_fromAngles  << endl;
// 		myfile << "cos_thetaPsi_SYRA            = " << cosThPsi <<  endl;
// 		myfile << " " << endl;
// // 		myfile << "phi_EPFL            = " << _phi << endl;
// 		myfile << "phi_EPFL            = " << phi_EPFL << endl;
// // 		myfile << "phi_EPFL_fromAngles = " << _phi_fromAngles << endl;
// 		myfile << "phi_CHECK_fromAngles= " << phi_EPFL_fromAngles << endl;
// 		myfile << "phi_SYRA            = " << phiKsPsi << endl;
// // 		myfile << "DIFF phi            = " << phiKsPsi - _phi <<endl;
// 		myfile << " " << endl;
// 		myfile << "cosThetaPsi_Z_EPFL            = " << cosThetaPsi_michal  << endl;
// 		myfile << "cosThetaPsi_Z_SYRA            = " << cosThPsi_Z <<  endl;

// 		myfile << " " << endl;
// 		myfile << "cosTheta_Z_EPFL              = " << cosThetaZ_michal  << endl;
// 		myfile << "cosTheta_Z_JOEL              = " << cos_Z_joel  << endl;
// // 		myfile << "cosTheta_Z_EPFL_fromAngles,  = " << cosThetaZ_fromAngles   << endl;
// 		myfile << "cosTheta_Z_SYRA              = " << cosThZ <<  endl;
// 		myfile << "diff cosTheta_Z              = " << acos(cosThZ) - acos(cosThetaZ_michal) <<  endl;

// 		myfile << " " << endl;
// 		myfile << "dphi_Z_EPFL           = " << dphi_michal  << endl;
// // 		myfile << "dphi_Z_EPFLfromAngles = " << dphi_fromAngles  << endl;
// // 		myfile << "dphi_Z_SYRA           = " << phiZPsi << endl;
// 		myfile << "dphi_Z_BIN            = " << phiZPsi_bin << endl;
// 		myfile << "diff phiZ angles      = " << phiZPsi_bin - dphi_michal<< endl;
// // 		myfile << "dphi_Z_BIN_EPFL       = " << phiZPsi_bin_EPFL << endl;
// // 		myfile << "dphi_Z_BIN_fromAngles = " << phiZPsi_bin_fromAngles << endl;
// 		myfile << " " << endl;
// 		myfile << " " << endl;
// 		myfile << " " << endl;
// 		myfile << " " << endl;
// 		myfile << " " << endl;


// // 		myfile << "m kpi = " << m23  << endl;
// // 		myfile << "m psi pi = " << m13  << endl;
// // // 		myfile << " " << endl;
// // // 		myfile << " " << endl;
// // 		myfile << "pPsi_fromAngles = " << pPsi_fromAngles.M() << endl;
// // 		myfile << "pKPi_fromAngles = " << pKPi_fromAngles.M() << endl;
// // 		myfile << "B_fromAngles = " << B_fromAngles.M() << endl;
// // 		myfile << "cosTheta_Z = " << cosThetaZ_michal  << endl;
// // 		myfile << "dphi_Z = " << dphi_michal  << endl;
// // 		myfile << "cosThetaPsi_Z = " << cosThetaPsi_michal  << endl;
// // 		myfile << "cosThetaZ_fromAngles,  = " << cosThetaZ_fromAngles   << endl;
// // 		myfile << "cosThetaPsi_fromAngles = " << cosThetaPsi_fromAngles  << endl;
// // 		myfile << "dphi_fromAngles = " << dphi_fromAngles  << endl;
// // 		myfile << "DIFF cos Z= " << cosThetaZ_fromAngles - cosThetaZ_michal <<endl;
// // 		myfile << "DIFF cos PsiZ= " << cosThetaPsi_fromAngles - cosThetaPsi_michal <<endl;
// // 		myfile << "DIFF phi Z= " << dphi_fromAngles - dphi_michal <<endl;

// // 		myfile << " " << endl;
// // 		myfile << " " << endl;
// // 		myfile << "======================== TOMASZ ======================== " << endl;
// // 		myfile << "  " << endl;
// // 		myfile << "cosThKs = " << cosThKs <<  endl;
// // 		myfile << "cosThPsi = " << cosThPsi <<  endl;
// // 		myfile << "phiKsPsi = " << phiKsPsi << endl;
// // 		myfile << "m_kpi = " << m_kpi <<  endl;
// // 		myfile << "m_psipi = " << m_psipi << endl;
// // 		myfile << "cosThZ = " << cosThZ <<  endl;
// // 		myfile << "cosThPsi_Z = " << cosThPsi_Z <<  endl;
// // 		myfile << "phiZPsi = " << phiZPsi << endl;
// // 		myfile << "  " << endl;
// // 		myfile << "DIFF phi Ks= " << phiKsPsi - _phi <<endl;
// // 		myfile << "DIFF phi Z= " << cosThPsi_Z - cosThetaZ_michal <<endl;
// // 		myfile << "  " << endl;
// // 		myfile << "  " << endl;
// // 		myfile << "  " << endl;


// 		/*
// 		cout << "sin_phi " << sin_phi << " cos_phi " << cos_phi << endl;
// 		cout << "    phi " << asin(sin_phi) << "    phi " << acos(cos_phi) << endl;
// 		cout << "phi " << _phi << " " << phi << " " << delta_phi <<" " << delta_phi2 <<" " << delta_phi_check <<" " << pion_ID <<endl;
// 		cout << "cos_thetaK " << cos_thetaK << " " << cosTheta2 << endl;
// 		cout << "cos_thetaL " << cos_thetaL << " " << cosTheta1 << endl;
// 		*/
// 	}

// 	tuple->Write();
// 	file->Close();
// 	myfile.close();

// }

// int main()
// {
// 	CalculateAngles t;
// 	t.Loop();
// 	return 1;

// }
