
#include "DPJpsiKaon.hh"
#include "DPBWResonanceShape.hh"
#include "DPNonresonant.hh"
#include "DPLassShape.hh"
#include "DPGLassShape.hh"
#include "DPBarrierFactor.hh"
#include "DPBarrierL0.hh"
#include "DPBarrierL1.hh"
#include "DPBarrierL2.hh"
#include "DPBarrierL3.hh"
#include "DPBarrierL4.hh"
#include "DPBarrierL5.hh"
#include "DPHelpers.hh"
#include "DPWignerFunctionGeneral.hh"

#include <iostream>

DPJpsiKaon::DPJpsiKaon(int fLB, int fLR, double fmB, double fmR,
                               double fgammaR, double fm1, double fm2,
		       double fRB, double fRR, double fmJpsi,
                               int spin,std::string fmShape,
                               double fa, double fr):
A0(0,0),
Aplus(0,0),
Aminus(0,0),
spinKaon(spin),
    m1(fm1),
    m2(fm2),
    LB(fLB),
    LR(fLR),
    mB(fmB),
    mR(fmR),
    gammaR(fgammaR),
    mJpsi(fmJpsi),
    RB(fRB),
    RR(fRR),
    a(fa),
    r(fr),
    mShape(fmShape)
{

    //std::cout << m1 << " " << m2 << " " << LB << " " << LR <<" " << gammaR <<  " " << mJpsi << " " << RB <<  " " << RR << " " << r << " " << a << " " << mShape << std::endl;
  if ( mShape == "LASS" )
  {
    std::cout<<"Using Lass parametrization\n";
    massShape = new DPLassShape(mR, gammaR, LR, m1, m2, RR, a, r);
  }
  else if ( mShape == "gLASS" )
  {
    std::cout<<"Using generalize Lass parametrization\n";
    massShape = new DPGLassShape(mR, gammaR, LR, m1, m2, RR, a, r);
  }
  else if ( mShape == "NR" )
  {
    std::cout<<"Using nonresonant component\n";
    massShape = new DPNonresonant(mR, gammaR, LR, m1, m2, RR);
  }
  else
  {
    std::cout<<"Using BW component\n";
    massShape = new DPBWResonanceShape(mR, gammaR, LR, m1, m2, RR);
  }
  switch (LB)
  {
    case 0: barrierB=new DPBarrierL0(RB);
            break;
    case 1: barrierB=new DPBarrierL1(RB);
            break;
    case 2: barrierB=new DPBarrierL2(RB);
            break;
    case 3: barrierB=new DPBarrierL3(RB);
            break;
    case 4: barrierB=new DPBarrierL4(RB);
            break;
    case 5: barrierB=new DPBarrierL5(RB);
            break;
    default: std::cout<<"WARNING DPJpsiKaon (LB): Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
             barrierB=new DPBarrierL0(RB);
             break;
  }
  switch (LR)
  {
    case 0: barrierR=new DPBarrierL0(RR);
            break;
    case 1: barrierR=new DPBarrierL1(RR);
            break;
    case 2: barrierR=new DPBarrierL2(RR);
            break;
    case 3: barrierR=new DPBarrierL3(RR);
            break;
    case 4: barrierR=new DPBarrierL4(RR);
            break;
    case 5: barrierR=new DPBarrierL5(RR);
            break;
    default: std::cout<<"WARNING DPJpsiKaon (LR): Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
             barrierR=new DPBarrierL0(RR);
             break;
  }
  switch(spinKaon)
  {
    case 0: wigner=new DPWignerFunctionJ0();
            break;
    case 1: wigner=new DPWignerFunctionJ1();
            break;
    case 2: wigner=new DPWignerFunctionJ2();
            break;
    case 3: wigner=new DPWignerFunctionGeneral(3);
            break;
    case 4: wigner=new DPWignerFunctionGeneral(4);
            break;
    case 5: wigner=new DPWignerFunctionGeneral(5);
            break;
  }
}

DPJpsiKaon::~DPJpsiKaon()
{
  //std::cout << "DPJpsiKaon dest" << massShape << std::endl;
  delete massShape;
  delete barrierR;
  delete barrierB;
  delete wigner;
}

DPJpsiKaon::DPJpsiKaon( const DPJpsiKaon& input ) : DPComponent( input ),
    A0(input.A0),
    Aplus(input.Aplus),
    Aminus(input.Aminus),
    spinKaon(input.spinKaon),
    m1(input.m1),
    m2(input.m2),
    LB(input.LB),
    LR(input.LR),
    mB(input.mB),
    mR(input.mR),
    gammaR(input.gammaR),
    mJpsi(input.mJpsi),
    RB(input.RB),
    RR(input.RR),
    mShape(input.mShape),
    a(input.a),
    r(input.r),
    massShape(NULL),
    barrierB(NULL),
    barrierR(NULL),
    wigner(NULL), wignerPsi(input.wignerPsi)
{
    //std::cout << "In DPJpsiKaon copy const " << mR << " " << gammaR << " " << LR << " " << LB << " " << RR << " " << RB << " " << mShape << std::endl;
    //std::cout << "In DPJpsiKaon copy const " << input.mR << " " << input.gammaR << " " << input.LR << " " << input.LB << " " << input.RR << " " << input.RB << " " << input.mShape << std::endl;
  if ( input.mShape == "LASS" )
  {
    massShape = new DPLassShape(input.mR, input.gammaR, input.LR, input.m1, input.m2, input.RR, input.a, input.r);
  }
  else if ( input.mShape == "gLASS" )
  {
    massShape = new DPGLassShape(input.mR, input.gammaR, input.LR, input.m1, input.m2, input.RR, input.a, input.r);
  }
  else if ( input.mShape == "NR" )
  {
    massShape = new DPNonresonant(input.mR, input.gammaR, input.LR, input.m1, input.m2, input.RR);
  }
  else
  {
    massShape = new DPBWResonanceShape(input.mR, input.gammaR, input.LR, input.m1, input.m2, input.RR);
  }
  switch (input.LB)
  {
      case 0: //std::cout << "Hi, I am about to construct DPBarrierL0" << std::endl;
            barrierB=new DPBarrierL0(input.RB);
            break;
    case 1: barrierB=new DPBarrierL1(input.RB);
            break;
    case 2: barrierB=new DPBarrierL2(input.RB);
            break;
    case 3: barrierB=new DPBarrierL3(input.RB);
            break;
    case 4: barrierB=new DPBarrierL4(input.RB);
            break;
    case 5: barrierB=new DPBarrierL5(input.RB);
            break;
    default: std::cout<<"WARNING DPJpsiKaon (LB): Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
             barrierB=new DPBarrierL0(input.RB);
             break;
  }
  switch (input.LR)
  {
    case 0: barrierR=new DPBarrierL0(input.RR);
            break;
    case 1: barrierR=new DPBarrierL1(input.RR);
            break;
    case 2: barrierR=new DPBarrierL2(input.RR);
            break;
    case 3: barrierR=new DPBarrierL3(input.RR);
            break;
    case 4: barrierR=new DPBarrierL4(input.RR);
            break;
    case 5: barrierR=new DPBarrierL5(input.RR);
            break;
    default: std::cout<<"WARNING DPJpsiKaon (LR): Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
             barrierR=new DPBarrierL0(input.RR);
             break;
  }
        if( input.wigner != NULL )
        {
                switch(input.spinKaon)
                {
                        case 0: wigner=new DPWignerFunctionJ0();
                        break;
                        case 1: wigner=new DPWignerFunctionJ1();
                        break;
                        case 2: wigner=new DPWignerFunctionJ2();
                        break;
                        case 3: wigner=new DPWignerFunctionGeneral(3);
                        break;
                        case 4: wigner=new DPWignerFunctionGeneral(4);
                        break;
                        case 5: wigner=new DPWignerFunctionGeneral(5);
                        break;
                }
        }
	//std::cout<<"DEBUG, used copy constructor\n";
	//std::cout<<"BarierFactors are "<<barrierB<<" "<<barrierR<<std::endl;
	//std::cout<<"Parameters are "<<m1 << " "<<mR<< " " << spinKaon<<" "<<A0<<std::endl;
}

TComplex DPJpsiKaon::amplitude(double m23, double cosTheta1,
			       double cosTheta2, double phi,
			       int twoLambda, int twoLambdaPsi, int pionID = -211)
{
  TComplex result(0,0);
  // Muon helicity
  if ( twoLambda != 2 && twoLambda != -2 ) // Not right option
  {
    return result;
  }

  // Psi helicity
  if ( twoLambdaPsi != 2 && twoLambdaPsi != -2 && twoLambdaPsi != 0 ) // Not right option
  {
    return result;
  }

  //
  if ( spinKaon==0 && twoLambdaPsi != 0 )
  {
    return result;
  }

  double m_min = m1 + m2;
  double m_max = mB - mJpsi;
  double m0_eff = m_min + (m_max - m_min)*(1+tanh( (mR - (m_min+m_max)/2.)/(m_max - m_min)))/2;

  //std::cout << "B" << mB << " Jpsi " << mJpsi << " m23 " << m23 << " m1 " << m1 << " m2 " << m2 << " mR " << mR << " m0_eff " << m0_eff << std::endl;
  double pB = DPHelpers::daughterMomentum(mB, mJpsi, m23);      // B, psi, K*
  double pR = DPHelpers::daughterMomentum(m23, m1, m2);         // K*, K, pi
  double pB0 = DPHelpers::daughterMomentum(mB, mJpsi, m0_eff);// B, psi, resonance
  double pR0 = DPHelpers::daughterMomentum(mR, m1, m2);   // resonance, K, pi // this is really mR here, not m0_eff

/*
 * There are different ways how this can be written. First two are
 * equivalent and differ just by multiplicative factor as denominators are
 * constants in both cases. Third way is one used in Belle Z+(4430) Dalitz
 * analysis, which is only one I know about using this. (MK)
 */
    //double orbitalFactor = TMath::Power(pB/pB0, LB)*
    //                       TMath::Power(pR/pR0, LR);
	//double orbitalFactor = TMath::Power(pB/mB, LB)*
	//	                     TMath::Power(pR/mR, LR);
  double orbitalFactor = TMath::Power(pB/mB, LB)*
                         TMath::Power(pR/m23, LR);

  double barrierFactor = barrierB->barrier( pB0, pB )*
                         barrierR->barrier( pR0, pR );

  TComplex massFactor = massShape->massShape(m23);
  if ( mShape == "NR" )
  {
    barrierFactor = 1.;
    orbitalFactor = pB/mB;
  }

  if (isnan(pR0)) std::cout << mShape << " " << mR << " " << pB << " " << pR <<  " " << pB0 << " " << pR0 << " " << orbitalFactor << " " << barrierFactor << "  " << massFactor << std::endl;

  // Angular part
  TComplex angular(0,0);
  angular=wignerPsi.function(cosTheta1,twoLambdaPsi/2,twoLambda/2)*
          wigner->function(cosTheta2,twoLambdaPsi/2,0)*
          TComplex::Exp(-0.5*twoLambdaPsi*TComplex::I()*phi);
  switch (twoLambdaPsi)
  {
    case 0: angular*=A0;
            break;
    case 2: angular*=Aplus;
            break;
    case -2: angular*=Aminus;
            break;
  }
  result = massFactor*barrierFactor*orbitalFactor*angular;

  return result;
}

void DPJpsiKaon::setHelicityAmplitudes(double magA0, double magAplus,
                     double magAminus, double phaseA0, double phaseAplus,
                     double phaseAminus)
{
  A0=TComplex(magA0*TMath::Cos(phaseA0),magA0*TMath::Sin(phaseA0));
  Aplus=TComplex(magAplus*TMath::Cos(phaseAplus),magAplus*TMath::Sin(phaseAplus));
  Aminus=TComplex(magAminus*TMath::Cos(phaseAminus),magAminus*TMath::Sin(phaseAminus));
}

void DPJpsiKaon::setResonanceParameters(double mass, double sigma)
{
        mR = mass;
        gammaR = sigma;
        massShape->setResonanceParameters( mass, sigma );
}

TComplex DPJpsiKaon::amplitudeProperVars(double m23, double cosTheta1,
					 double cosTheta2, double phi, int pionID,
                             int twoLambda, int twoLambdaPsi)
{
  return amplitude(m23,cosTheta1, cosTheta2, phi, twoLambda, twoLambdaPsi);
}
