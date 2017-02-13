
#include "DPZplusK.hh"
#include "DPBWResonanceShape.hh"
#include "DPBarrierFactor.hh"
#include "DPHelpers.hh"

#include <iostream>

DPZplusK::DPZplusK(int fLB, int fLR, double fmB, double fmR,
		   double fgammaR, double fm1, double fm2,
		   double fRB, double fRR, double fmJpsi,
		   int spin, int fresonanceIn):
	A0(1,0),
	Aplus(1,0),
	Aminus(1,0),
	spinZplus(spin),
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
    resonanceIn(fresonanceIn)
{
  if ( resonanceIn == 12 )
  {
    massShape = new DPBWResonanceShape(mR, gammaR, LR, m1, m2, RR);
  }
  else if ( resonanceIn == 13 )
  {
    massShape = new DPBWResonanceShape(mR, gammaR, LR, m1, mJpsi, RR);
  }
  else
  {
    massShape = new DPBWResonanceShape(mR, gammaR, LR, m2, mJpsi, RR);
  }
  double m_min = m1 + m2;
  double m_max = mB - mJpsi;
  double m0_eff = m_min + (m_max - m_min)*(1+tanh( (mR - (m_min+m_max)/2.)/(m_max - m_min)))/2;
  double pB0 = DPHelpers::daughterMomentum(mB, mJpsi, m0_eff);
  double pR0 = DPHelpers::daughterMomentum(mR, m1, m2);
  barrierB = new DPBarrierFactor(LB,RB,pB0);
  barrierR = new DPBarrierFactor(LR,RR,pR0);
	switch(spin)
	{
		case 0: wigner=new DPWignerFunctionJ0();
			break;
		case 1: wigner=new DPWignerFunctionJ1();
			break;
		case 2: wigner=new DPWignerFunctionJ2();
			break;
	}
}

DPZplusK::~DPZplusK()
{
	delete barrierR;
	delete barrierB;
	delete massShape;
	delete wigner;
}

DPZplusK::DPZplusK( const DPZplusK& input ) : DPComponent( input ),
    A0(input.A0),
	Aplus(input.Aplus),
	Aminus(input.Aminus),
	spinZplus(input.spinZplus),
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
    resonanceIn(input.resonanceIn),
    massShape(NULL),
    barrierB(NULL),
    barrierR(NULL),
	wigner(NULL), wignerPsi(input.wignerPsi)
{

  if ( input.resonanceIn == 12 )
  {
    massShape = new DPBWResonanceShape(input.mR, input.gammaR, input.LR, input.m1, input.m2, input.RR);
  }
  else if ( input.resonanceIn == 13 )
  {
    massShape = new DPBWResonanceShape(input.mR, input.gammaR, input.LR, input.m1, input.mJpsi, input.RR);
  }
  else
  {
    massShape = new DPBWResonanceShape(input.mR, input.gammaR, input.LR, input.m2, input.mJpsi, input.RR);
  }
  barrierB = new DPBarrierFactor(*input.barrierB);
  barrierR = new DPBarrierFactor(*input.barrierR);
	if( input.wigner != NULL )
	{
		switch(input.spinZplus)
		{
			case 0: wigner=new DPWignerFunctionJ0();
				break;
			case 1: wigner=new DPWignerFunctionJ1();
				break;
			case 2: wigner=new DPWignerFunctionJ2();
				break;
		}
	}
}

/*
 * Arguments of this function are proper for B0 --> psi K*, so first we
 * need to recalculate observables relevant for Z+ K- and then get
 * amplitude
 */
TComplex DPZplusK::amplitude(double m23, double cosTheta1,
			     double cosTheta2, double phi,
			     int twoLambda, int twoLambdaPsi, int pionID)
{
	TComplex result(0,0);


	// First get final state particles momenta
	TLorentzVector pMuPlus;
	TLorentzVector pMuMinus;
	TLorentzVector pPi;
	TLorentzVector pK;
	DPHelpers::calculateFinalStateMomenta(mB, m23, mJpsi,
					      cosTheta1,  cosTheta2, phi, pionID, 0.1056583715, 0.1056583715, m2, m1, pMuPlus, pMuMinus, pPi, pK);
//					      cosTheta1,  cosTheta2, phi, pionID, 0.105, 0.105, m2, m1, pMuPlus, pMuMinus, pPi, pK);
	TLorentzVector pB(0,0,0,mB);
	// Cos of the angle between psi reference axis (not needed in this place)
	//  double cosARefs=DPHelpers::referenceAxisCosAngle(pB, pMuPlus, pMuMinus, pPi, pK);
	// Proper Z+ variables
	double cosThetaZ;
	double cosThetaPsi;
	double dphi;
	DPHelpers::calculateZplusAngles(pB, pMuPlus, pMuMinus, pPi, pK,
					&cosThetaZ, &cosThetaPsi, &dphi, pionID);
	double m13=(pMuPlus+pMuMinus+pPi).M();

	// Call function which calculates amplitude using proper variables
	result=amplitudeProperVars(m13, cosThetaZ, cosThetaPsi, dphi, pionID, twoLambda,
			twoLambdaPsi);

	return result;
}

TComplex DPZplusK::amplitudeProperVars(double m13, double cosTheta1,
				       double cosTheta2, double phi, int pionID,
				       int twoLambda, int twoLambdaPsi)
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
	if ( spinZplus==0 && twoLambdaPsi != 0 )
	{
		return result;
	}

    double m_min = mJpsi + m2;
    double m_max = mB - m1;
    double m0_eff = m_min + (m_max - m_min)*(1+tanh( (mR - (m_min+m_max)/2.)/(m_max - m_min)))/2;
    //std::cout << mR << std::endl;
	double pB = DPHelpers::daughterMomentum(mB, m1, m13);
	double pR = DPHelpers::daughterMomentum(m13, mJpsi, m2);
	double pB0 = DPHelpers::daughterMomentum(mB, m1, m0_eff);
	double pR0 = DPHelpers::daughterMomentum(mR, mJpsi, m2);

/*
 * There are different ways how this can be written. First two are
 * equivalent and differ just by multiplicative factor as denominators are
 * constants in both cases. Third way is one used in Belle Z+(4430) Dalitz
 * analysis, which is only one I know about using this. (MK)
 */
	//double orbitalFactor = TMath::Power(pB/pB0, LB)*
	//	TMath::Power(pR/pR0, LR);
	//double orbitalFactor = TMath::Power(pB/mB, LB)*
	//	TMath::Power(pR/mR, LR);
  double orbitalFactor = TMath::Power(pB/mB, LB)*
      TMath::Power(pR/m13, LR);

	double barrierFactor = barrierB->barrier(pB)*
		barrierR->barrier(pR);
	std::complex<double> tmpmassfactor = massShape->massShape(m13);
	TComplex massFactor(tmpmassfactor.real(),tmpmassfactor.imag());


	// Angular part
	TComplex angular(0,0);
	angular=wignerPsi.function(cosTheta2,twoLambdaPsi/2,twoLambda/2)*
		wigner->function(cosTheta1,0,twoLambdaPsi/2)*
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

    //std::cout << "Z " << twoLambdaPsi << " " << twoLambda << " " << massFactor<< " " << barrierFactor<< " " << orbitalFactor<< " " << angular << std::endl;

	result = massFactor*barrierFactor*orbitalFactor*angular;

	return result;
}

void DPZplusK::setHelicityAmplitudes(double magA0, double magAplus,
		double magAminus, double phaseA0, double phaseAplus,
		double phaseAminus)
{
	A0=TComplex(magA0*TMath::Cos(phaseA0),magA0*TMath::Sin(phaseA0));
	Aplus=TComplex(magAplus*TMath::Cos(phaseAplus),magAplus*TMath::Sin(phaseAplus));
	Aminus=TComplex(magAminus*TMath::Cos(phaseAminus),magAminus*TMath::Sin(phaseAminus));
}

void DPZplusK::setResonanceParameters(double mass, double sigma)
{
    mR = mass;
    gammaR = sigma;
	massShape->setParameters( {mass, sigma} );
    //std::cout << mR << std::endl;
}
