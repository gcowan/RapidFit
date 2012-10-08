
#include "DPZplusK.hh"
#include "DPBWResonanceShape.hh"
#include "DPBarrierFactor.hh"
#include "DPBarrierL0.hh"
#include "DPBarrierL1.hh"
#include "DPBarrierL2.hh"
#include "DPBarrierL3.hh"
#include "DPHelpers.hh"

#include <iostream>

DPZplusK::DPZplusK(int fLB, int fLR, double fmB, double mR, 
		double gammaR, double m1, double m2, 
		double RB, double RR, double fmJpsi,
		int spin):
	A0(1,0),
	Aplus(1,0),
	Aminus(1,0),
	spinZplus(spin)
{
	this->LB=fLB;
	this->LR=fLR;
	this->mB=fmB;
	this->mR=mR;
	this->m1=m1;
	this->m2=m2;
	mJpsi=fmJpsi;
	massShape = new DPBWResonanceShape(mR, gammaR, LR, m1, m2, RR);
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
		default: std::cout<<"WARNING DPZplusK (LB): Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
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
		default: std::cout<<"WARNING DPZplusK (LR): Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
			 barrierR=new DPBarrierL0(RR);
			 break;
	}
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
	if( wigner != NULL ) delete wigner;
}

DPZplusK::DPZplusK( const DPZplusK& input ) : DPComponent( input ),
	mJpsi(input.mJpsi), m1(input.m1), m2(input.m2), pR0(input.pR0), A0(input.A0), Aplus(input.Aplus),
	Aminus(input.Aminus), spinZplus(input.spinZplus), wigner(NULL), wignerPsi(input.wignerPsi)
{  
	if( input.wigner != NULL )
	{
		switch(spinZplus)
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
		int twoLambda, int twoLambdaPsi)
{
	TComplex result(0,0);


	// First get final state particles momenta
	TLorentzVector pMuPlus;
	TLorentzVector pMuMinus;
	TLorentzVector pPi;
	TLorentzVector pK;
	DPHelpers::calculateFinalStateMomenta(this->mB, m23, mJpsi,
			cosTheta1,  cosTheta2, phi, 0.105, 0.105, m2, m1, pMuPlus, pMuMinus, pPi, pK);
	TLorentzVector pB(0,0,0,this->mB);
	// Cos of the angle between psi reference axis (not needed in this place)
	//  double cosARefs=DPHelpers::referenceAxisCosAngle(pB, pMuPlus, pMuMinus, pPi, pK);
	// Proper Z+ variables
	double cosThetaZ;
	double cosThetaPsi;
	double dphi;
	DPHelpers::calculateZplusAngles(pB, pMuPlus, pMuMinus, pPi, pK,
			&cosThetaZ, &cosThetaPsi, &dphi);
	double m13=(pMuPlus+pMuMinus+pPi).M();

	// Call function which calculates amplitude using proper variables
	result=amplitudeProperVars(m13, cosThetaZ, cosThetaPsi, dphi, twoLambda,
			twoLambdaPsi);

	return result;
}

TComplex DPZplusK::amplitudeProperVars(double m13, double cosTheta1, 
		double cosTheta2, double phi, 
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

	double pB = DPHelpers::daughterMomentum(this->mB, this->m1, m13);
	double pR = DPHelpers::daughterMomentum(m13, this->mJpsi, this->m2);

	double orbitalFactor = TMath::Power(pB/this->mB, this->LB)*
		TMath::Power(pR/m13, this->LR);

	double barrierFactor = barrierB->barrier( DPHelpers::daughterMomentum(this->mB,
				this->m1, this->mR), pB)*
		barrierR->barrier(DPHelpers::daughterMomentum(this->mR,
					this->mJpsi, this->m2), pR);

	TComplex massFactor = this->massShape->massShape(m13);


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
	massShape->setResonanceParameters( mass, sigma );
}
