
///	RapidFit Headers
#include "Bs2JpsiPhi_Angluar_Terms.h"
#include "Mathematics.h"
///	System Headers
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace::std;

/*!
 * These functions Assume that the they are provided with an input of vector<double>
 *
 * The format of the input is:
 *
 * Transversity:   [0] = cosTheta , [1] = cosPsi   , [2] = phi
 *
 * Helicity:       [0] = cosThetaK, [1] = cosThetaL, [2] = phi
 */

//------------------------------------------------------
// Angle factors for three angle PDFs  in transversity basis

//........ P Wave ..........


//........ for one angle tests ...................
double Bs2JpsiPhi_Angular_Terms::TangleFactorEven( vector<double> input )
{
	double cosTheta = input[0];
	double cosPsi = input[1];(void) cosPsi;
	double phi = input[2];(void) phi;

	return (1.0 + cosTheta*cosTheta) * Mathematics::Global_Frac()  * 4.0/3.0*Mathematics::Pi();
}

double Bs2JpsiPhi_Angular_Terms::TangleFactorOdd( vector<double> input )
{
	double cosTheta = input[0];
	double cosPsi = input[1];(void) cosPsi;
	double phi = input[2];(void) phi;

	return  sin( acos( cosTheta) ) * Mathematics::Global_Frac() * 8.0/3.0*Mathematics::Pi();
}


//...........................
double Bs2JpsiPhi_Angular_Terms::TangleFactorA0A0( vector<double> input )
{
	double cosTheta = input[0];
	double sinTheta = sin(acos(cosTheta));
	double cosPsi = input[1];
	double phi = input[2];
	double cosPhi = cos( phi );

	return 2.0 * cosPsi*cosPsi * (1.0 - sinTheta*sinTheta* cosPhi*cosPhi ) * Mathematics::Global_Frac();
}

//...........................
double Bs2JpsiPhi_Angular_Terms::TangleFactorAPAP( vector<double> input )
{
	double cosTheta = input[0];
	double sinTheta = sin(acos(cosTheta));
	double cosPsi = input[1];
	double sinPsi = sin(acos(cosPsi));
	double phi = input[2];
	double sinPhi = sin( phi );

	return sinPsi*sinPsi * (1.0 - sinTheta*sinTheta * sinPhi*sinPhi ) * Mathematics::Global_Frac();
}

//...........................
double Bs2JpsiPhi_Angular_Terms::TangleFactorATAT( vector<double> input )
{
	double cosTheta = input[0];
	double sinTheta = sin(acos(cosTheta));
	double cosPsi = input[1];
	double sinPsi = sin(acos(cosPsi));
	double phi = input[2];(void) phi;

	return sinPsi*sinPsi * sinTheta*sinTheta * Mathematics::Global_Frac();
}

//...........................
double Bs2JpsiPhi_Angular_Terms::TangleFactorImAPAT( vector<double> input )
{
	double cosTheta = input[0];
	double sin2Theta = sin(2.*acos(cosTheta));
	double cosPsi = input[1];
	double sinPsi = sin(acos(cosPsi));
	double phi = input[2];
	double sinPhi = sin(phi);

	return  -1. * sinPsi*sinPsi * sin2Theta * sinPhi * Mathematics::Global_Frac();
}

//...........................
double Bs2JpsiPhi_Angular_Terms::TangleFactorReA0AP( vector<double> input )
{
	double cosTheta = input[0];
	double sinTheta = sin(acos(cosTheta));
	double cosPsi = input[1];
	double sin2Psi = sin(2.*acos(cosPsi));
	double phi = input[2];
	double sin2Phi = sin(2.*phi);

	return Mathematics::_Over_SQRT_2() * sin2Psi * sinTheta*sinTheta * sin2Phi * Mathematics::Global_Frac();
}

//...........................
double Bs2JpsiPhi_Angular_Terms::TangleFactorImA0AT( vector<double> input )
{
	double cosTheta = input[0];
	double sin2Theta = sin(2.*acos(cosTheta));
	double cosPsi = input[1];
	double sin2Psi = sin(2.*acos(cosPsi));
	double phi = input[2];
	double cosPhi = cos(phi);

	return Mathematics::_Over_SQRT_2() * sin2Psi * sin2Theta * cosPhi * Mathematics::Global_Frac();
}

//......  S wave  ....

//.............................
double Bs2JpsiPhi_Angular_Terms::TangleFactorASAS( vector<double> input )
{
	double cosTheta = input[0];
	double sinTheta = sin(acos(cosTheta));
	double cosPsi = input[1];(void) cosPsi;
	double phi = input[2];
	double cosPhi = cos(phi);

	return  2.0*Mathematics::Third() * (1.0 - sinTheta*sinTheta * cosPhi*cosPhi ) * Mathematics::Global_Frac();
}

//...........................
double Bs2JpsiPhi_Angular_Terms::TangleFactorReASAP( vector<double> input )
{
	double cosTheta = input[0];
	double sinTheta = sin(acos(cosTheta));
	double cosPsi = input[1];
	double sinPsi = sin(acos(cosPsi));
	double phi = input[2];
	double sin2Phi = sin(2.*phi);

	// Corrected by PELC 130419 - it should have been sin^theta
    //return Mathematics::Root_6()*Mathematics::Third() * sinPsi * sinTheta * sin2Phi *  Mathematics::Global_Frac();
	return Mathematics::Root_6()*Mathematics::Third() * sinPsi * sinTheta*sinTheta * sin2Phi *  Mathematics::Global_Frac();
}

//...........................
// This appreas to be ifferent to the LHCB note, but on inspection it is not. It is the difference in sign of ImASAT <=> ImATAS
double Bs2JpsiPhi_Angular_Terms::TangleFactorImASAT( vector<double> input )
{
	double cosTheta = input[0];
	double sin2Theta = sin(2.*acos(cosTheta));
	double cosPsi = input[1];
	double sinPsi = sin(acos(cosPsi));
	double phi = input[2];
	double cosPhi = cos(phi);

	return  Mathematics::Root_6()*Mathematics::Third() *  sinPsi * sin2Theta *  cosPhi *  Mathematics::Global_Frac();
}


//...........................
double Bs2JpsiPhi_Angular_Terms::TangleFactorReASA0( vector<double> input )
{
	double cosTheta = input[0];
	double sinTheta = sin(acos(cosTheta));
	double cosPsi = input[1];
	double phi = input[2];
	double cosPhi = cos(phi);

	//There was a -1.0 here
	//Then I looked in the LHCb note being drafted and found this disagreed with it
	//For now i have made it +1.0 to agree with draft LHCb note.
	// Since then ive proved that this set of signs is consistent by my "pdfvalue < 0 test"
	return  4.0*Mathematics::Root_3()*Mathematics::Third() * cosPsi *  ( 1.0 - sinTheta*sinTheta * cosPhi*cosPhi  ) * Mathematics::Global_Frac();
}



//------------------------------------------------------
// Angle factors for three angle PDFs  in helicity basis

//........ for one angle tests ...................
double Bs2JpsiPhi_Angular_Terms::HangleFactorEven( vector<double> input )
{
	(void) input;
	cout<<"No Helicity One Angle Formula Yet"<<endl; exit(1); return 0;
}
double Bs2JpsiPhi_Angular_Terms::HangleFactorOdd( vector<double> input )
{
	(void) input;
	cout<<"No Helicity One Angle Formula Yet"<<endl; exit(1); return 0;
}


//........ P Wave ..........

// .................
double Bs2JpsiPhi_Angular_Terms::HangleFactorA0A0( vector<double> input )
{
	double costhetaK = input[0];
	double cossqthetaK = costhetaK*costhetaK;
	double costhetaL = input[1];
	double sinsqthetaL = 1.-costhetaL*costhetaL;
	double phi = input[2];(void) phi;

	return 4. * cossqthetaK * sinsqthetaL * Mathematics::Global_Frac() * 0.5;
}

//..................
double Bs2JpsiPhi_Angular_Terms::HangleFactorAPAP( vector<double> input )
{
	double costhetaK = input[0];
	double sinsqthetaK = 1.-costhetaK*costhetaK;
	double costhetaL = input[1];
	double cossqthetaL = costhetaL*costhetaL;
	double sinsqthetaL = 1.-cossqthetaL;
	double phi = input[2];
	double cos2Phi = cos(2.*phi);

	return (sinsqthetaK * (1.+ cossqthetaL) - sinsqthetaL * sinsqthetaK * cos2Phi ) * Mathematics::Global_Frac() * 0.5;
}

//..................
double Bs2JpsiPhi_Angular_Terms::HangleFactorATAT( vector<double> input )
{
	double costhetaK = input[0];
	double sinsqthetaK = 1.-costhetaK*costhetaK;
	double costhetaL = input[1];
	double cossqthetaL = costhetaL*costhetaL;
	double sinsqthetaL = 1.-cossqthetaL;
	double phi = input[2];
	double cos2Phi = cos(2.*phi);

	return (sinsqthetaK * (1.+ cossqthetaL) + sinsqthetaL * sinsqthetaK * cos2Phi ) * Mathematics::Global_Frac() * 0.5;
}

//..................
double Bs2JpsiPhi_Angular_Terms::HangleFactorImAPAT( vector<double> input )
{
	double costhetaK = input[0];
	double sinsqthetaK = 1.-costhetaK*costhetaK;
	double costhetaL = input[1];
	double sinsqthetaL = 1.-costhetaL*costhetaL;
	double phi = input[2];
	double sin2Phi = sin(2.*phi);

	return 2. * sinsqthetaK * sinsqthetaL * sin2Phi * Mathematics::Global_Frac()* 0.5;
}

//..................
double Bs2JpsiPhi_Angular_Terms::HangleFactorReA0AP( vector<double> input )
{
	double costhetaK = input[0];
	double sin2thetaK = sin(2.*acos(costhetaK));
	double costhetaL = input[1];
	double sin2thetaL = sin(2.*acos(costhetaL));
	double phi = input[2];
	double cosPhi = cos(phi);

	return -Mathematics::Root_2() * sin2thetaK * sin2thetaL * cosPhi * Mathematics::Global_Frac()* 0.5;
}

//..................
double Bs2JpsiPhi_Angular_Terms::HangleFactorImA0AT( vector<double> input )
{
	double costhetaK = input[0];
	double sin2thetaK = sin(2.*acos(costhetaK));
	double costhetaL = input[1];
	double sin2thetaL = sin(2.*acos(costhetaL));
	double phi = input[2];
	double sinPhi = sin(phi);

	return Mathematics::Root_2() * sin2thetaK * sin2thetaL * sinPhi * Mathematics::Global_Frac()* 0.5;
}

//......  S wave  ....

//.............................
double Bs2JpsiPhi_Angular_Terms::HangleFactorASAS( vector<double> input )
{
	double costhetaK = input[0];(void) costhetaK;
	double costhetaL = input[1];
	double sinsqthetaL = 1.-costhetaL*costhetaL;
	double phi = input[2];(void) phi;

	return 4.0*Mathematics::Third() * sinsqthetaL * Mathematics::Global_Frac() * 0.5;
}

//...........................
double Bs2JpsiPhi_Angular_Terms::HangleFactorReASAP( vector<double> input )
{
	double costhetaK = input[0];
	double sinthetaK = sin(acos(costhetaK));
	double costhetaL = input[1];
	double sin2thetaL = sin(2.*acos(costhetaL));
	double phi = input[2];
	double cosPhi = cos(phi);

	return -2.0*Mathematics::Third() * Mathematics::Root_6() * sin2thetaL * sinthetaK * cosPhi * Mathematics::Global_Frac() * 0.5;
}

//...........................
double Bs2JpsiPhi_Angular_Terms::HangleFactorImASAT( vector<double> input )
{
	double costhetaK = input[0];
	double sinthetaK = sin(acos(costhetaK));
	double costhetaL = input[1];
	double sin2thetaL = sin(2.*acos(costhetaL));
	double phi = input[2];
	double sinPhi = sin(phi);

	return 2.0*Mathematics::Third() * Mathematics::Root_6() * sin2thetaL * sinthetaK * sinPhi * Mathematics::Global_Frac() * 0.5;
}

//...........................
double Bs2JpsiPhi_Angular_Terms::HangleFactorReASA0( vector<double> input )
{
	double costhetaK = input[0];
	double costhetaL = input[1];
	double sinsqthetaL = 1.-costhetaL*costhetaL;
	double phi = input[2];(void) phi;

	return 8.0*Mathematics::Third() *Mathematics::Root_3() * sinsqthetaL * costhetaK * Mathematics::Global_Frac() * 0.5;
}

