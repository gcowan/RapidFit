// $Id: Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist_withAngDist.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution + non-trivial angular distribution
 *
 *  @author Pete Clarke
 *  @date 2009-10-04
 */

#ifndef Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist_H
#define Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist_H

#include "BasePDF.h"

class Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist : public BasePDF
{
	public:
		Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist();
		~Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);
		double buildPDFnumerator();
		double buildPDFdenominator();	
		double angularFactor();

		// Physics parameters
		string tauLL1Name;		// decay constant 1
		string tauLL2Name;		// decay constant 2
		string f_LL1Name;		// fraction
		string sigmaLL1Name;		// time res sigma 1
		string sigmaLL2Name;		// time res sigma 2
                string timeResLL1FracName;

		double tauLL1;
		double tauLL2;
		double f_LL1;
		double sigmaLL; // This is the member variable used in the "builder" functions 
		double sigmaLL1; // These are the physics parameters varied in the fit and passed from the XML;
		double sigmaLL2;
                double timeResLL1Frac;

		double tlow, thigh; // integration limits

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF. 
		string timeName;	// proper time
		string cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		string phiName;		// azimuthal angle of the mu+ in Jpsi frame
		string cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
		double time;
		double cosTheta;
		double phi;
		double cosPsi;
		
};

#endif
