// $Id: Bs2JpsiPhiLongLivedBkg_II.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Bs2JpsiPhiLongLivedBkg_II Bs2JpsiPhiLongLivedBkg_II.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution.
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef Bs2JpsiPhiLongLivedBkg_II_H
#define Bs2JpsiPhiLongLivedBkg_II_H

#include "BasePDF.h"

class Bs2JpsiPhiLongLivedBkg_II : public BasePDF
{
	public:
		Bs2JpsiPhiLongLivedBkg_II( PDFConfigurator* );
		~Bs2JpsiPhiLongLivedBkg_II();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);
		double buildPDFnumerator();
		double buildPDFdenominator();		

		// Physics parameters
		ObservableRef tauLL1Name;		// decay constant 1
		ObservableRef tauLL2Name;		// decay constant 2
		ObservableRef f_LL1Name;		// fraction
		ObservableRef sigmaLL1Name;		// time res sigma 1
		ObservableRef sigmaLL2Name;		// time res sigma 2
		ObservableRef timeResLL1FracName;

		ObservableRef timeName;
		ObservableRef constraint_timeName;

		double tauLL1;
		double tauLL2;
		double f_LL1;
		double sigmaLL; // This is the member variable used in the "builder" functions 
		double sigmaLL1; // These are the physics parameters varied in the fit and passed from the XML;
		double sigmaLL2;
		double timeResLL1Frac;

		double tlow, thigh; // integration limits

		double time;
};

#endif

