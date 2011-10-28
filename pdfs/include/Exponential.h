// $Id: Exponential.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Exponential Exponential.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution.
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef Exponential_H
#define Exponential_H

#include "BasePDF.h"

class Exponential : public BasePDF
{
	public:
		Exponential( PDFConfigurator* );
		~Exponential();

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

		// Physics parameters
		ObservableRef tauLL1Name;		// decay constant 1
		ObservableRef sigmaLL1Name;		// time res sigma 1
		ObservableRef sigmaLL2Name;		// time res sigma 2
                ObservableRef timeResLL1FracName;
		ObservableRef timeName;		// proper time

		double tauLL1;
		double sigmaLL; // This is the member variable used in the "builder" functions 
		double sigmaLL1; // These are the physics parameters varied in the fit and passed from the XML;
		double sigmaLL2;
                double timeResLL1Frac;

		double tlow, thigh; // integration limits

		// These contain the ObservableRefs that correspond
		// to the observable names that are used in the
		// PDF. 
		double time;
};

#endif
