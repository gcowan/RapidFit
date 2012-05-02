// $Id: LongLivedBkg.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class LongLivedBkg LongLivedBkg_withAngDist.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution + non-trivial angular distribution realised by a histogram
 *
 *  @author Ailsa Sparkes
 *  @date 2011-05-30
 */

#ifndef LongLivedBkg_H
#define LongLivedBkg_H

#include "TROOT.h"
#include "TFile.h"
#include "TH3D.h"

#include "BasePDF.h"
#include "SlicedAcceptance.h"

class LongLivedBkg : public BasePDF
{
	public:
		//LongLivedBkg();
		LongLivedBkg( PDFConfigurator* );
		//LongLivedBkg( const LongLivedBkg& );
		~LongLivedBkg();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);
		virtual double Norm(PhaseSpaceBoundary*);

	private:
		LongLivedBkg& operator=( const LongLivedBkg& );
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);
		double buildPDFnumerator();
		double buildPDFdenominator();

		// Physics parameters
		ObservableRef f_LL1Name;		// fraction of decay const 1
		ObservableRef tauLL1Name;            // decay constant 1
		ObservableRef tauLL2Name;            // decay constant 2

		//Detector parameters
		ObservableRef timeResLL1FracName; //fraction of timeres 1
		ObservableRef sigmaLL1Name;	// time res sigma 1
		ObservableRef sigmaLL2Name;	// time res sigma 2

		// These contain the strings that correspond
		ObservableRef timeName;		// proper time

		double tauLL1;
		double tauLL2;
		double f_LL1;
		double sigmaLL; 
		double sigmaLL1; 
		double sigmaLL2;
		double timeResLL1Frac;

		double tlow, thigh; // time integration limits

		double time;

		//Additions to deal with 3-D angular distribution via a histogram
		bool _useTimeAcceptance;
		SlicedAcceptance * timeAcc;
};

#endif
