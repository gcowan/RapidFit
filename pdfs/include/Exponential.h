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
#include "SlicedAcceptance.h"

class Exponential : public BasePDF
{
	public:
		Exponential( PDFConfigurator* );
		Exponential( const Exponential & );
		~Exponential();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);
		virtual vector<string> GetDoNotIntegrateList();

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation( PhaseSpaceBoundary* );
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);
		double buildPDFnumerator();
		double buildPDFdenominator();		

		// Physics parameters
		ObservableRef tauName;		// decay constant 1
		// Detector parameters
		ObservableRef eventResolutionName;                      // Scale to multiply all Gaussians with 
		ObservableRef resScale1Name;                     // Scale to multiply all Gaussians with 
		ObservableRef resScale2Name;                     // Scale to multiply all Gaussians with 
		ObservableRef resScale3Name;                     // Scale to multiply all Gaussians with 
		ObservableRef sigma1Name;		// time res sigma 1
		ObservableRef sigma2Name;		// time res sigma 2
		ObservableRef sigma3Name;		// time res sigma 2
		ObservableRef timeRes2FracName;
		ObservableRef timeRes3FracName;
		// Observable
		ObservableRef timeName;		// proper time

		double tau;
		
		// This stuff is all to do with resolution
		double sigma; 
		double sigma1; 
		double sigma2;
		double sigma3;
		double timeRes2Frac;
		double timeRes3Frac;
		double eventResolution;
		double resolutionScale1;
		double resolutionScale2;
		double resolutionScale3;
		bool _useEventResolution;
		inline bool useEventResolution() const {return _useEventResolution;}

		double tlow, thigh; // integration limits

		// These contain the ObservableRefs that correspond
		// to the observable names that are used in the
		// PDF. 
		double time;

		//Time acceptance 
		SlicedAcceptance * timeAcc;

		//Configurationparameters
		bool _useTimeAcceptance;
		bool _numericIntegralForce;
        	inline bool useTimeAcceptance() const { return _useTimeAcceptance;}

};

#endif
