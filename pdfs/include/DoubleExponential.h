// $Id: DoubleExponential.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class DoubleExponential DoubleExponential.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution.
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef DoubleExponential_H
#define DoubleExponential_H

#include "BasePDF.h"
#include "SlicedAcceptance.h"

class DoubleExponential : public BasePDF
{
	public:
		DoubleExponential( PDFConfigurator* );
		DoubleExponential( const DoubleExponential & );
		~DoubleExponential();

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
		ObservableRef tau1Name;		// decay constant 1
		ObservableRef tau2Name;		// decay constant 2
		ObservableRef fraction1Name;	// fraction of tau1
		// Detector parameters
		ObservableRef eventResolutionName;                      // Scale to multiply all Gaussians with 
		ObservableRef resScale1Name;                     // Scale to multiply all Gaussians with 
		ObservableRef resScale2Name;                     // Scale to multiply all Gaussians with 
		ObservableRef resScale3Name;                     // Scale to multiply all Gaussians with 
		ObservableRef timeOffsetName;  
		ObservableRef sigma1Name;		// time res sigma 1
		ObservableRef sigma2Name;		// time res sigma 2
		ObservableRef sigma3Name;		// time res sigma 2
		ObservableRef timeRes2FracName;
		ObservableRef timeRes3FracName;
		// Observable
		ObservableRef timeName;		// proper time

		double tau1;
		double tau2;
		double fraction1;
		
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
		double timeOffset;

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
		bool _useEventResolution;
		bool _usePunziSigmat;
		inline bool useEventResolution() const {return _useEventResolution;}
        	inline bool useTimeAcceptance() const { return _useTimeAcceptance;}
};

#endif
