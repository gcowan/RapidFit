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
		double Evaluate(DataPoint*);
		vector<string> GetDoNotIntegrateList();

		double EvaluateComponent( DataPoint*, ComponentRef* );

	protected:
		//Calculate the PDF normalisation
		double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		Exponential operator = ( const Exponential& );

		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);
		double buildPDFnumerator();
		double buildPDFdenominator();
		void prepareTimeInt();

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
		ObservableRef timeOffsetName;
		// Observable
		ObservableRef timeName;		// proper time

		double tau;
		double gamma;
		int sigmaNum;
		vector<PseudoObservable> _intexpIntObs_vec;
		DataPoint* _dataPoint;

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
		bool _useSteppedProjection;
		inline bool useTimeAcceptance() const { return _useTimeAcceptance;}
		inline bool useEventResolution() const {return _useEventResolution;}
};

#endif
