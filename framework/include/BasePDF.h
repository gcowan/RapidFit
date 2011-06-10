/**
        @class BasePDF

        Class that provides a general implementation of IPDF.
        Can inherit from this to make a PDF without worrying about the details.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef BASE_PDF_H
#define BASE_PDF_H

//	RapidFit Headers
#include "IPDF.h"
#include "ObservableRef.h"
#include "PDFConfigurator.h"
//	System Headers
#include <vector>

class BasePDF : public IPDF
{
	public:
		BasePDF();
		virtual ~BasePDF();

		//Indicate whether the function has been set up correctly
		virtual bool IsValid();

		//Set the function parameters
		virtual bool SetPhysicsParameters( ParameterSet* );

		//Return the integral of the function over the given boundary
		virtual double Integral( DataPoint*, PhaseSpaceBoundary* );

		//Return the function value at the given point
		virtual double Evaluate( DataPoint* );

		//Return the function value at the given point for use in numeric integral
		virtual double EvaluateForNumericIntegral( DataPoint* );
	
		//Return the components of the function value at the given point
		virtual vector<double> EvaluateComponents( DataPoint* );

		//Return a prototype data point
		virtual vector<string> GetPrototypeDataPoint();

		//Return a prototype set of physics parameters
		virtual vector<string> GetPrototypeParameterSet();

		//Return a list of parameters not to be integrated
		virtual vector<string> GetDoNotIntegrateList();

		//Update integral cache
		virtual void UpdateIntegralCache();

	protected:
		//Do the evaluation
		//virtual double Value(DataPoint*);

		//Do the integration
		virtual double Normalisation( PhaseSpaceBoundary* );
		virtual double Normalisation( DataPoint*, PhaseSpaceBoundary* );

		double cachedIntegral;
		bool cacheValid;
		ParameterSet allParameters;
		vector<string> allObservables;
		bool valid;
		ObservableRef observables;
};

#endif

