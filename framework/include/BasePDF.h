/**
        @class BasePDF

        Class that provides a general implementation of IPDF.
        Can inherit from this to make a PDF without worrying about the details.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef BASE_PDF_H
#define BASE_PDF_H

#include "IPDF.h"

class BasePDF : public IPDF
{
	public:
		BasePDF();
		~BasePDF();

		//Indicate whether the function has been set up correctly
		virtual bool IsValid();

		//Set the function parameters
		virtual bool SetPhysicsParameters( ParameterSet* );

		//Return the integral of the function over the given boundary
		virtual double Integral( DataPoint*, PhaseSpaceBoundary* );

		//Return the function value at the given point
		virtual double Evaluate( DataPoint* );

		//Return a prototype data point
		virtual vector<string> GetPrototypeDataPoint();

		//Return a prototype set of physics parameters
		virtual vector<string> GetPrototypeParameterSet();

                //Return a list of parameters not to be integrated
                virtual vector<string> GetDoNotIntegrateList();

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
};

#endif
