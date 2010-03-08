/**
        @class ProdPDF

        An implementation of IPDF for adding the values of two other IPDFs

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef PROD_PDF_H
#define PROD_PDF_H

#include "IPDF.h"

class ProdPDF : public IPDF
{
	public:
		ProdPDF();
		ProdPDF( IPDF*, IPDF* );
		~ProdPDF();

		void MakePrototypes();

		//Indicate whether the function has been set up correctly
		virtual bool IsValid();

		//Set the function parameters
		virtual bool SetPhysicsParameters(ParameterSet*);

		//Return the integral of the function over the given boundary
		virtual double Integral(DataPoint*, PhaseSpaceBoundary*);

		//Return the function value at the given point
		virtual double Evaluate(DataPoint*);

		//Return a prototype data point
		virtual vector<string> GetPrototypeDataPoint();

		//Return a prototype set of physics parameters
		virtual vector<string> GetPrototypeParameterSet();

		//Return a list of parameters not to be integrated
                virtual vector<string> GetDoNotIntegrateList();
	
		//Update cache
                virtual void UpdateIntegralCache();
	
	private:
		vector<string> prototypeDataPoint;
		vector<string> prototypeParameterSet;
		vector<string> doNotIntegrateList;
		IPDF * firstPDF;
		IPDF * secondPDF;
};

#endif
