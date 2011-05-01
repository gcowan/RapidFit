/**
        @class ProdPDF

        An implementation of IPDF for adding the values of two other IPDFs

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef PROD_PDF_H
#define PROD_PDF_H

//	RapidFit Headers
#include "IPDF.h"

class ProdPDF : public IPDF
{
	public:
		ProdPDF();
		ProdPDF( IPDF*, IPDF* );
		virtual ~ProdPDF();

		void MakePrototypes();

		//Indicate whether the function has been set up correctly
		virtual bool IsValid();

		//Set the function parameters
		virtual bool SetPhysicsParameters(ParameterSet*);

		//Return the integral of the function over the given boundary
		virtual double Integral(DataPoint*, PhaseSpaceBoundary*);

		//Return the function value at the given point
		virtual double Evaluate(DataPoint*);

		//Return the components of the function value at the given point
		virtual vector<double> EvaluateComponents(DataPoint*);
	
		//Return a prototype data point
		virtual vector<string> GetPrototypeDataPoint();

		//Return a prototype set of physics parameters
		virtual vector<string> GetPrototypeParameterSet();

		//Return a list of parameters not to be integrated
                virtual vector<string> GetDoNotIntegrateList();
	
		//Update cache
                virtual void UpdateIntegralCache();
	
	private:
		//	Uncopyable!
		ProdPDF ( const ProdPDF& );
		ProdPDF& operator = ( const ProdPDF& );

		void MakePrototypes( PhaseSpaceBoundary* );
		vector<string> prototypeDataPoint;
		vector<string> prototypeParameterSet;
		vector<string> doNotIntegrateList;
		IPDF * firstPDF;
		IPDF * secondPDF;
};

#endif
