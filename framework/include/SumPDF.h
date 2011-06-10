/**
        @class SumPDF

        An implementation of IPDF for adding the values of two other IPDFs

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef SUM_PDF_H
#define SUM_PDF_H

//	RapidFit Headers
//#include "IPDF.h"
#include "BasePDF.h"

class SumPDF : public BasePDF
{
	public:
		SumPDF();
		SumPDF( IPDF*, IPDF*, PhaseSpaceBoundary* );
		SumPDF( IPDF*, IPDF*, PhaseSpaceBoundary*, string );
		~SumPDF();

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
		SumPDF ( const SumPDF& );
		SumPDF& operator = ( const SumPDF& );
		void MakePrototypes( PhaseSpaceBoundary* );

		vector<string> prototypeDataPoint, prototypeParameterSet, doNotIntegrateList;
		IPDF * firstPDF;
	       	IPDF * secondPDF;
		double firstFraction, firstIntegralCorrection, secondIntegralCorrection;
		string fractionName;
};

#endif
