/**
        @class RooPDFWrapper

        The parent class for all wrappers for RooFit PDFs, implementing the basic methods.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef ROO_PDF_WRAPPER_H
#define ROO_PDF_WRAPPER_H

#include "IPDF.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"

class RooPDFWrapper : public IPDF
{
	public:
		RooPDFWrapper();
		~RooPDFWrapper();

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

	protected:
		bool SetUpWrapper( vector<string>, int, vector<string>, int );

		RooAbsPdf * wrappedPDF;
		vector< RooRealVar* > variables;
		vector< RooRealVar* > observables;
		bool valid;
		vector<string> prototypeParameterSet;
		vector<string> prototypeDataPoint;
		ParameterSet * pdfParameterSet;
};

#endif
