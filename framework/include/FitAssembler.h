/**
        @class FitAssembler

        The intention is for this class to formalise the process of assembling the components of a fit
	Ideally it will be a set of nested static methods, starting from more and more rudimentary components

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef FIT_ASSEMBLER_H
#define FIT_ASSEMBLER_H


#include "FitResult.h"
#include "IMinimiser.h"
#include "FitFunctionConfiguration.h"
#include "MinimiserConfiguration.h"
#include "PDFWithData.h"
#include "ScanParam.h"
#include "ToyStudyResult.h"
#include <vector>
#include <string>
#include <exception>

using namespace std;

class FitAssembler
{
	public:
		static FitResult * DoFit( IMinimiser*, FitFunction* );
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, PhysicsBottle* );
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* > );
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< IPDF* >, vector< IDataSet* >, vector< ConstraintFunction* > );
		static FitResult * DoSafeFit( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* > );

		//  New Interface to Scanning Code
		static vector<ToyStudyResult*> ContourScan( MinimiserConfiguration *, FitFunctionConfiguration *, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, OutputConfiguration*, string, string );
		static ToyStudyResult* SingleScan(  MinimiserConfiguration *, FitFunctionConfiguration *, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, OutputConfiguration*, string );


	private:
		static void DoScan( MinimiserConfiguration *, FitFunctionConfiguration *, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, ScanParam*, ToyStudyResult* );

		static void DoScan2D( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, pair<ScanParam*, ScanParam* >, vector<ToyStudyResult*>* );

		static void ShakeBottle( ParameterSet*, vector< PDFWithData* >, unsigned int );
};

#endif
