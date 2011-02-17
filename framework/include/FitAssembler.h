/**
        @class FitAssembler

        The intention is for this class to formalise the process of assembling the components of a fit
	Ideally it will be a set of nested static methods, starting from more and more rudimentary components

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef FIT_ASSEMBLER_H
#define FIT_ASSEMBLER_H

#include "LLscanResult.h"
#include "LLscanResult2D.h"
#include "FitResult.h"
#include "IMinimiser.h"
#include "FitFunctionConfiguration.h"
#include "MinimiserConfiguration.h"
#include "PDFWithData.h"
#include "ScanParam.h"
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
		static vector<LLscanResult*> DoCVScan( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, OutputConfiguration*, string );
		static vector<FitResult*> DoScan( MinimiserConfiguration *, FitFunctionConfiguration *, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, ScanParam* );
		static LLscanResult* DoLLScan( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, ScanParam* );
		static LLscanResult* DoLLScan( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, OutputConfiguration*, string );
		static LLscanResult2D * DoLLScan2D( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, OutputConfiguration*, string, string );
		static vector<vector<FitResult*> > DoScan2D( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, OutputConfiguration*, string, string, string="LLscan", string="LLscan" );
		static vector<FitResult* > Linearize( vector<vector<FitResult*> > Input_Results );
};

#endif
