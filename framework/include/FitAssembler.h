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
#include <vector>
#include <string>

using namespace std;

class FitAssembler
{
	public:
		static FitResult * DoFit( IMinimiser*, FitFunction* );
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, PhysicsBottle* );
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* > );
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< IPDF* >, vector< IDataSet* >, vector< ConstraintFunction* > );
		static LLscanResult * DoScan( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, string, int, double, double );
		static LLscanResult2D * DoScan2D( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, string, string, int, double, double, double, double );

};

#endif
