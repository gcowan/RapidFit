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
#include "FitFunction.h"
#include "PDFWithData.h"
#include <vector>

using namespace std;

class FitAssembler
{
	public:
		static FitResult * DoFit( IMinimiser*, FitFunction* );
		static FitResult * DoFit( string, FitFunction*, PhysicsBottle* );
		static FitResult * DoFit( string, FitFunction*, ParameterSet*, vector< PDFWithData* > );
};

#endif
