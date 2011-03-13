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
#include "XMLConfigReader.h"
#include <vector>
#include <string>
#include <exception>

using namespace std;

class FitAssembler
{
	public:
		//	New Interface to Fits, for a well behaved PDF it will not crash, it will always return a FitResult which makes it more stable for scanning and such
		//	Also provides way to control the output level during the fit (default = Everything)
		static FitResult * DoSafeFit( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, int=1);

		//	Provided for sWeightPrecalculator I would like to replace this with something closer to another DoSafeFit
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< IPDF* >, vector< IDataSet* >, vector< ConstraintFunction* > );


		//	New Interface to Scanning Code
		static vector<ToyStudyResult*> ContourScan( MinimiserConfiguration *, FitFunctionConfiguration *, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, OutputConfiguration*, string, string, int=-999 );
		static ToyStudyResult* SingleScan(  MinimiserConfiguration *, FitFunctionConfiguration *, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, OutputConfiguration*, string, int=-999 );

		//  Std Feldman-Cousins Code, relying on many already defined objects even though it pulls in the XMLConfigReader object
		static ToyStudyResult* FeldmanCousins( ToyStudyResult*, ToyStudyResult*, vector<unsigned short int>, unsigned short int, bool, OutputConfiguration*, MinimiserConfiguration*, FitFunctionConfiguration*, XMLConfigReader*, vector< PDFWithData* >, int=-999);

	private:
		static FitResult * DoFit( IMinimiser*, FitFunction* );
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, PhysicsBottle* );
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* > );

		static void DoScan( MinimiserConfiguration *, FitFunctionConfiguration *, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, ScanParam*, ToyStudyResult*, int );

		static void DoScan2D( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, pair<ScanParam*, ScanParam* >, vector<ToyStudyResult*>*, int );

//		static void ShakeBottle( ParameterSet*, vector< PDFWithData* >, unsigned int );
};

#endif
