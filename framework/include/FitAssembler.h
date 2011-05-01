/**
        @class FitAssembler

        The intention is for this class to formalise the process of assembling the components of a fit
	Ideally it will be a set of nested static methods, starting from more and more rudimentary components

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef FIT_ASSEMBLER_H
#define FIT_ASSEMBLER_H

//	RapidFit Headers
#include "FitResult.h"
#include "IMinimiser.h"
#include "FitFunctionConfiguration.h"
#include "MinimiserConfiguration.h"
#include "PDFWithData.h"
#include "ScanParam.h"
#include "ToyStudyResult.h"
#include "XMLConfigReader.h"
//	System Headers
#include <vector>
#include <string>
#include <exception>

using namespace std;

class FitAssembler
{
	public:
		//	New Interface to Fits, for a well behaved PDF it will not crash, it will always return a FitResult which makes it more stable for scanning and such
		//	Also provides way to control the output level during the fit (default = Everything)
		static FitResult * DoSafeFit( MinimiserConfiguration*, FitFunctionConfiguration*, const vector<ParameterSet*>, const vector< PDFWithData* >, const vector< ConstraintFunction* >, const int=1);

		//	Provided for sWeightPrecalculator I would like to replace this with something closer to another DoSafeFit
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, const vector<ParameterSet*>, const vector< IPDF* >, const vector< IDataSet* >, const vector< ConstraintFunction* > );


		//	New Interface to Scanning Code
		static vector<ToyStudyResult*> ContourScan( MinimiserConfiguration *, FitFunctionConfiguration *, const vector<ParameterSet*>, const vector< PDFWithData* >, const vector< ConstraintFunction* >, OutputConfiguration*, const string, const string, const int=-999 );
		static ToyStudyResult* SingleScan(  MinimiserConfiguration *, FitFunctionConfiguration *, const vector<ParameterSet*>, const vector< PDFWithData* >, const vector< ConstraintFunction* >, OutputConfiguration*, const string, const int=-999 );

		//  Std Feldman-Cousins Code, relying on many already defined objects even though it pulls in the XMLConfigReader object
		static ToyStudyResult* FeldmanCousins( ToyStudyResult*, ToyStudyResult*, const vector<unsigned int>, const unsigned int, const bool, OutputConfiguration*, MinimiserConfiguration*, FitFunctionConfiguration*, XMLConfigReader*, const vector< PDFWithData* >, const int=-999);

	private:
		static FitResult * DoFit( IMinimiser*, FitFunction* );
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, PhysicsBottle* );
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, vector<ParameterSet*>, const vector< PDFWithData* >, const vector< ConstraintFunction* > );

		static void DoScan( MinimiserConfiguration *, FitFunctionConfiguration *, const vector<ParameterSet*>, const vector< PDFWithData* >, const vector< ConstraintFunction* >, ScanParam*, ToyStudyResult*, const int );

		static void DoScan2D( MinimiserConfiguration*, FitFunctionConfiguration*, const vector<ParameterSet*>, const vector< PDFWithData* >, const vector< ConstraintFunction* >, const pair<ScanParam*, ScanParam* >, vector<ToyStudyResult*>*, const int );

//		static void ShakeBottle( ParameterSet*, vector< PDFWithData* >, unsigned int );
};

#endif
