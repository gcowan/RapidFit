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
#include "FitResultVector.h"
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
		static FitResult * DoSafeFit( MinimiserConfiguration*, FitFunctionConfiguration*, const vector<ParameterSet*>, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, const int=1);
		static FitResult * DoSafeFit_II( MinimiserConfiguration*, FitFunctionConfiguration*, const vector<ParameterSet*>, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, const int=1);

		//	This moves the 'safety' part of the DoFit into a dedicated function call which is much neater to work with
		static void SafeMinimise( IMinimiser* );


		static FitResult * Petes_DoSafeFit( MinimiserConfiguration*, FitFunctionConfiguration*, const vector<ParameterSet*>, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, const int=1);

		static FitResult * Robs_DoSafeFit( MinimiserConfiguration*, FitFunctionConfiguration*, const vector<ParameterSet*>, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, const int=1);

		static FitResult * DoSingleSafeFit( MinimiserConfiguration*, FitFunctionConfiguration*, const vector<ParameterSet*>, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, const int=1 );

		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, const vector<ParameterSet*>, const vector< IPDF* >, const vector< IDataSet* >,
				const vector< ConstraintFunction* > );


		//	Interface to Classic Scanning Code
		static vector<FitResultVector*> ContourScan( MinimiserConfiguration *, FitFunctionConfiguration *, const vector<ParameterSet*>, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, OutputConfiguration*, const string, const string, const int=-999 );

		static FitResultVector* SingleScan(  MinimiserConfiguration *, FitFunctionConfiguration *, const vector<ParameterSet*>, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, OutputConfiguration*, const string, const int=-999 );



		//  Std Feldman-Cousins Code, relying on many already defined objects even though it pulls in the XMLConfigReader object
		static FitResultVector* FeldmanCousins( FitResultVector*, FitResultVector*, const vector<unsigned int>, const unsigned int, const bool, OutputConfiguration*,
				MinimiserConfiguration*, FitFunctionConfiguration*, XMLConfigReader*, const vector< PDFWithData* >, const int=-999);

	private:
		//	Used for checking a ResultDataSet against the input Parameter to insert 'missing' parameters to avoid bugs downstream in further analysis in RapidFit
		static void CheckParameterSet( FitResult*, vector< ParameterSet* > );


		//	Last DoFit which initialises the minimiser and the FitFunction and does the first fit
		static FitResult * DoFit( IMinimiser*, FitFunction* );

		//	Second DoFit - Used for Constructing a FitFunction for passing to the 
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, PhysicsBottle* );

		//	First DoFit - Used For Constructing a PhysicsBottle to contain all of te 'Physics' input
		static FitResult * DoFit( MinimiserConfiguration*, FitFunctionConfiguration*, vector<ParameterSet*>, const vector< PDFWithData* >,
				const vector< ConstraintFunction* > );

		//	Classic DoScan engine
		static void DoScan( MinimiserConfiguration *, FitFunctionConfiguration *, const vector<ParameterSet*>, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, ScanParam*, FitResultVector*, const int );

		//	Classic DoScan2D engine
		static void DoScan2D( MinimiserConfiguration*, FitFunctionConfiguration*, const vector<ParameterSet*>, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, const pair<ScanParam*, ScanParam* >, vector<FitResultVector*>*, const int );

		//	This has likely been abandonded (it was intended to shake the physics bottle at a failed point in order to recall Minuit
		//	static void ShakeBottle( ParameterSet*, vector< PDFWithData* >, unsigned int );
};

#endif
