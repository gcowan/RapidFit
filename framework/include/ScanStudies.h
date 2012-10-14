
#pragma once
#ifndef ScanStudies_H
#define ScanStudies_H

//	RapidFit Headers
#include "FitResultVector.h"
#include "MinimiserConfiguration.h"
#include "ConstraintFunction.h"
#include "OutputConfiguration.h"
#include "ParameterSet.h"
#include "PDFWithData.h"
#include "FitFunctionConfiguration.h"
#include "ResultFormatter.h"
#include "DebugClass.h"
//	System Headers
#include <vector>
#include <string>

using namespace::std;

class ScanStudies
{

	public:

		//	New Interface to Scanning Code
		static vector<FitResultVector*> ContourScan( MinimiserConfiguration *, FitFunctionConfiguration *, ParameterSet*, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, OutputConfiguration*, const string, const string, const int=-999, DebugClass* debug=NULL, bool forceContinue=false );
		static FitResultVector* SingleScan(  MinimiserConfiguration *, FitFunctionConfiguration *, ParameterSet*, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, OutputConfiguration*, const string, const int=-999, DebugClass* debug=NULL, bool forceContinue=false );

	private:

		static void DoScan( MinimiserConfiguration *, FitFunctionConfiguration *, ParameterSet*, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, ScanParam*, FitResultVector*, const int, DebugClass* debug, bool forceContinue=false );

		static void DoScan2D( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, const pair<ScanParam*, ScanParam* >, vector<FitResultVector*>*, const int, DebugClass* debug, bool forceContinue=false );

};

#endif

