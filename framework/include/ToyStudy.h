/**
        @class ToyStudy

        A class that automates a whole toy study

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef TOY_STUDY_H
#define TOY_STUDY_H

//	RapidFit Headers
#include "IStudy.h"
#include "PDFWithData.h"
#include "ParameterSet.h"
#include "FitResultVector.h"
#include "FitFunctionConfiguration.h"
#include "MinimiserConfiguration.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class ToyStudy	:	public IStudy
{
	public:
		ToyStudy( string );
		ToyStudy( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, vector< ConstraintFunction* >, int );
		~ToyStudy();

		void DoWholeStudy( int = -999 );
		FitResultVector* GetStudyResult();

		void SetNumRepeats( int );			//	Set number of Repeats
		void SetCommandLineParams( vector<string> );	//	Set Command Line Physics Parameters

	private:
		//	Uncopyable!
		ToyStudy ( const ToyStudy& );
		ToyStudy& operator = ( const ToyStudy& );

		FitResult * GenerateAndMinimise();
};

#endif

