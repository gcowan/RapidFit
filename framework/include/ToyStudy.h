/**
        @class ToyStudy

        A class that automates a whole toy study

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef TOY_STUDY_H
#define TOY_STUDY_H

//	RapidFit Headers
#include "PDFWithData.h"
#include "ParameterSet.h"
#include "ToyStudyResult.h"
#include "FitFunctionConfiguration.h"
#include "MinimiserConfiguration.h"
//	System Headers
#include <string>
#include <vector>

using namespace std;

class ToyStudy
{
	public:
		ToyStudy();
		ToyStudy( string );
		ToyStudy( MinimiserConfiguration*, FitFunctionConfiguration*, vector<ParameterSet*>, vector< PDFWithData* >, vector< ConstraintFunction* >, int );
		~ToyStudy();

		ToyStudyResult * DoWholeStudy( bool=false );
		ToyStudyResult * GetToyStudyResult();

	private:
		//	Uncopyable!
		ToyStudy ( const ToyStudy& );
		ToyStudy& operator = ( const ToyStudy& );

		FitResult * GenerateAndMinimise();

		vector< PDFWithData* > pdfsAndData;
		vector< ParameterSet* > studyParameters;
		MinimiserConfiguration * theMinimiser;
		FitFunctionConfiguration * theFunction;
		ToyStudyResult * allResults;
		int numberStudies;
		vector< ConstraintFunction* > allConstraints;
};

#endif
