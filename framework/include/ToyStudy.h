/**
        @class ToyStudy

        A class that automates a whole toy study

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef TOY_STUDY_H
#define TOY_STUDY_H

#include "PDFWithData.h"
#include "ParameterSet.h"
#include "ToyStudyResult.h"
#include "FitFunctionConfiguration.h"
#include "MinimiserConfiguration.h"
#include <string>
#include <vector>

using namespace std;

class ToyStudy
{
	public:
		ToyStudy();
		ToyStudy( string );
		ToyStudy( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, vector< PDFWithData* >, int );
		~ToyStudy();

		ToyStudyResult * DoWholeStudy();
		ToyStudyResult * GetToyStudyResult();

	private:
		FitResult * GenerateAndMinimise();

		vector< PDFWithData* > pdfsAndData;
		ParameterSet * studyParameters;
		MinimiserConfiguration * theMinimiser;
		FitFunctionConfiguration * theFunction;
		ToyStudyResult * allResults;
		int numberStudies;
};

#endif
