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
#include "FitFunction.h"
#include <string>
#include <vector>

using namespace std;

class ToyStudy
{
	public:
		ToyStudy();
		ToyStudy( string );
		ToyStudy( string, FitFunction*, ParameterSet*, vector< PDFWithData* >, int );
		~ToyStudy();

		ToyStudyResult * DoWholeStudy();
		ToyStudyResult * GetToyStudyResult();

	private:
		FitResult * GenerateAndMinimise();

		vector< PDFWithData* > pdfsAndData;
		ParameterSet * studyParameters;
		string minimiserName;
		FitFunction * theFunction;
		ToyStudyResult * allResults;
		int numberStudies;
};

#endif
