/**
        @class ToyStudyResult

        The result of a toy study.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef TOY_STUDY_RESULT_H
#define TOY_STUDY_RESULT_H

#include "FitResult.h"
#include <vector>
#include <string>
#include "TStopwatch.h"

class ToyStudyResult
{
	public:
		ToyStudyResult();
		ToyStudyResult( vector<string> );
		~ToyStudyResult();

		bool AddFitResult( FitResult* );
		int NumberResults();
		FitResult * GetFitResult(int);
		void StartStopwatch();

		vector<string> GetAllNames();
		vector<double> GetAllRealTimes();
		vector<double> GetAllCPUTimes();
		vector<double> GetParameterValues(string);
		vector<double> GetParameterErrors(string);
		vector<double> GetParameterPulls(string);
		vector<double> GetFlatResult(int);
		TString GetFlatResultHeader();

	private:
		vector< FitResult* > allResults;
		vector<string> allNames;
		vector< vector<double> > allValues, allErrors, allPulls;
		vector<double> allRealTimes;
		vector<double> allCPUTimes;
		TStopwatch * clock;
};

#endif
