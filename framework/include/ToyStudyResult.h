/**
        @class ToyStudyResult

        The result of a toy study.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef TOY_STUDY_RESULT_H
#define TOY_STUDY_RESULT_H

#define LLSCAN_FIT_FAILURE_VALUE -9999

//	ROOT Headers
#include "TStopwatch.h"
//	RapidFit Headers
#include "FitResult.h"
//	System Headers
#include <vector>
#include <string>

class ToyStudyResult
{
	public:
		ToyStudyResult();
		ToyStudyResult( vector<ToyStudyResult*> );
		ToyStudyResult( vector<string> );
		~ToyStudyResult();

		bool AddFitResult( FitResult*, bool=true );
		int NumberResults();
		FitResult * GetFitResult(int);
		void StartStopwatch();

		vector<string> GetAllNames();
		vector<double> GetAllRealTimes();
		vector<double> GetAllCPUTimes();
		vector<double> GetAllMLL();
		vector<double> GetParameterValues(string);
		vector<double> GetParameterErrors(string);
		vector<double> GetParameterPulls(string);
		vector<double> GetFlatResult(int);
		TString GetFlatResultHeader();

		double GetCPUTime( int );
		double GetRealTime( int );
		void AddCPUTimes( vector<double> );
		void AddRealTimes( vector<double> );
		void AddCPUTime( double );
		void AddRealTime( double );
		void SetCPUTime( int, double );
		void SetRealTime( int, double );

	private:
		//	Uncopyable!
		ToyStudyResult ( const ToyStudyResult& );
		ToyStudyResult& operator = ( const ToyStudyResult& );

		vector< FitResult* > allResults;
		vector<string> allNames;
		vector< vector<double> > allValues, allErrors, allPulls, allGenValues;
		vector<double> allRealTimes;
		vector<double> allCPUTimes;
		TStopwatch * clock;
};

#endif
