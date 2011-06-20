/**
        @class FitResultVector

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

class FitResultVector
{
	public:
		FitResultVector();
		FitResultVector( vector<FitResultVector*> );
		FitResultVector( vector<string> );
		~FitResultVector();

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
		FitResultVector ( const FitResultVector& );
		FitResultVector& operator = ( const FitResultVector& );

		vector< FitResult* > allResults;
		vector<string> allNames;

		//	Was this wise storing a COPY of ALL of this data in the class?
		vector<vector<double> > allValues, allErrors, allPulls, allGenValues;
		vector<double> allRealTimes;
		vector<double> allCPUTimes;
		TStopwatch * clock;
};

#endif
