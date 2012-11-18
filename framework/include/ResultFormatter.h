/**
        @class ResultFormatter

        A collection of static methods for outputting RapidFit data objects

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef RESULT_FORMATTER_H
#define RESULT_FORMATTER_H

//	ROOT Headers
#include "TTree.h"
//	RapidFit Headers
#include "IDataSet.h"
#include "FitResult.h"
#include "FitResultVector.h"
//	System Headers
#include <string>
#include <vector>
#include <sstream>

using namespace::std;

class ResultFormatter
{
	public:
		static void MakeRootDataFile( string, vector<IDataSet*> );
		static void DebugOutputFitResult( FitResult* );

		static void LatexDocHeader( stringstream& latex );
		static void LatexDocFooter( stringstream& latex );

		static void TableFooter( stringstream& latex );
		static void TableFooterLandscape( stringstream& latex );
		static void TableHeader( stringstream& latex, unsigned int colnum );
		static void TableHeaderLandscape( stringstream& latex, unsigned int colnum );

		static void AddBranch( TTree* inputTree, const string& BranchName, const vector<int>& IntData );
		static void AddBranch( TTree* inputTree, const string& BranchName, const vector<double>& DoubleData );

		static void WriteOutputLatex( FitResult* OutputData );

		static void LatexSimpleFitResultTable( FitResult * OutputData, stringstream& latex );
		static void LatexFullFitResultTable( FitResult * OutputData, stringstream& latex );
		static void LatexMinimumFitResultTable( FitResult * OutputData, stringstream& latex );
		static void LatexCovMatrix( FitResult * OutputData, stringstream& latex );

		static void LatexOutputFitResult( FitResult* );
		static void ReviewOutput( FitResult * OutputData );
		static void LatexOutputCovarianceMatrix( FitResult* );
		static void PlotFitContours( FitResult*, string );
		static bool IsParameterFree( FitResult*, string );

		//MakePullPlots chooses the appropriate method based on the first string argument
		static void MakePullPlots( string, string, FitResultVector* );
		static void FlatNTuplePullPlots( string, FitResultVector* );
		static void WriteFlatNtuple( string , FitResultVector* );
		static void SeparateParameterPullPlots( string, FitResultVector* );

		static vector<TString> get_branch_names( TTree* );
};

#endif

