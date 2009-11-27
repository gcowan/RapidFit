/**
        @class ResultFormatter

        A collection of static methods for outputting RapidFit data objects

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef RESULT_FORMATTER_H
#define RESULT_FORMATTER_H

#include "IDataSet.h"
#include <string>
#include "FitResult.h"
#include "ToyStudyResult.h"

using namespace std;

class ResultFormatter
{
	public:
		static void MakeRootDataFile( string, IDataSet* );
		static void DebugOutputFitResult( FitResult* );
		static void LatexOutputFitResult( FitResult* );
		static void LatexOutputCovarianceMatrix( FitResult* );
		static void PlotFitContours( FitResult*, string );
		static string FindAndReplaceString( string );
		static double GetElementFromCovarianceMatrix( vector<double>, int, int);		
		static bool IsParameterFree( FitResult*, string );
		static void MakePullPlots( string, ToyStudyResult* );
};

#endif
