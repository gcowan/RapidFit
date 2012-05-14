/**
        @class StatisticsFunctions

        A collection of static methods for statistics calculations

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef STATISTICS_FUNCTIONS_H
#define STATISTICS_FUNCTIONS_H

//	RapidFit Headers
#include "PhaseSpaceBoundary.h"
#include "IPDF.h"
#include "IDataSet.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class StatisticsFunctions
{
	public:
		static double Mean( vector<double> );
		static double Variance( vector<double> );
		static int OptimumBinNumber( vector<double> );
		static double Maximum( vector<double> );
		static double Minimum( vector<double> );
		static vector< vector<double> > DiscreteCombinations( vector<string>*, const PhaseSpaceBoundary*, vector<string>&, vector<string>&, vector< vector<double> >& );
		static void DoDontIntegrateLists( IPDF*, const PhaseSpaceBoundary*, const vector<string>*, vector<string>&, vector<string>& );
		static vector<DataPoint*> DataAverage( IDataSet*, vector< vector<double> >, vector< vector<double> >, vector<string>, vector<string>, vector<string>&, vector<double>& );
};

#endif

