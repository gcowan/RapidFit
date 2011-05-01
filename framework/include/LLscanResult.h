/**
        @class LLscanResult

        Container for all results from a LL scan of one parameter

        @author Pete Clarke 
	@date 2010-11-22
*/


#ifndef LLSCAN_RESULT_H
#define LLSCAN_RESULT_H

//	ROOT Headers
#include "TGraph.h"
//	System Headers
#include <string>
#include <vector>

#define LLSCAN_FIT_FAILURE_VALUE -9999.

using namespace std;

class LLscanResult
{
	public:
		LLscanResult();
		LLscanResult( string _parameterName, vector<double> _parameterValues, vector<double> _llvalues  );
		~LLscanResult();
		void print() ;
	
		vector<double> GetParameterValues();
		vector<double> GetLLvalues();
		vector<double> GetRawLLvalues();
	
		TGraph * GetGraph();
	
	private:
	
		string parameterName ;
		double llmin;
		vector<double> parameterValues;
		vector<double> llvalues;
		vector<double> llvalues_offset;
	
};

#endif
