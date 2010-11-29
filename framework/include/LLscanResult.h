/**
        @class LLscanResult

        Container for all results from a LL scan of one parameter

        @author Pete Clarke 
	@date 2010-11-22
*/


#ifndef LLSCAN_RESULT_H
#define LLSCAN_RESULT_H

#include <string> 
#include <vector> 
#include "TGraph.h" 

using namespace std;

class LLscanResult
{
	public:
		LLscanResult();
		LLscanResult( string _parameterName, double _centralParameterValue, double _parameterError, double _llmin, vector<double> _parameterValues, vector<double> _llvalues  );
		~LLscanResult();
		void print() ;
	
		vector<double> GetParameterValues();
		vector<double>GetLLvalues();
	
	TGraph * GetGraph();
	
	private:
	
		string parameterName ;
		double centralParameterValue; 
		double parameterError ;
		double llmin;
		vector<double> parameterValues;
		vector<double> llvalues;
	
};

#endif
