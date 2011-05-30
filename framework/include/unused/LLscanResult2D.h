/**
        @class LLscanResult2D

        Container for all results from a LL scan of two parameters

        @author Pete Clarke 
	@date 2010-12-22
*/


#ifndef LLSCAN_RESULT2D_H
#define LLSCAN_RESULT2D_H

//	ROOT Headers
#include "TGraph2D.h" 
#include "TH2D.h"
//	RapidFit Headers
#include "LLscanResult.h" 
//	System Headers
#include <string> 
#include <vector> 

using namespace std;

class LLscanResult2D
{
	public:
		LLscanResult2D();
		LLscanResult2D( string _parameterName, vector<double> _parameterValues, string _parameterName2, vector<double> _parameterValues2, vector<LLscanResult*> _llscans  );
		~LLscanResult2D();
		void print() ;
	
		TH2D * GetTH2D();
		
	private:
	
		string parameterName ;
		vector<double> parameterValues;
		string parameterName2 ;
		vector<double> parameterValues2;
		vector<LLscanResult*> listOfLLscans;
};

#endif
