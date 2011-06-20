/**
        @class ResultParameterSet

        A set of physics parameters after fitting

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef RESULT_PARAMETER_SET_H
#define RESULT_PARAMETER_SET_H

//	RapidFit Headers
#include "ResultParameter.h"
#include "ParameterSet.h"
//	System Headers
#include <vector>
#include <string>

using namespace std;

class ResultParameterSet
{
	public:
		ResultParameterSet();
		ResultParameterSet( vector<string> );
		~ResultParameterSet();

		vector<string> GetAllNames();
		ResultParameter * GetResultParameter( int );
		ResultParameter * GetResultParameter(string);
		bool SetResultParameter( string, ResultParameter* );
		bool SetResultParameter( string, double, double, double, double, double, string, string );
		bool ForceNewResultParameter( string, double, double, double, double, double, string, string );
		ParameterSet* GetDummyParameterSet();

	private:
		vector<ResultParameter> allParameters;
		vector<string> allNames;
};

#endif
