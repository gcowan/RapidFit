/**
        @class ResultParameter

        A physics parameter after it has been fitted

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef RESULT_PARAMETER_H
#define RESULT_PARAMETER_H

//	RapidFit Headers
#include "PhysicsParameter.h"
//	System Headers
#include <string>

using namespace std;

class ResultParameter
{
	public:
		ResultParameter();
		ResultParameter( string, double, double, double, double, double, string, string );
		~ResultParameter();

		double GetValue();
		double GetOriginalValue();
		double GetError();
		double GetPull();
		double GetMinimum();
		double GetMaximum();
		string GetType();
		string GetUnit();
		void ForceOriginalValue( double );
		void ForcePullValue( double );
		PhysicsParameter* GetDummyPhysicsParameter();

	private:
		string name;
		double value;
		double originalValue;
		double error;
		double minimum;
		double maximum;
		string type;
		string unit;
};

#endif
