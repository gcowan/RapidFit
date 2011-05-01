/**
        @class LLscanList

        Information needed to perform an LLscan.

        @author Pete Clarke
	@date 2009-10-02
*/

#ifndef RESULT_PARAMETER_H
#define RESULT_PARAMETER_H

//	System Headers
#include <string>

using namespace std;

class LLscanList
{
	public:
		LLscanList();
		LLscanList( string name, double lo, double hi, int npoints );
		~LLscanList();

		double Name();
		double LoLim();
		double HiLim();
		double Npoints();

	private:
		double value;
		double originalValue;
		double error;
		double minimum;
		double maximum;
		string type;
		string unit;
};

#endif
