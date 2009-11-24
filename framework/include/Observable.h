/**
        @class Observable

        The measured value of a particular variable for an event

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include <string>

using namespace std;

class Observable
{
	public:
		Observable();
		Observable( string, double, double, string );
		~Observable();

		double GetValue();
		void SetValue(double);

		double GetError();
		void SetError(double);

		string GetUnit();

	private:
		double value;
		double error;
		string unit;
};

#endif
