/**
  @class ExternalConstraint

  A class that holds experimentally dervied constraints on fit parameters

  @author Benjamin M Wynne bwynne@cern.ch
  @date 21-01-10
  */

#ifndef EXTERNAL_CONSTRAINT_H
#define EXTERNAL_CONSTRAINT_H

#include <string>

using namespace std;

class ExternalConstraint
{
	public:
		ExternalConstraint();
		ExternalConstraint( string, double, double );
		~ExternalConstraint();

		string GetName();
		double GetValue();
		double GetError();

	private:
		string name;
		double value, error;
};

#endif
