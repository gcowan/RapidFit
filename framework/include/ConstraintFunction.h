/**
  @class ConstraintFunction

  Where external, experimental constraints on PhysicsParameters are calculated

  @author Benjamin M Wynne bwynne@cern.ch
  @date 21-01-10
  */

#ifndef CONSTRAINT_FUNCTION_H
#define CONSTRAINT_FUNCTION_H

#include "ParameterSet.h"
#include "ExternalConstraint.h"
#include <vector>

using namespace std;

class ConstraintFunction
{
	public:
		ConstraintFunction();
		ConstraintFunction( vector< ExternalConstraint* > );
		~ConstraintFunction();

		double Evaluate( ParameterSet* );

	private:
		vector< ExternalConstraint* > allConstraints;
};

#endif
