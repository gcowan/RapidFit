/**
  @class ConstraintFunction

  Where external, experimental constraints on PhysicsParameters are calculated

  @author Benjamin M Wynne bwynne@cern.ch
  @date 21-01-10
  */

#include "ConstraintFunction.h"
#include "StringProcessing.h"
#include <iostream>

using namespace std;

//Default constructor
ConstraintFunction::ConstraintFunction()
{
}

//Constructor with correct arguments
ConstraintFunction::ConstraintFunction( vector< ExternalConstraint* > NewConstraints) : allConstraints(NewConstraints)
{
}

//Destructor
ConstraintFunction::~ConstraintFunction()
{
}

//Perform the constraint calculation
double ConstraintFunction::Evaluate( ParameterSet * NewParameters )
{
	vector<string> parameterNames = NewParameters->GetAllNames();
	double constraintValue = 0.0;

	//Loop over all ExternalConstraints
	for ( int constraintIndex = 0; constraintIndex < allConstraints.size(); constraintIndex++ )
	{
		string name = allConstraints[constraintIndex]->GetName();
		if ( name == "GammaL" )
		{
			//Insert code here
		}
		else if ( StringProcessing::VectorContains( &parameterNames, &name ) >= 0 )
		{
			//Do standard gaussian constraint calculation
			double parameterValue = NewParameters->GetPhysicsParameter(name)->GetValue();
			double gaussSqrt = ( parameterValue - allConstraints[constraintIndex]->GetValue() ) / allConstraints[constraintIndex]->GetError();
			constraintValue += gaussSqrt * gaussSqrt;
		}
	}

	return constraintValue;
}
