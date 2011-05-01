/**
  @class ExternalConstraint

  A class that holds experimentally dervied constraints on fit parameters

  @author Benjamin M Wynne bwynne@cern.ch
  @date 21-01-10
  */

//	RapidFit Headers
#include "ExternalConstraint.h"

//Default constructor
ExternalConstraint::ExternalConstraint(): name(), value(), error()
{
}

//Constructor with correct arguments
ExternalConstraint::ExternalConstraint( string NewName, double NewValue, double NewError ) : name(NewName), value(NewValue), error(NewError)
{
}

//Destructor
ExternalConstraint::~ExternalConstraint()
{
}

//Return the information held
string ExternalConstraint::GetName()
{
	return name;
}
double ExternalConstraint::GetValue()
{
	return value;
}
double ExternalConstraint::GetError()
{
	return error;
}
