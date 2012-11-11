/**
  @class ExternalConstraint

  A class that holds experimentally dervied constraints on fit parameters

  @author Benjamin M Wynne bwynne@cern.ch
  @date 21-01-10
  */

//	RapidFit Headers
#include "ExternalConstraint.h"
///	System Headers
#include <string>
#include <sstream>
#include <iostream>

using namespace::std;

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
string ExternalConstraint::GetName() const
{
	return name;
}
double ExternalConstraint::GetValue() const
{
	return value;
}

double ExternalConstraint::GetError() const
{
	return error;
}

void ExternalConstraint::Print() const
{
	cout << "External Constrint: " << name << "\tValue: " << value << "\tError: " << error << endl;
}

string ExternalConstraint::XML() const
{
	stringstream xml;

	xml << "<ToFit>" << endl;
	xml << "\t<ConstraintFunction>" << endl;
	xml << "\t\t<ExternalConstraint>" << endl;
	xml << "\t\t\t<Name>" << name << "</Name>" << endl;
	xml << "\t\t\t<Value>" << value << "</Value>" << endl;
	xml << "\t\t\t<Error>" << error << "</Error>" << endl;
	xml << "\t\t</ExternalConstraint>" << endl;
	xml << "\t</ConstraintFunction>" << endl;
	xml << "</ToFit>" << endl;

	return xml.str();
}

