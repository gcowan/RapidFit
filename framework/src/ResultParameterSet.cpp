/**
        @class ResultParameterSet

        A set of physics parameters after fitting

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "ResultParameterSet.h"
#include <iostream>

//Default constructor
ResultParameterSet::ResultParameterSet()
{
}

//Constructor with correct arguments
ResultParameterSet::ResultParameterSet( vector<string> NewNames )
{
	//Populate the map
	for (int nameIndex = 0; nameIndex < NewNames.size(); nameIndex++)
	{
		allParameters.push_back( ResultParameter() );
	}

	allNames = NewNames;
}

//Destructor
ResultParameterSet::~ResultParameterSet()
{
}

//Retrieve names of all parameters stored
vector<string> ResultParameterSet::GetAllNames()
{
	return allNames;
}

//Retrieve a physics parameter by its name
ResultParameter * ResultParameterSet::GetResultParameter(string Name)
{
	//Check if the name is stored in the map
	for ( int nameIndex = 0; nameIndex < allNames.size(); nameIndex++)
	{
		if ( allNames[nameIndex] == Name )
		{
			return &allParameters[nameIndex];
		}
	}

	//If the parameter is not found, return an error
	return new ResultParameter( Name, 0.0, 0.0, 0.0, 0.0, 0.0, "Error", "NameNotFoundError");
}

//Set a physics parameter by name
bool ResultParameterSet::SetResultParameter( string Name, ResultParameter * NewResultParameter )
{
	//Check if the name is stored in the map
	for ( int nameIndex = 0; nameIndex < allNames.size(); nameIndex++)
	{
		if ( allNames[nameIndex] == Name )
		{
			//Delete old parameter before overwriting the pointer
			//delete allParameters[nameIndex];
			allParameters[nameIndex] = *NewResultParameter;
			return true;
		}
	}

	//If the parameter is not found, return false
	return false;
}

//Initialise a physics parameter
bool ResultParameterSet::SetResultParameter( string Name, double Value, double OriginalValue, double Error, double Minimum, double Maximum, string Type, string Unit )
{
	ResultParameter * newParameter = new ResultParameter( Name, Value, OriginalValue, Error, Minimum, Maximum, Type, Unit );
	bool returnValue = SetResultParameter( Name, newParameter );
	delete newParameter;
	return returnValue;
}

//Initialise a physics parameter
bool ResultParameterSet::ForceNewResultParameter( string Name, double Value, double OriginalValue, double Error, double Minimum, double Maximum, string Type, string Unit )
{
	ResultParameter * newParameter = new ResultParameter( Name, Value, OriginalValue, Error, Minimum, Maximum, Type, Unit );
	bool found=false;
	for( short int i=0; i< allNames.size(); i++ )
	{
		if( allNames[i] == Name )  found = true;
	}
	if( !found )
	{
		allNames.push_back( Name );
		allParameters.push_back( *newParameter );
	}
	delete newParameter;
	return !found;
}