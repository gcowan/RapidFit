/**
        @class ResultParameterSet

        A set of physics parameters after fitting

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "ResultParameterSet.h"
#include "ParameterSet.h"
//	System Headers
#include <iostream>

//Default constructor
ResultParameterSet::ResultParameterSet() : allParameters(), allNames()
{
}

//Constructor with correct arguments
ResultParameterSet::ResultParameterSet( vector<string> NewNames ) : allParameters(), allNames(NewNames)
{
	//Populate the map
	for (unsigned short int nameIndex = 0; nameIndex < NewNames.size(); ++nameIndex)
	{
		allParameters.push_back( ResultParameter(NewNames[nameIndex],0,-9999.,0.,0.,0.,0.,"type","unit") );
	}
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

ResultParameter * ResultParameterSet::GetResultParameter( int number )
{
	if( number < (int)allNames.size() )
	{
		return &allParameters[ (unsigned)number ];
	} else {
		return new ResultParameter( "DummyResult", 0.0, 0.0, 0.0, 0.0, 0.0, 0., "Error", "NameNotFoundError");
	}
	return NULL;
}

//Retrieve a physics parameter by its name
ResultParameter * ResultParameterSet::GetResultParameter( string Name )
{
	//Check if the name is stored in the map
	for ( unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex)
	{
		if ( allNames[nameIndex] == Name )
		{
			return &allParameters[nameIndex];
		}
	}

	//If the parameter is not found, return an error
	return new ResultParameter( Name, 0.0, 0.0, 0.0, 0.0, 0.0, 0., "Error", "NameNotFoundError");
}

//Set a physics parameter by name
bool ResultParameterSet::SetResultParameter( string Name, ResultParameter * NewResultParameter )
{
	//Check if the name is stored in the map
	for ( unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex)
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
bool ResultParameterSet::SetResultParameter( string Name, double Value, double OriginalValue, double Error, double Minimum, double Maximum, double Step, string Type, string Unit )
{
	ResultParameter * newParameter = new ResultParameter( Name, Value, OriginalValue, Error, Minimum, Maximum, Step, Type, Unit );
	bool returnValue = SetResultParameter( Name, newParameter );
	delete newParameter;
	return returnValue;
}

//Initialise a physics parameter
bool ResultParameterSet::ForceNewResultParameter( string Name, double Value, double OriginalValue, double Error, double Minimum, double Maximum, double Step, string Type, string Unit )
{
	ResultParameter * newParameter = new ResultParameter( Name, Value, OriginalValue, Error, Minimum, Maximum, Step, Type, Unit );
	bool found=false;
	for( unsigned short int i=0; i< allNames.size(); ++i )
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

ParameterSet* ResultParameterSet::GetDummyParameterSet()
{
	ParameterSet* new_Set = new ParameterSet( allNames );
	for( unsigned short int i=0; i < allParameters.size(); ++i )
	{
		new_Set->SetPhysicsParameter( allNames[i], GetResultParameter( allNames[i] )->GetDummyPhysicsParameter() );
	}
	return new_Set;
}
