/**
        @class ParameterSet

        A collection of physics parameters

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "StringProcessing.h"
#include "ParameterSet.h"
#include <iostream>

//Default constructor
ParameterSet::ParameterSet()
{
}

//Constructor with correct arguments
ParameterSet::ParameterSet( vector<string> NewNames )
{
	//Populate the map
	for (int nameIndex = 0; nameIndex < NewNames.size(); nameIndex++)
	{
		allParameters.push_back( PhysicsParameter() );
	}

	allNames = NewNames;
}

//Destructor
ParameterSet::~ParameterSet()
{
}

//Retrieve names of all parameters stored
vector<string> ParameterSet::GetAllNames()
{
	return allNames;
}

//Retrieve a physics parameter by its name
PhysicsParameter * ParameterSet::GetPhysicsParameter(string Name)
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if ( nameIndex == -1 )
	{
		return new PhysicsParameter( Name, 0.0, 0.0, 0.0, "Error", "NameNotFoundError");
	}
	else
	{
		return &allParameters[nameIndex];
	}
}

//Set a physics parameter by name
bool ParameterSet::SetPhysicsParameter( string Name, PhysicsParameter * NewPhysicsParameter )
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
        if ( nameIndex == -1 )
	{
		return false;
	}
	else
	{
		allParameters[nameIndex] = *NewPhysicsParameter;
                return true;
	}
}

//Initialise a physics parameter
bool ParameterSet::SetPhysicsParameter( string Name, double Value, double Minimum, double Maximum, string Type, string Unit )
{
	PhysicsParameter * newParameter = new PhysicsParameter( Name, Value, Minimum, Maximum, Type, Unit );
	bool returnValue = SetPhysicsParameter( Name, newParameter );
	delete newParameter;
	return returnValue;
}

//Set all physics parameters
bool ParameterSet::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	for ( int nameIndex = 0; nameIndex < allNames.size(); nameIndex++)
	{
		PhysicsParameter * inputParameter = NewParameterSet->GetPhysicsParameter( allNames[nameIndex] );
		if ( inputParameter->GetUnit() == "NameNotFoundError" )
		{
			//Fail if a required parameter is missing
			cerr << "Parameter \"" << allNames[nameIndex] << "\" expected but not found" << endl;
			return false;
		}
		else
		{
			allParameters[nameIndex] = *inputParameter;
		}
	}

	return true;
}

//Not very pleasant in OO terms, and unsafe. Quick however.
bool ParameterSet::SetPhysicsParameters( double * NewValues )
{
	for ( int parameterIndex = 0; parameterIndex < allParameters.size(); parameterIndex++ )
	{
		allParameters[parameterIndex].SetValue( NewValues[parameterIndex] );
	}
	return true;
}
bool ParameterSet::SetPhysicsParameters( vector<double> NewValues )
{
	if ( NewValues.size() == allParameters.size() )
	{
		for ( int parameterIndex = 0; parameterIndex < allParameters.size(); parameterIndex++ )
		{
			allParameters[parameterIndex].SetValue( NewValues[parameterIndex] );
		}
		return true;
	}
	else
	{
		cerr << "Parameter number mismatch: " << NewValues.size() << " vs " << allParameters.size();
		return false;
	}
}
