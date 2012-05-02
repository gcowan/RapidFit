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
#include <sstream>
using namespace::std;

//Constructor with correct arguments
ResultParameterSet::ResultParameterSet( vector<string> NewNames ) : allParameters(), allNames(NewNames)
{
	//Populate the map
	for (unsigned short int nameIndex = 0; nameIndex < NewNames.size(); ++nameIndex)
	{
		allParameters.push_back( new ResultParameter( NewNames[nameIndex], 0., -9999., 0., 0., 0., "type", "unit" ) );
	}
}

ResultParameterSet::ResultParameterSet( const ResultParameterSet& input ) : allParameters(), allNames( input.allNames )
{
	for( vector<ResultParameter*>::const_iterator par_i = input.allParameters.begin(); par_i != input.allParameters.end(); ++par_i )
	{
		allParameters.push_back( new ResultParameter(**par_i) );
	}
}

//Destructor
ResultParameterSet::~ResultParameterSet()
{
	while( !allParameters.empty() )
	{
		if( allParameters.back() != NULL ) delete allParameters.back();
		allParameters.pop_back();
	}
}

//Retrieve names of all parameters stored
vector<string> ResultParameterSet::GetAllNames() const
{
	return allNames;
}

ResultParameter * ResultParameterSet::GetResultParameter( int number ) const
{
	if( number < (int)allNames.size() )
	{
		return allParameters[ (unsigned)number ];
	} else {
		return new ResultParameter( "DummyResult", 0.0, 0.0, 0.0, 0.0, 0.0, "Error", "NameNotFoundError" );
	}
	return NULL;
}

//Retrieve a physics parameter by its name
ResultParameter * ResultParameterSet::GetResultParameter( string Name ) const
{
	//Check if the name is stored in the map
	for ( unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex)
	{
		if ( allNames[nameIndex] == Name )
		{
			return allParameters[nameIndex];
		}
	}

	//If the parameter is not found, return an error
	return new ResultParameter( Name, 0.0, 0.0, 0.0, 0.0, 0.0, "Error", "NameNotFoundError");
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
			if( allParameters[nameIndex] != NULL ) delete allParameters[nameIndex];
			allParameters[nameIndex] = NewResultParameter;
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
	//delete newParameter;
	return returnValue;
}

//Initialise a physics parameter
bool ResultParameterSet::ForceNewResultParameter( string Name, double Value, double OriginalValue, double Error, double Minimum, double Maximum, string Type, string Unit )
{
	ResultParameter * newParameter = new ResultParameter( Name, Value, OriginalValue, Error, Minimum, Maximum, Type, Unit );
	bool found=false;
	for( unsigned short int i=0; i< allNames.size(); ++i )
	{
		if( allNames[i] == Name )  found = true;
	}
	if( !found )
	{
		allNames.push_back( Name );
		allParameters.push_back( newParameter );
	}
	//delete newParameter;
	return !found;
}

ParameterSet* ResultParameterSet::GetDummyParameterSet() const
{
	ParameterSet* new_Set = new ParameterSet( allNames );
	for( unsigned short int i=0; i < allParameters.size(); ++i )
	{
		new_Set->SetPhysicsParameter( allNames[i], GetResultParameter( allNames[i] )->GetDummyPhysicsParameter() );
	}
	return new_Set;
}

void ResultParameterSet::Print() const
{
	cout << "ResultParameterSet:" << endl;
	for( vector<ResultParameter*>::const_iterator par_i = allParameters.begin(); par_i != allParameters.end(); ++par_i )
	{
		(*par_i)->Print();
	}
}

string ResultParameterSet::XML( const bool fit ) const
{
	stringstream xml;
	xml << "<ParameterSet>" << endl;
	for( vector<ResultParameter*>::const_iterator par_i = allParameters.begin(); par_i != allParameters.end(); ++par_i )
	{
		if( fit == true )
		{
			xml << (*par_i)->FitXML() << endl;
		}
		else
		{
			xml << (*par_i)->ToyXML() << endl;
		}
	}
	xml << "</ParameterSet>" << endl;
	return xml.str();
}

string ResultParameterSet::FitXML() const
{
	return this->XML( true );
}

string ResultParameterSet::ToyXML() const
{
	return this->XML( false );
}

