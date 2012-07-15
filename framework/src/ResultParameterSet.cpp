/**
  @class ResultParameterSet

  A set of physics parameters after fitting

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

///	RapidFit Headers
#include "ResultParameterSet.h"
#include "ParameterSet.h"
#include "StringProcessing.h"
#include "RapidFitMatrix.h"
///	System Headers
#include <iostream>
#include <sstream>
#include <cmath>

using namespace::std;

//Constructor with correct arguments
ResultParameterSet::ResultParameterSet( vector<string> NewNames ) : allParameters(), allNames()
{
	//Populate the map
	for (unsigned short int nameIndex = 0; nameIndex < NewNames.size(); ++nameIndex)
	{
		allParameters.push_back( new ResultParameter( NewNames[nameIndex], 0., -9999., 0., 0., 0., "type", "unit" ) );
	}

	vector<string> duplicates;
	allNames = StringProcessing::RemoveDuplicates( NewNames, duplicates );
	if( allNames.size() != NewNames.size() )
	{
		cerr << "WARNING: Cannot Generate a ResultParameterSet with 2 Occurances of the same name" << endl;
		for( vector<string>::iterator str_i = duplicates.begin(); str_i != duplicates.end(); ++str_i )
		{
			cout << *str_i << endl;
		}
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

//Retrieve names of all parameters stored that are fixed in the PDF
vector<string> ResultParameterSet::GetAllFixedNames() const
{
	vector<string> Fixed_List;
	for(unsigned short int i=0; i<allNames.size(); ++i )
	{
		if( this->GetResultParameter( allNames[i] )->GetType() == "Fixed" )  Fixed_List.push_back( allNames[i] );
	}
	return Fixed_List;
}

//Retrieve names of all parameters stored that are floated in the pdf
vector<string> ResultParameterSet::GetAllFloatNames() const
{
	vector<string> Not_Fixed_List;
	for(unsigned short int i=0; i<allNames.size(); ++i )
	{
		if( this->GetResultParameter( allNames[i] )->GetType() != "Fixed" )  Not_Fixed_List.push_back( allNames[i] );
	}
	return Not_Fixed_List;
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

ResultParameter* ResultParameterSet::GetResultParameter( const ObservableRef& object ) const
{
	if( object.GetIndex() < 0 )
	{
		object.SetIndex( StringProcessing::VectorContains( &allNames, object.NameRef() ) );
		if( object.GetIndex() >= 0 )
		{
			return allParameters[ (unsigned) object.GetIndex() ];
		}
	}
	else
	{
		//      This has to be here to ensure that badly constructed parameters don't cause headaches!
		return allParameters[ (unsigned) object.GetIndex() ];
	}
	cerr << "ResultParameter " << object.Name().c_str() << " not found(2)" << endl;
	return NULL;
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
			allParameters[nameIndex] = new ResultParameter( *NewResultParameter );
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
bool ResultParameterSet::ForceNewResultParameter( ResultParameter* Input )
{
	string Name = Input->GetName();
	bool found=false;
	for( unsigned short int i=0; i< allNames.size(); ++i )
	{
		if( allNames[i] == Name )  found = true;
	}
	if( !found )
	{
		allNames.push_back( Name );
		allParameters.push_back( new ResultParameter( *Input ) );
	}
	return !found;
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

void ResultParameterSet::ApplyCovarianceMatrix( RapidFitMatrix* Input )
{
	vector<string> appliedParmaters = Input->theseParameters;
	cout << "Applying ResultParameterSet:  ";
	int i=0;
	for( vector<string>::iterator name_i = appliedParmaters.begin(); name_i != appliedParmaters.end(); ++name_i, ++i )
	{
		cout << *name_i << "\t";
		ResultParameter* thisParam = this->GetResultParameter( *name_i );
		double thisErrorSq = (*(Input->thisMatrix))(i,i);
		cout << i << "\t" << thisErrorSq << "\tsqrt=" << sqrt(thisErrorSq) << endl;
		thisParam->SetError( sqrt( thisErrorSq ) );
	}
	cout << endl;
	return;
}


