/**
  @class ParameterSet

  A collection of physics parameters

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

//	RapidFit Headers
#include "StringProcessing.h"
#include "ObservableRef.h"
#include "ParameterSet.h"
//	System Headers
#include <iostream>
#include <stdlib.h>
#include <sstream>

using namespace::std;

ParameterSet::ParameterSet( vector<ParameterSet*> input ) : allParameters(), allNames(), trusted_set(), trusted( false )
{
	for( vector<ParameterSet*>::iterator set_i = input.begin(); set_i != input.end(); ++set_i )
	{
		vector<string> input_names = (*set_i)->GetAllNames();

		for( vector<string>::iterator param_i = input_names.begin(); param_i != input_names.end(); ++param_i )
		{
			string temp_par = *param_i;
			if( StringProcessing::VectorContains( &allNames, &temp_par ) == -1 )
			{
				allNames.push_back( temp_par );
				allParameters.push_back( new PhysicsParameter( *((*set_i)->GetPhysicsParameter( temp_par )) ) );
				//this->AddPhysicsParameter( (*set_i)->GetPhysicsParameter( temp_par ) );
			}
			else
			{
				cerr << "Physics Parameter: " << *param_i << " already defined." << endl;
				cerr << "Using last definition of this Parameter." << endl;
				(*set_i)->GetPhysicsParameter( temp_par )->Print();
				//this->AddPhysicsParameter( ((*set_i)->GetPhysicsParameter( temp_par )) );
				allNames.push_back( temp_par );
				allParameters.push_back( new PhysicsParameter( *((*set_i)->GetPhysicsParameter( temp_par )) ) );
			}
		}
	}
}

ParameterSet::ParameterSet( const ParameterSet& input ) : allParameters(), allNames(input.allNames)
{
	for( vector<PhysicsParameter*>::const_iterator param_i = input.allParameters.begin(); param_i != input.allParameters.end(); ++param_i )
	{
		allParameters.push_back( new PhysicsParameter( *(*param_i) ) );
	}
}

ParameterSet& ParameterSet::operator= ( const ParameterSet& input )
{
	if( this != &input )
	{
		while( !allParameters.empty() )
		{
			if( allParameters.back() != NULL ) delete allParameters.back();
			allParameters.pop_back();
		}
		for(vector<PhysicsParameter*>::const_iterator param_i = input.allParameters.begin(); param_i != input.allParameters.end(); ++param_i )
		{
			this->allParameters.push_back( new PhysicsParameter( *(*param_i) ) );
		}
		this->allNames = input.allNames;
	}
	return *this;
}

//Constructor with correct arguments
ParameterSet::ParameterSet( vector<string> NewNames ) : allParameters(), allNames(NewNames)
{
	//Populate the map
	for( unsigned int nameIndex = 0; nameIndex < NewNames.size(); ++nameIndex )
	{
		allParameters.push_back( new PhysicsParameter( NewNames[nameIndex] ) );
	}
}

//Destructor
ParameterSet::~ParameterSet()
{
	while( !trusted_set.empty() )
	{
		if( trusted_set.back() != NULL ) delete trusted_set.back();
		trusted_set.pop_back();
	}
	while( !allParameters.empty() )
	{
		if( allParameters.back() != NULL ) delete allParameters.back();
		allParameters.pop_back();
	}
}

//Retrieve names of all parameters stored
vector<string> ParameterSet::GetAllNames() const
{
	return allNames;
}

//Retrieve names of all parameters stored that are fixed in the PDF
vector<string> ParameterSet::GetAllFixedNames() const
{
	vector<string> Fixed_List;
	for(unsigned short int i=0; i<allNames.size(); ++i )
	{
		if( ParameterSet::GetPhysicsParameter( allNames[i] )->GetType() == "Fixed" )  Fixed_List.push_back( allNames[i] );
	}
	return Fixed_List;
}

//Retrieve names of all parameters stored that are floated in the pdf
vector<string> ParameterSet::GetAllFloatNames() const
{
	vector<string> Not_Fixed_List;
	for(unsigned short int i=0; i<allNames.size(); ++i )
	{
		if( ParameterSet::GetPhysicsParameter( allNames[i] )->GetType() != "Fixed" )  Not_Fixed_List.push_back( allNames[i] );
	}
	return Not_Fixed_List;
}

//Retrieve a physics parameter it's cached index, or find it and save it's index for reference
PhysicsParameter * ParameterSet::GetPhysicsParameter( pair<string,int>* wanted_param ) const
{
	if( wanted_param->second != -1 )
	{
		return allParameters[unsigned(wanted_param->second)];
	} else {
		wanted_param->second = StringProcessing::VectorContains( &allNames, &(wanted_param->first) );
	}
	if( wanted_param->second == -1 ){
		cerr << "PhysicsParameter " << wanted_param->first << " not found" <<endl;
		exit(-3);
	}else{
		return allParameters[unsigned(wanted_param->second)];}
	throw(-20);
}


//Retrieve a physics parameter by its name
PhysicsParameter * ParameterSet::GetPhysicsParameter( const string Name ) const
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if ( nameIndex == -1 )
	{
		cerr << "PhysicsParameter " << Name << " not found(1)" << endl;
		throw(-20);
		//return new PhysicsParameter( Name, 0.0, 0.0, 0.0, "Error", "NameNotFoundError");
	}
	else
	{
		//	This has to be here to ensure that badly constructed parameters don't cause headaches!
		if( allParameters[unsigned(nameIndex)]->GetName().empty() || allParameters[unsigned(nameIndex)]->GetName() == "" )
		{
			allParameters[unsigned(nameIndex)]->SetName( Name );
		}
		return allParameters[unsigned(nameIndex)];
	}
}

PhysicsParameter* ParameterSet::GetPhysicsParameter( const ObservableRef& object ) const
{
	if( object.GetIndex() < 0 )
	{
		object.SetIndex( StringProcessing::VectorContains( &allNames, object.NameRef() ) );
		if( object.GetIndex() >= 0 )
		{
			//      This has to be here to ensure that badly constructed parameters don't cause headaches!
			if( allParameters[ (unsigned) object.GetIndex() ]->GetName().empty() || allParameters[ (unsigned) object.GetIndex() ]->GetName() == "" )
			{
				allParameters[ (unsigned) object.GetIndex() ]->SetName( allNames[(unsigned) object.GetIndex()] );
			}
			return allParameters[ (unsigned) object.GetIndex() ];
		}
	}
	else
	{
		//      This has to be here to ensure that badly constructed parameters don't cause headaches!
		if( allParameters[ (unsigned) object.GetIndex() ]->GetName().empty() || allParameters[ (unsigned) object.GetIndex() ]->GetName() == "" )
		{
			allParameters[ (unsigned) object.GetIndex() ]->SetName( allNames[(unsigned) object.GetIndex()]  );
		}
		return allParameters[ (unsigned) object.GetIndex() ];
	}
	cerr << "PhysicsParameter " << object.Name().c_str() << " not found(2)" << endl;
	throw(-20);
}

//Set a physics parameter by name
bool ParameterSet::SetPhysicsParameter( string Name, PhysicsParameter * NewPhysicsParameter )
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if ( nameIndex == -1 )
	{
		cerr << "PhysicsParameter " << Name << " not found(3)" << endl;
		//exit(1);
		return false;
	}
	else
	{
		//	Copy the new parameter into the old one so no need to delete anything
		allParameters[unsigned(nameIndex)] = NewPhysicsParameter;
		return true;
	}
}

//Initialise a physics parameter
bool ParameterSet::SetPhysicsParameter( string Name, double Value, double Minimum, double Maximum, double Step, string Type, string Unit )
{
	PhysicsParameter * newParameter = new PhysicsParameter( Name, Value, Minimum, Maximum, Step, Type, Unit );
	bool returnValue = SetPhysicsParameter( Name, newParameter );
	delete newParameter;
	return returnValue;
}

//Set all physics parameters
bool ParameterSet::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	if( NewParameterSet->GetTrusted() )
	{
		if( trusted_set.empty() )
		{
			for( unsigned int i=0; i < allNames.size(); ++i )
			{
				ObservableRef* new_ref = new ObservableRef( allNames[i] );
				allParameters[i]->SetValue( NewParameterSet->GetPhysicsParameter( *new_ref )->GetValue() );
				trusted_set.push_back( new_ref );
			}
		}
		else
		{
			for( unsigned int i=0; i< trusted_set.size(); ++i )
			{
				allParameters[i]->SetValue( NewParameterSet->GetPhysicsParameter( *(trusted_set[i]) )->GetValue() );
			}
		}
	}
	else
	{
		for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); nameIndex++)
		{
			PhysicsParameter* inputParameter = new PhysicsParameter( *(NewParameterSet->GetPhysicsParameter( allNames[nameIndex] )) );
			if ( inputParameter->GetUnit() == "NameNotFoundError" )
			{
				//Fail if a required parameter is missing
				cerr << "Parameter \"" << allNames[nameIndex] << "\" expected but not found" << endl;
				return false;
			}
			else
			{
				if( allParameters[nameIndex] != NULL ) delete allParameters[nameIndex];
				allParameters[nameIndex] = inputParameter;
			}
		}
	}

	return true;
}

bool ParameterSet::AddPhysicsParameter( PhysicsParameter* NewParameter, bool replace )
{
	string name = NewParameter->GetName();
	int paramIndex = StringProcessing::VectorContains( &allNames, &name );
	if( paramIndex == -1 )
	{
		allNames.push_back( name );
		allParameters.push_back( new PhysicsParameter( *(NewParameter) ) );
	}
	else
	{
		if( replace )
		{
			if( allParameters[(unsigned)paramIndex] != NULL ) delete allParameters[(unsigned)paramIndex];
			allParameters[(unsigned)paramIndex] = new PhysicsParameter( *(NewParameter) );

		}
	}
	return true;
}

//Set all physics parameters
bool ParameterSet::AddPhysicsParameters( ParameterSet * NewParameterSet, bool replace )
{
	for (unsigned short int nameIndex = 0; nameIndex < NewParameterSet->GetAllNames().size(); nameIndex++)
	{
		//PhysicsParameter * inputParameter = NewParameterSet->GetPhysicsParameter( NewParameterSet->GetAllNames()[nameIndex] );
		int paramIndex = StringProcessing::VectorContains( &allNames, &(NewParameterSet->GetAllNames()[nameIndex]) );
		if( paramIndex == -1 )
		{
			allNames.push_back( NewParameterSet->GetAllNames()[nameIndex] );
			allParameters.push_back( new PhysicsParameter( *(NewParameterSet->GetPhysicsParameter( NewParameterSet->GetAllNames()[nameIndex] )) ) );
		}
		else
		{
			if( replace )
			{
				if( allParameters[(unsigned)paramIndex] != NULL ) delete allParameters[(unsigned)paramIndex];
				allParameters[(unsigned)paramIndex] = new PhysicsParameter( *(NewParameterSet->GetPhysicsParameter( NewParameterSet->GetAllNames()[nameIndex] )) );
			}
		}
	}
	return true;
}

//Not very pleasant in OO terms, and unsafe. Quick however.
bool ParameterSet::SetPhysicsParameters( double * NewValues )
{
	unsigned int temp=0;
	for( vector<PhysicsParameter*>::iterator parameterIndex = allParameters.begin(); parameterIndex != allParameters.end(); ++parameterIndex, ++temp )
	{
		(*parameterIndex)->SetValue( NewValues[temp] );
	}
	return true;
}

double* ParameterSet::GetPhysicsParameters() const
{
	unsigned int temp = 0;
	double* returnable = new double[allParameters.size()];
	for( vector<PhysicsParameter*>::const_iterator parameterIndex = allParameters.begin(); parameterIndex != allParameters.end(); ++parameterIndex, ++temp )
	{
		returnable[temp] = (*parameterIndex)->GetValue();
	}
	return returnable;
}

bool ParameterSet::SetPhysicsParameters( vector<double> NewValues )
{
	if( NewValues.size() == allParameters.size() )
	{
		vector<double>::iterator temp = NewValues.begin();
		for( vector<PhysicsParameter*>::iterator parameterIndex = allParameters.begin(); parameterIndex != allParameters.end(); ++parameterIndex, ++temp )
		{
			(*parameterIndex)->SetValue( *temp );
		}
		return true;
	}
	else
	{
		cerr << "Parameter number mismatch: " << NewValues.size() << " vs " << allParameters.size();
		return false;
	}
}

//General Print method for a dataset
void ParameterSet::Print() const
{
	for(unsigned short int i=0; i< allParameters.size(); ++i)
	{
		cout << "Parameter name: " << allNames[i] <<  endl;
		allParameters[i]->Print();
	}
}

bool ParameterSet::GetTrusted() const
{
	return trusted;
}

void ParameterSet::SetTrusted( bool input )
{
	trusted = input;
}

string ParameterSet::XML() const
{
	stringstream xml;
	xml << "<ParameterSet>" << endl;
	xml << endl;
	for( vector<PhysicsParameter*>::const_iterator parameterIndex = allParameters.begin(); parameterIndex != allParameters.end(); ++parameterIndex )
	{
		string param_str = (*parameterIndex)->XML();
		xml << param_str << endl;
		//xml << endl;
	}
	xml << "</ParameterSet>" << endl;
	return xml.str();
}


