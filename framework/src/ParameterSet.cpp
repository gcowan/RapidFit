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
#include <cstdlib>
#include <sstream>

using namespace::std;

ParameterSet::ParameterSet( vector<ParameterSet*> input ) :
	allParameters(), allNames(), uniqueID(0), allInternalNames(), allForeignNames()
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
	uniqueID = reinterpret_cast<size_t>(this);
	for( unsigned int i=0; i< allNames.size(); ++i )
	{
		allInternalNames.push_back( ObservableRef(allNames[i]) );
		allForeignNames.push_back( ObservableRef(allNames[i]) );
	}
}

ParameterSet::ParameterSet( const ParameterSet& input ) :
	allParameters(), allNames(input.allNames), uniqueID(0), allInternalNames(input.allInternalNames), allForeignNames(input.allForeignNames)
{
	vector<PhysicsParameter*>::const_iterator param_i = input.allParameters.begin();
	for( ; param_i != input.allParameters.end(); ++param_i )
	{
		if( *param_i != NULL ) allParameters.push_back( new PhysicsParameter( *(*param_i) ) );
		else allParameters.push_back( NULL );
	}

	uniqueID = reinterpret_cast<size_t>(this)+1;
}

ParameterSet& ParameterSet::operator= ( const ParameterSet& input )
{
	if( this != &input )
	{
		while( !this->allParameters.empty() )
		{
			if( this->allParameters.back() != NULL ) delete this->allParameters.back();
			this->allParameters.pop_back();
		}
		for(vector<PhysicsParameter*>::const_iterator param_i = input.allParameters.begin(); param_i != input.allParameters.end(); ++param_i )
		{
			this->allParameters.push_back( new PhysicsParameter( *(*param_i) ) );
		}
		this->allNames = input.allNames;
		this->allInternalNames = input.allInternalNames;
		this->allForeignNames = input.allForeignNames;
		this->uniqueID = reinterpret_cast<size_t>(this)+1;
	}
	return *this;
}

bool ParameterSet::operator== ( const ParameterSet& input_rhs )
{
	bool names_match=true;
	vector<string> input_names_lhs = this->allNames;
	vector<string> input_names_rhs = input_rhs.allNames;
	//	Simple logical check for all wanted objects
	if( input_names_rhs.size() < input_names_lhs.size() ) return false;

	//	Check to see if wanted objects exist
	for( unsigned int i=0; i< input_names_lhs.size(); ++i )
	{
		//int this_index = StringProcessing::VectorContains( &input_names_rhs, &(input_names_lhs[i]) );
		input_rhs.GetPhysicsParameter( allInternalNames[i] );
		int this_index = allInternalNames[i].GetIndex();
		if( this_index == -1 )
		{
			names_match = false;
			break;
		}
	}

	if( !names_match )
	{
		return false;
	}
	else
	{
		//	Finally if needed compare the numerical values
		bool changed = false;
		for( unsigned int i=0; i< input_names_lhs.size(); ++i )
		{
			double lhs_val = this->GetPhysicsParameter( allInternalNames[i] )->GetValue();
			double rhs_val = input_rhs.GetPhysicsParameter( allForeignNames[i] )->GetValue();
			if( lhs_val != rhs_val )
			{
				changed = true;
				break;
			}
		}
		return !changed;
		//if( !changed ) return true;
		//else return false;
	}
}

//Constructor with correct arguments
ParameterSet::ParameterSet( vector<string> NewNames ) : allParameters(), allNames(), uniqueID(0), allInternalNames(), allForeignNames()
{
	vector<string> duplicates;
	allNames = StringProcessing::RemoveDuplicates( NewNames, duplicates );
	if( allNames.size() != NewNames.size() )
	{
		cerr << "WARNING: Cannot Generate a ParameterSet with 2 Occurances of the same name" << endl;
		for( vector<string>::iterator str_i = duplicates.begin(); str_i != duplicates.end(); ++str_i )
		{
			cout << *str_i << endl;
		}
		cerr << "This is harmless, but you will now have some merged PhysicsParameter(s)" << endl;
	}
	//Populate the map
	for( unsigned int nameIndex = 0; nameIndex < NewNames.size(); ++nameIndex )
	{
		allParameters.push_back( new PhysicsParameter( NewNames[nameIndex] ) );
	}
	uniqueID = reinterpret_cast<size_t>(this)+1;
	for( unsigned int i=0; i< allNames.size(); ++i )
	{
		allInternalNames.push_back( ObservableRef(allNames[i]) );
		allForeignNames.push_back( ObservableRef(allNames[i]) );
	}
}

//Destructor
ParameterSet::~ParameterSet()
{
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

//Retrieve a physics parameter by its name
PhysicsParameter * ParameterSet::GetPhysicsParameter( const string Name ) const
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	//cout << Name << "\t" << nameIndex << endl;
	if( nameIndex == -1 )
	{
		cerr << "PhysicsParameter " << Name << " not found(1)" << endl;
		throw(-20);
		//return new PhysicsParameter( Name, 0.0, 0.0, 0.0, "Error", "NameNotFoundError");
	}
	else
	{
		if( allParameters[(unsigned)nameIndex] == NULL ) allParameters[(unsigned)nameIndex] = new PhysicsParameter( allNames[(unsigned)nameIndex] );
		//	This has to be here to ensure that badly constructed parameters don't cause headaches!
		if( allParameters[(unsigned)nameIndex]->GetName().empty() || allParameters[(unsigned)nameIndex]->GetName() == "" )
		{
			allParameters[(unsigned)nameIndex]->SetName( Name );
		}
		return allParameters[(unsigned)nameIndex];
	}
}

PhysicsParameter* ParameterSet::GetPhysicsParameter( const ObservableRef& object ) const
{
	if( object.GetExternalID() != uniqueID ) object.SetIndex(-1);

	if( object.GetIndex() < 0 )
	{
		object.SetIndex( StringProcessing::VectorContains( &allNames, object.NameRef() ) );
		object.SetExternalID( uniqueID );
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
	object.Print();
	cerr << "ParameterSet: PhysicsParameter " << object.Name().c_str() << " not found(2)" << endl;
	cerr << "ParamererSet: The Likely Cause of this is that your PDF is NOT adveritising that it requires: " << object.Name().c_str() << " please check and remedy this!" << endl << endl;
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
		if( allParameters[unsigned(nameIndex)] != NULL ) delete allParameters[unsigned(nameIndex)];
		//	Copy the new parameter into the old one so no need to delete anything
		allParameters[unsigned(nameIndex)] = new PhysicsParameter(*NewPhysicsParameter);
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
bool ParameterSet::SetPhysicsParameters( const ParameterSet * NewParameterSet )
{
	vector<string> input_names = NewParameterSet->GetAllNames();
	for (unsigned short int nameIndex = 0; nameIndex < input_names.size(); nameIndex++)
	{
		string thisName = input_names[nameIndex];
		int lookup = StringProcessing::VectorContains( &allNames, &thisName );
		if( lookup == -1 )
		{
			//Fail if a required parameter is missing
			//cerr << "SetParmaters: Parameter \"" << thisName << "\" expected but not found" << endl;
		}
		else
		{
			if( allParameters[(unsigned)lookup] != NULL ) delete allParameters[(unsigned)lookup];
			//cout << thisName << endl;
			allParameters[(unsigned)lookup] = new PhysicsParameter( *NewParameterSet->GetPhysicsParameter(thisName) );
		}
	}
	return true;
}

bool ParameterSet::AddPhysicsParameter( const PhysicsParameter* NewParameter, bool replace )
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

	vector<ObservableRef> emptystring, emptystring2;
	allInternalNames.swap( emptystring );	allForeignNames.swap( emptystring2 );
	for( unsigned int i=0; i< allNames.size(); ++i )
	{
		allInternalNames.push_back( allNames[i] );
		allForeignNames.push_back( allNames[i] );
	}
	++uniqueID;
	return true;
}

//Set all physics parameters
bool ParameterSet::AddPhysicsParameters( const ParameterSet * NewParameterSet, bool replace )
{
	for (unsigned short int nameIndex = 0; nameIndex < NewParameterSet->GetAllNames().size(); nameIndex++)
	{
		string thisName = NewParameterSet->GetAllNames()[nameIndex];
		//PhysicsParameter * inputParameter = NewParameterSet->GetPhysicsParameter( NewParameterSet->GetAllNames()[nameIndex] );
		int paramIndex = StringProcessing::VectorContains( &allNames, &thisName );
		if( paramIndex == -1 )
		{
			allNames.push_back( thisName );
			PhysicsParameter* temp = new PhysicsParameter( *(NewParameterSet->GetPhysicsParameter( thisName )) );
			allParameters.push_back( temp );
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

	vector<ObservableRef> emptystring, emptystring2;
	allInternalNames.swap( emptystring );       allForeignNames.swap( emptystring2 );
	for( unsigned int i=0; i< allNames.size(); ++i )
	{
		allInternalNames.push_back( ObservableRef(allNames[i]) );
		allForeignNames.push_back( ObservableRef(allNames[i]) );
	}
	++uniqueID;
	return true;
}

//Not very pleasant in OO terms, and unsafe. Quick however.
void ParameterSet::UpdatePhysicsParameters( double* NewValues, int npar )
{
	if( npar <= -1 )
	{
		unsigned int temp=0;
		for( vector<PhysicsParameter*>::iterator parameterIndex = allParameters.begin(); parameterIndex != allParameters.end(); ++parameterIndex, ++temp )
		{
			(*parameterIndex)->SetValue( NewValues[temp] );
		}
	}
	else
	{
		vector<PhysicsParameter*>::iterator parameterIndex = allParameters.begin();
		for( unsigned int i=0; i< (unsigned) npar; ++i, ++parameterIndex )
		{
			(*parameterIndex)->SetValue( NewValues[i] );
		}
	}
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

void ParameterSet::SetPhysicsParameters( vector<double> NewValues )
{
	if( NewValues.size() == allParameters.size() )
	{
		vector<double>::iterator temp = NewValues.begin();
		for( vector<PhysicsParameter*>::iterator parameterIndex = allParameters.begin(); parameterIndex != allParameters.end(); ++parameterIndex, ++temp )
		{
			(*parameterIndex)->SetValue( *temp );
		}
	}
	else
	{
		cerr << "Parameter number mismatch: " << NewValues.size() << " vs " << allParameters.size();
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

size_t ParameterSet::GetUniqueID() const
{
	return uniqueID;
}

void ParameterSet::SetUniqueID( size_t input )
{
	uniqueID = input;
}

void ParameterSet::FloatedFirst()
{
	++uniqueID;
	vector<string> floated = this->GetAllFloatNames();
	vector<string> fixed = this->GetAllFixedNames();

	vector<PhysicsParameter*> sorted_parameters;
	vector<string> sorted_names;
	for( vector<string>::iterator name_i = floated.begin(); name_i != floated.end(); ++name_i )
	{
		sorted_parameters.push_back( this->GetPhysicsParameter( *name_i ) );
		sorted_names.push_back( *name_i );
	}
	for( vector<string>::iterator name_i = fixed.begin(); name_i != fixed.end(); ++name_i )
	{
		sorted_parameters.push_back( this->GetPhysicsParameter( *name_i ) );
		sorted_names.push_back( *name_i );
	}

	allParameters = sorted_parameters;
	allNames = sorted_names;

	vector<ObservableRef> emptyNames, emptyNames2;
	allInternalNames.swap( emptyNames );       allForeignNames.swap( emptyNames2 );
	for( unsigned int i=0; i< allNames.size(); ++i )
	{
		allInternalNames.push_back( allNames[i] );
		allForeignNames.push_back( allNames[i] );
	}
}

