
//	RapidFit Headers
#include "ObservableRef.h"
//	System Headers
#include <string>
#include <vector>
#include <iostream>

using namespace::std;

ObservableRef::ObservableRef( string ObsName ) : Observable_Name( ObsName ), Observable_Index(-1), Observable_Names(), Observable_Refs()
{
}

ObservableRef& ObservableRef::operator= ( const ObservableRef& input )
{
	if( this != &input )
	{
		this->Observable_Name = input.Observable_Name;
		this->Observable_Index = input.Observable_Index;
		this->Observable_Names = vector<string>( input.Observable_Names );
		this->Observable_Refs = vector<ObservableRef>();
		for( vector<ObservableRef>::const_iterator obs_i = input.Observable_Refs.begin(); obs_i != input.Observable_Refs.end(); ++obs_i )
		{
			this->Observable_Refs.push_back( ObservableRef( *obs_i ) );
		}
	}

	return *this;
}

ObservableRef::ObservableRef( const ObservableRef& input ) :
	Observable_Name( input.Observable_Name ), Observable_Index( input.Observable_Index ), Observable_Names( input.Observable_Names ), Observable_Refs()
{
	for( vector<ObservableRef>::const_iterator obs_i = input.Observable_Refs.begin(); obs_i != input.Observable_Refs.end(); ++obs_i )
	{
		Observable_Refs.push_back( ObservableRef( *obs_i ) );
	}
}

ObservableRef::ObservableRef( vector<string> ObsList ) : Observable_Name(), Observable_Index(-2), Observable_Names(ObsList), Observable_Refs()
{
	for( vector<string>::iterator name_i=ObsList.begin(); name_i!=ObsList.end(); ++name_i )
	{
		Observable_Refs.push_back( ObservableRef( *name_i ) );
	}
}

ObservableRef::~ObservableRef()
{
}

string ObservableRef::Name() const
{
	return Observable_Name;
}

const string* ObservableRef::NameRef() const
{
	return &Observable_Name;
}

void ObservableRef::SetIndex( const int Index ) const
{
	Observable_Index = Index;
}

int ObservableRef::GetIndex() const
{
	return Observable_Index;
}

size_t ObservableRef::size() const
{
	return Observable_Refs.size();
}

void ObservableRef::push_back( string input_name )
{
	Observable_Refs.push_back( ObservableRef( input_name ) );
	Observable_Names.push_back( input_name );
}

