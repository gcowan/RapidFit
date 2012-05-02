
//	RapidFit Headers
#include "ObservableRef.h"
//	System Headers
#include <string>
#include <vector>
#include <iostream>

using namespace::std;

ObservableRef::ObservableRef() : Observable_Name( "unknown" ), Observable_Index(-1), Observable_Names(), Observable_Refs()
{
}

ObservableRef::ObservableRef( string ObsName ) : Observable_Name( ObsName ), Observable_Index(-1), Observable_Names(), Observable_Refs()
{
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

string* ObservableRef::NameRef() const
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

