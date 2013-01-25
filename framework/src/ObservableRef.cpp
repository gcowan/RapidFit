
//	RapidFit Headers
#include "ObservableRef.h"
//	System Headers
#include <string>
#include <vector>
#include <iostream>

using namespace::std;

ObservableRef::ObservableRef( string ObsName ) : Observable_Name( ObsName ), Observable_Index(-1), externID(0)
{
}

ObservableRef& ObservableRef::operator= ( const ObservableRef& input )
{
	if( this != &input )
	{
		this->Observable_Name = input.Observable_Name;
		this->Observable_Index = input.Observable_Index;
		this->externID = 0;
	}

	return *this;
}

ObservableRef::ObservableRef( const ObservableRef& input ) :
	Observable_Name( input.Observable_Name ), Observable_Index( input.Observable_Index ), externID(0)
{
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

void ObservableRef::Print() const
{
	cout << "ObservableRef:" << endl;
	cout << "Name:\t" << Observable_Name << endl;
	cout << "Index:\t" << Observable_Index << endl;
	cout << "UniqueID:\t" << externID << endl;
	return;
}

void ObservableRef::SetExternalID( size_t input ) const
{
	externID = input;
}

size_t ObservableRef::GetExternalID() const
{
	return externID;
}

