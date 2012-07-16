
///	RapidFit Headers
#include "PseudoObservable.h"
#include "Observable.h"
#include "ObservableRef.h"
///	System Headers
#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using namespace::std;

PseudoObservable::PseudoObservable() :
	internal_function(NULL), Dependencies(), Value(0.), valid(false), Name(" "),
	Index(-1), internal_Input()
{
}

PseudoObservable::PseudoObservable( string InputName ) : 
	internal_function(NULL), Dependencies(), Value(0.), valid(false), Name(InputName),
	Index(-1), internal_Input()
{
}

PseudoObservable::PseudoObservable( const PseudoObservable& Input ) :
	internal_function( Input.internal_function ), Dependencies(Input.Dependencies),
	Value( Input.Value), valid( false ), Name( Input.Name ), Index( Input.Index ),
	internal_Input( Input.internal_Input )
{
}

PseudoObservable& PseudoObservable::operator= ( const PseudoObservable& input )
{
	if( this != &input )
	{
		this->internal_function = input.internal_function;
		this->Dependencies = input.Dependencies;
		this->Value = input.Value;
		this->valid = input.valid;
		this->Name = input.Name;
		this->Index = input.Index;
		this->internal_Input = input.internal_Input;
	}
	return *this;
}

PseudoObservable::~PseudoObservable()
{
}

void PseudoObservable::AddFunction( double (*pseudoRelation)(vector<double>) )
{
	internal_function = pseudoRelation;
}

void PseudoObservable::AddDependency( ObservableRef Input )
{
	Dependencies.push_back( ObservableRef( Input ) );
}

void PseudoObservable::AddDependencies( vector<ObservableRef> Input )
{
	for( vector<ObservableRef>::iterator ref_i = Input.begin(); ref_i != Input.end(); ++ref_i )
	{
		Dependencies.push_back( ObservableRef( (*ref_i) ) );
	}
}

void PseudoObservable::SetInput( const vector<double> Input )
{
	internal_Input = Input;
}

void PseudoObservable::SetValid( const bool Input ) const
{
	valid = Input;
}

bool PseudoObservable::GetValid() const
{
	return valid;
}

vector<ObservableRef>* PseudoObservable::GetDependencies()
{
	return &Dependencies;
}

double PseudoObservable::GetValue() const
{
	if( internal_function == NULL ) return 0.;
	else
	{
		if( valid ) return Value;
		else
		{
			Value = internal_function( internal_Input );
			return Value;
		}
	}
}

string PseudoObservable::GetName() const
{
	return Name;
}

double PseudoObservable::GetPseudoObservable()
{
	if( internal_Input.empty() ) return 0.;
	if( valid == true )
	{
		double returnable=0.;
		if( Value == 0. )
		{
			valid = false;	//	Valid=false here to trigger the re-calculation
			returnable = this->GetValue();
			valid = true;
		}
		returnable = this->GetValue();
		return returnable;
	}
	else
	{
		double returnable = this->GetValue();
		return returnable;
	}
}

int PseudoObservable::GetIndex() const
{
	return Index;
}

void PseudoObservable::SetIndex( const int Input ) const
{
	Index = Input;
}

void PseudoObservable::Print() const
{
	return;
}

bool PseudoObservable::GetValid( const vector<double> Input ) const
{
	if( Input.size() != internal_Input.size() )
	{
		valid = false;
		return false;
	}
	else
	{
		bool check_ok = true;
		for( unsigned int i=0; i< internal_Input.size(); ++i )
		{
			if( internal_Input[i] != Input[i] )
			{
				check_ok = false;
				break;
			}
		}
		valid = check_ok;
		return check_ok;
	}
}

vector<double> PseudoObservable::GetCacheInput() const
{
	return internal_Input;
}

