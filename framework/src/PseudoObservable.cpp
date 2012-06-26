
///	RapidFit Headers
#include "PseudoObservable.h"
#include "Observable.h"
#include "ObservableRef.h"
///	System Headers
#include <vector>
#include <string>
#include <iostream>

using namespace::std;

PseudoObservable::PseudoObservable() :
	internal_function(NULL), internal_error_function(NULL), Dependencies(), Value(0.), Error(0.), valid(false), Name("undefinded"), Unit("unitless"),
	Index(-1), internal_Input(), internalObservable(NULL), internal_function_extra(NULL)
{
}

PseudoObservable::PseudoObservable( string InputName ) : 
	internal_function(NULL), internal_error_function(NULL), Dependencies(), Value(0.), Error(0.), valid(false), Name(InputName), Unit("unitless"),
	Index(-1), internal_Input(), internalObservable(NULL), internal_function_extra(NULL)
{
}

PseudoObservable::PseudoObservable( const PseudoObservable& Input ) :
	internal_function( Input.internal_function ), internal_error_function( Input.internal_error_function ), Dependencies(Input.Dependencies),
	Value( Input.Value), Error( Input.Error ), valid( Input.valid ), Name( Input.Name ), Unit( Input.Unit ), Index( Input.Index ),
	internal_Input( Input.internal_Input ), internalObservable(NULL), internal_function_extra(NULL)
{
	if( Input.internalObservable != NULL ) internalObservable = new Observable( *(Input.internalObservable) );
}

PseudoObservable& PseudoObservable::operator= ( const PseudoObservable& input )
{
	if( this != &input )
	{
		this->internal_function = input.internal_function;
		this->internal_error_function = input.internal_error_function;
		this->Dependencies = input.Dependencies;
		this->Value = input.Value;
		this->Error = input.Error;
		this->valid = input.valid;
		this->Name = input.Name;
		this->Unit = input.Unit;
		this->Index = input.Index;
		this->internal_Input = input.internal_Input;
		this->internalObservable = (input.internalObservable==NULL)?NULL:new Observable( *(input.internalObservable) );
		this->internal_function_extra = input.internal_function_extra;
	}
	return *this;
}

PseudoObservable::~PseudoObservable()
{
	if( internalObservable != NULL ) delete internalObservable;
}

void PseudoObservable::AddFunction( double (*pseudoRelation)(vector<double>) )
{
	internal_function = pseudoRelation;
}

void PseudoObservable::AddErrorFunction( double (*pseudoRelation)(vector<double>) )
{
	internal_error_function = pseudoRelation;
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

void PseudoObservable::SetInput( vector<double> Input )
{
	internal_Input = Input;
}

void PseudoObservable::SetValid( bool Input )
{
	valid = Input;
}

bool PseudoObservable::GetValid()
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
			Error = GetError();
			valid = true;
			return Value;
		}
	}
}

double PseudoObservable::GetError() const
{
	if( internal_error_function == NULL ) return 0.;
	else
	{
		if( valid ) return Value;
		else
		{
			Error = internal_error_function( internal_Input );
			Value = GetValue();
			valid = true;
			return Error;
		}
	}
}

string PseudoObservable::GetName() const
{
	return Name;
}

void PseudoObservable::SetUnit( string Input )
{
	Unit = Input;
}

string PseudoObservable::GetUnit() const
{
	return Unit;
}

Observable* PseudoObservable::GetPseudoObservable()
{
	if( internal_Input.empty() ) return NULL;
	if( valid == true )
	{
		if( internalObservable == NULL )
		{
			valid = false;
			internalObservable = new Observable( this->GetName(), this->GetValue(), this->GetError(), this->GetUnit() );
		}
		return internalObservable;
	}
	else
	{
		if( internalObservable != NULL ) delete internalObservable;

		internalObservable = new Observable( this->GetName(), this->GetValue(), this->GetError(), this->GetUnit() );
		valid = true;
		return internalObservable;
	}
}

int PseudoObservable::GetIndex() const
{
	return Index;
}

void PseudoObservable::SetIndex( int Input )
{
	Index = Input;
}

void PseudoObservable::Print() const
{
	return;
}


