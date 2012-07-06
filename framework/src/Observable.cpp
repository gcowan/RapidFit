/**
  @class Observable

  The measured value of a particular variable for an event

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

//	RapidFit Headers
#include "Observable.h"
//	System Headers
#include <iostream>
#include <vector>
#include <string>

using namespace::std;

//Constructor with correct argument
Observable::Observable( string Name, double NewValue, double NewError, string NewUnit )
	: name(Name), value(NewValue), error(NewError), unit(NewUnit), bin_num(-1), acceptance(-1.)
{
	if (unit == "")
	{
		cerr << "Observable \"" << Name << "\" has no unit! What kind of physicist are you?" << endl;
	}
}

Observable::Observable( string Name ) : name(Name), value(0.), error(0.), unit("Uninitialised"), bin_num(-1), acceptance(-1.)
{
}

//Destructor
Observable::~Observable()
{
}

//Get and set the value
double Observable::GetValue() const
{
	return value;
}

//Get and set the error on the value
double Observable::GetError() const
{
	return error;
}

//Get the unit
string Observable::GetUnit() const
{
	return unit;
}

void Observable::SetBinNumber( int input ) const
{
	bin_num = input;
}

int Observable::GetBinNumber() const
{
	return bin_num;
}

void Observable::SetAcceptance( double input ) const
{
	acceptance = input;
}

double Observable::GetAcceptance() const
{
	return acceptance;
}

string Observable::GetName() const
{
	return name;
}

void Observable::Print() const
{
	cout << "Name: " << name;
	cout << "\tValue: " << value;
	cout << "\tError: " << error;
	cout << "\tUnit: " << unit ;
	if( bin_num != -1 ) cout << "\tBinNum: " << bin_num << "\tAcceptance: " << acceptance;
	cout << endl;
}

void Observable::SetObservable( Observable* input )
{
	name = input->GetName();
	value = input->GetValue();
	error = input->GetError();
	unit = input->GetUnit();
}



