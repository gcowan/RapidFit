/**
        @class Observable

        The measured value of a particular variable for an event

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "Observable.h"
#include <iostream>

//Default constructor
Observable::Observable() : value(0.0), error(0.0), unit("Uninitialised")
{
}

//Constructor with correct argument
Observable::Observable( string Name, double NewValue, double NewError, string NewUnit )
	: value(NewValue), error(NewError), unit(NewUnit)
{
	if (unit == "")
	{
		cerr << "Observable \"" << Name << "\" has no unit! What kind of physicist are you?" << endl;
	}
}

//Destructor
Observable::~Observable()
{
}

//Get and set the value
double Observable::GetValue()
{
	return value;
}
void Observable::SetValue(double NewValue)
{
	value = NewValue;
}

//Get and set the error on the value
double Observable::GetError()
{
	return error;
}
void Observable::SetError(double NewError)
{
	error = NewError;
}

//Get the unit
string Observable::GetUnit()
{
	return unit;
}
