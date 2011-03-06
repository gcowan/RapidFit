/**
        @class ObservableDiscreteConstraint

        A constraint defining an observable that can take any of a list of discrete values.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "ObservableDiscreteConstraint.h"
#include <iostream>
#include <math.h>

#define DOUBLE_TOLERANCE 1E-6

//Default constructor
ObservableDiscreteConstraint::ObservableDiscreteConstraint() :  unit("Uninitialised")
{
}

//Constructor with correct argument
ObservableDiscreteConstraint::ObservableDiscreteConstraint( string Name, vector<double> NewValues, string NewUnit ) : allValues(NewValues), unit(NewUnit)
{
	if (unit == "")
	{
		cerr << "Single bound \"" << Name << "\" has no unit! What kind of physicist are you?" << endl;
	}
}

//Destructor
ObservableDiscreteConstraint::~ObservableDiscreteConstraint()
{
}

//There is no minimum or maximum, so return an error
double ObservableDiscreteConstraint::GetMinimum()
{
	cerr << "Minimum requested, but constraint is discrete" << endl;	
	return 0.0;
}
double ObservableDiscreteConstraint::GetMaximum()
{
	cerr << "Maximum requested, but constraint is discrete" << endl;
	return 0.0;
}

//Get and set all values
vector<double> ObservableDiscreteConstraint::GetValues()
{
	return allValues;
}
void ObservableDiscreteConstraint::SetValues( vector<double> NewValues )
{
	allValues = NewValues;
}

//Add a value to the constraint
void ObservableDiscreteConstraint::AddValue( double NewValue )
{
	allValues.push_back(NewValue);
}

//Get the unit
string ObservableDiscreteConstraint::GetUnit()
{
	return unit;
}

//Check whether an observable fits with this constraint
bool ObservableDiscreteConstraint::CheckObservable( Observable * TestObservable )
{
	//Check the units are the same
	if ( TestObservable->GetUnit() != unit )
	{
		cerr << "Unit mismatch: boundary expects \"" << unit << "\" but observable is \"" << TestObservable->GetUnit() << "\"" << endl;
	}

	double value = TestObservable->GetValue();

	//Check if the observable value is one of those listed
	vector<double>::iterator valueIterator;
	for ( valueIterator = allValues.begin(); valueIterator != allValues.end(); valueIterator++ )
	{
		if ( ( fabs( value - *valueIterator ) < DOUBLE_TOLERANCE ) )
		{
			return true;
		}
	}

	//The value was not found
	return false;
}

//Create an observable within this constraint, without specifying a random number generator
Observable * ObservableDiscreteConstraint::CreateObservable()
{
	TRandom3 * random = new TRandom3(0);
	Observable * returnObservable = CreateObservable(random);
	delete random;
	return returnObservable;
}

//Create an observable within this constraint, using the specified random number generator
Observable * ObservableDiscreteConstraint::CreateObservable( TRandom3 * RandomNumberGenerator )
{
	int randomIndex = (int)floor( int(allValues.size()) * RandomNumberGenerator->Rndm() );
	return new Observable( "Unknown", allValues[randomIndex], 0.0, unit );
}

bool ObservableDiscreteConstraint::IsDiscrete()
{
	return true;
}
