/**
        @class ObservableContinuousConstraint

        A constraint defining an observable that can take any value between a given maximum and minimum.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "ObservableContinuousConstraint.h"
#include <iostream>
#define DOUBLE_TOLERANCE 1E-6

//Default constructor
ObservableContinuousConstraint::ObservableContinuousConstraint() : minimum(0.0), maximum(0.0), unit("")
{
}

//Constructor with correct argument
ObservableContinuousConstraint::ObservableContinuousConstraint( string Name, double NewMinimum, double NewMaximum, string NewUnit ) : minimum(NewMinimum),
	maximum(NewMaximum), unit(NewUnit)
{
	if (maximum < minimum)
	{
		cerr << "Single bound \"" << Name << "\" has maximum less than minimum: values swapped" << endl;
		minimum = NewMaximum;
		maximum = NewMinimum;
	}

	if (unit == "")
	{
		cerr << "Single bound \"" << Name << "\" has no unit! What kind of physicist are you?" << endl;
	}
}

//Destructor
ObservableContinuousConstraint::~ObservableContinuousConstraint()
{
}

//Get and set the minimum
double ObservableContinuousConstraint::GetMinimum()
{       
	return minimum;
}
void ObservableContinuousConstraint::SetMinimum(double NewMinimum)
{
	if ( NewMinimum < maximum )
	{
		minimum = NewMinimum;
	}
	else
	{
		cerr << "Attempted to set parameter minimum to " << NewMinimum << " which is greater than maximum ( " << maximum << " ) : ignored" << endl;
	}
}

//Get and set the maximum
double ObservableContinuousConstraint::GetMaximum()
{       
	return maximum;
}
void ObservableContinuousConstraint::SetMaximum(double NewMaximum)
{
	if ( NewMaximum > minimum )
	{
		maximum = NewMaximum;
	}
	else
	{
		cerr << "Attempted to set parameter maximum to " << NewMaximum << " which is less than minimum ( " << minimum << " ) : ignored" << endl;
	}
}

//Return an error
vector<double> ObservableContinuousConstraint::GetValues()
{
	cerr << "Values of constraint requested, but constraint is continuous" << endl;
	vector<double> emptyVector;
	return emptyVector;
}

//Set max and min at same time: avoids annoying situations
void ObservableContinuousConstraint::SetLimits(double NewMaximum, double NewMinimum)
{
	if ( NewMaximum > NewMinimum )
	{
		maximum = NewMaximum;
		minimum = NewMinimum;
	}
	else
	{
		cerr << "Attempted to set new bound maximum lower than new minimum ( " << NewMaximum << " < " << NewMinimum << " ) : inverted" << endl;
		maximum = NewMinimum;
		minimum = NewMaximum;
	}
}

//Get the unit
string ObservableContinuousConstraint::GetUnit()
{
	return unit;
}

//Check whether an observable fits with this constraint
bool ObservableContinuousConstraint::CheckObservable( Observable * TestObservable )
{
	//Check the units are the same
	if ( TestObservable->GetUnit() != unit )
	{
		cerr << "Unit mismatch: boundary expects \"" << unit << "\" but observable is \"" << TestObservable->GetUnit() << "\"" << endl;
	}

	double value = TestObservable->GetValue();
	return (value < (maximum + DOUBLE_TOLERANCE) ) && ( value > (minimum + DOUBLE_TOLERANCE) );
}

//Create an observable within this constraint, without specifying a random number generator
Observable * ObservableContinuousConstraint::CreateObservable()
{
	TRandom3 * random = new TRandom3(0);
	Observable * returnObservable = CreateObservable(random);
	delete random;
	return returnObservable;
}

//Create an observable within this constraint, using the specified random number generator
Observable * ObservableContinuousConstraint::CreateObservable( TRandom3 * RandomNumberGenerator )
{
	double value = minimum + ( ( maximum - minimum ) * RandomNumberGenerator->Rndm() );
	return new Observable( "Unknown", value, 0.0, unit );
}

bool ObservableContinuousConstraint::IsDiscrete()
{
	return false;
}
