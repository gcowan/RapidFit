/**
  @class PhysicsParameter

  A parameter to be adjusted by the fitter, with a starting value and limits
  Modified by Pete Clarke december 2010 to add blinding 

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

#include "PhysicsParameter.h"
#include <iostream>

//Default constructor
PhysicsParameter::PhysicsParameter() : value(0.0), originalValue(0.0), minimum(0.0), maximum(0.0), type("Uninitialised"), unit("Uninitialised"), toBeBlinded(false), blindOffset(0.0)
{
}

//Constructor with correct argument
PhysicsParameter::PhysicsParameter( string Name, double NewValue, double NewMinimum, double NewMaximum, string NewType, string NewUnit )
	: value(NewValue), minimum(NewMinimum), maximum(NewMaximum), type(NewType), unit(NewUnit), toBeBlinded(false), blindOffset(0.0)
{
	if ( maximum < minimum )
	{
		cerr << "Physics parameter \"" << Name << "\" has maximum less than minimum: values swapped" << endl;
		minimum = NewMaximum;
		maximum = NewMinimum;
	}

	if ( value < minimum )
	{
		cerr << "Physics parameter \"" << Name << "\" has value less than minimum: value set to minimum" << endl;
		value = minimum;
	}

	if ( value > maximum )
	{
		cerr << "Physics parameter \"" << Name << "\" has value greater than maximum: value set to maximum" << endl;
		value = maximum;
	}

	if ( unit == "" )
	{
		cerr << "Physics parameter \"" << Name << "\" has no unit! What kind of physicist are you?" << endl;
		unit = "Not specified";
	}

	originalValue = value;
}

//Constructor for unbounded parameter
PhysicsParameter::PhysicsParameter( string Name, double NewValue, string NewType, string NewUnit ) : value(NewValue), type(NewType), unit(NewUnit),
	maximum(0.0), minimum(0.0), toBeBlinded(false), blindOffset(0.0)
{
	//You could define a fixed parameter with no maximum or minimum, but it must be unbounded if not fixed.
	if ( type != "Fixed" )
	{
		type = "Unbounded";
	}

	if ( unit == "" )
	{
		cerr << "Physics parameter \"" << Name << "\" has no unit! What kind of physicist are you?" << endl;
	}

	originalValue = value;
}

//Destructor
PhysicsParameter::~PhysicsParameter()
{
}

//Get and set the value
double PhysicsParameter::GetValue()
{
	if( toBeBlinded ) 
	{
		return value + blindOffset;
	}
	else return value ;
}
void PhysicsParameter::SetValue(double NewValue)
{
	value = NewValue;
}

//Get and set the minimum
double PhysicsParameter::GetMinimum()
{
	if ( type == "Unbounded" )
	{
		cerr << "Minimum of unbounded parameter requested" << endl;
	}

	return minimum;
}
void PhysicsParameter::SetMinimum(double NewMinimum)
{
	if ( type == "Unbounded" )
	{
		cerr << "Tried to change the minimum of an unbounded parameter" << endl;
	}
	else
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
}

//Get and set the maximum
double PhysicsParameter::GetMaximum()
{
	if ( type == "Unbounded" )
	{
		cerr << "Maximum of unbounded parameter requested" << endl;
	}

	return maximum;
}
void PhysicsParameter::SetMaximum(double NewMaximum)
{
	if ( type == "Unbounded" )
	{
		cerr << "Tried to change the maximum of an unbounded parameter" << endl;
	}
	else
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
}

//Set max and min at same time: avoids annoying situations
void PhysicsParameter::SetLimits(double NewMaximum, double NewMinimum)
{
	if ( type == "Unbounded" )
	{
		cerr << "Attempted to set the limits for an unbounded parameter" << endl;
	}
	else
	{
		if ( NewMaximum >= NewMinimum )
		{
			maximum = NewMaximum;
			minimum = NewMinimum;
		}
		else
		{
			cerr << "Attempted to set new parameter maximum lower than new minimum ( " << NewMaximum << " < " << NewMinimum << " ) : inverted" << endl;
			maximum = NewMinimum;
			minimum = NewMaximum;
		}

		if (value < minimum)
		{
			cerr << "New minimum greater than parameter value ( " << value << " ) : changing value to minimum ( " << minimum << " )" << endl;
			value = minimum;
		}

		if (value > maximum)
		{
			cerr << "New maximum less than parameter value ( " << value << " ) : changing value to maximum ( " << maximum << " )" << endl;
			value = maximum;
		}
	}
}

//Get and set the type
string PhysicsParameter::GetType()
{
	return type;
}
void PhysicsParameter::SetType(string NewType)
{
	type = NewType;
}

//Get the original value
double PhysicsParameter::GetOriginalValue()
{
	return originalValue;
}

//Get the unit
string PhysicsParameter::GetUnit()
{
	return unit;
}

//Set blinding offset
void PhysicsParameter::SetBlindOffset( double offset )
{
	blindOffset = offset ;
	toBeBlinded = true ;
	return ;
}

//Set blinding on or off
void PhysicsParameter::SetBlinding( bool state )
{
	toBeBlinded = state ;
	return ;
}

//General print
void PhysicsParameter::print()
{
	cout << "   value       " << value << endl ;
	cout << "   blindOffset " << blindOffset << endl ;
	cout << "   minimum     " << minimum << endl ;
	cout << "   maximum     " << maximum << endl ;
	cout << "   type        " << type << endl ;
	cout << "   unit        " << unit << endl ;
}
