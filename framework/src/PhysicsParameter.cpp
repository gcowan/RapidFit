/*!
 *
 * @class PhysicsParameter
 *
 * A parameter to be adjusted by the fitter, with a starting value and limits
 * Modified by Pete Clarke december 2010 to add blinding 
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

//	RapidFit Headers
#include "PhysicsParameter.h"
//	System Headers
#include <iostream>
#include <sstream>

using namespace::std;

const double default_val = -9999.;

//Default constructor
PhysicsParameter::PhysicsParameter( string Name ) :
	name(Name), value(default_val), originalValue(default_val), minimum(default_val), maximum(default_val), stepSize(default_val),
	type("Uninitialised"), unit("Uninitialised"), toBeBlinded(false), blindOffset(default_val), blindString("uninitialized"), blindScale(-999.)
{
}

//Constructor with correct argument
PhysicsParameter::PhysicsParameter( string Name, double NewValue, double NewMinimum, double NewMaximum, double StepSize, string NewType, string NewUnit )
	: name(Name), value(NewValue), originalValue(NewValue), minimum(NewMinimum), maximum(NewMaximum), stepSize(StepSize),
	type(NewType), unit(NewUnit), toBeBlinded(false), blindOffset(0.0), blindString("uninitialized"), blindScale(-999.)
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
PhysicsParameter::PhysicsParameter( string Name, double NewValue, double StepSize, string NewType, string NewUnit ) :
	value(NewValue), originalValue(NewValue), minimum(0.0), maximum(0.0), stepSize(StepSize), type(NewType), unit(NewUnit), toBeBlinded(false), blindOffset(0.0), blindString("uninitialized"), blindScale(-999.), name(Name)
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

void PhysicsParameter::SetName( string Input )
{
	name = Input;
}

string PhysicsParameter::GetName() const
{
	return name;
}

//............ Get ans Set methods became complex since adding blinding  .......

//Get the unblinded value. 
double PhysicsParameter::GetValue() const
{
	return this->GetTrueValue();
}

// Set the blinded value 
// {This sets a different thing to what GetValue gets for historic reasons.}
void PhysicsParameter::SetValue(double NewValue)
{
	this->SetBlindedValue( NewValue );
}

//Get the blinded value
double PhysicsParameter::GetBlindedValue() const
{
	return value;
}
//Set the blinded value
void PhysicsParameter::SetBlindedValue(double NewValue)
{
	value = NewValue;
}

//Get the true value
double PhysicsParameter::GetTrueValue() const
{
	double new_value=-9999;
	if( toBeBlinded ) 
	{
		new_value = value + blindOffset;   
	}
	else new_value = value ;
	return new_value;
}

//Set the true value
void PhysicsParameter::SetTrueValue(double NewValue)
{
	if( toBeBlinded ) 
	{
		value = NewValue - blindOffset;   
	}
	else value = NewValue;	
}

//.....................
//Get and set the minimum
double PhysicsParameter::GetMinimum() const
{
	//if ( type == "Unbounded" )	//	We have defined sensible behaviour in this instance, this is not an error
	//{
	//	cerr << "Minimum of unbounded parameter requested" << endl;
	//}

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
double PhysicsParameter::GetMaximum() const
{
	//if ( type == "Unbounded" )	//	We have defined sensible behviour in this instance, this is not an error
	//{
	//	cerr << "Maximum of unbounded parameter requested" << endl;
	//}

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
string PhysicsParameter::GetType() const
{
	return type;
}

void PhysicsParameter::SetType(string NewType)
{
	if( NewType == "Fixed" ) originalValue = this->GetBlindedValue();
	type = NewType;
}

//Get the original value
double PhysicsParameter::GetOriginalValue() const
{
	return originalValue;
}

void PhysicsParameter::ForceOriginalValue( double new_original_value )
{
	originalValue = new_original_value;
}

//Get the unit
string PhysicsParameter::GetUnit() const
{
	return unit;
}

//Set blinding offset
void PhysicsParameter::SetBlindOffset( double offset )
{
	blindOffset = offset;
	toBeBlinded = true;
	return;
}

//Set blinding on or off
void PhysicsParameter::SetBlinding( bool state )
{
	toBeBlinded = state;
	return;
}

//General print
void PhysicsParameter::Print() const
{
        cout << "   value       " << value << endl;
        //cout << "   blindOffset " << blindOffset << endl;
        cout << "   blinded?    " << string( toBeBlinded ? "Yes" : "No" ) << endl;
        cout << "   originalVal " << originalValue << endl;
        cout << "   minimum     " << minimum << endl;
        cout << "   maximum     " << maximum << endl;
        cout << "   stepsize    " << stepSize << endl;
        cout << "   type        " << type << endl;
        cout << "   unit        " << unit << endl;
}

double PhysicsParameter::GetStepSize() const
{
	return stepSize;
}

void PhysicsParameter::SetStepSize( double newStep )
{
	stepSize = newStep;
}

string PhysicsParameter::XML() const
{
	stringstream xml;
	xml << "\t<PhysicsParameter>" << endl;
	xml << "\t\t<Name>" << name << "</Name>" << endl;
	xml << "\t\t<Value>" << value << "</Value>" << endl;
	if( toBeBlinded )
	{
		xml << "\t\t<BlindString>" << blindString << "</BlindString>" << endl;
		xml << "\t\t<BlindScale>" << blindScale << "</BlindScale>" << endl;
	}
	if( type != "Fixed" )
	{
		xml << "\t\t<Minimum>" << minimum << "</Minimum>" << endl;
		xml << "\t\t<Maximum>" << maximum << "</Maximum>" << endl;
	}
	xml << "\t\t<Type>";
	if( type != "Uninitialised" ) xml << type << "</Type>" << endl;
	else xml << "Free" << "</Type>" << endl;
	xml << "\t\t<Unit>";
	if( unit != "Uninitialised" ) xml << unit << "</Unit>" << endl;
	else xml << "someUnit" << "</Unit>" << endl;
	xml << "\t</PhysicsParameter>" << endl;
	return xml.str();
}

void PhysicsParameter::SetBlindingInfo( string input_str, double input_val )
{
	blindString = input_str;
	blindScale = input_val;
}

