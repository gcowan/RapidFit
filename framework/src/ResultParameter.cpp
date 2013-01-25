/**
  @class ResultParameter

  A physics parameter after it has been fitted

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

//	RapidFit Headers
#include "ResultParameter.h"
#include "PhysicsParameter.h"
#include "ExternalConstraint.h"
//	System Headers
#include <iostream>
#include <math.h>
#include <float.h>
#include <sstream>

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

using namespace::std;

ResultParameter::ResultParameter( const PhysicsParameter* Input )
	: name("Undefined"), value(0.0), originalValue(-9999.), error(0.0), minimum(0.0), maximum(0.0), stepSize(0.0), type("Uninitialised"), unit("Uninitialised"), ScanStatus(false),
	assym_err(false), err_hi(0.), err_low(0.)
{
	name = Input->GetName();
	value = Input->GetValue();
	originalValue = Input->GetOriginalValue();
	error = Input->GetStepSize();
	if( error < 0 ) error = 0.;
	minimum = Input->GetMinimum();
	maximum = Input->GetMaximum();
	stepSize = Input->GetStepSize();
	type = Input->GetType();
	unit = Input->GetUnit();
	if (maximum < minimum)
	{
		cerr << "Result parameter \"" << name << "\" has maximum less than minimum: values swapped" << endl;
		double temp = minimum;
		minimum = maximum;
		maximum = temp;
	}

	if ( (value < minimum) && ( fabs( minimum ) > DOUBLE_TOLERANCE ) )
	{
		cerr << "Result parameter \"" << name << "\" has value less than minimum" << endl;
	}

	if ( (value > maximum) && ( fabs( maximum ) > DOUBLE_TOLERANCE ) )
	{
		cerr << "Result parameter \"" << name << "\" has value greater than maximum" << endl;
	}

	if (unit == "")
	{
		cerr << "Result parameter \"" << name << "\" has no unit! What kind of physicist are you?" << endl;
	}
}

//Constructor with correct argument
ResultParameter::ResultParameter( string Name, double NewValue, double NewOriginalValue, double NewError, double NewMinimum, double NewMaximum, string NewType, string NewUnit ) :
	name(Name), value(NewValue), originalValue(NewOriginalValue), error(NewError), minimum(NewMinimum), maximum(NewMaximum), stepSize(NewError), type(NewType), unit(NewUnit), ScanStatus(false),
	assym_err(false), err_hi(0.), err_low(0.)
{
	if (maximum < minimum)
	{
		cerr << "Result parameter \"" << Name << "\" has maximum less than minimum: values swapped" << endl;
		minimum = NewMaximum;
		maximum = NewMinimum;
	}

	if ( (value < minimum) && ( fabs( minimum ) > DOUBLE_TOLERANCE ) ) 
	{
		cerr << "Result parameter \"" << Name << "\" has value less than minimum" << endl;
	}

	if ( (value > maximum) && ( fabs( maximum ) > DOUBLE_TOLERANCE ) )
	{
		cerr << "Result parameter \"" << Name << "\" has value greater than maximum" << endl;
	}

	if (unit == "")
	{
		cerr << "Result parameter \"" << Name << "\" has no unit! What kind of physicist are you?" << endl;
	}
}

//Destructor
ResultParameter::~ResultParameter()
{
}

bool ResultParameter::GetAssym() const
{
	return assym_err;
}

void ResultParameter::SetAssymErrors( const double err_plus, const double err_minus, const double err_para )
{
	assym_err=true;
	error = err_para;
	err_low = err_minus;
	err_hi = err_plus;
}

double ResultParameter::GetErrLow() const
{
	return err_low;
}

double ResultParameter::GetErrHi() const
{
	return err_hi;
}

PhysicsParameter* ResultParameter::GetDummyPhysicsParameter() const
{
	if( (type == "Fixed") || ( fabs(minimum - maximum) < 1E-6 ) )
	{
		return new PhysicsParameter( name, value, stepSize, type, unit );
	}
	else
	{
		return new PhysicsParameter( name, value, minimum, maximum, stepSize, type, unit );
	}
}

string ResultParameter::GetName() const
{
	return name;
}

//Get the value
double ResultParameter::GetValue() const
{
	return value;
}

//Get the original value
double ResultParameter::GetOriginalValue() const
{
	return originalValue;
}

void ResultParameter::ForceOriginalValue( double new_value )
{
	originalValue = new_value;
}

void ResultParameter::ForcePullValue( double new_value )
{
	originalValue = new_value;
}

//Get the error
double ResultParameter::GetError() const
{
	return error;
}

void ResultParameter::SetError( double input )
{
	error = input;
}

//Get the pull
double ResultParameter::GetPull() const
{
	if( error < 1E-6 ) return 0;
	else return (value - originalValue) / error;
}

//Get the minimum
double ResultParameter::GetMinimum() const
{       
	return minimum;
}

//Get the maximum
double ResultParameter::GetMaximum() const
{       
	return maximum;
}

//Get the type
string ResultParameter::GetType() const
{
	return type;
}

void ResultParameter::ForceType( string input )
{
	type = input;
}

//Get the unit
string ResultParameter::GetUnit() const
{
	return unit;
}

double ResultParameter::GetStepSize() const
{
	return stepSize;
}

bool ResultParameter::GetScanStatus() const
{
	return ScanStatus;
}

void ResultParameter::SetScanStatus( bool input )
{
	ScanStatus = input;
}

void ResultParameter::Print() const
{
	cout << "ResultParameter:" << endl;
	cout << "\tName:\t" << name << endl;
	cout << "\tValue:\t" << value << endl;
	cout << "\tError:\t" << error << endl;
	cout << "\tGen:\t" << originalValue << endl;
	cout << "\tPull:\t" << this->GetPull() << endl;
	cout << "\tType:\t" << type << endl;
	cout << "\tUnit:\t" << unit << endl;
	if( ScanStatus ) cout << "\tControlled" << endl;
	else cout << "\tNot Controlled" << endl;
}

string ResultParameter::XML( const bool fit ) const
{
	stringstream xml;
	xml << "\t<PhysicsParameter>" << endl;
	xml << "\t\t<Name>" << name << "</Name>" << endl;
	xml << "\t\t<Value>" << value << "</Value>" << endl;
	xml << "\t\t<StepSize>" << error << "</StepSize>" << endl;
	if( fit == true || type=="Fixed" )
	{
		xml << "\t\t<Type>Fixed</Type>" << endl;
	}
	else
	{
		xml << "\t\t<Type>Free</Type>" << endl;
		if( fabs(minimum-maximum)>1E-6 )
		{
			xml << "\t\t<Maximum>" << maximum << "</Maximum>" << endl;
			xml << "\t\t<Minimum>" << minimum << "</Minimum>" << endl;
		}
	}
	xml << "\t\t<Unit>" << unit << "</Unit>" << endl;
	xml << "\t</PhysicsParameter>" << endl;
	return xml.str();
}

string ResultParameter::FitXML() const
{
	return this->XML( true );
}

string ResultParameter::ToyXML() const
{
	return this->XML( false );
}

string ResultParameter::ConstraintXML() const
{
	ExternalConstraint* newConstraint = new ExternalConstraint( this->GetName(), this->GetValue(), this->GetError() );
	string output = newConstraint->XML();
	delete newConstraint;
	return output;
}

