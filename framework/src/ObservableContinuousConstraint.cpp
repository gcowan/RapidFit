/*!
 * @class ObservableContinuousConstraint
 *
 * A constraint defining an observable that can take any value between a given maximum and minimum.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

//	RapidFit Headers
#include "ObservableContinuousConstraint.h"
//	System Headers
#include <iostream>
#include <float.h>
#include <sstream>
#include <cstdlib>

using namespace::std;

//#define CONTINUOUS_TOLERANCE DBL_MIN
#define CONTINUOUS_TOLERANCE 1E-9

//Constructor with correct argument
ObservableContinuousConstraint::ObservableContinuousConstraint( string Name, double NewMinimum, double NewMaximum, string NewUnit, string TF1 ) : 
	name(Name), minimum(NewMinimum), maximum(NewMaximum), unit(NewUnit), tf1( TF1 )
{
	if(maximum < minimum)
	{
		cerr << "Single bound \"" << Name << "\" has maximum less than minimum: values swapped" << endl;
		minimum = NewMaximum;
		maximum = NewMinimum;
	}

	if(unit == "")
	{
		cerr << "Single bound \"" << Name << "\" has no unit! What kind of physicist are you?" << endl;
	}

	if( tf1 == "" || tf1.empty() ) tf1 = name;
}

ObservableContinuousConstraint::ObservableContinuousConstraint( const IConstraint* input ) : name(), minimum(), maximum(), unit(), tf1()
{
	if( !input->IsDiscrete() )
	{
		name = input->GetName();
		minimum = input->GetMinimum();
		maximum = input->GetMaximum();
		unit = input->GetUnit();
		tf1 = input->GetTF1();
	}
	else
	{
		cerr << "Trying to Construct a Discrete Constraint from a Continuouse one" << endl;
		exit(-962);
	}
}

//Destructor
ObservableContinuousConstraint::~ObservableContinuousConstraint()
{
}

string ObservableContinuousConstraint::GetName() const
{
	return name;
}

//Get and set the minimum
double ObservableContinuousConstraint::GetMinimum() const
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
double ObservableContinuousConstraint::GetMaximum() const
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
vector<double> ObservableContinuousConstraint::GetValues() const
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
string ObservableContinuousConstraint::GetUnit() const
{
	return unit;
}

//Check whether an observable fits with this constraint
bool ObservableContinuousConstraint::CheckObservable( Observable * TestObservable ) const
{
	//Check the units are the same
	if ( TestObservable->GetUnit() != unit )
	{
		cerr << "Unit mismatch: boundary expects \"" << unit << "\" but observable is \"" << TestObservable->GetUnit() << "\"" << endl;
	}

	double value = TestObservable->GetValue();
	return (value < (maximum + CONTINUOUS_TOLERANCE) ) && ( value > (minimum - CONTINUOUS_TOLERANCE)  );
}

//Create an observable within this constraint, without specifying a random number generator
Observable * ObservableContinuousConstraint::CreateObservable() const
{
	TRandom3 * random = new TRandom3(0);
	Observable * returnObservable = this->CreateObservable(random);
	delete random;
	return returnObservable;
}

//Create an observable within this constraint, using the specified random number generator
Observable * ObservableContinuousConstraint::CreateObservable( TRandom3 * RandomNumberGenerator ) const
{
	double value = minimum + ( ( maximum - minimum ) * RandomNumberGenerator->Rndm() );
	return new Observable( name, value/*, 0.0*/, unit );
}

bool ObservableContinuousConstraint::IsDiscrete() const
{
	return false;
}

void ObservableContinuousConstraint::Print() const
{
	cout << "Maximum: " << maximum << "\tMinimum: " << minimum << "\tUnit: " << unit << endl;
}

Observable* ObservableContinuousConstraint::GetMidRangeValue() const
{
	double diff=(maximum-minimum)/2.;
	return new Observable( name, minimum+diff, unit ); 
}

string ObservableContinuousConstraint::GetTF1() const
{
	return tf1;
}

void ObservableContinuousConstraint::SetTF1( const string input )
{
	tf1 = input;
}

string ObservableContinuousConstraint::XML() const
{
	stringstream xml;

	xml << "\t<Observable>" << endl;
	xml << "\t\t<Name>" << name << "</Name>" << endl;
	xml << "\t\t<Minimum>" << minimum << "</Minimum>" << endl;
	xml << "\t\t<Maximum>" << maximum << "</Maximum>" << endl;
	xml << "\t\t<Unit>" << unit << "</Unit>" << endl;
	if( !tf1.empty() ) xml << "\t\t<TF1>" << tf1 << "</TF1>" << endl;
	xml << "\t</Observable>" << endl;

	return xml.str();
}

