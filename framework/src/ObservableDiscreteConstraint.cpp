/*!
 * @class ObservableDiscreteConstraint
 *
 * A constraint defining an observable that can take any of a list of discrete values.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
*/

///	RapidFit Headers
#include "ObservableDiscreteConstraint.h"
///	System Headers
#include <iostream>
#include <math.h>
#include <float.h>
#include <sstream>
#include <cstdlib>

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

using namespace::std;

//Constructor with correct argument
ObservableDiscreteConstraint::ObservableDiscreteConstraint( string Name, vector<double> NewValues, string NewUnit, string TF1 )
: name(Name), allValues(NewValues), unit(NewUnit), tf1(TF1)
{
	if (unit == "")
	{
		cerr << "Single bound \"" << Name << "\" has no unit! What kind of physicist are you?" << endl;
	}

	if( tf1 == "" || tf1.empty() ) tf1 = name;
}

ObservableDiscreteConstraint::ObservableDiscreteConstraint( const IConstraint* input ) : name(), allValues(), unit(), tf1()
{
	if( input->IsDiscrete() )
	{
		name=input->GetName();
		allValues=input->GetValues();
		unit=input->GetUnit();
		tf1=input->GetTF1();
	}
	else
	{
		cerr << "Trying to Construct a Discrete Constraint from a Continuouse one" << endl;
		exit(-962);
	}
}

//Destructor
ObservableDiscreteConstraint::~ObservableDiscreteConstraint() 
{
}

string ObservableDiscreteConstraint::GetName() const
{
	return name;
}

//There is no minimum or maximum, so return an error
double ObservableDiscreteConstraint::GetMinimum() const
{
	cerr << "Minimum requested, but constraint is discrete" << endl;	
	return 0.0;
}
double ObservableDiscreteConstraint::GetMaximum() const
{
	cerr << "Maximum requested, but constraint is discrete" << endl;
	return 0.0;
}

//Get and set all values
vector<double> ObservableDiscreteConstraint::GetValues() const
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
string ObservableDiscreteConstraint::GetUnit() const
{
	return unit;
}

//Check whether an observable fits with this constraint
bool ObservableDiscreteConstraint::CheckObservable( Observable * TestObservable ) const
{
	//Check the units are the same
	if ( TestObservable->GetUnit() != unit )
	{
		cerr << "Unit mismatch: boundary expects \"" << unit << "\" but observable is \"" << TestObservable->GetUnit() << "\"" << endl;
	}

	double value = TestObservable->GetValue();

	//Check if the observable value is one of those listed
	for ( unsigned int i=0; i< allValues.size(); ++i )
	{
		if ( ( fabs( value - allValues[i] ) < DOUBLE_TOLERANCE ) )
		{
			return true;
		}
	}

	//The value was not found
	return false;
}

//Create an observable within this constraint, without specifying a random number generator
Observable * ObservableDiscreteConstraint::CreateObservable() const
{
	TRandom3 * random = new TRandom3(0);
	Observable * returnObservable = this->CreateObservable(random);
	delete random;
	return returnObservable;
}

//Create an observable within this constraint, using the specified random number generator
Observable * ObservableDiscreteConstraint::CreateObservable( TRandom3 * RandomNumberGenerator ) const
{
	int randomIndex = (int)floor( int(allValues.size()) * RandomNumberGenerator->Rndm() );
	return new Observable( name, allValues[unsigned(randomIndex)]/*, 0.0*/, unit );
}

bool ObservableDiscreteConstraint::IsDiscrete() const
{
	return true;
}

void ObservableDiscreteConstraint::Print() const
{
	cout << "Value: " << allValues[0];
	for( unsigned int i=1; i< allValues.size(); ++i )
	{
		cout << ",  " << allValues[i] ;
	}
	cout << "\tUnit: " << unit << endl;
}

//	This is intentionally BAD BY DESIGN
//	IT ONLY DOES WHAT IS WRITTEN ON THE TIN FOR TAG +1, -1 DISCRETE CONSTRAINT
//	IF THERE IS EVER A MORE COMPLEX USE CASE THEN THIS MUST BE RE-WRITTEN!!!!
Observable* ObservableDiscreteConstraint::GetMidRangeValue() const
{
	double middle_val = 0.;
	for( vector<double>::const_iterator val_i = allValues.begin(); val_i != allValues.end(); ++val_i )
	{
		middle_val+=*val_i;
	}
	return new Observable( name, middle_val, unit );
}

string ObservableDiscreteConstraint::GetTF1() const
{
	return tf1;
}

void ObservableDiscreteConstraint::SetTF1( const string input )
{
	tf1 = input;
}

string ObservableDiscreteConstraint::XML() const
{                           
	stringstream xml;   

	xml << "\t<Observable>" << endl;
	xml << "\t\t<Name>" << name << "</Name>" << endl;
	for( vector<double>::const_iterator val_i = allValues.begin(); val_i != allValues.end(); ++val_i )
	{
		xml << "\t\t<Value>" << *val_i << "</Value>" << endl;
	}
	xml << "\t\t<Unit>" << unit << "</Unit>" << endl;
	if( !tf1.empty() ) xml << "\t\t<TF1>" << tf1 << "</TF1>" << endl;
	xml << "\t</Observable>" << endl;                                                                                                                                                                                                    

	return xml.str();   
}

