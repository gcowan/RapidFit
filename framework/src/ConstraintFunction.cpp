/*!
 * @class ConstraintFunction
 *
 * Where external, experimental constraints on PhysicsParameters are calculated
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern,ch
 */

//	RapidFit Headers
#include "ConstraintFunction.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <vector>
#include <sstream>

using namespace::std;

//Constructor with correct arguments
ConstraintFunction::ConstraintFunction( const vector< IConstraintFunction* > NewConstraints ) : Found_Position(), allConstraints()
{
	for( vector<IConstraintFunction*>::const_iterator const_i = NewConstraints.begin(); const_i != NewConstraints.end(); ++const_i )
	{
		if( (*const_i)->isExternalConst() == true )
		{
			allConstraints.push_back( (IConstraintFunction*) new ExternalConstraint( *((ExternalConstraint*)(*const_i)) ) );
		}
		else if( (*const_i)->isExternalMatrix() == true )
		{
			allConstraints.push_back( (IConstraintFunction*)new ExternalConstMatrix( *(ExternalConstMatrix*)(*const_i) ) );
		}
		cout << "Adding Constraint: " << allConstraints.back()->GetName() << endl;
	}
}

ConstraintFunction::ConstraintFunction( const ConstraintFunction& input ) : Found_Position( input.Found_Position ), allConstraints()
{
	for( vector<IConstraintFunction*>::const_iterator const_i = input.allConstraints.begin(); const_i != input.allConstraints.end(); ++const_i )
	{
		if( (*const_i)->isExternalConst() ) allConstraints.push_back( (IConstraintFunction*)new ExternalConstraint( *(ExternalConstraint*)(*const_i) ) );
		else if( (*const_i)->isExternalMatrix() ) allConstraints.push_back( (IConstraintFunction*)new ExternalConstMatrix( *(ExternalConstMatrix*)(*const_i) ) );
	}
}


//Destructor
ConstraintFunction::~ConstraintFunction()
{
	while( !allConstraints.empty() )
	{
		if( allConstraints.back() != NULL ) delete allConstraints.back();
		allConstraints.pop_back();
	}
}

//Perform the constraint calculation
double ConstraintFunction::Evaluate( const ParameterSet * NewParameters )
{
	vector<string> parameterNames = NewParameters->GetAllNames();
	double constraintValue = 0.0;

	if( Found_Position.size() != allConstraints.size() )
	{
		for( unsigned int constraintIndex = 0; constraintIndex < allConstraints.size(); ++constraintIndex )
		{
			string name = allConstraints[constraintIndex]->GetName();
			Found_Position.push_back( StringProcessing::VectorContains( &parameterNames, &name ) );
		}
	}

	unsigned int num_i=0;
	vector< vector<IConstraintFunction*>::iterator > to_be_removed;
	//Loop over all ExternalConstraints
	for( vector<IConstraintFunction*>::iterator constraintIndex = allConstraints.begin(); constraintIndex != allConstraints.end(); ++constraintIndex, ++num_i )
	{
		bool canApply = (*constraintIndex)->CanApply( NewParameters );
		if( !canApply )
		{
			to_be_removed.push_back( constraintIndex );
			continue;
		}
		else
		{
			(*constraintIndex)->SetPhysicsParameters( NewParameters );
			constraintValue += (*constraintIndex)->GetChi2();
		}
	}

	for( unsigned int i=0; i< to_be_removed.size(); ++i )
	{
		allConstraints.erase( to_be_removed[i] );
	}

	return 0.5 * constraintValue;
}

void ConstraintFunction::Print() const
{
	cout << "ConstraintFunction:" << endl;
	cout << "This is to be coded up when needed" << endl;
}

string ConstraintFunction::XML() const
{
	stringstream xml;

	for( vector<IConstraintFunction*>::const_iterator const_i = allConstraints.begin(); const_i != allConstraints.end(); ++const_i )
	{
		xml << (*const_i)->XML() << endl;
		xml << endl;
	}

	return xml.str();
}

vector<string> ConstraintFunction::ConstrainedParameter() const
{
	vector<string> constnames;
	for( unsigned int i=0; i< allConstraints.size(); ++i )
	{
		constnames = StringProcessing::CombineUniques( constnames, vector<string>(1,allConstraints[i]->GetName()) );
	}
	return constnames;
}

