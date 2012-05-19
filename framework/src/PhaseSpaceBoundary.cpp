/**
  @class PhaseSpaceBoundary

  A collection of constraints on observables, defining the phase space in which a data point exists

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

//	RapidFit Headers
#include "StringProcessing.h"
#include "StatisticsFunctions.h"
#include "PhaseSpaceBoundary.h"
#include "ObservableContinuousConstraint.h"
#include "ObservableDiscreteConstraint.h"
#include "ObservableRef.h"
//	System Headers
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <cmath>

#define DOUBLE_TOLERANCE_PHASE 1E-6


//Constructor with correct arguments
PhaseSpaceBoundary::PhaseSpaceBoundary( vector<string> NewNames ) :
	allConstraints(), allNames(NewNames), DiscreteCombinationNumber(-1)
{
	allConstraints.reserve(NewNames.size());
	//Populate the map
	for (unsigned int nameIndex = 0; nameIndex < NewNames.size(); ++nameIndex)
	{
		allConstraints.push_back(NULL);
	}
}

PhaseSpaceBoundary::PhaseSpaceBoundary( const PhaseSpaceBoundary& NewBoundary ) :
	allConstraints(), allNames( NewBoundary.allNames ), DiscreteCombinationNumber(NewBoundary.DiscreteCombinationNumber)
{
	allConstraints.reserve( NewBoundary.allConstraints.size() );
	for( unsigned int i=0; i< NewBoundary.allConstraints.size(); ++i )
	{
		if( NewBoundary.allConstraints[i]->IsDiscrete() )
		{
			allConstraints.push_back( new ObservableDiscreteConstraint( string(NewBoundary.allNames[i]),
						vector<double>(NewBoundary.allConstraints[i]->GetValues()), string(NewBoundary.allConstraints[i]->GetUnit()),
						string(NewBoundary.allConstraints[i]->GetTF1()) ) );
		}
		else
		{
			allConstraints.push_back( new ObservableContinuousConstraint( string(NewBoundary.allNames[i]),
						double(NewBoundary.allConstraints[i]->GetMinimum()), double(NewBoundary.allConstraints[i]->GetMaximum()),
						string(NewBoundary.allConstraints[i]->GetUnit()), string(NewBoundary.allConstraints[i]->GetTF1()) ) );
		}
	}
}

//Destructor
PhaseSpaceBoundary::~PhaseSpaceBoundary()
{
}

//Return the names of all bounds stored
vector<string> PhaseSpaceBoundary::GetAllNames() const
{
	return allNames;
}

vector<string> PhaseSpaceBoundary::GetDiscreteNames() const
{
	vector<string> disc_names;

	vector<string>::const_iterator name_i = allNames.begin();
	vector<IConstraint*>::const_iterator constr_i = allConstraints.begin();
	for( ; name_i != allNames.end(); ++name_i, ++constr_i )
	{
		if( (*constr_i)->IsDiscrete() == true ) disc_names.push_back( (*name_i) );
	}
	return disc_names;
}

vector<string> PhaseSpaceBoundary::GetContinuousNames() const             
{
	vector<string> cont_names;

	vector<string>::const_iterator name_i = allNames.begin();
	vector<IConstraint*>::const_iterator constr_i = allConstraints.begin();
	for( ; name_i != allNames.end(); ++name_i, ++constr_i )
	{
		if( (*constr_i)->IsDiscrete() == false ) cont_names.push_back( (*name_i) );
	}
	return cont_names;
}

//Retrieve a constraint by its name
IConstraint * PhaseSpaceBoundary::GetConstraint( pair<string,int>* wanted_constraint )
{
	if( wanted_constraint->second != -1 )
	{
		return allConstraints[unsigned(wanted_constraint->second)];
	} else {
		wanted_constraint->second = StringProcessing::VectorContains( &allNames, &(wanted_constraint->first) );
	}
	if( wanted_constraint->second == -1 ){
		cerr << "PhysicsParameter " << wanted_constraint->first << " not found" <<endl;
	}else{
		return allConstraints[unsigned(wanted_constraint->second)];}
	exit(-1);
}

IConstraint * PhaseSpaceBoundary::GetConstraint( ObservableRef& object ) const
{
	if( object.GetIndex() < 0 ) {
		object.SetIndex( StringProcessing::VectorContains( &allNames, object.NameRef()) );
		if( object.GetIndex() >= 0 ) return allConstraints[ (unsigned) object.GetIndex() ];
	} else {
		return allConstraints[ (unsigned) object.GetIndex() ];
	}
	cerr << "Observable name " << object.Name().c_str() << " not found (PhaseSpaceBoundary)" << endl;
	throw(-20);
}

//Retrieve a constraint by its name
IConstraint * PhaseSpaceBoundary::GetConstraint(string Name) const
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if ( nameIndex == -1 )
	{
		cerr << "Constraint on " << Name << " not found" << endl;
		exit(1);
		//return new ObservableContinuousConstraint( Name, 0.0, 0.0, "NameNotFoundError" );
	}
	else
	{
		return allConstraints[unsigned(nameIndex)];
	}
}

//Set a constraint by name
bool PhaseSpaceBoundary::SetConstraint( string Name, IConstraint * NewConstraint )
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if ( nameIndex == -1 )
	{
		cerr << "Constraint on " << Name << " not found" << endl;
		exit(1);
		//return false;
	}
	else
	{
		//Delete the old constraint before overwriting pointer
		if( allConstraints[unsigned(nameIndex)] != NULL ) delete allConstraints[unsigned(nameIndex)];
		allConstraints[unsigned(nameIndex)] = NewConstraint;
		return true;
	}
}

//Initialise bound
bool PhaseSpaceBoundary::SetConstraint( string Name, double Minimum, double Maximum, string Unit )
{
	ObservableContinuousConstraint * newConstraint = new ObservableContinuousConstraint( Name, Minimum, Maximum, Unit );
	bool returnValue = SetConstraint( Name, newConstraint );
	return returnValue;
}

bool PhaseSpaceBoundary::SetConstraint( string Name, vector<double> Values, string Unit )
{
	ObservableDiscreteConstraint * newConstraint = new ObservableDiscreteConstraint( Name, Values, Unit );
	bool returnValue = SetConstraint( Name, newConstraint );
	return returnValue;
}

void PhaseSpaceBoundary::AddConstraint( string Name, IConstraint* NewConstraint )
{
	if( StringProcessing::VectorContains( &allNames, &Name ) == -1 )
	{
		allNames.push_back( Name );
		allConstraints.push_back( NewConstraint );
	}
	else
	{
		SetConstraint( Name, NewConstraint );
	}
}

//Returns whether a point is within the boundary
bool PhaseSpaceBoundary::IsPointInBoundary( DataPoint * TestDataPoint )
{
	for (unsigned int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex )
	{
		//Check if test Observable exists in the DataPoint
		Observable * testObservable = TestDataPoint->GetObservable( allNames[nameIndex] );
		if ( testObservable->GetUnit() == "NameNotFoundError" )
		{
			cerr << "Observable \"" << allNames[nameIndex] << "\" expected but not found" << endl;
			return false;
		}
		else
		{
			//Check if the Observable fits
			if ( !allConstraints[nameIndex]->CheckObservable(testObservable) )
			{
				return false;
			}
			else if( allConstraints[nameIndex]->IsDiscrete() )
			{
				vector<double> temp_vec = allConstraints[nameIndex]->GetValues();
				for( unsigned int i=0; i< temp_vec.size(); ++i )
				{
					if( fabs( testObservable->GetValue() - temp_vec[i] ) < DOUBLE_TOLERANCE_PHASE )
					{
						Observable* temp = new Observable( testObservable->GetName(), temp_vec[i], testObservable->GetError(), testObservable->GetUnit() );
						TestDataPoint->SetObservable( testObservable->GetName(), temp );
						delete temp;
					}
				}
			}
		}
	}

	//Point is within the boundary
	return true;
}

void PhaseSpaceBoundary::Print() const
{
	cout << "PhaseSpaceBoundary:" << endl;
	for( unsigned int i=0; i< allConstraints.size(); ++i )
	{
		cout << "Name: " << allNames[i] << "\t";
		allConstraints[i]->Print();
	}
	cout << endl;
}


DataPoint* PhaseSpaceBoundary::GetMidPoint() const
{
	DataPoint* returnable_point = new DataPoint( allNames );
	vector<string>::const_iterator name_i = allNames.begin();
	for( vector<IConstraint*>::const_iterator const_i = allConstraints.begin(); const_i != allConstraints.end(); ++const_i, ++name_i )
	{
		Observable* thisMiddleObservable = (*const_i)->GetMidRangeValue();
		returnable_point->AddObservable( *name_i, thisMiddleObservable );
		delete thisMiddleObservable;
	}
	return returnable_point;
}


string PhaseSpaceBoundary::DiscreteDescription( const DataPoint* input ) const
{
	stringstream description;

	vector<string> desc_observables = this->GetDiscreteNames();

	for( vector<string>::iterator obs_i = desc_observables.begin(); obs_i != desc_observables.end(); ++obs_i )
	{
		double value = input->GetObservable( *obs_i )->GetValue();
		description << *obs_i << " : " << value << "\t";
	}

	description << endl;
	return description.str();
}


//Return a list of data points
//Each should take the data average value of each continuous observable
//Each should represent one combination of possible discrete values
vector<DataPoint*> PhaseSpaceBoundary::GetDiscreteCombinations() const
{
	//Calculate all possible combinations of discrete observables
	vector<string> allNames = this->GetAllNames();
	vector<vector<double> > discreteValues;
	vector<string> discreteNames, continuousNames;
	vector<vector<double> > discreteCombinations = StatisticsFunctions::DiscreteCombinations( &allNames, this, discreteNames, continuousNames, discreteValues );

	(void) continuousNames;

	//Create the data points to return
	vector<DataPoint*> newDataPoints;

	DataPoint* tempPoint = this->GetMidPoint();

	for( unsigned int combinationIndex = 0; combinationIndex < discreteCombinations.size(); ++combinationIndex )
	{
		DataPoint* templateDataPoint = new DataPoint( *tempPoint );

		//Output the discrete values for this combination
		for( unsigned int discreteIndex = 0; discreteIndex < discreteNames.size(); ++discreteIndex )
		{
			//Set the data point
			Observable* oldValue = templateDataPoint->GetObservable( discreteNames[discreteIndex] );

			Observable* newValue = new Observable( oldValue->GetName(), discreteCombinations[combinationIndex][discreteIndex],
					0, oldValue->GetUnit() );

			templateDataPoint->SetObservable( discreteNames[discreteIndex], newValue );
			delete newValue;
		}

		newDataPoints.push_back( templateDataPoint );
	}

	delete tempPoint;

	return newDataPoints;
}

string PhaseSpaceBoundary::XML() const
{
	stringstream xml;

	xml << "\t<PhaseSpaceBoundary>" << endl;
	vector<string>::const_iterator name_i = allNames.begin();
	vector<IConstraint*>::const_iterator const_i = allConstraints.begin();
	for( ; name_i != allNames.end(); ++name_i, ++const_i )
	{
		xml << endl;
		if( (*const_i) == NULL )
		{
			xml << "\t<Observable>" << endl;
			xml << "\t\t<Name>" << *name_i << "</Name>" << endl;
			xml << "\t\t<Minimum>" << "minimum" << "</Minimum>" << endl;
			xml << "\t\t<Maximum>" << "maximum" << "</Maximum>" << endl;
			xml << "\t\t<Value>" << "Value1" << "</Value>" << endl;
			xml << "\t\t<Unit>" << "someUnit" << "</Unit>" << endl;
			xml << "\t</Observable>" << endl;
		}
		else
		{
			xml << (*const_i)->XML();
		}
	}
	xml << endl;
	xml << "\t</PhaseSpaceBoundary>" << endl;

	return xml.str();
}

unsigned int PhaseSpaceBoundary::GetDiscreteIndex( DataPoint* Input ) const
{
	int thisIndex = Input->GetDiscreteIndex();
	if( thisIndex != -1 ) return (unsigned)thisIndex;

	if( this->GetDiscreteNames().empty() ) return 0;

	vector<DataPoint*> allCombinations = this->GetDiscreteCombinations();

	vector<string> allDiscreteNames = this->GetDiscreteNames();

	vector<string>::iterator name_i = allDiscreteNames.begin();
	vector<DataPoint*>::iterator comb_i = allCombinations.begin();

	for( int index=0 ; comb_i != allCombinations.end(); ++comb_i, ++index )
	{
		bool match = true;
		for( ; name_i != allDiscreteNames.end(); ++name_i )
		{
			double wanted_val = (*comb_i)->GetObservable( *name_i )->GetValue();
			double this_val = Input->GetObservable( *name_i )->GetValue();
			if( fabs( wanted_val - this_val ) > 1E-6 )
			{
				match = false;
				break;
			}
		}
		if( match )
		{
			thisIndex = index;
			break;
		}
	}

	while( !allCombinations.empty() )
	{
		if( allCombinations.back() != NULL ) delete allCombinations.back();
		allCombinations.pop_back();
	}

	if( thisIndex == -1 )
	{
		cerr << "This DataPoint does not Lie within this PhaseSpace." << endl;
		cerr << "This is a SERIOUS MISCONFIGURATION!!! Exiting :(" << endl;
		Input->Print();
		this->Print();
		exit(-99865);
	}

	Input->SetDiscreteIndex( thisIndex );
	return (unsigned)thisIndex;
}


int PhaseSpaceBoundary::GetNumberCombinations() const
{
	if( DiscreteCombinationNumber != -1 ) return DiscreteCombinationNumber;

	vector<DataPoint*> thisMany = this->GetDiscreteCombinations();

	if( thisMany.empty() || (thisMany.size() == 0) )
	{
		DiscreteCombinationNumber = 1;
	}
	else
	{
		DiscreteCombinationNumber = (int)thisMany.size();

		while( !thisMany.empty() )
		{
			if( thisMany.back() != NULL ) delete thisMany.back();
			thisMany.pop_back();
		}
	}

	return DiscreteCombinationNumber;
}

