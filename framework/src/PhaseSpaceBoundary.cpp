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
PhaseSpaceBoundary::PhaseSpaceBoundary( const vector<string> NewNames ) :
	allConstraints(), allNames(), DiscreteCombinationNumber(-1), uniqueID(0)
{
	allConstraints.reserve(NewNames.size());
	//Populate the map
	for (unsigned int nameIndex = 0; nameIndex < NewNames.size(); ++nameIndex)
	{
		allConstraints.push_back(NULL);
	}
	vector<string> duplicates;
	allNames = StringProcessing::RemoveDuplicates( NewNames, duplicates );
	if( allNames.size() != NewNames.size() )
	{
		cerr << "WARNING: PhaseSpace Cannot be created with multiple occurances of the same Observable!" << endl;
		for( vector<string>::iterator str_i = duplicates.begin(); str_i != duplicates.end(); ++str_i )
		{
			cout << *str_i << endl;
		}
	}

	uniqueID = reinterpret_cast<size_t>(this);
}

PhaseSpaceBoundary::PhaseSpaceBoundary( const PhaseSpaceBoundary& NewBoundary ) :
	allConstraints(), allNames( NewBoundary.allNames ), DiscreteCombinationNumber(NewBoundary.DiscreteCombinationNumber), uniqueID(0)
{
	for( unsigned int i=0; i< allNames.size(); ++i )
	{
		if( NewBoundary.allConstraints[i] != NULL )
		{
			if( NewBoundary.allConstraints[i]->IsDiscrete() )
			{
				allConstraints.push_back( new ObservableDiscreteConstraint( *((ObservableDiscreteConstraint*)NewBoundary.allConstraints[i]) ));
			}
			else
			{
				allConstraints.push_back( new ObservableContinuousConstraint( *((ObservableContinuousConstraint*)NewBoundary.allConstraints[i]) ));
			}
		}
		else
		{
			allConstraints.push_back(NULL);
		}
	}

	uniqueID = reinterpret_cast<size_t>(this)+1;
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

pair< vector<ObservableRef>, vector<double> > PhaseSpaceBoundary::GetDiscreteInfo( DataPoint* input ) const
{
	pair<vector<ObservableRef>, vector<double> > stored_pair;

	vector<ObservableRef> discreteInfoNames;
	for( unsigned int i=0; i< this->GetDiscreteNames().size(); ++i ) discreteInfoNames.push_back( ObservableRef(this->GetDiscreteNames()[i]) );
	vector<double> discreteInfoValues;
	for( unsigned int i=0; i< discreteInfoNames.size(); ++i )
	{
		ObservableDiscreteConstraint* thisDisc = (ObservableDiscreteConstraint*) this->GetConstraint( discreteInfoNames[i] );
		vector<double> possibleValues = thisDisc->GetValues();
		double thisValue = input->GetObservable( discreteInfoNames[i] )->GetValue();
		for( unsigned int j=0; j< possibleValues.size(); ++j )
		{
			if( fabs( possibleValues[j] - thisValue ) < 1E-5 )
			{
				discreteInfoValues.push_back( possibleValues[j] );
				break;
			}
		}
	}

	if( discreteInfoNames.size() == discreteInfoValues.size() )
	{
		stored_pair = make_pair( discreteInfoNames, discreteInfoValues );
	}
	else
	{
		stored_pair = make_pair( vector<ObservableRef>(), vector<double>() );
	}

	return stored_pair;
}

void PhaseSpaceBoundary::RemoveConstraint( string Name )
{
	int lookup = StringProcessing::VectorContains( &allNames, &Name );

	if( lookup != -1 )
	{
		vector<string>::iterator name_i = allNames.begin();
		vector<string>::iterator remove_name;
		vector<IConstraint*>::iterator const_i = allConstraints.begin();
		vector<IConstraint*>::iterator remove_const;
		for( ; name_i != allNames.end(); ++name_i, ++const_i )
		{
			if( Name == *name_i )
			{
				remove_name = name_i;
				remove_const = const_i;
				break;
			}
		}
		allNames.erase( remove_name );
		if( *remove_const != NULL ) delete *remove_const;
		allConstraints.erase( remove_const );
		++uniqueID;
	}
}

IConstraint * PhaseSpaceBoundary::GetConstraint( ObservableRef& object ) const
{
	if( object.GetExternalID() != uniqueID ) object.SetIndex(-1);

	if( object.GetIndex() < 0 ) {
		object.SetIndex( StringProcessing::VectorContains( &allNames, object.NameRef()) );
		object.SetExternalID( uniqueID );
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
		cerr << "Constraint on " << Name << " not found(1)" << endl;
		throw(3823);
		return NULL;//new ObservableContinuousConstraint( Name, 0.0, 0.0, "NameNotFoundError" );
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
		cerr << "Constraint on " << Name << " not found(2)" << endl;
		cerr << "If you wish to Add a NEW constraint, and know it may not yet exist, use:" << endl;
		cerr << "PhaseSpaceBoundary::AddConstraint()" << endl;
		cerr << endl;
		throw(81258);
		return false;
	}
	else
	{
		//Delete the old constraint before overwriting pointer
		if( allConstraints[unsigned(nameIndex)] != NULL ) delete allConstraints[unsigned(nameIndex)];
		if( NewConstraint->IsDiscrete() )
		{
			allConstraints[unsigned(nameIndex)] = new ObservableDiscreteConstraint( *((ObservableDiscreteConstraint*)NewConstraint) );
		}
		else
		{
			allConstraints[unsigned(nameIndex)] = new ObservableContinuousConstraint( *((ObservableContinuousConstraint*)NewConstraint) );
		}
		return true;
	}
}

//Initialise bound
bool PhaseSpaceBoundary::SetConstraint( string Name, double Minimum, double Maximum, string Unit )
{
	ObservableContinuousConstraint * newConstraint = new ObservableContinuousConstraint( Name, Minimum, Maximum, Unit );
	bool returnValue = this->SetConstraint( Name, newConstraint );
	delete newConstraint;
	return returnValue;
}

bool PhaseSpaceBoundary::SetConstraint( string Name, vector<double> Values, string Unit )
{
	ObservableDiscreteConstraint * newConstraint = new ObservableDiscreteConstraint( Name, Values, Unit );
	bool returnValue = this->SetConstraint( Name, newConstraint );
	delete newConstraint;
	return returnValue;
}

void PhaseSpaceBoundary::AddConstraint( string Name, IConstraint* NewConstraint, bool overwrite )
{
	int lookup = StringProcessing::VectorContains( &allNames, &Name );
	if( lookup == -1 )
	{
		allNames.push_back( Name );
		allConstraints.push_back( NULL );
		this->SetConstraint( Name, NewConstraint );
	}
	else if( allConstraints[(unsigned)lookup] == NULL )
	{
		this->SetConstraint( Name, NewConstraint );
	}
	else
	{
		if( overwrite ) this->SetConstraint( Name, NewConstraint );
	}
}

//Returns whether a point is within the boundary
bool PhaseSpaceBoundary::IsPointInBoundary( DataPoint * TestDataPoint, bool silence )
{
	for (unsigned int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex )
	{
		//Check if test Observable exists in the DataPoint
		Observable * testObservable = TestDataPoint->GetObservable( allNames[nameIndex], silence );
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
						Observable* temp = new Observable( testObservable->GetName(), temp_vec[i], testObservable->GetUnit() );
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
	for( unsigned int i=0; i< allNames.size(); ++i )
	{
		cout << "Name: " << allNames[i] << "\t";
		if( allConstraints[i] != NULL ) allConstraints[i]->Print();
		else cout << "NULL" << endl;
	}
	cout << endl;
}


DataPoint* PhaseSpaceBoundary::GetMidPoint() const
{
	DataPoint* returnable_point = new DataPoint( allNames );
	vector<string>::const_iterator name_i = allNames.begin();
	for( vector<IConstraint*>::const_iterator const_i = allConstraints.begin(); const_i != allConstraints.end(); ++const_i, ++name_i )
	{
		Observable* thisMiddleObservable=NULL;
		if( (*const_i) !=NULL )
		{
			thisMiddleObservable = (*const_i)->GetMidRangeValue();
		}
		else
		{
			delete returnable_point;
			return NULL;
		}
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
	vector<string> thisAllNames = this->GetAllNames();
	vector<vector<double> > discreteValues;
	vector<string> discreteNames, continuousNames;
	vector<vector<double> > discreteCombinations = StatisticsFunctions::DiscreteCombinations( &thisAllNames, this, discreteNames, continuousNames, discreteValues );

	(void) continuousNames;

	//Create the data points to return
	vector<DataPoint*> newDataPoints;

	DataPoint* tempPoint = this->GetMidPoint();
	if( tempPoint == NULL ) return newDataPoints;

	for( unsigned int combinationIndex = 0; combinationIndex < discreteCombinations.size(); ++combinationIndex )
	{
		DataPoint* templateDataPoint = new DataPoint( *tempPoint );

		//Output the discrete values for this combination
		for( unsigned int discreteIndex = 0; discreteIndex < discreteNames.size(); ++discreteIndex )
		{
			//Set the data point
			Observable* oldValue = templateDataPoint->GetObservable( discreteNames[discreteIndex] );

			Observable* newValue = new Observable( oldValue->GetName(), discreteCombinations[combinationIndex][discreteIndex], oldValue->GetUnit() );

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

unsigned int PhaseSpaceBoundary::GetDiscreteIndex( DataPoint* Input, const bool silence ) const
{
	(void) silence;
	//	Exit on simple case
	int thisIndex = Input->GetDiscreteIndex();
	size_t thisID = Input->GetDiscreteIndexID();
	if( thisIndex != -1 && thisID == uniqueID ) return (unsigned)thisIndex;

	if( this->GetDiscreteNames().empty() || (this->GetDiscreteNames().size() == 1) )
	{
		Input->SetDiscreteIndex( 0 );
		return 0;
	}

	/*
	if( this->GetConstraint(this->GetDiscreteNames()[0])->GetValues().size() > 1 )
	{
		cout << "This DataPoint is:" << endl;
		Input->Print();
		cout << this->GetConstraint(this->GetDiscreteNames()[0])->GetValues().size() << endl;
	}
	*/

	//	Get all possible discrete combination datapoints and the names of all discrete observables
	vector<DataPoint*> allCombinations = this->GetDiscreteCombinations();


	vector<string> allDiscreteNames = this->GetDiscreteNames();

	//	Construct array of ObservableRef objects to pick out Discrete Observables
	vector<string>::iterator name_i = allDiscreteNames.begin();
	vector<string>::iterator end_name_i = allDiscreteNames.end();
	vector<ObservableRef> allDiscreteObs;
	for( ; name_i != end_name_i; ++name_i ) allDiscreteObs.push_back( ObservableRef( *name_i ) );

	//	Initialize iterators
	vector<ObservableRef>::iterator start_Obsname_i = allDiscreteObs.begin();
	vector<ObservableRef>::iterator Obsname_i = start_Obsname_i;
	vector<ObservableRef>::iterator end_Obsname_i = allDiscreteObs.end();
	//vector<DataPoint*>::iterator start_comb_i = allCombinations.begin();
	vector<DataPoint*>::iterator comb_i = allCombinations.begin();
	vector<DataPoint*>::iterator end_comb_i = allCombinations.end();

	//	Construct objects
	bool match=true;
	double wanted_val=0., this_val=0.;
	string thisName;
	//	Loop over all possible combinations
	for( int index=0 ; comb_i != end_comb_i; ++comb_i, ++index )
	{
		//	Check if this datapoint is the same as this combination
		match = true;
		for( Obsname_i = start_Obsname_i; Obsname_i != end_Obsname_i; ++Obsname_i )
		{
			wanted_val = (*comb_i)->GetObservable( *Obsname_i, true )->GetValue();
			this_val = Input->GetObservable( *Obsname_i, true )->GetValue();
			//if( this->GetConstraint(this->GetDiscreteNames()[0])->GetValues().size() > 1 )
			//{
			//	cout << "wanted: " << wanted_val << endl;
			//	cout << "this_val: " << this_val << endl;
			//}

			//	Stop checking once we have at least one discrtete observable different
			if( fabs( wanted_val - this_val ) > 1E-6 )
			{
				match = false;
				thisName=string(*Obsname_i);
				break;
			}
		}
		//	If this combination matches set the index and leave
		if( match )
		{
			thisIndex = index;
			thisName="end";
			break;
		}
	}

	//if( this->GetConstraint(this->GetDiscreteNames()[0])->GetValues().size() > 1 )
	//{
	//	cout << "This Index is: " << thisIndex << endl;
	//	//exit(-520);
	//}

	//	Destory temporary objects
	while( !allCombinations.empty() )
	{
		if( allCombinations.back() != NULL ) delete allCombinations.back();
		allCombinations.pop_back();
	}

	//	Check for error
	if( thisIndex == -1 )
	{
		cerr << thisName << ":" << endl;
		cerr << "wanted_val: " << wanted_val << endl;
		cerr << "this_val: " << this_val << endl;
		cerr << "This DataPoint does not Lie within this PhaseSpace." << endl;
		cerr << "This is a SERIOUS MISCONFIGURATION!!! Exiting :(" << endl;
		Input->Print();
		this->Print();
		//exit(-99865);
		DebugClass::SegFault();
	}

	//	Store and return the lookup of the DataPoint's index
	Input->SetDiscreteIndex( thisIndex );
	Input->SetDiscreteIndexID( uniqueID );
	return (unsigned)thisIndex;
}


int PhaseSpaceBoundary::GetNumberCombinations() const
{
	if( DiscreteCombinationNumber != -1 ) return DiscreteCombinationNumber;

	vector<DataPoint*> thisMany = this->GetDiscreteCombinations();

	if( thisMany.empty() )
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

void PhaseSpaceBoundary::CheckPhaseSpace( IPDF* toCheck ) const
{
	cout << endl;
	vector<string> needed = toCheck->GetPrototypeDataPoint();

	for( unsigned int i=0; i< needed.size(); ++i )
	{
		int lookup = StringProcessing::VectorContains( &allNames, &(needed[i]) );
		if( lookup == -1 )
		{
			cout << endl;
			cout << "ERROR: Missing Constraint on: " << needed[i] << endl;
			cout << endl;
		}
		else
		{
			cout << "Constraint Found:\t";
			allConstraints[i]->Print();
		}
	}
}

size_t PhaseSpaceBoundary::GetID() const
{
	return uniqueID;
}

