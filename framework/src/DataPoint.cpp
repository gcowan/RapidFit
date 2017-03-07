/**
  @class DataPoint

  Holds all observables for a given event

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

#include "TRandom.h"

//	RapidFit Headers
#include "DataPoint.h"
#include "ObservableRef.h"
#include "StringProcessing.h"
#include "DebugClass.h"
//	System Headers
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <limits>

using namespace::std;

//	Required for Sorting
DataPoint::DataPoint() : allObservables(), allNames(), myPhaseSpaceBoundary(NULL), thisDiscreteIndex(-1),
	WeightValue(1.), storedID(0), initialNLL( numeric_limits<double>::quiet_NaN() ), PerEventData(), nameIndex(), DiscreteIndexMap()
{
}

//Constructor with correct arguments
DataPoint::DataPoint( vector<string> NewNames ) : allObservables(), allNames(), myPhaseSpaceBoundary(NULL),
	thisDiscreteIndex(-1), WeightValue(1.), storedID(0), initialNLL( numeric_limits<double>::quiet_NaN() ),
	PerEventData(), nameIndex(), DiscreteIndexMap()
{
	allObservables.reserve( NewNames.size() );
	//Populate the map
	for( nameIndex = 0; nameIndex < (int)NewNames.size(); ++nameIndex )
	{
		allObservables.push_back( Observable(NewNames[(unsigned)nameIndex]) );
	}
	vector<string> duplicates;
	allNames = StringProcessing::RemoveDuplicates( NewNames, duplicates );
	if( allNames.size() != NewNames.size() )
	{
		cerr << "WARNING: Cannot Generate a DataPoint with 2 Occurances of the same Observable" << endl;
		for( vector<string>::iterator str_i = duplicates.begin(); str_i != duplicates.end(); ++str_i )
		{
			cout << *str_i << endl;
		}
		cerr << "This is harmless, but you will now have some merged Observable(s)" << endl;
	}
}

//Assignment operator
DataPoint& DataPoint::operator=( const DataPoint& NewPoint )
{
	if( &NewPoint != this )
	{
		this->allNames = NewPoint.allNames;
		this->myPhaseSpaceBoundary = NewPoint.myPhaseSpaceBoundary;
		this->thisDiscreteIndex = NewPoint.thisDiscreteIndex;
		this->WeightValue = NewPoint.WeightValue;
		this->storedID = NewPoint.storedID;
		this->initialNLL = NewPoint.initialNLL;
		this->PerEventData = NewPoint.PerEventData;
		this->nameIndex = NewPoint.nameIndex;
		for( unsigned int i=0; i< NewPoint.allObservables.size(); ++i )
		{
			this->allObservables.push_back( Observable( (NewPoint.allObservables[i]) ) );
		}
		this->DiscreteIndexMap = NewPoint.DiscreteIndexMap;
	}
	return *(this);
}

DataPoint::DataPoint( const DataPoint& input ) :
	allObservables(), allNames(input.allNames), myPhaseSpaceBoundary(input.myPhaseSpaceBoundary),
	thisDiscreteIndex(input.thisDiscreteIndex), WeightValue(input.WeightValue), storedID(input.storedID),
	initialNLL( input.initialNLL ), PerEventData(input.PerEventData), nameIndex(), DiscreteIndexMap(input.DiscreteIndexMap)
{
	for( unsigned int i=0; i< input.allObservables.size(); ++i )
	{
		allObservables.push_back( Observable( (input.allObservables[i])) );
	}
}

//Destructor
DataPoint::~DataPoint()
{
}

//Retrieve names of all observables stored
vector<string> DataPoint::GetAllNames() const
{
	return allNames;
}

void DataPoint::RemoveObservable( const string input )
{
	vector<string>::iterator name_i = allNames.begin();
	vector<Observable>::iterator obs_i = allObservables.begin();

	vector<string>::iterator name_to_remove;
	vector<Observable>::iterator observable_to_remove;

	for( ; name_i != allNames.end(); ++name_i, ++obs_i )
	{
		if( (*name_i) == input )
		{
			name_to_remove = name_i;
			observable_to_remove = obs_i;
			break;
		}
	}

	allNames.erase( name_to_remove );
	allObservables.erase( observable_to_remove );
}

Observable* DataPoint::GetObservable( unsigned int wanted )
{
	return &(allObservables[ wanted ]);
}

//Retrieve an observable by its name
//	!!!THIS IS VERY, VERY, VERY WASTEFUL FOR LARGE DATASETS!!!
Observable* DataPoint::GetObservable(string const Name, const bool silence )
{
	//Check if the name is stored in the map
	nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if( nameIndex == -1 )
	{
		if( !silence ) cerr << "Observable name " << Name << " not found (2)" << endl;
		throw(-1543);
	}
	else
	{
		return &(allObservables[(unsigned)nameIndex]);
	}
}

Observable* DataPoint::GetObservable( const ObservableRef& object, const bool silence )
{
	if( object.GetIndex() < 0 )
	{
		object.SetIndex( StringProcessing::VectorContains( &allNames, object.NameRef()) );
		if( object.GetIndex() >= 0 ) return &(allObservables[ (unsigned) object.GetIndex() ]);
	}
	else
	{
		return &(allObservables[ (unsigned) object.GetIndex() ]);
	}
	if( !silence ) cerr << "Observable name " << object.Name().c_str() << " not found (3)" << endl;
	throw(-20);
}

//Set an observable by name
bool DataPoint::SetObservable( string Name, Observable * NewObservable )
{
	//Check if the name is stored in the map
	nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if ( nameIndex == -1 )
	{
		cerr << "Observable name " << Name << " not found (4)" << endl;
		throw(438);
	}
	else
	{
		allObservables[(unsigned)nameIndex].SetObservable(NewObservable);
		return true;
	}
}

bool DataPoint::SetObservable( ObservableRef& Name, Observable * NewObservable )
{
	//Check if the name is stored in the map
	if( Name.GetIndex() == -1 )
	{
		string thisName = Name.Name();
		nameIndex = StringProcessing::VectorContains( &allNames, &thisName );
		if( nameIndex == -1 )
		{
			cerr << "Observable name " << thisName << " not found (5)" << endl;
			throw(4389);
		}
		else
		{
			Name.SetIndex( nameIndex );
			allObservables[(unsigned)nameIndex].SetObservable(NewObservable);
			return true;
		}
		//return false;
	}
	else
	{
		allObservables[(unsigned)Name.GetIndex()].SetObservable(NewObservable);
		return true;
	}
}

void DataPoint::AddObservable( string Name, Observable* NewObservable )
{
	if( StringProcessing::VectorContains( &allNames, &Name ) == -1 )
	{
		allNames.push_back( Name );
		allObservables.push_back( Observable(*NewObservable) );
	}
	else
	{
		this->SetObservable( Name, NewObservable );
	}
}

void DataPoint::AddObservable( string Name, double Value, string Unit, bool trusted, int thisnameIndex )
{
	Observable *tempObservable = new Observable( Name, Value, Unit );
	if( trusted )
	{
		allObservables[(unsigned)thisnameIndex].SetObservable( tempObservable );
	}
	else
	{
		this->AddObservable( Name, tempObservable );
	}
	delete tempObservable;
}

//Initialise observable
bool DataPoint::SetObservable( string Name, double Value, string Unit, bool trusted, int thisnameIndex )
{
	Observable * temporaryObservable = new Observable( Name, Value, Unit );
	bool returnValue=false;
	if( trusted )
	{
		returnValue=true;
		allObservables[(unsigned)thisnameIndex].SetObservable( temporaryObservable );
	}
	else
	{
		returnValue = SetObservable( Name, temporaryObservable );
	}
	delete temporaryObservable;
	return returnValue;
}

//	Used for Sorting DataPoints
bool DataPoint::operator() ( pair<DataPoint,ObservableRef> first, pair<DataPoint,ObservableRef> second )
{
	double param_val_1 = first.first.GetObservable( first.second )->GetValue();
	double param_val_2 = second.first.GetObservable( second.second )->GetValue();
	return (param_val_1 < param_val_2 );
}

void DataPoint::Print() const
{
	cout << "DataPoint:" << endl;
	for( unsigned int i=0; i< allObservables.size(); ++i )
	{
		//cout << "Obs: " << allObservables[i] << " ";
		allObservables[i].Print();
	}

	cout << endl;
}

PhaseSpaceBoundary* DataPoint::GetPhaseSpaceBoundary() const
{
	return myPhaseSpaceBoundary;
}

void DataPoint::SetPhaseSpaceBoundary( PhaseSpaceBoundary* input )
{
	myPhaseSpaceBoundary = input;
}

int DataPoint::GetDiscreteIndex() const
{
	return thisDiscreteIndex;
}

void DataPoint::SetDiscreteIndex( int Input )
{
	thisDiscreteIndex = Input;
}

double DataPoint::GetEventWeight() const
{
	return WeightValue;
}

void DataPoint::SetEventWeight( const double Input )
{
	WeightValue = Input;
}

void DataPoint::SetDiscreteIndexID( size_t thisID )
{
	storedID = thisID;
}

size_t DataPoint::GetDiscreteIndexID() const
{
	return storedID;
}

void DataPoint::SetInitialNLL( const double input )
{
	initialNLL = input;
}

double DataPoint::GetInitialNLL() const
{
	return initialNLL;
}

vector<double> DataPoint::GetPerEventData() const
{
	return PerEventData;
}

void DataPoint::SetPerEventData( const vector<double> input )
{
	PerEventData = input;
}

void DataPoint::ClearPerEventData()
{
	PerEventData.clear();
}

void DataPoint::SetDiscreteIndexIDMap( size_t thisID, int index )
{
	DiscreteIndexMap.insert( pair<size_t,int>(thisID, index) );
}

bool DataPoint::FindDiscreteIndexID( size_t thisID )
{
	return (DiscreteIndexMap.find( thisID ) != DiscreteIndexMap.end() );
}

int DataPoint::GetDiscreteIndexMap( size_t thisID )
{
	if( this->FindDiscreteIndexID( thisID ) )
	{
		return DiscreteIndexMap.at( thisID );
	}
	else
	{
		return -1;
	}
}

void DataPoint::ClearDiscreteIndexMap()
{
	DiscreteIndexMap.clear();
}

