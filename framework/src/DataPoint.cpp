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
#include "PseudoObservable.h"
//	System Headers
#include <iostream>
#include <stdlib.h>
#include <iomanip>

using namespace::std;

//	Required for Sorting
DataPoint::DataPoint() : allObservables(), allNames(), allPseudoNames(), allPseudoObservables(), myPhaseSpaceBoundary(NULL), thisDiscreteIndex(-1), WeightValue(1.), storedID(0)
{
}

//Constructor with correct arguments
DataPoint::DataPoint( vector<string> NewNames ) : allObservables(), allNames(), allPseudoNames(), allPseudoObservables(), myPhaseSpaceBoundary(NULL), thisDiscreteIndex(-1), WeightValue(1.), storedID(0)
{
	allObservables.reserve( NewNames.size() );
	//Populate the map
	for (unsigned short int nameIndex = 0; nameIndex < NewNames.size(); nameIndex++)
	{
		allObservables.push_back( new Observable(NewNames[nameIndex]) );
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

DataPoint::DataPoint( const DataPoint& input ) :
	allObservables(), allNames(input.allNames), allPseudoNames(input.allPseudoNames), allPseudoObservables(input.allPseudoObservables), myPhaseSpaceBoundary(input.myPhaseSpaceBoundary),
	thisDiscreteIndex(input.thisDiscreteIndex), WeightValue(input.WeightValue), storedID(input.storedID)
{
	for( unsigned int i=0; i< input.allObservables.size(); ++i )
	{
		allObservables.push_back( new Observable(*(input.allObservables[i])) );
	}
}

//Destructor
DataPoint::~DataPoint()
{
	while( !allObservables.empty() )
	{
		if( allObservables.back() != NULL ) { delete allObservables.back(); }
		allObservables.pop_back();
	}
}

//Retrieve names of all observables stored
vector<string> DataPoint::GetAllNames() const
{
	return allNames;
}

void DataPoint::RemoveObservable( const string input )
{
	vector<string>::iterator name_i = allNames.begin();
	vector<Observable*>::iterator obs_i = allObservables.begin();

	vector<string>::iterator name_to_remove;
	vector<Observable*>::iterator observable_to_remove;

	for( ; name_i != allNames.end(); ++name_i, ++obs_i )
	{
		if( (*name_i) == input )
		{
			name_to_remove = name_i;
			observable_to_remove = obs_i;
			break;
		}
	}

	if( *observable_to_remove != NULL ) delete *observable_to_remove;

	allNames.erase( name_to_remove );
	allObservables.erase( observable_to_remove );
}

Observable* DataPoint::GetObservable( unsigned int wanted ) const
{
	return allObservables[ wanted ];
}

//Retrieve an observable by its name
//	!!!THIS IS VERY, VERY, VERY WASTEFUL FOR LARGE DATASETS!!!
Observable* DataPoint::GetObservable(string const Name, const bool silence ) const
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if( nameIndex == -1 )
	{
		if( !silence ) cerr << "Observable name " << Name << " not found (2)" << endl;
		//this->Print();
		throw(-1543);
		//return new Observable( Name, 0.0, 0.0, "NameNotFoundError");
	}
	else
	{
		return allObservables[unsigned(nameIndex)];
	}
}

Observable* DataPoint::GetObservable( const ObservableRef& object, const bool silence ) const
{
	if( object.GetIndex() < 0 )
	{
		object.SetIndex( StringProcessing::VectorContains( &allNames, object.NameRef()) );
		if( object.GetIndex() >= 0 ) return allObservables[ (unsigned) object.GetIndex() ];
	}
	else
	{
		return allObservables[ (unsigned) object.GetIndex() ];
	}
	if( !silence ) cerr << "Observable name " << object.Name().c_str() << " not found (3)" << endl;
	throw(-20);
}

//Set an observable by name
bool DataPoint::SetObservable( string Name, Observable * NewObservable )
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if ( nameIndex == -1 )
	{
		cerr << "Observable name " << Name << " not found (4)" << endl;
		throw(438);
		//return false;
	}
	else
	{
		allObservables[unsigned(nameIndex)]->SetObservable(NewObservable);
		return true;
	}
}

void DataPoint::AddObservable( string Name, Observable* NewObservable )
{
	if( StringProcessing::VectorContains( &allNames, &Name ) == -1 )
	{
		allNames.push_back( Name );
		allObservables.push_back( new Observable(*NewObservable) );
	}
	else
	{
		this->SetObservable( Name, NewObservable );
	}
}

void DataPoint::AddObservable( string Name, double Value, string Unit, bool trusted, int nameIndex )
{
	Observable *tempObservable = new Observable( Name, Value, Unit );
	if( trusted )
	{
		allObservables[unsigned(nameIndex)]->SetObservable( tempObservable );
	}
	else
	{
		this->AddObservable( Name, tempObservable );
	}
}

//Initialise observable
bool DataPoint::SetObservable( string Name, double Value, string Unit, bool trusted, int nameIndex )
{
	Observable * temporaryObservable = new Observable( Name, Value, Unit );
	bool returnValue=false;
	if( trusted )
	{
		returnValue=true;
		allObservables[unsigned(nameIndex)]->SetObservable( temporaryObservable );
	}
	else
	{
		returnValue = SetObservable( Name, temporaryObservable );
	}
	delete temporaryObservable;
	return returnValue;
}

//	Used for Sorting DataPoints
bool DataPoint::operator() ( pair<DataPoint*,ObservableRef> first, pair<DataPoint*,ObservableRef> second )
{
	double param_val_1 = first.first->GetObservable( first.second )->GetValue();
	double param_val_2 = second.first->GetObservable( second.second )->GetValue();
	return (param_val_1 < param_val_2 );
}

double DataPoint::GetPseudoObservable( PseudoObservable& Input, vector<double> Values )
{
	if(Input.GetIndex()<0)
	{
		string name = Input.GetName();
		int lookup = StringProcessing::VectorContains( &allPseudoNames, &name );
		if( lookup == -1 )
		{
			allPseudoObservables.push_back( &Input );
			allPseudoNames.push_back( name );
			Input.SetIndex( (int)allPseudoObservables.size()-1 );
			allPseudoObservables.back()->SetIndex( (int)allPseudoObservables.size()-1 );
		}
		else
		{
			Input.SetIndex( lookup );
			allPseudoObservables[ (unsigned)Input.GetIndex() ] = &Input;
		}
	}
	else
	{
		if( Input.GetIndex()>(int)(allPseudoObservables.size()-1) )
		{
			string name = Input.GetName();
			int lookup = StringProcessing::VectorContains( &allPseudoNames, &name );
			if( lookup == -1 )
			{
				allPseudoObservables.push_back( &Input);
				allPseudoNames.push_back( name );
				Input.SetIndex( (int)allPseudoObservables.size()-1 );
				allPseudoObservables.back()->SetIndex( (int)allPseudoObservables.size()-1 );
			}
			else
			{
				Input.SetIndex( lookup );
				allPseudoObservables[ (unsigned)Input.GetIndex() ] = &Input;
			}
		}
		else
		{
			Input = *(allPseudoObservables[ (unsigned)Input.GetIndex() ]);
		}
	}

	PseudoObservable* thisObservable = (allPseudoObservables[ (unsigned)Input.GetIndex() ]);

	double outputObservable = 0.;

	if( Values.empty() )
	{
		if( !thisObservable->GetValid() )
		{
			vector<ObservableRef>* deps = thisObservable->GetDependencies();
			vector<double> input;   input.resize( deps->size() );
			unsigned int i=0;
			for( vector<ObservableRef>::iterator dep_i = deps->begin(); dep_i != deps->end(); ++dep_i, ++i )
			{
				input[ i ] = ( this->GetObservable( *(dep_i) )->GetValue() );
			}
			thisObservable->SetInput( input );
		}
	}
	else
	{
		if( !thisObservable->GetValid( Values ) )
		{
			thisObservable->SetInput( Values );
		}
	}

	outputObservable = thisObservable->GetPseudoObservable();

	return outputObservable;
}

void DataPoint::Print() const
{
	cout << "DataPoint:" << endl;
	for( unsigned int i=0; i< allObservables.size(); ++i )
	{
		allObservables[i]->Print();
	}

        if( !allPseudoObservables.empty() )
        {
                for( unsigned int i=0; i< allPseudoObservables.size(); ++i )
                {
                        allPseudoObservables[i]->Print();
                }
        }
	cout << endl;
}

void DataPoint::ClearPseudoObservable()
{
	vector<string> name1;
	allPseudoNames.swap( name1 );
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

