/**
  @class DataPoint

  Holds all observables for a given event

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

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
DataPoint::DataPoint() : allObservables(), allNames(), allPseudoNames(), allPseudoNames2(), allPseudoObservables(), allPseudoObservables2(), myPhaseSpaceBoundary(NULL), thisDiscreteIndex(-1), WeightValue(1.)
{
}

//Constructor with correct arguments
DataPoint::DataPoint( vector<string> NewNames ) : allObservables(), allNames(), allPseudoNames(), allPseudoNames2(), allPseudoObservables(), allPseudoObservables2(), myPhaseSpaceBoundary(NULL), thisDiscreteIndex(-1), WeightValue(1.)
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
	}
}

DataPoint::DataPoint( const DataPoint& input ) :
	allObservables(), allNames(input.allNames), allPseudoNames(input.allPseudoNames), allPseudoObservables(), myPhaseSpaceBoundary(input.myPhaseSpaceBoundary),
	thisDiscreteIndex(input.thisDiscreteIndex), allPseudoNames2(input.allPseudoNames2), allPseudoObservables2(), WeightValue(input.WeightValue)
{
	for( unsigned int i=0; i< input.allObservables.size(); ++i )
	{
		allObservables.push_back( new Observable(*(input.allObservables[i])) );
	}
	for( unsigned int i=0; i< input.allPseudoObservables.size(); ++i )
	{
		allPseudoObservables.push_back( PseudoObservable( (input.allPseudoObservables[i]) ) );
	}
	for( unsigned int i=0; i< input.allPseudoObservables2.size(); ++i )
	{
		allPseudoObservables2.push_back( PseudoObservable( (input.allPseudoObservables2[i]) ) );
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

void DataPoint::AddObservable( string Name, double Value, double Error, string Unit, bool trusted, int nameIndex )
{
	Observable *tempObservable = new Observable( Name, Value, Error, Unit );
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
bool DataPoint::SetObservable( string Name, double Value, double Error, string Unit, bool trusted, int nameIndex )
{
	Observable * temporaryObservable = new Observable( Name, Value, Error, Unit );
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

//#pragma GCC diagnostic ignored "-Wfloat-equal"
double DataPoint::GetPseudoObservable( PseudoObservable& Input )
{
	if(Input.GetIndex()<0)
	{
		string name = Input.GetName();
		int lookup = StringProcessing::VectorContains( &allPseudoNames, &name );
		if( lookup == -1 )
		{
			allPseudoObservables.push_back( PseudoObservable( Input ) );
			allPseudoNames.push_back( name );
			allPseudoObservables.back().SetIndex( (int)allPseudoObservables.size()-1 );
			Input.SetIndex( (int)allPseudoObservables.size()-1 );
			//cout << "<0 " << "creating" << " at " << allPseudoObservables.size()-1 << endl;
		}
		else
		{
			Input.SetIndex( lookup );
			allPseudoObservables[ (unsigned)Input.GetIndex() ] = PseudoObservable( Input );
			//cout << "found at " << lookup << endl;
		}
	}
	else
	{
		if( Input.GetIndex()>(int)(allPseudoObservables.size()-1) )
		{
			//int var = (int)allPseudoObservables.size();
			allPseudoObservables.resize( (unsigned)Input.GetIndex()+1 );
			allPseudoObservables[(unsigned)Input.GetIndex()] = PseudoObservable( Input );
			//cout << "resize to " << Input.GetIndex()+1 << " from " << var << endl;
		}
		else
		{
			allPseudoObservables[(unsigned)Input.GetIndex()] = PseudoObservable( Input );
			//cout << "define " << Input.GetIndex() << " of " << allPseudoObservables.size()-1 << endl;
			//cout << "requested at (2) " << (unsigned)Input.GetIndex() << " of " << allPseudoObservables.size()-1 << endl;
		}
	}

	PseudoObservable thisObservable = allPseudoObservables[ (unsigned)Input.GetIndex() ];

	double outputObservable = 0.;

	if( thisObservable.GetValid() )
	{
		outputObservable = thisObservable.GetPseudoObservable();
	}

	if( !(thisObservable.GetValid()) ) //|| outputObservable == 0. )
	{
		vector<ObservableRef>* deps = thisObservable.GetDependencies();
		vector<double> input;
		for( vector<ObservableRef>::iterator dep_i = deps->begin(); dep_i != deps->end(); ++dep_i )
		{
			input.push_back( this->GetObservable( *(dep_i) )->GetValue() );
		}
		thisObservable.SetInput( input );

		outputObservable = thisObservable.GetPseudoObservable();
	}

	thisObservable.SetValid( Input.GetValid() );

	return outputObservable;
}
//#pragma GCC diagnostic pop

double DataPoint::GetPseudoObservable( PseudoObservable& Input, vector<double> Values )
{
	if(Input.GetIndex()<0)
	{
		string name = Input.GetName();
		int lookup = StringProcessing::VectorContains( &allPseudoNames2, &name );
		if( lookup == -1 )
		{
			allPseudoObservables2.push_back( PseudoObservable( Input ) );
			allPseudoNames2.push_back( name );
			allPseudoObservables2.back().SetIndex( (int)allPseudoObservables2.size()-1 );
			Input.SetIndex( (int)allPseudoObservables2.size()-1 );
			//cout << "<0 " << "creating" << " at " << allPseudoObservables2.size()-1 << endl;
		}
		else
		{
			Input.SetIndex( lookup );
			allPseudoObservables2[ (unsigned)Input.GetIndex() ] = PseudoObservable( Input );
			//cout << "found at " << lookup << endl;
		}
	}
	else
	{
		if( Input.GetIndex()>(int)(allPseudoObservables2.size()-1) )
		{
			//int var = (int)allPseudoObservables2.size();
			allPseudoObservables2.resize( (unsigned)Input.GetIndex()+1 );
			allPseudoObservables2[ (unsigned)Input.GetIndex() ] = PseudoObservable( Input );
			//cout << "resize to " << Input.GetIndex()+1 << " from " << var << endl;
		}
		else
		{
			allPseudoObservables2[ (unsigned)Input.GetIndex() ] = PseudoObservable( Input );
			//cout << "define " << Input.GetIndex() << " of " << allPseudoObservables2.size()-1 << endl;
			//cout << "requested at (3) " << (unsigned)Input.GetIndex() << " of " << allPseudoObservables2.size()-1 << endl;
		}
	}

	PseudoObservable thisObservable = allPseudoObservables2[ (unsigned)Input.GetIndex() ];

	double outputObservable = 0.;

	if( thisObservable.GetValid( Values ) )
	{
		outputObservable = thisObservable.GetPseudoObservable();
	}
	else
	{
		thisObservable.SetInput( Values );

		outputObservable = thisObservable.GetPseudoObservable();
	}

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
			allPseudoObservables[i].Print();
		}
	}
        if( !allPseudoObservables2.empty() )
        {
                for( unsigned int i=0; i< allPseudoObservables2.size(); ++i )
                {
                        allPseudoObservables2[i].Print();
                }
        }
	cout << endl;
}

void DataPoint::ClearPsuedoObservable()
{
	allPseudoObservables.clear();
	allPseudoObservables2.clear();
	allPseudoNames.clear();
	allPseudoNames2.clear();
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

