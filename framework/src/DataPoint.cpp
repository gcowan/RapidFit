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
DataPoint::DataPoint() : allObservables(), allNames(), allPseudoNames(), allPseudoObservables(), myPhaseSpaceBoundary(NULL), thisDiscreteIndex(-1)
{
}

//Constructor with correct arguments
DataPoint::DataPoint( vector<string> NewNames ) : allObservables(), allNames(), allPseudoNames(), allPseudoObservables(), myPhaseSpaceBoundary(NULL), thisDiscreteIndex(-1)
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
	thisDiscreteIndex(input.thisDiscreteIndex)
	, allPseudoNames2(input.allPseudoNames2), allPseudoObservables2()
{
	for( unsigned int i=0; i< input.allObservables.size(); ++i )
	{
		allObservables.push_back( new Observable(*(input.allObservables[i])) );
	}
	for( unsigned int i=0; i< input.allPseudoObservables.size(); ++i )
	{
		allPseudoObservables.push_back( new PseudoObservable(*(input.allPseudoObservables[i])) );
	}
	for( unsigned int i=0; i< input.allPseudoObservables2.size(); ++i )
	{
		allPseudoObservables2.push_back( new PseudoObservable(*(input.allPseudoObservables2[i])) );
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
	while( !allPseudoObservables.empty() )
	{
		if( allPseudoObservables.back() != NULL ) { delete allPseudoObservables.back(); }
		allPseudoObservables.pop_back();
	}
	while( !allPseudoObservables2.empty() )
	{
		if( allPseudoObservables2.back() != NULL ) { delete allPseudoObservables2.back(); }
		allPseudoObservables2.pop_back();
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

//Retrieve an observable it's cached index, or find it and save it's index for reference
Observable* DataPoint::GetObservable( pair<string,int>* wanted_param ) const
{
	if( wanted_param->second != -1 )
	{
		return allObservables[(unsigned)wanted_param->second];
	} else {
		wanted_param->second = StringProcessing::VectorContains( &allNames, &(wanted_param->first) );
	}
	if( wanted_param->second == -1 ){
		cerr << "Observable name " << wanted_param->first << " not found (1)" << endl;
	}else{
		return allObservables[unsigned(wanted_param->second)];
	}
	exit(-1);
}

//Retrieve an observable by its name
//	!!!THIS IS VERY, VERY, VERY WASTEFUL FOR LARGE DATASETS!!!
Observable* DataPoint::GetObservable(string const Name) const
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if( nameIndex == -1 )
	{
		cerr << "Observable name " << Name << " not found (2)" << endl;
		this->Print();
		exit(1);
		//return new Observable( Name, 0.0, 0.0, "NameNotFoundError");
	}
	else
	{
		return allObservables[unsigned(nameIndex)];
	}
}

Observable* DataPoint::GetObservable( const ObservableRef& object ) const
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
	cerr << "Observable name " << object.Name().c_str() << " not found (3)" << endl;
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
bool DataPoint::operator() ( pair<DataPoint*,pair<string,int> > first, pair<DataPoint*,pair<string,int> > second )
{
	double param_val_1 = first.first->GetObservable( &first.second )->GetValue();
	double param_val_2 = second.first->GetObservable( &second.second )->GetValue();
	return (param_val_1 < param_val_2 );
}

double DataPoint::GetPseudoObservable( PseudoObservable& Input )
{
	if(Input.GetIndex()<0)
	{
		string name = Input.GetName();
		int lookup = StringProcessing::VectorContains( &allPseudoNames, &name );
		if( lookup == -1 )
		{
			allPseudoObservables.push_back( new PseudoObservable( Input ) );
			allPseudoNames.push_back( name );
			allPseudoObservables.back()->SetIndex( (int)allPseudoObservables.size()-1 );
			Input.SetIndex( (int)allPseudoObservables.size()-1 );
			//cout << "<0 " << "creating" << " at " << allPseudoObservables.size()-1 << endl;
		}
		else
		{
			Input.SetIndex( lookup );
			if( allPseudoObservables[ (unsigned)Input.GetIndex() ] == NULL )
			{
				allPseudoObservables[ (unsigned)Input.GetIndex() ] = new PseudoObservable( Input );
			}
			//cout << "found at " << lookup << endl;
		}
	}
	else
	{
		if( Input.GetIndex()>(int)(allPseudoObservables.size()-1) )
		{
			//int var = (int)allPseudoObservables.size();
			allPseudoObservables.resize( Input.GetIndex()+1, NULL );
			if( allPseudoObservables[(unsigned)Input.GetIndex()] != NULL ) delete allPseudoObservables[(unsigned)Input.GetIndex()];
			allPseudoObservables[(unsigned)Input.GetIndex()] = new PseudoObservable( Input );
			//cout << "resize to " << Input.GetIndex()+1 << " from " << var << endl;
		}
		else
		{
			if( allPseudoObservables[(unsigned)Input.GetIndex()] == NULL )
			{
				allPseudoObservables[(unsigned)Input.GetIndex()] = new PseudoObservable( Input );
				//cout << "define " << Input.GetIndex() << " of " << allPseudoObservables.size()-1 << endl;
			}
			//cout << "requested at (2) " << (unsigned)Input.GetIndex() << " of " << allPseudoObservables.size()-1 << endl;
		}
	}

	PseudoObservable* thisObservable = allPseudoObservables[ (unsigned)Input.GetIndex() ];

	double outputObservable = 0.;

	if( thisObservable->GetValid() )
	{
		outputObservable = thisObservable->GetPseudoObservable();
	}

//#pragma GCC diagnostic ignored "-Wfloat-equal"
	if( !thisObservable->GetValid() ) //|| outputObservable == 0. )
//#pragma GCC diagnostic pop
	{
		vector<ObservableRef>* deps = thisObservable->GetDependencies();
		vector<double> input;
		for( vector<ObservableRef>::iterator dep_i = deps->begin(); dep_i != deps->end(); ++dep_i )
		{
			input.push_back( this->GetObservable( *(dep_i) )->GetValue() );
		}
		thisObservable->SetInput( input );

		outputObservable = thisObservable->GetPseudoObservable();
	}

	thisObservable->SetValid( Input.GetValid() );

	return outputObservable;
}

double DataPoint::GetPseudoObservable( PseudoObservable& Input, vector<double> Values )
{
	if(Input.GetIndex()<0)
	{
		string name = Input.GetName();
		int lookup = StringProcessing::VectorContains( &allPseudoNames2, &name );
		if( lookup == -1 )
		{
			allPseudoObservables2.push_back( new PseudoObservable( Input ) );
			allPseudoNames2.push_back( name );
			allPseudoObservables2.back()->SetIndex( (int)allPseudoObservables2.size()-1 );
			Input.SetIndex( (int)allPseudoObservables2.size()-1 );
			//cout << "<0 " << "creating" << " at " << allPseudoObservables2.size()-1 << endl;
		}
		else
		{
			Input.SetIndex( lookup );
			if( allPseudoObservables2[ (unsigned)Input.GetIndex() ] == NULL )
			{
				allPseudoObservables2[ (unsigned)Input.GetIndex() ] = new PseudoObservable( Input );
			}
			//cout << "found at " << lookup << endl;
		}
	}
	else
	{
		if( Input.GetIndex()>(int)(allPseudoObservables2.size()-1) )
		{
			//int var = (int)allPseudoObservables2.size();
			allPseudoObservables2.resize( Input.GetIndex()+1, NULL );
			allPseudoObservables2[(unsigned)Input.GetIndex()] = new PseudoObservable( Input );
			//cout << "resize to " << Input.GetIndex()+1 << " from " << var << endl;
		}
		else
		{
			if( allPseudoObservables2[(unsigned)Input.GetIndex()] == NULL )
			{
				allPseudoObservables2[(unsigned)Input.GetIndex()] = new PseudoObservable( Input );
				//cout << "define " << Input.GetIndex() << " of " << allPseudoObservables2.size()-1 << endl;
			}
			//cout << "requested at (3) " << (unsigned)Input.GetIndex() << " of " << allPseudoObservables2.size()-1 << endl;
		}
	}

	PseudoObservable* thisObservable = allPseudoObservables2[ (unsigned)Input.GetIndex() ];

	double outputObservable = 0.;

	if( thisObservable->GetValid( Values ) )
	{
		outputObservable = thisObservable->GetPseudoObservable();
	}
	else
	{
		thisObservable->SetInput( Values );

		outputObservable = thisObservable->GetPseudoObservable();
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
			allPseudoObservables[i]->Print();
		}
	}
	cout << endl;
}

void DataPoint::ClearPsuedoObservable()
{
	while( !allPseudoObservables.empty() )
	{
		if( allPseudoObservables.back() != NULL )
		{
			delete allPseudoObservables.back();
		}
		allPseudoObservables.pop_back();
	}
	while( !allPseudoObservables2.empty() )                          
	{                                                               
		if( allPseudoObservables2.back() != NULL )               
		{                                                       
			delete allPseudoObservables2.back();             
		}                                                       
		allPseudoObservables2.pop_back();                        
	}
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

