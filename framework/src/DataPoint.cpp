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
//	System Headers
#include <iostream>
#include <stdlib.h>

//Default constructor
DataPoint::DataPoint() : allObservables(), allNames(), allPseudoNames(), allPseudoObservables()
{
}

//Constructor with correct arguments
DataPoint::DataPoint( vector<string> NewNames ) : allObservables(), allNames(NewNames), allPseudoNames(), allPseudoObservables()
{
	allObservables.reserve( NewNames.size() );
	//Populate the map
	for (unsigned short int nameIndex = 0; nameIndex < NewNames.size(); nameIndex++)
	{
		allObservables.push_back( Observable() );
	}
}

//Destructor
DataPoint::~DataPoint()
{
}

//Retrieve names of all observables stored
vector<string> DataPoint::GetAllNames()
{
	return allNames;
}

//Retrieve an observable it's cached index, or find it and save it's index for reference
Observable * DataPoint::GetObservable( pair<string,int>* wanted_param )
{
	if( wanted_param->second != -1 )
	{
		return &allObservables[(unsigned)wanted_param->second];
	} else {
		wanted_param->second = StringProcessing::VectorContains( &allNames, &(wanted_param->first) );
	}
	if( wanted_param->second == -1 ){
		cerr << "Observable name " << wanted_param->first << " not found" <<endl;
	}else{
		return &allObservables[unsigned(wanted_param->second)];}
	exit(-1);
}

//Retrieve an observable by its name
//	!!!THIS IS VERY, VERY, VERY WASTEFUL FOR LARGE DATASETS!!!
Observable * DataPoint::GetObservable(string const Name)
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if ( nameIndex == -1 )
	{
		cerr << "Observable name " << Name << " not found" << endl;
		exit(1);
		//return new Observable( Name, 0.0, 0.0, "NameNotFoundError");
	}
	else
	{
		return &allObservables[unsigned(nameIndex)];
	}
}

Observable * DataPoint::GetObservable( ObservableRef& object )
{
	if( object.GetIndex() < 0 ) {
		object.SetIndex( StringProcessing::VectorContains( &allNames, object.NameRef()) );
		if( object.GetIndex() >= 0 ) return &allObservables[ (unsigned) object.GetIndex() ];
	} else {
		return &allObservables[ (unsigned) object.GetIndex() ];
	}
	cerr << "Observable name " << object.Name().c_str() << " not found" << endl;
	throw(-20);
}

//Set an observable by name
bool DataPoint::SetObservable( string Name, Observable * NewObservable )
{
	//Check if the name is stored in the map
	int nameIndex = StringProcessing::VectorContains( &allNames, &Name );
	if ( nameIndex == -1 )
	{
		cerr << "Observable name " << Name << " not found" << endl;
		exit(1);
		//return false;
	}
	else
	{
		allObservables[unsigned(nameIndex)] = *NewObservable;
		return true;
	}
}

//Initialise observable
bool DataPoint::SetObservable( const string Name, const double Value, const double Error, const string Unit, const bool trusted, const int nameIndex )
{
	Observable * temporaryObservable = new Observable( Name, Value, Error, Unit );
	bool returnValue=false;
	if( trusted )
	{
		returnValue=true;
		allObservables[unsigned(nameIndex)] = *temporaryObservable;
	} else {
		returnValue = SetObservable( Name, temporaryObservable );
	}
	delete temporaryObservable;
	return returnValue;
}

//	Used for Sorting DataPoints
bool DataPoint::operator() ( pair<DataPoint,pair<string,int> > first, pair<DataPoint,pair<string,int> > second )
{
	double param_val_1 = first.first.GetObservable( &first.second )->GetValue();
	double param_val_2 = second.first.GetObservable( &second.second )->GetValue();
	return (param_val_1 < param_val_2 );
}

Observable* DataPoint::GetPseudoObservable( ObservableRef& final_observable, string dependencies, double (*pseudoRelation)(vector<double>) )
{
	if( ( final_observable.GetIndex() >=0 ) && ( ( (int)allPseudoObservables.size() -1 ) >= final_observable.GetIndex() ) )
	{
		return &( allPseudoObservables[ final_observable.GetIndex() ] );
	}
	else
	{
		//	Test if all observables before this one have been created
		//	They should have been, if they haven't then this _will_ lead to logic errors
		if( ( (int) allPseudoObservables.size() -2 ) != ( final_observable.GetIndex() - 1 ) )
		{
			cerr << endl << "\tWarning: Internal Logic Error in the handling of PsuedoObservables, resorting to look-up method" << endl << endl;
			final_observable.SetIndex( -1 );
		}

		vector<string> split_deps = StringProcessing::SplitString( dependencies, ':' );

		vector<double> input_from_deps;
		for( unsigned int i=0; i< split_deps.size(); ++i )
		{
			input_from_deps.push_back( this->GetObservable( split_deps[i] )->GetValue() );
		}
						//		Name		Value			     Error   Unit
		Observable* new_pseudo = new Observable( final_observable, pseudoRelation( input_from_deps ), 0, "none" );

		allPseudoObservables.push_back( *new_pseudo );
		allPseudoNames.push_back( final_observable );

		if( final_observable.GetIndex() == -1 )
		{
			final_observable.SetIndex( (int)allPseudoObservables.size() -1 );
		}

		return &( allPseudoObservables.back() );
	}
	return NULL;
}

//	Same function but for pair<double,double> (*pseudoRelation)(vector<double>)  which provides Value and error
//	This is unlikely to be a generic function for more than 1 pdf but is here for extensibility
Observable* DataPoint::GetPseudoObservable( ObservableRef& final_observable, string dependencies, pair<double,double> (*pseudoRelation)(vector<double>) )
{
	if( ( final_observable.GetIndex() >=0 ) && ( ( (int)allPseudoObservables.size() -1 ) >= final_observable.GetIndex() ) )
	{
		return &( allPseudoObservables[ final_observable.GetIndex() ] );
	}
	else
	{
		//      Test if all observables before this one have been created
		//      They should have been, if they haven't then this _will_ lead to logic errors
		if( ( (int) allPseudoObservables.size() -2 ) != ( final_observable.GetIndex() - 1 ) )
		{
			cerr << endl << "\tWarning: Internal Logic Error in the handling of PsuedoObservables, resorting to look-up method" << endl << endl;
			final_observable.SetIndex( -1 );
		}

		vector<string> split_deps = StringProcessing::SplitString( dependencies, ':' );

		vector<double> input_from_deps;
		for( unsigned int i=0; i< split_deps.size(); ++i )
		{
			input_from_deps.push_back( this->GetObservable( split_deps[i] )->GetValue() );
		}

		pair<double,double> relation_result = pseudoRelation( input_from_deps );

						//              Name            Value                Error                Unit
		Observable* new_pseudo = new Observable( final_observable, relation_result.first, relation_result.second, "none" );

		allPseudoObservables.push_back( *new_pseudo );
		allPseudoNames.push_back( final_observable );

		if( final_observable.GetIndex() == -1 )
		{
			final_observable.SetIndex( (int)allPseudoObservables.size() -1 );
		}

		return &( allPseudoObservables.back() );
	}
	return NULL;
}

