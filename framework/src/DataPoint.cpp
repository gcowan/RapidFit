/**
  @class DataPoint

  Holds all observables for a given event

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

//	RapidFit Headers
#include "DataPoint.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <stdlib.h>

//Default constructor
DataPoint::DataPoint() : allObservables(), allNames()
{
}

//Constructor with correct arguments
DataPoint::DataPoint( vector<string> NewNames ) : allObservables(), allNames(NewNames)
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

Observable const * DataPoint::GetSafeObservable( string const Name ) const
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
