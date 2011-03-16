/**
  @class DataPoint

  Holds all observables for a given event

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */


#include "DataPoint.h"
#include "StringProcessing.h"
#include <iostream>
#include <stdlib.h>

//Default constructor
DataPoint::DataPoint()
{
}

//Constructor with correct arguments
DataPoint::DataPoint( vector<string> NewNames )
{
	allObservables.reserve( NewNames.size() );
	//Populate the map
	for (unsigned short int nameIndex = 0; nameIndex < NewNames.size(); nameIndex++)
	{
		allObservables.push_back( Observable() );
	}

	allNames = NewNames;
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
		return &allObservables[wanted_param->second];
	} else {
		wanted_param->second = StringProcessing::VectorContains( &allNames, &(wanted_param->first) );
	}
	if( wanted_param->second == -1 ){
		cerr << "Observable name " << wanted_param->first << " not found" <<endl;
	}else{
		return &allObservables[wanted_param->second];}
	exit(-1);
}

//Retrieve an observable by its name
//	!!!THIS IS VERY, VERY, VERY WASTEFUL FOR LARGE DATASETS!!!
Observable * DataPoint::GetObservable(string Name)
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
		return &allObservables[nameIndex];
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
		allObservables[nameIndex] = *NewObservable;
		return true;
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
		allObservables[nameIndex] = *temporaryObservable;
	} else {
		returnValue = SetObservable( Name, temporaryObservable );
	}
	delete temporaryObservable;
	return returnValue;
}
