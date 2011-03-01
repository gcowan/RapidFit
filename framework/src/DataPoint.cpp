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

//Retrieve an observable by its name
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
bool DataPoint::SetObservable( string Name, double Value, double Error, string Unit )
{
	Observable * temporaryObservable = new Observable( Name, Value, Error, Unit );
	bool returnValue = SetObservable( Name, temporaryObservable );
	delete temporaryObservable;
	return returnValue;
}
