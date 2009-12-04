// $Id: DataFileLoader.cpp,v 1.1 2009/11/10 10:35:46 gcowan Exp $
/**
  @class DataFileLoader

  Holds methods for loading given file formats into an IDataSet

  @author Benjamin M Wynne bwynne@cern.ch
  @author Greig A Cowan greig.cowan@cern.ch
  @date 2009-10-02
 */

#include "DataFileLoader.h"
#include "StringProcessing.h"
#include <stdlib.h>
#include "RootFileDataSet.h"
#include "MemoryDataSet.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "TFile.h"
#include "TBranch.h"
#include "TLeaf.h"
#include <stdlib.h>

//Default constructor
DataFileLoader::DataFileLoader()
{
}

/*
   DataFileLoader::DataFileLoader( string fileName, string tuplePath, PhaseSpaceBoundary * dataBoundary, long numberEventsToRead )
   {
//Find the file type, and treat appropriately
vector<string> splitFileName = StringProcessing::SplitString( fileName, '.' );
string fileNameExtension = splitFileName[ splitFileName.size() - 1 ];
if ( fileNameExtension == "root" )
{
// Make a MemoryDataSet from a root file
data = new MemoryDataSet( dataBoundary );
LoadRootFile( fileName, tuplePath, numberEventsToRead );
}
}
 */

//Constructor with correct arguments
DataFileLoader::DataFileLoader( string fileName, PhaseSpaceBoundary * dataBoundary, long numberEventsToRead )
{
	//Find the file type, and treat appropriately
	vector<string> splitFileName = StringProcessing::SplitString( fileName, '.' );
	string fileNameExtension = splitFileName[ splitFileName.size() - 1 ];
	if ( fileNameExtension == "root")
	{
		//Make a RootFileDataSet from a root file
		//data = new RootFileDataSet( fileName, dataBoundary );
		data = new MemoryDataSet( dataBoundary );
		LoadRootFileIntoMemory( fileName, "dataNTuple", numberEventsToRead );
	}
	else if ( fileNameExtension == "csv" )
	{
		//Load a csv file into memory (or make a root file data set if too big?)
		cerr << "CSV file not implemented yet" << endl;
		exit(1);
	}
	else if ( fileNameExtension == "txt" )
	{
		//Load a ascii text file into memory (or make a root file data set if too big?)
		data = new MemoryDataSet( dataBoundary );
		LoadAsciiFileIntoMemory( fileName, numberEventsToRead );
	}
	else
	{
		cerr << "Unrecognised file type: " << fileName << endl;
		exit(1);
	}
}

//Destructor
DataFileLoader::~DataFileLoader()
{
}

//Return the new DataSet
IDataSet * DataFileLoader::GetDataSet()
{
	return data;
}

void DataFileLoader::LoadRootFileIntoMemory( string fileName, string ntuplePath, long numberEventsToRead )
{
	std::vector<string> observableNames = (data->GetBoundary())->GetAllNames();
	int numberOfObservables = observableNames.size();

	TFile * inputFile = new TFile( fileName.c_str(), "UPDATE" );
	TNtuple * ntuple = (TNtuple*)inputFile->Get( ntuplePath.c_str() );
	int totalNumberOfEvents = ntuple->GetEntries();

	// Check that ntuple is compatible with the PhaseSpaceBoundary that is defined
	// in the XML file.
	if ( !CheckTNtupleWithBoundary( ntuple, data->GetBoundary() ) )
	{
		cerr << "NTuple is incompatible with boundary" << endl;
		return;
	}

	// Need to set the branch addresses
	vector<TBranch *> branches( numberOfObservables );
	vector<Float_t> observableValues( numberOfObservables );
	for ( int obsIndex = 0; obsIndex < numberOfObservables; obsIndex++ )
	{
		ntuple->SetBranchAddress(observableNames[obsIndex].c_str(), &(observableValues[obsIndex]), &(branches[obsIndex]));
	}

	// Now populate the dataset
	int numberOfDataPoints = 0;
	while ( numberOfDataPoints < numberEventsToRead && numberOfDataPoints < totalNumberOfEvents )
	{
		DataPoint point( observableNames );
		ntuple->GetEntry( numberOfDataPoints );

		for ( int obsIndex = 0; obsIndex < numberOfObservables; obsIndex++ )
		{
			string name = observableNames[obsIndex];
			string unit = data->GetBoundary()->GetConstraint( name )->GetUnit();
			point.SetObservable( name, observableValues[obsIndex], 0.0, unit);
		}
		numberOfDataPoints += 1;
		data->AddDataPoint( &point );
	}

	inputFile->Close();
	delete inputFile;
	cout << "Read " << numberOfDataPoints << " events from ROOT file: " << fileName << endl;
}

bool DataFileLoader::CheckTNtupleWithBoundary( TNtuple * TestTuple, PhaseSpaceBoundary * TestBoundary )
{
	bool compatible = true;

	vector<string> allNames = TestBoundary->GetAllNames();
	for ( int observableIndex = 0; observableIndex < allNames.size(); observableIndex++ )
	{
		bool found = false;
		//Find the requested observable name in the ntuple
		TIter observableIterator( TestTuple->GetListOfLeaves() );
		TLeaf * observableLeaf;
		while ( observableLeaf = (TLeaf*)observableIterator() )
		{
			string leafName = observableLeaf->GetName();
			if ( leafName == allNames[observableIndex] )
			{
				found = true;
				break;
			}
		}
		//If the name isn't present, the ntuple is incompatible
		if ( !found )
		{
			compatible = false;
			break;
		}
	}
	return compatible;
}

void DataFileLoader::LoadAsciiFileIntoMemory( string fileName, long numberEventsToRead )
{
	std::vector<string> observableNamesInFile;
	std::map<string, string> observableNamesToUnits;
	std::vector<string> observableNames = (data->GetBoundary())->GetAllNames();

	//Make a map of the observable names to their units
	for (int i = 0; i < observableNames.size(); i++ )
	{
		string name = observableNames[i];
		string unit = data->GetBoundary()->GetConstraint( name )->GetUnit();
		observableNamesToUnits[ name ] =  unit;
	}
	bool readFirstLine = false;
	int numberOfDataPoints = 0;

	//Open the file
	std::ifstream file_to_read( fileName.c_str() );
	if ( file_to_read.is_open() )
	{
		//Stop reading at EOF or when sufficient events loaded
		while( ( !file_to_read.eof() ) && (numberOfDataPoints < numberEventsToRead) )
		{
			//Read a line, and split on white space
			string line;
			getline(file_to_read, line);
			vector<string> splitVec;
			splitVec = StringProcessing::SplitString( line, ' ' );

			if (readFirstLine)
			{
				//Make a data point with the values on the line
				DataPoint point(observableNames);
				for (int i = 0; i < splitVec.size(); i++)
				{
					string name = observableNamesInFile[i];
					string unit = observableNamesToUnits[ name ];
					point.SetObservable( name, strtod( splitVec[i].c_str(), NULL ), 0.0, unit );
				}
				numberOfDataPoints += 1;
				data->AddDataPoint( &point );
			}
			else
			{
				//Load the obsevable names
				for (int i = 0; i < splitVec.size(); i++)
				{
					observableNamesInFile.push_back( splitVec[i] );
				}

				//Check that all required observable names are specified
				for ( int observableIndex = 0; observableIndex < observableNames.size(); observableIndex++ )
				{
					if ( StringProcessing::VectorContains( &observableNamesInFile, &( observableNames[observableIndex] ) ) == -1 )
					{
						cerr << "Observable " << observableNames[observableIndex] << " not provided by file " << fileName << endl;
						exit(1);
					}
				}

				readFirstLine = true;
			}
		}
		file_to_read.close();
		cout << "Read " << numberOfDataPoints << " events from file: " << fileName << endl;
	}
	else
	{
		cerr << "Failed to open input data file: " << fileName << endl;
		exit(1);
	}
}
