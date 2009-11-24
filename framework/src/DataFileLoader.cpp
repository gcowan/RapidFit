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
#include "InvalidObject.h"
#include "RootFileDataSet.h"
#include "MemoryDataSet.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "TFile.h"
#include "TBranch.h"
#include "TLeaf.h"


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
		data = new InvalidObject( "Unrecognised file type: " + fileName );
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
	return;
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

  std::ifstream file_to_read( fileName.c_str() );

  std::vector<string> observableNamesInFile;
  std::map<string, string> observableNamesToUnits;
  std::vector<string> observableNames = (data->GetBoundary())->GetAllNames();
  
  // Get a map of the observable names to their units
  for (int i = 0; i < observableNames.size(); i++ ){
        string name = observableNames[i];
        string unit = data->GetBoundary()->GetConstraint( name )->GetUnit();
        observableNamesToUnits[ name ] =  unit;
  }
  bool readFirstLine = false;
  int numberOfDataPoints = 0;

  if ( file_to_read.is_open() ) {
    while( ( !file_to_read.eof() ) && (numberOfDataPoints < numberEventsToRead) ) {
      string line;
      getline(file_to_read, line);
      vector<string> splitVec;

      // Need to do something here to allow for any comment statement using any amount of
      // white space between the column names.
      boost::split( splitVec, line, boost::is_any_of(" ") );

      if ( !readFirstLine ){
        if ( splitVec[0] != "#" ){
          cerr << "Ascii data set must start with a header containing the names of the" << endl;
          cerr << "observables that are used in the PDF." << endl;
        } else {
          for (int i = 1; i < splitVec.size(); i++){
            observableNamesInFile.push_back( splitVec[i] );
          }
        }
        readFirstLine = true;
      }

      DataPoint point( observableNames );
      try {
        for (int i = 0; i < splitVec.size(); i++){
          string name = observableNamesInFile[i];
          string unit = observableNamesToUnits[ name ];
	  // This is a little bit nasty as it assumes that anything
	  // with no units "" is not allowed.
          if (unit == "") continue;
          // Also, what happens if the observable isn't a double? // BEN_SAYS: there's no provision for this in the whole data model. I don't think we need it.
          point.SetObservable( name, boost::lexical_cast<double>(splitVec[i]), 0.0, unit);
        }
        numberOfDataPoints += 1;
      }
      catch(boost::bad_lexical_cast&){
        cout << "Input data incorrectly formatted" << endl;
      }
      //cout << "about to add data point" <<endl;
      data->AddDataPoint( &point );
    }
    file_to_read.close();
    cout << "Read " << numberOfDataPoints << " events from file: " << fileName << endl;
  }
  return;
}
