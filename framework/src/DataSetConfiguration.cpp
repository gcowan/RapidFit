/**
  @class DataSetConfiguration

  A class for holding the data to create a data set, and creating that data set when requested

  @author Benjamin M Wynne bwynne@cern.ch
  @author Greig A Cowan greig.cowan@cern.ch
  @data 2009-12-16
 */

#include "StringProcessing.h"
#include "MemoryDataSet.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "TFile.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TEventList.h"
#include "DataSetConfiguration.h"
#include "ClassLookUp.h"
#include <stdlib.h>

//Default constructor
DataSetConfiguration::DataSetConfiguration()
{
}

//Constructor with correct argument
DataSetConfiguration::DataSetConfiguration( string DataSource, long DataNumber, string cut, vector<string> DataArguments, vector<string> DataArgumentNames ) : source(DataSource),
	numberEvents(DataNumber),
	cutString(cut),
	arguments(DataArguments), argumentNames(DataArgumentNames), separateGeneratePDF(false), parametersAreSet(false)
{
}

//Constructor with separate data generation PDF
DataSetConfiguration::DataSetConfiguration( string DataSource, long DataNumber, string cut, vector<string> DataArguments, vector<string> DataArgumentNames, IPDF * DataPDF ) : source(DataSource),
	numberEvents(DataNumber),
	cutString(cut),
	arguments(DataArguments), argumentNames(DataArgumentNames), generatePDF(DataPDF), separateGeneratePDF(true), parametersAreSet(false)
{
}

//Destructor
DataSetConfiguration::~DataSetConfiguration()
{
}

//Set the parameters of the generation PDF
bool DataSetConfiguration::SetPhysicsParameters( ParameterSet * InputParameters )
{
	if (separateGeneratePDF)
	{
		if ( generatePDF->SetPhysicsParameters(InputParameters) )
		{
			parametersAreSet = true;
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		parametersAreSet = true;
		return true;
	}
}

//Create the DataSet
IDataSet * DataSetConfiguration::MakeDataSet( PhaseSpaceBoundary * DataBoundary, IPDF * FitPDF )
{
	//Some kind of decision about what kind of data set to use?
	IDataSet * newDataSet;

	if ( source == "File" )
	{
		//Load data from file
		newDataSet = LoadDataFile( arguments, argumentNames, DataBoundary, numberEvents );
	}
	else
	{
		//PDF parameters must be set before data can be generated
		if (parametersAreSet)
		{
			//Assume it's an accept/reject generator, or some child of it
			IDataGenerator * dataGenerator;
			if (separateGeneratePDF)
			{
				dataGenerator = ClassLookUp::LookUpDataGenerator( source, DataBoundary, generatePDF );
			}
			else
			{
				dataGenerator = ClassLookUp::LookUpDataGenerator( source, DataBoundary, FitPDF );
			}
			dataGenerator->GenerateData(numberEvents);
			newDataSet = dataGenerator->GetDataSet();
			delete dataGenerator;
		}
		else
		{
			cerr << "PDF parameters must be set before data can be generated" << endl;
			exit(1);
		}
	}

	return newDataSet;
}

//Constructor with correct arguments
IDataSet * DataSetConfiguration::LoadDataFile( vector<string> Arguments, vector<string> ArgumentNames, PhaseSpaceBoundary * DataBoundary, long NumberEventsToRead )
{
	//Find file name
	string searchName = "FileName";
	int fileNameIndex = StringProcessing::VectorContains( &ArgumentNames, &searchName );
	string fileName = "NotFound";
	if ( fileNameIndex >= 0 )
	{
		fileName = Arguments[fileNameIndex];
	}
	else
	{
		cerr << "FileName argument not found" << endl;
		exit(1);
	}

	//Find nTuple path if specified
	searchName = "NTuplePath";
	int nTuplePathIndex = StringProcessing::VectorContains( &ArgumentNames, &searchName );
	string nTuplePath = "NotFound";
	if ( nTuplePathIndex >= 0 )
	{
		nTuplePath = Arguments[nTuplePathIndex];
	}

	//Find the file type, and treat appropriately
	vector<string> splitFileName = StringProcessing::SplitString( fileName, '.' );
	string fileNameExtension = splitFileName[ splitFileName.size() - 1 ];
	if ( fileNameExtension == "root")
	{
		//Make a RootFileDataSet from a root file
		//data = new RootFileDataSet( fileName, dataBoundary );
		return LoadRootFileIntoMemory( fileName, nTuplePath, NumberEventsToRead, DataBoundary );
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
		return LoadAsciiFileIntoMemory( fileName, NumberEventsToRead, DataBoundary );
	}
	else
	{
		cerr << "Unrecognised file type: " << fileName << endl;
		exit(1);
	}
}

IDataSet * DataSetConfiguration::LoadRootFileIntoMemory( string fileName, string ntuplePath, long numberEventsToRead, PhaseSpaceBoundary * DataBoundary )
{
	MemoryDataSet * data = new MemoryDataSet(DataBoundary);
	std::vector<string> observableNames = DataBoundary->GetAllNames();
	int numberOfObservables = observableNames.size();

	TFile * inputFile = new TFile( fileName.c_str(), "READ" );
	TNtuple * ntuple = (TNtuple*)inputFile->Get( ntuplePath.c_str() );
	if ( ntuple == NULL )
	{
		cerr << "Ntuple not found. Are you specifying the correct file and ntuple path?" << endl;
		exit(1);
	}
	int totalNumberOfEvents = ntuple->GetEntries();
	int numberOfEventsAfterCut = ntuple->Draw(">>evtList", cutString.c_str()); // apply the cut which automatically creates the evtList object
	if ( numberOfEventsAfterCut == -1 )
	{
		cerr << "Please check the cut string you are using!" << endl;
		exit(1);
	}
	TEventList * evtList = (TEventList*)gDirectory->Get("evtList"); // ROOT is weird

	cout << "Total number of events in file: " << totalNumberOfEvents << endl;
	cout << "You have applied this cut to the data: '" << cutString << "'" << endl;
	cout << "Total number of events after cut: " << numberOfEventsAfterCut << endl;

	// Check that ntuple is compatible with the PhaseSpaceBoundary that is defined
	// in the XML file.
	if ( !CheckTNtupleWithBoundary( ntuple, data->GetBoundary() ) )
	{
		cerr << "NTuple is incompatible with boundary" << endl;
		exit(1);
	}

	// Need to set the branch addresses
	vector<TBranch *> branches( numberOfObservables );
	vector<Float_t> observableValues( numberOfObservables );
	for ( int obsIndex = 0; obsIndex < numberOfObservables; obsIndex++ )
	{
		ntuple->SetBranchAddress(observableNames[obsIndex].c_str(), &(observableValues[obsIndex]), &(branches[obsIndex]));
	}

	// Now populate the dataset
	int numberOfDataPointsAdded = 0;
	int numberOfDataPointsRead = 0;
	while ( numberOfDataPointsAdded < numberEventsToRead && numberOfDataPointsAdded < numberOfEventsAfterCut )
	{
		if ( numberOfDataPointsRead > numberOfEventsAfterCut) break;
		DataPoint point( observableNames );
		ntuple->GetEntry( evtList->GetEntry(numberOfDataPointsRead) );

		for ( int obsIndex = 0; obsIndex < numberOfObservables; obsIndex++ )
		{
			string name = observableNames[obsIndex];
			string unit = data->GetBoundary()->GetConstraint( name )->GetUnit();
			point.SetObservable( name, observableValues[obsIndex], 0.0, unit);
		}
		bool dataPointAdded = data->AddDataPoint( &point );
		if (dataPointAdded) numberOfDataPointsAdded += 1;
		numberOfDataPointsRead += 1;
	}

	inputFile->Close();
	delete inputFile;
	cout << "Added " << numberOfDataPointsAdded << " events from ROOT file: " << fileName << " which are consistent with the PhaseSpaceBoundary" << endl;
	time_t timeNow;
        time(&timeNow);
        cout << "Time: " << ctime( &timeNow ) << endl;
	return data;
}

bool DataSetConfiguration::CheckTNtupleWithBoundary( TNtuple * TestTuple, PhaseSpaceBoundary * TestBoundary )
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

IDataSet * DataSetConfiguration::LoadAsciiFileIntoMemory( string fileName, long numberEventsToRead, PhaseSpaceBoundary * DataBoundary )
{
	MemoryDataSet * data = new MemoryDataSet(DataBoundary);
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

			//Check for EOF
			if ( line == "" && file_to_read.eof() )
			{
				cout << "End of data file reached after loading only " << numberOfDataPoints << " data points" << endl;
				file_to_read.close();
				return data;
			}

			//Split line on white space
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
		return data;
	}
	else
	{
		cerr << "Failed to open input data file: " << fileName << endl;
		exit(1);
	}
}
