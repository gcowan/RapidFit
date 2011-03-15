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
#include "TCanvas.h"
#include "TString.h"

//Default constructor
DataSetConfiguration::DataSetConfiguration()
{
}

//Constructor with correct argument
DataSetConfiguration::DataSetConfiguration( string DataSource, long DataNumber, string cut, vector<string> DataArguments, vector<string> DataArgumentNames ) : source(DataSource), cutString(cut), numberEvents(DataNumber), arguments(DataArguments), argumentNames(DataArgumentNames), separateGeneratePDF(false), parametersAreSet(false)
{
}

//Constructor with separate data generation PDF
DataSetConfiguration::DataSetConfiguration( string DataSource, long DataNumber, string cut, vector<string> DataArguments, vector<string> DataArgumentNames, IPDF * DataPDF ) : source(DataSource), cutString(cut), numberEvents(DataNumber), arguments(DataArguments), argumentNames(DataArgumentNames), generatePDF(DataPDF), separateGeneratePDF(true), parametersAreSet(false)
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
IDataSet * DataSetConfiguration::MakeDataSet( PhaseSpaceBoundary * DataBoundary, IPDF * FitPDF, int real_numberEvents)
{
	//Some kind of decision about what kind of data set to use?
	IDataSet * newDataSet;
	if( real_numberEvents!=-1 ) numberEvents = real_numberEvents;

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
			dataGenerator->GenerateData( int(numberEvents) );
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
	int numberOfObservables = int(observableNames.size());

	TFile * inputFile = new TFile( fileName.c_str(), "READ" );
	TNtuple * ntuple = (TNtuple*)inputFile->Get( ntuplePath.c_str() );
	if ( ntuple == NULL )
	{
		cerr << "Ntuple not found. Are you specifying the correct file and ntuple path?" << endl;
		exit(1);
	}
	int totalNumberOfEvents = int(ntuple->GetEntries());
	ntuple->SetEstimate(ntuple->GetEntries());  // Fix the size of the array of doubles to be created (There will never be more than this many)
	int numberOfEventsAfterCut = int(ntuple->Draw(">>evtList", cutString.c_str())); // apply the cut which automatically creates the evtList object
	if ( numberOfEventsAfterCut == -1 )
	{
		cerr << "Please check the cut string you are using!" << endl;
		exit(1);
	}
//	TEventList * evtList = (TEventList*)gDirectory->Get("evtList"); // ROOT is weird

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

	//  Container for all of the data read in from a root file
	vector<vector<Double_t> > real_data_array;
	real_data_array.reserve( numberOfObservables );


	//time_t timeNow1;
	//time(&timeNow1);
	//cout << "Time Before Draw Read: " << ctime( &timeNow1 ) << endl;

	//  Instead of Reading in in Row form as with events in NTuples
	//  (Convenient for many columns and few Rows)
	//  Read in the data in Column form
	//
	//  Fundamental logic behind this is reading in 200k events from MC takes forever but plotting these variables in TBrowser is click of a button fast
	//  Use the internal Root Structures to return an array of Double_t objects from the plot (No matter what the Branch Type)
	//  Plus as Tuple structures are disk objects this no longer calls many thousand of reads from a very large file to cache from disk
	//  This now simply reads all the data in a few fast requests which removes the slowest part of startup
	//  This could probably be done with branches and loops but this is gauranteed to be quick as it's used behind the scenes by TBrowser
	//  (why this has to be so hard with root I will never know...)

	//  Hold the Data in a temp object
	vector<Double_t *> data_array;
	TString PlotString("");
	for ( int obsIndex = 0; obsIndex < numberOfObservables; obsIndex+=3 )
	{
		data_array.reserve(3);
		//  Construct a Plot String to use the TTree->Draw Method
		PlotString = observableNames[obsIndex];
		for( short int i=1; ((obsIndex+i)<numberOfObservables)&&((i<3)); ++i )
		{
			PlotString.Append(":");
			PlotString.Append(observableNames[obsIndex+i]);
		}
		//cout << "PlotString: " << PlotString << endl;

		//  Draw 3 observables at a time in some large plot
		//  use the 'goff' option to turn graphical output (and annoying text output from default co/de-structors) off
		//  (it doesn't matter what this looks like and we can throw it away from here)
		ntuple->Draw( PlotString , cutString.c_str(), "goff" );
		
		//  Store pointers to the objects for ease of access
		data_array.push_back( ntuple->GetV1() );
		//  GetV1 returns a pointer to an array of doubles for the first ('X') axis in the plot we've just drawn
		if( (obsIndex+1) < numberOfObservables) {  data_array.push_back( ntuple->GetV2() ); }
		if( (obsIndex+2) < numberOfObservables) {  data_array.push_back( ntuple->GetV3() ); }
		
		//  As You can only plot 3 at a time using this mechanism and Root will overide the results between drawing the plots
		//  Save the actual data to somewhere protected in our memory
		for(unsigned short int i=0; i < data_array.size(); ++i )
		{
			vector<Double_t> temp_vector;
			temp_vector.reserve( numberOfEventsAfterCut );
			for(int j=0; j < numberOfEventsAfterCut; ++j )
			{
				temp_vector.push_back( data_array[i][j] );
			}
			real_data_array.push_back( temp_vector );
		}
		//	Cleanup
		while( !data_array.empty() ) { data_array.pop_back(); }
	}
        //time_t timeNow2;
	//time(&timeNow2);
	//cout << "Time After Draw Read: " << ctime( &timeNow2 ) << endl;


        // Now populate the dataset
        int numberOfDataPointsAdded = 0;
        int numberOfDataPointsRead = 0;
	//  Now we have all of the data stored in memory in real_data_array which has a 1<->1 with observableName
	//  Create and store data points for each event as before and throw away events outside of the PhaseSpace 
	if ( int(numberEventsToRead) < int(numberOfEventsAfterCut*0.75) )	data->ReserveDataSpace( int(numberEventsToRead) );
	else	data->ReserveDataSpace( int(numberOfEventsAfterCut*0.75) );
	for( ; (numberOfDataPointsRead < numberOfEventsAfterCut) && (numberOfDataPointsAdded < numberEventsToRead) ; ++numberOfDataPointsRead )
	{
		DataPoint point( observableNames );
		for(int obsIndex = 0; obsIndex < numberOfObservables; ++obsIndex )
		{
			string name = observableNames[obsIndex];
			string unit = data->GetBoundary()->GetConstraint( name )->GetUnit();
			point.SetObservable( name, real_data_array[obsIndex][numberOfDataPointsRead], 0.0, unit, true, obsIndex);
		}
		bool dataPointAdded = data->AddDataPoint( &point );
		if (dataPointAdded) ++numberOfDataPointsAdded;
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
	for (unsigned int observableIndex = 0; observableIndex < allNames.size(); ++observableIndex )
	{
		bool found = false;
		//Find the requested observable name in the ntuple
		TIter observableIterator( TestTuple->GetListOfLeaves() );
		TLeaf * observableLeaf;
		while( ( observableLeaf = (TLeaf*)observableIterator() ) )
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
	for (unsigned int i = 0; i < observableNames.size(); ++i )
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
				for (unsigned int i = 0; i < splitVec.size(); ++i)
				{
					string name = observableNamesInFile[i];
					string unit = observableNamesToUnits[ name ];
					point.SetObservable( name, strtod( splitVec[i].c_str(), NULL ), 0.0, unit );
				}
				++numberOfDataPoints;
				data->AddDataPoint( &point );
			}
			else
			{
				//Load the obsevable names
				for (unsigned int i = 0; i < splitVec.size(); ++i)
				{
					observableNamesInFile.push_back( splitVec[i] );
				}

				//Check that all required observable names are specified
				for (unsigned int observableIndex = 0; observableIndex < observableNames.size(); ++observableIndex )
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
