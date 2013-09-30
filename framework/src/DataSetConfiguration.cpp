/*!
 * @class DataSetConfiguration
 *
 * @brief A class for holding the data to create a data set, and creating that data set when requested
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Greig A Cowan greig.cowan@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

///	ROOT Headers
#include "TFile.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TEventList.h"
#include "TCanvas.h"
#include "TString.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TTreeFormula.h"
///	RapidFit Headers
#include "StringProcessing.h"
#include "MemoryDataSet.h"
#include "DataSetConfiguration.h"
#include "ClassLookUp.h"
#include "ResultFormatter.h"
///	System Headers
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <sstream>

using namespace::std;

//Constructor with correct argument
DataSetConfiguration::DataSetConfiguration( string DataSource, long DataNumber, string cut, vector<string> DataArguments, vector<string> DataArgumentNames, int starting_entry, PhaseSpaceBoundary* Boundary ) :
	source(DataSource), cutString(cut), numberEvents(DataNumber), arguments(DataArguments), argumentNames(DataArgumentNames),
	generatePDF(NULL), separateGeneratePDF(false), parametersAreSet(false), Start_Entry(starting_entry), DEBUG_DATA(false), internalBoundary(NULL),
	internalRef(NULL), debug( new DebugClass(false) ), fileName("undefined")
{
	if( Boundary != NULL )
	{
		internalBoundary = new PhaseSpaceBoundary(*Boundary);
		//internalBoundary->Print();
	}
}

//Constructor with separate data generation PDF
DataSetConfiguration::DataSetConfiguration( string DataSource, long DataNumber, string cut, vector<string> DataArguments, vector<string> DataArgumentNames, IPDF * DataPDF, PhaseSpaceBoundary* Boundary ) :
	source(DataSource), cutString(cut), numberEvents(DataNumber), arguments(DataArguments), argumentNames(DataArgumentNames),
	generatePDF( ClassLookUp::CopyPDF(DataPDF) ), separateGeneratePDF(true), parametersAreSet(false), Start_Entry(0), DEBUG_DATA(false), internalBoundary(NULL),
	internalRef(NULL), debug( new DebugClass(false) ), fileName("undefined")
{
	if( Boundary != NULL )
	{
		internalBoundary = new PhaseSpaceBoundary(*Boundary);
		//internalBoundary->Print();
	}
}

DataSetConfiguration::DataSetConfiguration ( const DataSetConfiguration& input ) :
	source(input.source), cutString(input.cutString), numberEvents(input.numberEvents), arguments(input.arguments), argumentNames(input.argumentNames),
	generatePDF( (input.generatePDF==NULL)?NULL:ClassLookUp::CopyPDF(input.generatePDF) ), separateGeneratePDF(input.separateGeneratePDF), parametersAreSet(input.parametersAreSet),
	Start_Entry(input.Start_Entry), DEBUG_DATA(input.DEBUG_DATA), internalBoundary(NULL), internalRef(NULL), debug( new DebugClass(*input.debug)), fileName( input.fileName )
{
	if( input.internalBoundary != NULL )
	{
		internalBoundary = new PhaseSpaceBoundary( *input.internalBoundary );
		//internalBoundary->Print();
	}
}

//Destructor
DataSetConfiguration::~DataSetConfiguration()
{
	if( generatePDF != NULL ) delete generatePDF;
	if( internalBoundary != NULL ) delete internalBoundary;
	if( debug != NULL ) delete debug;
}

IPDF* DataSetConfiguration::GetGenerationPDF() const
{
	return generatePDF;
}

//Set the parameters of the generation PDF
bool DataSetConfiguration::SetPhysicsParameters( ParameterSet* InputParameters )
{
	if( separateGeneratePDF )
	{
		generatePDF->UpdatePhysicsParameters( InputParameters );
		parametersAreSet = true;
		return true;
	}
	else
	{
		parametersAreSet = true;
		return true;
	}
}

void DataSetConfiguration::SetDebug( bool Flag )
{
	DEBUG_DATA = Flag;
}

bool DataSetConfiguration::SetSource( string NewSource )
{
	source = NewSource;
	return true;
}

string DataSetConfiguration::GetSource() const
{
	return source;
}

//Create the DataSet
IDataSet * DataSetConfiguration::MakeDataSet( PhaseSpaceBoundary * DataBoundary, IPDF * FitPDF, int real_numberEvents )
{
	if(DataBoundary != NULL )
	{
		if( internalBoundary != NULL ) delete internalBoundary;
		internalBoundary = new PhaseSpaceBoundary( *DataBoundary );
	}

	//Some kind of decision about what kind of data set to use?
	IDataSet * newDataSet = NULL;
	if( real_numberEvents!=-1 ) numberEvents = real_numberEvents;

	if( source == "File" )
	{
		//Load data from file
		newDataSet = this->LoadDataFile( arguments, argumentNames, internalBoundary, numberEvents );
	}
	else if( source == "FoamFile" )
	{
		newDataSet = this->FoamFile( arguments, argumentNames, internalBoundary, (int)numberEvents, FitPDF );
	}
	else
	{
		newDataSet = this->LoadGeneratorDataset( source, internalBoundary, (int)numberEvents, FitPDF );
	}

	internalRef = newDataSet;
	return newDataSet;
}

IDataSet* DataSetConfiguration::LoadGeneratorDataset( string Source, PhaseSpaceBoundary* InternalBoundary, int NumberEvents, IPDF* FitPDF )
{
	//Assume it's an accept/reject generator, or some child of it
	IDataGenerator * dataGenerator = NULL;
	if (separateGeneratePDF)
	{
		dataGenerator = ClassLookUp::LookUpDataGenerator( Source, InternalBoundary, generatePDF );
	}
	else
	{
		dataGenerator = ClassLookUp::LookUpDataGenerator( Source, InternalBoundary, FitPDF );
	}
	if( dataGenerator == NULL )
	{
		cerr << "Generator NOT found!" << endl;
		exit(-9864);
	}
	dataGenerator->GenerateData( (int)NumberEvents );
	IDataSet* newDataSet = dataGenerator->GetDataSet();
	delete dataGenerator;

	return newDataSet;
}

IDataSet* DataSetConfiguration::FoamFile( vector<string> Arguments, vector<string> ArgumentNames, PhaseSpaceBoundary* InternalBoundary, int NumberEvents, IPDF* FitPDF )
{
	string searchName = "FileName";
	int fileNameIndex = StringProcessing::VectorContains( &ArgumentNames, &searchName );
	if( fileNameIndex == -1 )
	{
		fileName = "tempFile.root";
		ArgumentNames.push_back( "FileName" );
		Arguments.push_back( fileName );
	}
	else
	{
		fileName = Arguments[(unsigned)fileNameIndex];
	}

	string Source="Foam";

	IDataSet* FoamDataSet = this->LoadGeneratorDataset( Source, InternalBoundary, NumberEvents, FitPDF );

	vector<IDataSet*> TotalFoamDataSet(1, FoamDataSet);

	remove( fileName.c_str() );

	ResultFormatter::MakeRootDataFile( fileName, TotalFoamDataSet );

	delete FoamDataSet;

	IDataSet* FileDataSet = this->LoadDataFile( Arguments, ArgumentNames, InternalBoundary, NumberEvents );

	remove( fileName.c_str() );

	TotalFoamDataSet[0] = FileDataSet;

	if( fileNameIndex != -1 )
	{
		ResultFormatter::MakeRootDataFile( fileName, TotalFoamDataSet );
	}

	return FileDataSet;
}

string DataSetConfiguration::getNtuplePath( string this_fileName )
{
	TFile* temp_file = TFile::Open( this_fileName.c_str(), "READ" );

	if( temp_file == NULL )
	{
		cerr << "Can't open File: " << this_fileName << endl << endl;
		exit(8547);
	}

	TString current_path = gDirectory->GetPath();

	string tuple_path;

	vector<pair<string, string> > found_names;

	get_object_list( current_path, &found_names, TTree::Class(), new TString(current_path) );

	tuple_path = found_names[0].first + "/" + found_names[0].second;

	temp_file->Close();

	return tuple_path;
}

//  From the root current path look for all keys (objects in root file) and loop over them
//  For each one that is actually an object of inherit_type store it's name and the number of events it has
//  Slightly stupid and ineffient as it always has to decend into the file stuture
//  This shouldn't be a huge problem as I never anticipate finding files with hundereds of files
void DataSetConfiguration::get_object_list( TString current_path, vector<pair<string,string> > *found_names, TClass* inherit_type, TString *relative_path )
{
	//  goto current search path
	gDirectory->cd( current_path );
	TDirectory *current_path_directory = gDirectory;
	TString current_path_str( current_path );

	//  get all keys present in current path
	TIter nextkey( current_path_directory->GetListOfKeys() );
	TKey *key, *oldkey=0;

	//cout << current_path << endl;
	while ( ( key = (TKey*) nextkey() ) )
	{
		//  Loop over all keys (objects) in the current path
		//  Leave on the last key

		if ( oldkey && !strcmp( oldkey->GetName(), key->GetName() ) ) continue;
		//  Should already be here?
		//gDirectory->cd( current_path );
		//  Read in the object
		TObject *obj = key->ReadObj();

		//  Get the name of the object in current directory
		string obj_name( (char*) obj->GetName() );
		TString temp( obj_name );

		//  If this is a Tuple object we want to store some info about it
		if ( obj->IsA()->InheritsFrom( inherit_type ) )
		{
			//  Store the full name of the object and it's full path to determine if we've already found it
			string full_name( relative_path->Data() );
			full_name.append( obj_name );

			vector<string> all_found_names;

			for( vector<pair<string,string> >::iterator found_i=found_names->begin(); found_i!=found_names->end(); ++found_i )
			{
				string this_full_name = found_i->first;
				this_full_name.append( string("/") );
				this_full_name.append( found_i->second );
				all_found_names.push_back( this_full_name );
			}

			//  If this is an object of the wanted type we foud ignore it
			if( StringProcessing::VectorContains( &all_found_names, &full_name ) == -1 )
			{
				//  If this is a new class object save it's name and event_number
				//relative_path->Append( obj_name );
				//found_names->push_back( (*relative_path).Data() );
				pair<string,string> path_and_name;
				path_and_name.first = string( (*relative_path).Data() );
				path_and_name.second = obj_name;
				found_names->push_back( path_and_name );
				continue;
			}
		} else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {

			//  If this is a directory we want to decend down into it
			current_path_str.Append( obj_name );
			current_path_str.Append( "/" );
			relative_path->Append( obj_name );
			relative_path->Append( "/" );
			get_object_list( current_path_str, found_names, inherit_type, relative_path );
		}
	}
}


//Constructor with correct arguments
IDataSet * DataSetConfiguration::LoadDataFile( vector<string> Arguments, vector<string> ArgumentNames, PhaseSpaceBoundary * DataBoundary, long NumberEventsToRead )
{
	//Find file name
	string searchName = "FileName";
	int fileNameIndex = StringProcessing::VectorContains( &ArgumentNames, &searchName );
	if ( fileNameIndex >= 0 )
	{
		fileName = Arguments[unsigned(fileNameIndex)];
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
		nTuplePath = Arguments[unsigned(nTuplePathIndex)];
	}
	else
	{
		nTuplePath = this->getNtuplePath( fileName );
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

IDataSet * DataSetConfiguration::LoadRootFileIntoMemory( string this_fileName, string ntuplePath, long numberEventsToRead, PhaseSpaceBoundary * DataBoundary )
{
	MemoryDataSet * data = new MemoryDataSet(DataBoundary);
	vector<string> observableNames = DataBoundary->GetAllNames();
	int numberOfObservables = int(observableNames.size());

	TFile * inputFile = new TFile( this_fileName.c_str(), "READ" );
	TTree * ntuple = (TTree*)inputFile->Get( ntuplePath.c_str() );
	if( ntuple == NULL )
	{
		cerr << "Ntuple not found. Are you specifying the correct file and ntuple path?" << endl;
		exit(2374);
	}
	if( Start_Entry != 0 ) cout << "Starting From Entry: " << Start_Entry<< " in the ntuple." << endl;
	int totalNumberOfEvents = int(ntuple->GetEntries());
	if( totalNumberOfEvents <= 0 )
	{
		cerr << "\t\tInvalid number of Events! exiting" << endl << endl;
		exit(-2374);
	}
	ntuple->SetEstimate(ntuple->GetEntries());  // Fix the size of the array of doubles to be created (There will never be more than this many)
	int numberOfEventsAfterCut = int(ntuple->Draw(">>evtList", cutString.c_str(),"goff",ntuple->GetEntries(),Start_Entry)); // apply the cut which automatically creates the evtList object
	if ( numberOfEventsAfterCut == -1 )
	{
		cerr << "Please check the cut string you are using!" << endl;
		exit(97823);
	}
	//	TEventList * evtList = (TEventList*)gDirectory->Get("evtList"); // ROOT is weird

	cout << "Total number of events in file: " << totalNumberOfEvents << endl;
	cout << "You have applied this cut to the data: '" << cutString << "'" << endl;
	cout << "Total number of events after cut: " << numberOfEventsAfterCut << endl;

	//  Container for all of the data read in from a root file
	vector<vector<Double_t> > real_data_array;
	//real_data_array.reserve( unsigned(numberOfObservables) );


	//time_t timeNow1;
	//time(&timeNow1);
	//cout << "Time Before Draw Read: " << ctime( &timeNow1 ) << endl;

	//  Instead of Reading in in Row form as with events in NTuples
	//  (Convenient for many columns and few Rows)
	//  Read in the data in Column form
	//

	//	This has been rewritten at the behest of ROOT BEING CRAP
	//
	//	However...
	//
	//	This does allow the oportunity to 'check the data as it is read into RapidFit :D

	//  Hold the Data in a temp object
	vector<Double_t > temp_vector;
	TString PlotString("");
	TString Name("");
	TString Plot_Options = "";
	if( DEBUG_DATA )
	{
		Plot_Options = "";
	} else {
		Plot_Options = "goff";
	}

	//DataBoundary->Print();

	TCanvas* bob = new TCanvas( Name, Name, 1680, 1050 );

	TString FormulaName="Fomula_";
	for( int obsIndex = 0; obsIndex < numberOfObservables; ++obsIndex )
	{
		//	If we want to debug the selection construct a canvas
		Name="Canvas_Name_"; Name+=obsIndex;

		//  Construct a Plot String to use the TTree->Draw Method

		ObservableRef refName = ObservableRef( observableNames[unsigned(obsIndex)] );

		IConstraint* this_const = DataBoundary->GetConstraint( refName );
		if( this_const != NULL )
		{
			PlotString = "("+this_const->GetTF1()+")";
		}
		else
		{
			cerr << "CANNOT FIND CONSTRAINT: " << string(refName) << endl;
			exit(-8734);
		}

		TString thisFormulaName = FormulaName; thisFormulaName+=obsIndex;
		TTreeFormula* tempFormula = new TTreeFormula( thisFormulaName, PlotString, ntuple );

		//	We may one day come across a situation where we need to evaluate a 0-dim Observable (but your likely being a lazy developer)
		//	I have had to make the choice here to check for 0-dim as the Init code is too difficult/under-documented for me to parse
		//	If you have 0-dim it is LIKELY your formula is invalid
		//	If you have a NULL internal pointer to a Tree your forumla IS invalid!
		//	ROOT did NOT code up a perfect way of checking for these things
		if( (tempFormula->GetTree() == NULL) || (tempFormula->GetNdim() == 0 ) )
		{
			cerr << "Error evaluating: " << PlotString << " within your NTuple " << ntuple->GetName() << " titled: " << ntuple->GetTitle() << endl;
			cerr << "This was for file: " << this_fileName << endl;
			cerr << "Cannot process DataSet. exiting!" << endl << endl;
			exit(-765);
		}

		//  Draw upto 3 observables at a time in some large plot
		//  use the 'goff' option to turn graphical output (and annoying text output from default co/de-structors) off
		//  (it doesn't matter what this looks like and we can throw it away from here)

		ntuple->Draw( PlotString , cutString.c_str(), Plot_Options, Long64_t(ntuple->GetEntries()), Long64_t(Start_Entry) );

		//  Save the actual data to somewhere protected in our memory
		Double_t* values = ntuple->GetV1();

		//	Use object recasting to make a copy of the data
		for(int j=0; j < numberOfEventsAfterCut; ++j )	temp_vector.push_back( (double)(values[j]) );


		if( DEBUG_DATA )
		{
			bob->Update();
			bob->Print(TString("Observable_"+observableNames[(unsigned)obsIndex]+"_selected.png"));
		}

		//	Store the data
		real_data_array.push_back( vector<double>(temp_vector) );

		//	Cleanup
		vector<Double_t> empty;
		temp_vector.swap(empty);

		delete tempFormula;
	}
	delete bob;
	//time_t timeNow2;
	//time(&timeNow2);
	//cout << "Time After Draw Read: " << ctime( &timeNow2 ) << endl;


	// Now populate the dataset
	int numberOfDataPointsAdded = 0;
	int numberOfDataPointsRead = 0;

	//  Now we have all of the data stored in memory in real_data_array which has a 1<->1 with observableName
	//  Create and store data points for each event as before and throw away events outside of the PhaseSpace
	for( ; (numberOfDataPointsRead < numberOfEventsAfterCut) && (numberOfDataPointsAdded < numberEventsToRead) ; ++numberOfDataPointsRead )
	{
		DataPoint* point = new DataPoint( observableNames );
		for(int obsIndex = 0; obsIndex < numberOfObservables; ++obsIndex )
		{
			string name = observableNames[unsigned(obsIndex)];
			string unit = data->GetBoundary()->GetConstraint( name )->GetUnit();
			point->SetObservable( name, real_data_array[unsigned(obsIndex)][unsigned(numberOfDataPointsRead)], unit, true, obsIndex);
			//cout << real_data_array[unsigned(obsIndex)][unsigned(numberOfDataPointsRead)] << endl;
		}
		bool dataPointAdded = data->AddDataPoint( point );
		if (dataPointAdded) ++numberOfDataPointsAdded;
	}

	if( DEBUG_DATA )
	{
		data->Print();
	}

	inputFile->Close();
	delete inputFile;
	cout << "Added " << numberOfDataPointsAdded << " events from ROOT file: " << this_fileName << " which are consistent with the PhaseSpaceBoundary" << endl;
	time_t timeNow;
	time(&timeNow);
	cout << "Time: " << ctime( &timeNow );
	while( !real_data_array.empty() ){ while( !real_data_array.back().empty() ){real_data_array.back().pop_back(); }; real_data_array.pop_back(); }
	return data;
}

IDataSet * DataSetConfiguration::LoadAsciiFileIntoMemory( string this_fileName, long numberEventsToRead, PhaseSpaceBoundary * DataBoundary )
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
	std::ifstream file_to_read( this_fileName.c_str() );
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
					point.SetObservable( name, strtod( splitVec[i].c_str(), NULL ), unit );
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
						cerr << "Observable " << observableNames[observableIndex] << " not provided by file " << this_fileName << endl;
						exit(1);
					}
				}

				readFirstLine = true;
			}
		}
		file_to_read.close();
		cout << "Read " << numberOfDataPoints << " events from file: " << this_fileName << endl;
		return data;
	}
	else
	{
		cerr << "Failed to open input data file: " << this_fileName << endl;
		exit(1);
	}
}

string DataSetConfiguration::XML() const
{
	stringstream xml;

	xml << "<DataSet>" << endl;
	if( source == "Foam" )
	{
		xml << generatePDF->XML();
	}
	xml << "\t" << "<Source>" << source << "</Source>" << endl;
	if( source == "File" || source == "FoamFile" )
	{
		if( !cutString.empty() )
		{
			xml << "\t" << "<CutString>"  << cutString << "</CutString>" << endl;
		}
		if( !fileName.empty() )
		{
			xml << "\t" << "<FileName>" << fileName << "</FileName>" << endl;
		}
	}
	xml << "\t<NumberEvents>";
	if( internalRef != NULL ) xml << internalRef->GetDataNumber() << "</NumberEvents>" << endl;
	else xml << 10000 << "</NumberEvents>" << endl;
	if( source != "File" && generatePDF != NULL ) xml << generatePDF->XML() << endl;
	if( internalBoundary != NULL ) xml << internalBoundary->XML() << endl;
	xml << "</DataSet>" << endl;

	return xml.str();
}

void DataSetConfiguration::SetDebug( DebugClass* input_debug )
{
	if( debug != NULL ) delete debug;
	debug = new DebugClass( *input_debug );
}

PhaseSpaceBoundary* DataSetConfiguration::GetPhaseSpace() const
{
	return internalBoundary;
}

int DataSetConfiguration::GetDataSetSize() const
{
	return (int)numberEvents;
}

