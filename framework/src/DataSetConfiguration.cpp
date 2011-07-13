/**
  @class DataSetConfiguration

  A class for holding the data to create a data set, and creating that data set when requested

  @author Benjamin M Wynne bwynne@cern.ch
  @author Greig A Cowan greig.cowan@cern.ch
  @data 2009-12-16
 */

//	ROOT Headers
#include "TFile.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TEventList.h"
#include "TCanvas.h"
#include "TString.h"
#include "TTree.h"
//	RapidFit Headers
#include "StringProcessing.h"
#include "MemoryDataSet.h"
#include "DataSetConfiguration.h"
#include "ClassLookUp.h"
#include "ResultFormatter.h"
//	System Headers
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <stdlib.h>

//Default constructor
DataSetConfiguration::DataSetConfiguration() : source(), cutString(), numberEvents(), arguments(), argumentNames(), generatePDF(), separateGeneratePDF(), parametersAreSet(), Start_Entry(0), DEBUG_DATA(false)
{
}

//Constructor with correct argument
DataSetConfiguration::DataSetConfiguration( string DataSource, long DataNumber, string cut, vector<string> DataArguments, vector<string> DataArgumentNames, int starting_entry ) : source(DataSource), cutString(cut), numberEvents(DataNumber), arguments(DataArguments), argumentNames(DataArgumentNames), generatePDF(NULL), separateGeneratePDF(false), parametersAreSet(false), Start_Entry(starting_entry), DEBUG_DATA(false)
{
}

//Constructor with separate data generation PDF
DataSetConfiguration::DataSetConfiguration( string DataSource, long DataNumber, string cut, vector<string> DataArguments, vector<string> DataArgumentNames, IPDF * DataPDF ) : source(DataSource), cutString(cut), numberEvents(DataNumber), arguments(DataArguments), argumentNames(DataArgumentNames), generatePDF(DataPDF), separateGeneratePDF(true), parametersAreSet(false), Start_Entry(0), DEBUG_DATA(false)
{
}

//Destructor
DataSetConfiguration::~DataSetConfiguration()
{
}

IPDF* DataSetConfiguration::GetGenerationPDF()
{
	return generatePDF;
}

//Set the parameters of the generation PDF
bool DataSetConfiguration::SetPhysicsParameters( vector<ParameterSet*> InputParameters )
{
	if (separateGeneratePDF)
	{
		if ( generatePDF->SetPhysicsParameters(InputParameters.back()) )
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

void DataSetConfiguration::SetDebug( bool Flag )
{
	DEBUG_DATA = Flag;
}

bool DataSetConfiguration::SetSource( string NewSource )
{
	source = NewSource;
	return true;
}

string DataSetConfiguration::GetSource()
{
	return source;
}

//Create the DataSet
IDataSet * DataSetConfiguration::MakeDataSet( PhaseSpaceBoundary * DataBoundary, IPDF * FitPDF, int real_numberEvents )
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
	vector<string> observableNames = DataBoundary->GetAllNames();
	int numberOfObservables = int(observableNames.size());

	TFile * inputFile = new TFile( fileName.c_str(), "READ" );
	TNtuple * ntuple = (TNtuple*)inputFile->Get( ntuplePath.c_str() );
	if ( ntuple == NULL )
	{
		cerr << "Ntuple not found. Are you specifying the correct file and ntuple path?" << endl;
		exit(2374);
	}
	if( Start_Entry != 0 ) cout << "Starting From Entry: " << Start_Entry<< " in the ntuple." << endl;
	int totalNumberOfEvents = int(ntuple->GetEntries());
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

	// Check that ntuple is compatible with the PhaseSpaceBoundary that is defined
	// in the XML file.
	if ( !CheckTNtupleWithBoundary( ntuple, data->GetBoundary() ) )
	{
		cerr << "NTuple is incompatible with boundary" << endl;
		cerr << "Does this file contain all of the Parameters in your PhaseSpace?" << endl;
		exit(1);
	}

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

	TCanvas* bob = new TCanvas( Name, Name, 1680, 1050 );

	for ( int obsIndex = 0; obsIndex < numberOfObservables; ++obsIndex )
	{
		//	If we want to debug the selection construct a canvas
		Name="Canvas_Name_"; Name+=obsIndex;

		//  Construct a Plot String to use the TTree->Draw Method
		PlotString = observableNames[unsigned(obsIndex)];

		//  Draw upto 3 observables at a time in some large plot
		//  use the 'goff' option to turn graphical output (and annoying text output from default co/de-structors) off
		//  (it doesn't matter what this looks like and we can throw it away from here)

		ntuple->Draw( PlotString , cutString.c_str(), Plot_Options, Long64_t(ntuple->GetEntries()), Long64_t(Start_Entry) );

		//  Save the actual data to somewhere protected in our memory

		//	Use object recasting to make a copy of the data
		for(int j=0; j < numberOfEventsAfterCut; ++j )
		{
			temp_vector.push_back( double(ntuple->GetV1()[j]) );
		}

                if( DEBUG_DATA )
		{
			bob->Update();
			bob->Print(TString("Observable_"+observableNames[obsIndex]+"_selected.png"));
		}

		//	Store the data
		real_data_array.push_back( temp_vector );

		//	Cleanup
		while( !temp_vector.empty() ) { temp_vector.pop_back(); }
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
	if ( int(numberEventsToRead) < int(numberOfEventsAfterCut) )	data->ReserveDataSpace( int(numberEventsToRead) );
	else	data->ReserveDataSpace( int(numberOfEventsAfterCut) );
	for( ; (numberOfDataPointsRead < numberOfEventsAfterCut) && (numberOfDataPointsAdded < numberEventsToRead) ; ++numberOfDataPointsRead )
	{
		DataPoint point( observableNames );
		for(int obsIndex = 0; obsIndex < numberOfObservables; ++obsIndex )
		{
			string name = observableNames[unsigned(obsIndex)];
			string unit = data->GetBoundary()->GetConstraint( name )->GetUnit();
			point.SetObservable( name, real_data_array[unsigned(obsIndex)][unsigned(numberOfDataPointsRead)], 0.0, unit, true, obsIndex);
			//cout << real_data_array[unsigned(obsIndex)][unsigned(numberOfDataPointsRead)] << endl;
		}
		bool dataPointAdded = data->AddDataPoint( &point );
		if (dataPointAdded) ++numberOfDataPointsAdded;
	}
	
	inputFile->Close();
	delete inputFile;
	cout << "Added " << numberOfDataPointsAdded << " events from ROOT file: " << fileName << " which are consistent with the PhaseSpaceBoundary" << endl;
	time_t timeNow;
        time(&timeNow);
        cout << "Time: " << ctime( &timeNow );
	while( !real_data_array.empty() ){ while( !real_data_array.back().empty() ){real_data_array.back().pop_back(); }; real_data_array.pop_back(); }
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
		TIter observableIterator( TestTuple->GetListOfBranches() );
		TBranch * observableBranch;
		while( ( observableBranch = (TBranch*)observableIterator() ) )
		{
			string BranchName = observableBranch->GetName();
			if ( BranchName == allNames[observableIndex] )
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

FitResultVector* DataSetConfiguration::LoadFitResult( TString input_file, ParameterSet* ParametersFromXML )
{
	TFile* input_LL = new TFile( input_file, "READ" );
	//	FitResult file format is ALWAYS for now until the end of time defined as:
	//	TTree Element 0, 2, 4, ....... 0+2n (int)n are ALWAYS IDENTICAL AND THE GLOBAL FIT RESULT!!!
	//	FitResults 1, 3, 5, 7 ....... 0+n (int)n   ARE THE INPUT FITRESULTS

	TTree* input_tree = (TTree*) gDirectory->Get( "RapidFitResult" );

	vector<TString> all_branches = ResultFormatter::get_branch_names( input_tree );

	vector<TString> all_parameter_values = StringProcessing::GetStringContaining( all_branches, TString("_value") );

	vector<TString> all_parameters = StringProcessing::StripStrings( all_parameter_values, TString("_value") );

	vector<Float_t> all_values; all_values.resize( all_parameters.size() );
	vector<Float_t> all_generated; all_generated.resize( all_parameters.size() );
	vector<Float_t> all_errors; all_errors.resize( all_parameters.size() );
	vector<Float_t> all_min; all_min.resize( all_parameters.size() );
	vector<Float_t> all_max; all_max.resize( all_parameters.size() );
	vector<Float_t> all_steps; all_steps.resize( all_parameters.size() );
	Float_t Fit_Status=-1;
	Float_t NLL=-9999;
	Float_t CPUTime=-1;
	Float_t RealTime=-1;

	TString value, gen, error, min_str, max_str, step_str;

	for( unsigned int i=0; i< all_parameter_values.size(); ++i )
	{
		//cout << all_parameters[i] << endl;

		value = all_parameters[i] + "_value";
		gen = all_parameters[i] + "_gen";
		error = all_parameters[i] + "_error";
		min_str = all_parameters[i] + "_min";
		max_str = all_parameters[i] + "_max";
		step_str = all_parameters[i] + "_step";

		input_tree->SetBranchAddress( value, &all_values[i] );
		input_tree->SetBranchAddress( gen, &all_generated[i] );
		input_tree->SetBranchAddress( error, &all_errors[i] );
		input_tree->SetBranchAddress( min_str, &all_min[i] );
		input_tree->SetBranchAddress( max_str, &all_max[i] );
		input_tree->SetBranchAddress( "NLL", &NLL );
		input_tree->SetBranchAddress( "Fit_Status", &Fit_Status );
	}

	ResultParameterSet* global_result = new ResultParameterSet( StringProcessing::Convert( all_parameters ) );
	FitResult* global_fitresult = NULL;
	FitResultVector* global_fitresultvector = NULL;

	for( unsigned int j=0; j< all_parameters.size(); ++j )
	{
		ResultParameter* new_param = new ResultParameter( all_parameters[j].Data(), (double)all_values[j], (double)all_generated[j], (double)all_errors[j],
								(double)all_min[j], (double)all_max[j], (double)all_steps[j],
							ParametersFromXML->GetPhysicsParameter( all_parameters[j].Data() )->GetType(),
							ParametersFromXML->GetPhysicsParameter( all_parameters[j].Data() )->GetUnit() );
		global_result->SetResultParameter( all_parameters[j].Data(), new_param );
		//	Don't worry, the ResultParameter is actually copied into a memory persistant standard class
		delete new_param;
	}

	global_fitresult = new FitResult( (double)NLL, global_result, (int)Fit_Status );
	global_fitresultvector = new FitResultVector( StringProcessing::Convert( all_parameters ) );

	global_fitresultvector->AddFitResult( global_fitresult, false );
	global_fitresultvector->AddCPUTime( (double)CPUTime );
	global_fitresultvector->AddRealTime( (double)RealTime );

	int number_fits = (int)input_tree->GetEntries();
	int number_grid_points = number_fits/2;

	if( number_fits != 2*number_grid_points )
	{
		cerr << "Badly Defined RapidFit File!" << endl;
		exit(-89);
	}

	vector<ResultParameterSet*> all_sets_in_file;
	vector<FitResult*> all_sets_in_file_fitresult;
	vector<FitResultVector*> all_sets_in_file_fitresultvector;

	for( int i=0; i< number_grid_points; ++i )
	{
		input_tree->GetEntry( i*2+1 );

		ResultParameterSet* new_set = new ResultParameterSet( StringProcessing::Convert( all_parameters ) );

		for( unsigned int j=0; j< all_parameters.size(); ++j )
		{
			ResultParameter* new_param = new ResultParameter( all_parameters[j].Data(), (double)all_values[j], (double)all_generated[j], (double)all_errors[j],
									 (double)all_min[j], (double)all_max[j], (double)all_steps[j],
							ParametersFromXML->GetPhysicsParameter( all_parameters[j].Data() )->GetType(),
							ParametersFromXML->GetPhysicsParameter( all_parameters[j].Data() )->GetUnit() ); 
			new_set->SetResultParameter( all_parameters[j].Data(), new_param );
			delete new_param;
		}
		all_sets_in_file.push_back( new_set );

		all_sets_in_file_fitresult.push_back( new FitResult( (double)NLL, all_sets_in_file.back(), (int)Fit_Status ) );
		all_sets_in_file_fitresultvector.push_back( new FitResultVector( StringProcessing::Convert( all_parameters ) ) );

		all_sets_in_file_fitresultvector.back()->AddFitResult( all_sets_in_file_fitresult.back(), false );
		all_sets_in_file_fitresultvector.back()->AddCPUTime( (double)CPUTime );
		all_sets_in_file_fitresultvector.back()->AddRealTime( (double)RealTime );
	}

	//for( unsigned int i=0; i< all_sets_in_file.size(); ++i )
	//{
	//	cout << "SET:\t" << i << endl;
	//	vector<string> all_stored_params = all_sets_in_file[i]->GetAllNames();
	//	for( unsigned j=0; j< all_stored_params.size(); ++j )
	//	{
	//		cout << "Name:\t" << all_stored_params[j].c_str() << "\tValue:\t" << (double)all_sets_in_file[i]->GetResultParameter( all_stored_params[j] )->GetValue() << "\tMax:\t" << (double)all_sets_in_file[i]->GetResultParameter( all_stored_params[j] )->GetMaximum() << "\tMin:\t" << (double)all_sets_in_file[i]->GetResultParameter( all_stored_params[j] )->GetMinimum() << "\tType:\t" << all_sets_in_file[i]->GetResultParameter( all_stored_params[j] )->GetType() << "\tUnit:\t" << all_sets_in_file[i]->GetResultParameter( all_stored_params[j] )->GetUnit() << endl;
	//	}
	//	cout << endl;
	//}

	vector<FitResultVector*> all_vectors;
	all_vectors.push_back( global_fitresultvector );
	for( unsigned int i=0; i< all_sets_in_file_fitresultvector.size(); ++i )
	{
		all_vectors.push_back( all_sets_in_file_fitresultvector[i] );
	}

	FitResultVector* all_Results = new FitResultVector( all_vectors );

	input_LL->Close();	
	return all_Results;
}

