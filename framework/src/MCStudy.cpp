//	MCStudy Class to sequentially step through a data file and perform fits on subsets

//	RapidFit Headers
#include "MCStudy.h"
#include "XMLConfigReader.h"
#include "FitResultVector.h"
#include "FitAssembler.h"
#include "PDFWithData.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

//	Destructor
MCStudy::~MCStudy()
{
}

//	XML Constructor
MCStudy::MCStudy( XMLConfigReader* new_xmlConfig ) :
IStudy(), events_to_step_over( new_xmlConfig->GetAllDataSetSizes() ), CommandLineParams(), StartingEntries( new_xmlConfig->GetAllStartEntries() ), PhysParams(), input_xml(new_xmlConfig)
{
	xmlConfig = new_xmlConfig;
	numberStudies = xmlConfig->GetNumberRepeats();
	theMinimiser = xmlConfig->GetMinimiserConfiguration();
	theFunction = xmlConfig->GetFitFunctionConfiguration();
	studyParameters = xmlConfig->GetFitParameters( CommandLineParams );
	allConstraints = xmlConfig->GetConstraints();
	delete_objects = true;
}

//	XML Constructor with events_to_step_over, num_repeats & starting_entry defined
MCStudy::MCStudy( XMLConfigReader* new_xmlConfig, vector<int> new_events_to_step_over, int new_num_repeats, vector<int> starting_entries ) :
IStudy(), events_to_step_over( new_events_to_step_over ), CommandLineParams(), StartingEntries( starting_entries ), PhysParams(), input_xml(new_xmlConfig)
{
	xmlConfig = new_xmlConfig;
	numberStudies = new_num_repeats;
	theMinimiser = xmlConfig->GetMinimiserConfiguration();
	theFunction = xmlConfig->GetFitFunctionConfiguration();
	studyParameters = xmlConfig->GetFitParameters( CommandLineParams );
	allConstraints = xmlConfig->GetConstraints();
	delete_objects = true;
}

//	Set the number of repeats
void MCStudy::SetNumRepeats( int new_num_repeats )
{
	numberStudies = new_num_repeats;
}

//	Set the number of repeats
void MCStudy::SetNumEvents( vector<int> new_events_to_step_over )
{
	events_to_step_over = new_events_to_step_over;
}

//	Set the number of repeats
void MCStudy::SetStartingEntry( vector<int> new_starting_entries )
{
	StartingEntries = new_starting_entries;
}

void MCStudy::SetPhysParams( vector<string> new_Phys_params )
{
	PhysParams = new_Phys_params;
}

//	Set the command line physics parameters
void MCStudy::SetCommandLineParams( vector<string> new_CommandLineParams )
{
	CommandLineParams = new_CommandLineParams;
}

//	Perform the Study
void MCStudy::DoWholeStudy( int OutputLevel )
{
	//		THIS WILL NEED TO BE CHANGED IF WE MOVE TO MULTIPLE PARAMETERSETS IN THE XML STRUCTURE
	allResults = new FitResultVector( xmlConfig->GetFitParameters( CommandLineParams )->GetAllNames() );

	if( ( StartingEntries.size() != 1) && ( StartingEntries.size() != xmlConfig->GetAllStartEntries().size() ) )
	{
		cerr << "The number of Starting Entries passed is neither 1 nor same as the number present in the XML" << endl;
		cerr << "I dont' know what to do!" << endl;	return;
	}
	if( ( events_to_step_over.size() != 1) && ( events_to_step_over.size() != xmlConfig->GetAllDataSetSizes().size()))
	{
		cerr << "The number of Step_Sizes passed is neither 1 nor the same as the number present in the XML" <<endl;
		cerr << "I don't know what to do!" << endl;	return;
	}
	if( events_to_step_over.size() != StartingEntries.size() )
	{
		cerr << "Sorry don't quite know how this happened!" << endl; return;
	}

	for( int i=0; i < numberStudies; ++i )
	{
		vector<int> local_Start_entries;
		vector<int>::iterator start_i=StartingEntries.begin();
		vector<int>::iterator entries_i=events_to_step_over.begin();
		for( ; ( start_i!=StartingEntries.end() && entries_i != events_to_step_over.end() );  )
		{
			//	Written such that this always addresses the correct element
			local_Start_entries.push_back( (*start_i) + (i)* (*entries_i) );

			cout << (*start_i) << "\t" <<(*entries_i) << "\t" << local_Start_entries.back() << endl;

			//	We either have 1 step/start for the whole xml or we have a unique step/start for each PDFWithData from the XML
			if( events_to_step_over.size() != 1 )	{	++entries_i;	} else { entries_i = events_to_step_over.end(); }
			if( StartingEntries.size() != 1 )	{	++start_i;	} else { start_i = StartingEntries.end(); }
		}

		if( StartingEntries.size() == 1 )
		{
			for( unsigned int j=0; j < xmlConfig->GetAllDataSetSizes().size(); ++j )
			{
				local_Start_entries.push_back( local_Start_entries[0] );
			}
		}

		vector<PDFWithData*> local_PDF_w_Data = xmlConfig->GetPDFsAndData( local_Start_entries );

		allResults->StartStopwatch();

		FitResult * newResult = FitAssembler::DoSafeFit( theMinimiser, theFunction, studyParameters, local_PDF_w_Data, allConstraints, OutputLevel );

		allResults->AddFitResult( newResult );

	}

}

FitResultVector* MCStudy::GetStudyResult()
{
	return allResults;
}

