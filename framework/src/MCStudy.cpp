//	MCStudy Class to sequentially step through a data file and perform fits on subsets

//	RapidFit Headers
#include "MCStudy.h"
#include "XMLConfigReader.h"
#include "ToyStudyResult.h"
#include "FitAssembler.h"
#include "PDFWithData.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

//	Default
MCStudy::MCStudy() : input_xml(NULL), events_to_step_over(), num_repeats(), CommandLineParams(), StartingEntries(), ALL_Fit_Results(NULL), PhysParams()
{}

//	Destructor
MCStudy::~MCStudy()
{
	if( ALL_Fit_Results != NULL ) delete ALL_Fit_Results;
}

//	XML Constructor
MCStudy::MCStudy( XMLConfigReader* new_input_xml ) : input_xml( new_input_xml ), events_to_step_over( new_input_xml->GetAllDataSetSizes() ), num_repeats( new_input_xml->GetNumberRepeats() ), CommandLineParams(), StartingEntries( new_input_xml->GetAllStartEntries() ), ALL_Fit_Results(NULL), PhysParams()
{}

//	XML Constructor with events_to_step_over, num_repeats & starting_entry defined
MCStudy::MCStudy( XMLConfigReader* new_input_xml, vector<int> new_events_to_step_over, int new_num_repeats, vector<int> starting_entries ) : input_xml( new_input_xml ), events_to_step_over( new_events_to_step_over ), num_repeats( new_num_repeats ), CommandLineParams(), StartingEntries( starting_entries ), ALL_Fit_Results(NULL), PhysParams()
{}

//	Set the number of repeats
void MCStudy::SetNumRepeats( int new_num_repeats )
{
	num_repeats = new_num_repeats;
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
void MCStudy::DoWholeStudy()
{
	//		THIS WILL NEED TO BE CHANGED IF WE MOVE TO MULTIPLE PARAMETERSETS IN THE XML STRUCTURE
	ALL_Fit_Results = new ToyStudyResult( input_xml->GetFitParameters( CommandLineParams )[0]->GetAllNames() );

	if( ( StartingEntries.size() != 1) && ( StartingEntries.size() != input_xml->GetAllStartEntries().size() ) )
	{
		cerr << "The number of Starting Entries passed is neither 1 nor same as the number present in the XML" << endl;
		cerr << "I dont' know what to do!" << endl;	return;
	}
	if( ( events_to_step_over.size() != 1) && ( events_to_step_over.size() != input_xml->GetAllDataSetSizes().size()))
	{
		cerr << "The number of Step_Sizes passed is neither 1 nor the same as the number present in the XML" <<endl;
		cerr << "I don't know what to do!" << endl;	return;
	}
	if( events_to_step_over.size() != StartingEntries.size() )
	{
		cerr << "Sorry don't quite know how this happened!" << endl; return;
	}

	for( int i=0; i < num_repeats; ++i )
	{
		vector<int> local_Start_entries;
		vector<int>::iterator start_i=StartingEntries.begin();
		vector<int>::iterator entries_i=events_to_step_over.begin();
		for( ; ( start_i!=StartingEntries.end() && entries_i != events_to_step_over.end() );  )
		{
			//	Written such that this always addresses the correct element
			local_Start_entries.push_back( (*start_i) + (i+1)* (*entries_i) );

			cout << (*start_i) << "\t" <<(*entries_i) << "\t" << local_Start_entries.back() << endl;

			//	We either have 1 step/start for the whole xml or we have a unique step/start for each PDFWithData from the XML
			if( events_to_step_over.size() != 1 )	{	++entries_i;	}
			if( StartingEntries.size() != 1 )	{	++start_i;	}
		}

		if( StartingEntries.size() == 1 )
		{
			for( unsigned int i=0; i < input_xml->GetAllDataSetSizes().size(); ++i )
			{
				local_Start_entries.push_back( local_Start_entries[0] );
			}
		}

		vector<PDFWithData*> local_PDF_w_Data = input_xml->GetPDFsAndData( local_Start_entries );

		ALL_Fit_Results->StartStopwatch();

		FitResult * newResult = FitAssembler::DoSafeFit( input_xml->GetMinimiserConfiguration(), input_xml->GetFitFunctionConfiguration(), input_xml->GetFitParameters( CommandLineParams ), local_PDF_w_Data, input_xml->GetConstraints() );

		ALL_Fit_Results->AddFitResult( newResult );

	}

}


ToyStudyResult* MCStudy::GetStudyResult()
{
	return ALL_Fit_Results;
}


