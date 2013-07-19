//	This is the 3rd complete re-write of the plotting code used in the analysis of 2DLL and FC plots from the Edinburgh RapidFit Fitter
//	The reasoning behind this re-write is complex but the complexity is driven by shortcomings in the root framework,
//	Whilst the speed and actual plots are a credit to the things that root does well

//	ROOT Headers
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TString.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph2D.h"
#include "TPolyMarker3D.h"
#include "TList.h"
#include "TRandom3.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TString.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TMultiGraph.h"
#include "TMatrixDSym.h"
//	RapidFit Headers
#include "EdStyle.h"
//	RapidFit Utils Headers
#include "Toy_Study.h"
#include "RapidLL.h"
#include "Rapid2DLL.h"
#include "DoFCAnalysis.h"
#include "StringOperations.h"
#include "ROOT_File_Processing.h"
#include "TTree_Processing.h"
#include "Histo_Processing.h"
#include "Component_Projections.h"
#include "RapidFit_Output_File.h"
#include "CorrMatrix.h"
//	System Headers
#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace::std;

struct study_to_plot
{
	vector<string> control_parameters;
	TTree* tree_to_plot;
	bool found_toys;
	TObject* produced_plot;
};

void title()
{
	cout << endl;
	cout << "\t8888888b.                    d8b      888      8888888b.  888          888    " << endl;
	cout << "\t888   Y88b                   Y8P      888      888   Y88b 888          888    " << endl;
	cout << "\t888    888                            888      888    888 888          888    " << endl;
	cout << "\t888   d88P  8888b.  88888b.  888  .d88888      888   d88P 888  .d88b.  888888 " << endl;
	cout << "\t8888888P\"      \"88b 888 \"88b 888 d88\" 888      8888888P\"  888 d88\"\"88b 888    " << endl;
	cout << "\t888 T88b   .d888888 888  888 888 888  888      888        888 888  888 888    " << endl;
	cout << "\t888  T88b  888  888 888 d88P 888 Y88b 888      888        888 Y88..88P Y88b.  " << endl;
	cout << "\t888   T88b \"Y888888 88888P\"  888  \"Y88888      888        888  \"Y88P\"   \"Y888 " << endl;
	cout << "\t                    888                                                       " << endl;
	cout << "\t                    888                                                       " << endl;
	cout << "\t                    888                                                       " << endl;
	cout << endl << endl;
	return;
};

void usage()
{
	cout << "Usage:" << endl;
	cout << "\t" << "RapidPlot" << "\tSomeFile.root" << "\tSomeFile2.root\t..." << endl;
	cout << endl << "\teg:\t" << "RapidPlot" << "\tFile1.root\tFile2.root"<<endl;
	cout << endl;
	cout << "Any options passed starting with the '-' option are treated as runtime arguments not filenames.... sorry I don't like files starting with '-' but can't you think of a nicer filename?" << endl;
	cout << endl;
	cout << "For more optional commands and advanced usage try running with --help" << endl;
	cout << endl;
	return;
};

void helpFunction()
{
	cout << endl;
	cout << "Welcome to RapidPlot :)" << endl << endl;
	cout << "RapidPlot has been designed to perform post-processing of data from the RapidFit fitter." << endl;
	cout << endl << "RapidPlot currently knows how to handle:" << endl;
	cout << endl;
	cout << "RapidPlot now also knows how to summarise the content of a RapidFit Fit output, latest versions only" << endl;
	cout << "In order to produce one of these summaries you need to run with: --Summarise to give a summary of the contents of a file" << endl;
	cout << endl;
	cout << "Toy Studies" << endl;
	cout << "1DLL scans" << endl;
	cout << "2DLL scans" << endl;
	cout << "FC scans" << endl;
	cout << "Plotting of Correlation Matricies" << endl;
	cout << endl;
	cout << "Each of these has additional command line options which alter there behaviour" << endl;
	cout << "The additional Commands are:" << endl;
	ToyStudyAnalysis::Help();
	RapidLL::Help();
	Rapid2DLL::Help();
	//FeldmanCousinsAnalysis::Help();		Not yet Implemented
	CorrMatrix::Help();
	cout << endl;
	cout << "Did you Know that the full xml and runtime arguments that you ran your fit with are stored in the global .root output?!" << endl;
	cout << endl;
	cout << "To regenerate a copy of the XML you ran with run with:" << endl << endl;
	cout << "RapidPlot Global_Fit_Output.root --RestoreXML" << endl;
	cout << endl;
	cout << "or you can chose the output XML name with:" << endl;
	cout << "RapidPlot Global_Fit_Output.root --RestoreXML --SaveAs MyOutput.xml" << endl;
	cout << endl;
	cout << "To how this fit was executed run with:" << endl << endl;
	cout << "RapidPlot Global_Fit_Output.root --Summarise" << endl;
	cout << endl;
}

void Summarise( vector<string> input_filenames, vector<string> other_params )
{
	(void) other_params;

	for( unsigned int i=0; i< input_filenames.size(); ++i )
	{
		TTree* runtimeArgs = ROOT_File_Processing::GetTree( input_filenames[i], "RuntimeArgs" );

		vector<string>* thisArgs = new vector<string>();

		runtimeArgs->SetBranchAddress( "RuntimeArgs", &thisArgs );

		runtimeArgs->GetEntry( 0 );

		cout << "File:\t" << input_filenames[i] << "\tgenerated with the following args:" << endl;
		cout << endl;
		cout << "\t";
		for( unsigned int j=0; j< thisArgs->size(); ++j )
		{
			cout << (*thisArgs)[j] << " ";
		}
		cout << endl << endl;
	}
}

void RestoreXML( vector<string> input_filenames, vector<string> other_params )
{
	string SaveAs="--SaveAs";
	string SaveAs2="-SaveAs";

	TString XMLName = "XMLFile_";

	for( unsigned int i=0; i< other_params.size(); ++i )
	{
		string thisString = other_params[i];
		size_t found = thisString.find(SaveAs);
		size_t found2 = thisString.find(SaveAs2);

		if( found2 != string::npos )
		{
			found = found2;
		}

		if( found != string::npos )
		{
			string newName = StringOperations::SplitString( thisString, ':' )[1];
			string xmlext=".xml";
			size_t foundext = newName.find(xmlext);
			if( foundext != string::npos ) newName = StringOperations::SplitString( thisString, '.' )[0];

			XMLName = newName.c_str();
			XMLName.Append("_");
		}
	}

	for( unsigned int i=0; i< input_filenames.size(); ++i )
	{
		TString thisXMLName = XMLName;

		thisXMLName+=i; thisXMLName.Append("_");

		thisXMLName.Append( StringOperations::TimeString() );

		thisXMLName.Append( ".xml" );

		TTree* runtimeXML = ROOT_File_Processing::GetTree( input_filenames[i], "FittingXML" );

		vector<string>* thisXML = new vector<string>();

		runtimeXML->SetBranchAddress( "FittingXML", &thisXML );

		runtimeXML->GetEntry(0);

		cout << "Saving output XML to file:\t" << thisXMLName << endl;
		cout << endl;

		stringstream full_xml;

		for( unsigned int j=0; j< thisXML->size(); ++j )
		{
			full_xml << (*thisXML)[j] << endl;
		}

		ofstream output_xmlFile;
		output_xmlFile.open( thisXMLName.Data() );

		output_xmlFile << full_xml.str();

		output_xmlFile.close();
	}
}

int main( int argc, char* argv[] )
{
	//	Mute ROOT
	gErrorIgnoreLevel = kFatal;
	title();
	usage();

	//	If we were provided no input files, exit 'cleanly'
	if( argc < 2 ) exit(-42);

	//	Use UUID based seed from ROOT, just used primarily for unique identification of ROOT objects
	TRandom3* rand_gen = new TRandom3();

	//	Setup the Canvas and such
	EdStyle* RapidFit_Style = new EdStyle();
	RapidFit_Style->SetStyle();

	vector<string> input_filenames;
	vector<string> other_params;
	for( int i=1; i< argc; ++i )
	{
		if( argv[i][0] == '-' )	other_params.push_back( argv[i] );
		else	input_filenames.push_back( argv[i] );
	}

	string helpOption="--help";
	string helpOption2="-help";
	if( StringOperations::VectorContains( &other_params, &helpOption ) != -1 || StringOperations::VectorContains( &other_params, &helpOption2 ) != -1 )
	{
		helpFunction();
		exit(0);
	}

	string SummeriseText="--Summarise";
	string SummeriseText2="-Summarise";
	if( StringOperations::VectorContains( &other_params, &SummeriseText ) != -1 || StringOperations::VectorContains( &other_params, &SummeriseText2 ) != -1 )
	{
		Summarise( input_filenames, other_params );
		exit(0);
	}

	string restoreXML="--RestoreXML";
	string restoreXML2="-RestoreXML";
	if( StringOperations::VectorContains( &other_params, &restoreXML ) != -1 || StringOperations::VectorContains( &other_params, &restoreXML2 ) != -1 )
	{
		RestoreXML( input_filenames, other_params );
		exit(0);
	}

	vector<TTree*> input_trees = ROOT_File_Processing::GetMultipleTrees( input_filenames, RapidFitOutputTupleName );

	string CorrMatrixName="corr_matrix";
	vector<TTree*> corr_trees = ROOT_File_Processing::GetMultipleTrees( input_filenames, CorrMatrixName );

	if( !corr_trees.empty() )
	{
		vector<string> argv_str;
		for( unsigned int i=0; i< (unsigned) argc; ++i ) argv_str.push_back( argv[i] );
		string CorrMatrix = "--CorrMatrix";
		if( StringOperations::VectorContains( &argv_str , &CorrMatrix ) != -1 )
		{
			CorrMatrix::Analyse( corr_trees, other_params );
		}
	}

	int good_files=0;
	for( unsigned int i=0; i< input_trees.size(); ++i )
	{
		if( input_trees[i] != NULL ) ++good_files;
	}

	vector<pair<string,string> > Directories_in_file;
	TFile* ProjFile=NULL;
	if( good_files==0 )
	{
		ProjFile = ROOT_File_Processing::OpenFile( argv[1] );
		(void) ProjFile;
		TString top_dir = gDirectory->GetPath();

		ROOT_File_Processing::get_TDirectory_list( top_dir, &Directories_in_file );
		if( Directories_in_file.size() != 0 ) ++good_files;
	}

	if( good_files == 0 )
	{
		cerr << "\n\tNo usable files found as arguments, exiting...\n" << endl;
		exit(-99);
	}


	//	Construct study_to_plot objects for each input tree found

	vector<struct study_to_plot*> Studies_to_Plot;

	for( unsigned int file_i=0; file_i < input_trees.size(); ++file_i )
	{
		vector<string> controlled_parameters_scan = RapidFit_Output_File::get_control_parameters( input_trees[file_i] );

		vector<string> controlled_parameters;
		for( vector<string>::iterator index_i = controlled_parameters_scan.begin(); index_i != controlled_parameters_scan.end(); ++index_i )
		{
			controlled_parameters.push_back( StringOperations::RemoveSuffix( *index_i, "_scan" ).Data() );
		}
		struct study_to_plot* this_study = new study_to_plot;

		if( controlled_parameters.empty() )
		{
			cout << "No controlled parameters." << endl;
			this_study->tree_to_plot = input_trees[file_i];
			Studies_to_Plot.push_back( this_study );
			continue;
		}


		bool has_toys = RapidFit_Output_File::HasToys( input_trees[file_i], controlled_parameters, rand_gen );

		this_study->control_parameters = controlled_parameters;
		this_study->tree_to_plot = input_trees[file_i];
		this_study->found_toys = has_toys;

		Studies_to_Plot.push_back( this_study );
	}

	//	Perform analysis

	//	Construct Objects which must be used to store the result from mutliple scans for merging the output
	TMultiGraph* GraphsToOverlay = new TMultiGraph( "top_level_overlay", "top_level_overlay" );

	TString Here = TString( gSystem->pwd() );
	for( vector<struct study_to_plot*>::iterator study_i =  Studies_to_Plot.begin(); study_i != Studies_to_Plot.end(); ++study_i )
	{
		if( (*study_i)->control_parameters.empty() == true )
		{
			cout << "\n\tFound 0 Scanned Parameters but found fit data, performing Toy Study Analysis.\n" << endl;
			ToyStudyAnalysis::Toy_Study( (*study_i)->tree_to_plot, rand_gen, other_params );
		}
		else
		{
			if( (*study_i)->found_toys == true )
			{
				cout << "\n\tFound Toys in a file used to perform a scan, performing a FelmanCousins Analysis.\n" << endl;
				FeldmanCousinsAnalysis::DoFCAnalysis( (*study_i)->tree_to_plot, (*study_i)->control_parameters, rand_gen, other_params );
			}
			else
			{
				if( (*study_i)->control_parameters.size() == 1 )
				{
					cout << "\n\tPlotting the Results in LLscan format.\n" << endl;
					TGraph* thisGraph = RapidLL::PlotRapidLL( (*study_i)->control_parameters[0], (*study_i)->tree_to_plot, rand_gen, other_params );
					GraphsToOverlay->Add( thisGraph );
				}
				else if( (*study_i)->control_parameters.size() == 2 )
				{
					//Rapid2DLL::PlotRapid2DLL( (*study_i)->control_parameters[0], (*study_i)->control_parameters[1], (*study_i)->tree_to_plot, rand_gen, other_params );
					Rapid2DLL::PlotRapidFit2DLL( (*study_i)->control_parameters[0], (*study_i)->control_parameters[1], (*study_i)->tree_to_plot, rand_gen, other_params );
				}
			}
		}
		gSystem->cd( Here );
	}

	//	Merge multiple outputs
	if( GraphsToOverlay->GetListOfGraphs() != NULL )
	{
		if( GraphsToOverlay->GetListOfGraphs()->Capacity() > 1 )
		{
			RapidLL::OverlayMutliplePlots( GraphsToOverlay );
		}
	}
	//	2D case possible to be written but as of this time 2012/10 There is no useful call to write this so it has not been addressed

	cout << endl;
	cout << "Any segfaults beyond here are likely problems with ROOT..." << endl;
	cout << "Goodbye from RapidPlot :D" << endl;
	cout << endl;

	return 0;
}

