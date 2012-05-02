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
//	System Headers
#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>
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
	return;
};

int main( int argc, char* argv[] )
{
	title();
	usage();

	//	If we were provided no input files, exit 'cleanly'
	if( argc < 2 ) exit(-42);

	//	Setup the theme while we remember
	EdStyle* thisStyle = new EdStyle();
	thisStyle->SetStyle();

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

	vector<TTree*> input_trees = ROOT_File_Processing::GetMultipleTrees( input_filenames, RapidFitOutputTupleName );

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
					RapidLL::PlotRapidLL( (*study_i)->control_parameters[0], (*study_i)->tree_to_plot, rand_gen, other_params );
				}
				else if( (*study_i)->control_parameters.size() == 2 )
				{
					//Rapid2DLL::PlotRapid2DLL( (*study_i)->control_parameters[0], (*study_i)->control_parameters[1], (*study_i)->tree_to_plot, rand_gen, other_params );
					Rapid2DLL::PlotRapidFit2DLL( (*study_i)->control_parameters[0], (*study_i)->control_parameters[1], (*study_i)->tree_to_plot, rand_gen, other_params );
				}
			}
		}
	}

	//	Merge multiple outputs

	//	To be written, but this will likely call member functions in each Analysis space to add multiple plots which are returned from each file analysed.

	return 0;
}

