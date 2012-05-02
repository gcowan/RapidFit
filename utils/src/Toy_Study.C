
#include "TString.h"
#include "TTree.h"
#include "TH1.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TPad.h"
#include "TROOT.h"
#include "TPaveStats.h"
#include "TList.h"

#include "EdStyle.h"

#include "Histo_Processing.h"
#include "TTree_Processing.h"
#include "StringOperations.h"
#include "RapidFit_Output_File.h"
#include "Toy_Study.h"
#include "Mathematics.h"

#include <vector>
#include <string>

using namespace::std;

int ToyStudyAnalysis::Toy_Study( TTree* input_tree, TRandom3* rand_gen, vector<string> OtherOptions )
{
	(void) OtherOptions;

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(111);

	TString Output_Dir;

	Output_Dir = "Output" ;

	gSystem->mkdir( Output_Dir );

	//      Get a list of all branches in allresults
	vector<TString> all_parameters = TTree_Processing::get_branch_names( input_tree );
	//      Get a list of all branches in allresults with '_value' in their name
	vector<TString> all_parameter_values = StringOperations::filter_names( all_parameters, string(value_suffix.Data()) );
	vector<TString> all_parameter_errors = StringOperations::filter_names( all_parameters, string(error_suffix.Data()) );
	vector<TString> all_parameter_pulls = StringOperations::filter_names( all_parameters, string(pull_suffix.Data()) );

	vector<TString> to_be_removed;
	for( vector<TString>::iterator j=all_parameter_values.begin(); j!= all_parameter_values.end(); ++j )
	{
		TString temp = StringOperations::RemoveSuffix( *j, value_suffix );
		temp+=error_suffix;
		//	If I have a value, but NOT an error then it was a fixed parameter
		if( StringOperations::VectorContains( all_parameter_errors, temp ) == -1 )
		{
			to_be_removed.push_back( *j );
		}
	}

	for( vector<TString>::iterator j= to_be_removed.begin(); j!= to_be_removed.end(); ++j )
	{
		int position = StringOperations::VectorContains( all_parameter_values, *j );
		if( position != -1 ) all_parameter_values.erase( all_parameter_values.begin() + position );
	}

	vector<vector<TString> > all; all.push_back( all_parameter_values ); all.push_back( all_parameter_errors ); all.push_back( all_parameter_pulls );

	vector<TString> all_parameter_plots = concatonnate( all );

	cout << all_parameter_plots.size() << " Plots to Draw" << endl;

	for( unsigned int j=0; j< all_parameter_plots.size(); ++j )
	{
		input_tree->Draw( all_parameter_plots[j], "(Fit_Status==3)", "goff" );

		TH1* input_histo = (TH1*) input_tree->GetHistogram();
		TString Histo_Name( all_parameter_plots[j] );
		Histo_Name+=rand_gen->Rndm();
		input_histo->SetName( Histo_Name );

		Histogram_Processing::OptimallyRebin( input_histo );

		TString fit_type = Histogram_Processing::Best_Fit_Function( input_histo );

		Histogram_Processing::Silent_Fit( input_histo, fit_type );

		TString Canvas_Name("Canvas");	Canvas_Name+=rand_gen->Rndm();
		TCanvas* c1 = new TCanvas( Canvas_Name, Canvas_Name, 1680, 1050 );
		Histogram_Processing::Silent_Draw( c1, input_histo );
		TAxis* x_axis = input_histo->GetXaxis();
		x_axis->SetTitle( EdStyle::GetParamRootName( all_parameter_plots[j] ) );
		input_histo->SetTitle( "Toy Study " + EdStyle::GetParamRootName( all_parameter_plots[j] ) + " distribution" );

		c1->Update();

		Histogram_Processing::Silent_Print( c1, Output_Dir+"/"+all_parameter_plots[j]+".png" );

		TPaveStats* thisStats = (TPaveStats*)input_histo->GetListOfFunctions()->FindObject("stats");
		thisStats->SetFillStyle(3023);
		thisStats->SetTextColor(1);

		thisStats->Draw();

		c1->Update();

		Histogram_Processing::Silent_Print( c1, Output_Dir+"/"+all_parameter_plots[j]+"_c_thru.png" );
	}

	return 0;
}

