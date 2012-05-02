
#include "TString.h"
#include "TTree.h"
#include "TH1.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "EdStyle.h"

#include "Histo_Processing.h"
#include "NTuple_Processing.h"
#include "TString_Processing.h"

#include <vector>

int Toy_Study( vector<TTree*> input_trees, TRandom3* rand_gen, vector<TString> Output_Dirs )
{
        //      Strings that are universally defined
        TString Double_Tolerance="0.000001";
        TString Copy_Option = "fast";

        gStyle->SetOptStat(0);
        gStyle->SetOptFit(111);

	TString Output_Dir;

	for( unsigned int i=0; i< input_trees.size(); ++i )
	{
		Output_Dir = Output_Dirs[i];

        	if( is_empty( Output_Dir ) )
        	{
        	        Output_Dir = "Output" ;
        	}

		gSystem->mkdir( Output_Dir );

		//      Get a list of all branches in allresults
		vector<TString> all_parameters = get_branch_names( input_trees[0] );
		//      Get a list of all branches in allresults with '_value' in their name
		vector<TString> all_parameter_values = filter_names( all_parameters, string(value_suffix.Data()) );
		vector<TString> all_parameter_errors = filter_names( all_parameters, string(error_suffix.Data()) );
		vector<TString> all_parameter_pulls = filter_names( all_parameters, string(pull_suffix.Data()) );

		vector<TString> to_be_removed;
		for( vector<TString>::iterator j=all_parameter_values.begin(); j!= all_parameter_values.end(); ++j )
		{
			TString temp = RemoveSuffix( *j, value_suffix );
			temp+=error_suffix;
			//	If I have a value, but NOT an error then it was a fixed parameter
			if( VectorContains( all_parameter_errors, temp ) == -1 )
			{
				to_be_removed.push_back( *j );
			}
		}

		for( vector<TString>::iterator j= to_be_removed.begin(); j!= to_be_removed.end(); ++j )
		{
			int position = VectorContains( all_parameter_values, *j );
			if( position != -1 ) all_parameter_values.erase( all_parameter_values.begin() + position );
		}

		vector<vector<TString> > all; all.push_back( all_parameter_values ); all.push_back( all_parameter_errors ); all.push_back( all_parameter_pulls );

		vector<TString> all_parameter_plots = concatonnate( all );

		for( unsigned int j=0; j< all_parameter_plots.size(); ++j )
		{
			input_trees[i]->Draw( all_parameter_plots[j], "Fit_Status==3", "goff" );

			TH1* input_histo = (TH1*) input_trees[i]->GetHistogram();
			TString Histo_Name( all_parameter_plots[j] );
			Histo_Name+=rand_gen->Rndm();
			input_histo->SetName( Histo_Name );

			OptimallyRebin( input_histo );

			TString fit_type = Best_Fit_Function( input_histo );

			Silent_Fit( input_histo, fit_type );

			TString Canvas_Name("Canvas");	Canvas_Name+=rand_gen->Rndm();
			TCanvas* c1 = new TCanvas( Canvas_Name, Canvas_Name, 1680, 1050 );
			Silent_Draw( c1, input_histo );
			TAxis* x_axis = input_histo->GetXaxis();
			x_axis->SetTitle( EdStyle::GetParamRootName( all_parameter_plots[j] ) );
			input_histo->SetTitle( "Toy Study " + EdStyle::GetParamRootName( all_parameter_plots[j] ) + " distribution" );
			Silent_Print( c1, Output_Dir+"/"+all_parameter_plots[j]+".png" );
		}
	}

	return 0;
}

