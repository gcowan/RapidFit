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
//	RapidFit Headers
#include "EdStyle.h"
//	RapidFit Utils Headers
#include "TString_Processing.h"
#include "NTuple_Processing.h"
#include "Histo_Processing.h"
#include "Rapid2DLL.h"
#include "DoFCAnalysis.h"
//	System Headers
#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace::std;

int Rapid2DLL( vector<TString> controlled_parameter1, vector<TString> controlled_parameter2, vector<TTree*> input_trees, TRandom3* rand_gen, vector<string> other_params )
{
	vector<string> controls; controls.push_back( controlled_parameter1[0].Data() ); controls.push_back( controlled_parameter2[0].Data() );

	bool hasToys = HasToys( input_trees[0], controls, rand_gen );

	//if( hasToys == true )
	//{
	//	cout << "I HAVE TOYS" << endl;
	//	DoFCAnalysis( input_trees[0], controls, rand_gen );
	//	exit(-9);
	//}
	//else
	//{
	//	cout << "I DO NOT HAVE TOYS" << endl;
	//}

	//      Strings that are universally defined
	TString Double_Tolerance="0.000001";
	TString Fit_Status = "Fit_Status";
	TString NLL = "NLL";
	TString notgen = "-9999.";
	TString error_suffix = "_error";
	TString value_suffix = "_value";
	TString pull_suffix = "_pull";
	TString gen_suffix = "_gen";
	TString ScanStatus_suffix = "_scan";
	TString Copy_Option = "fast";

	if( input_trees.size() > 1 )
	{
		cerr << "\nBad number of files to process, may include multiple files at some time later, but not yet, exiting.'\n" << endl;
		exit(-45);
	}

	//	Setup the Canvas and such
	EdStyle* RapidFit_Style = new EdStyle();
	RapidFit_Style->SetStyle();

	//	By Design
	TString param1string = controlled_parameter1[0];
	TString param2string = controlled_parameter2[0];

	TString outputdir( param1string+"_"+param2string );

	//	Make OutputDir
	gSystem->mkdir( outputdir );

	//	Open the File
	TTree* allresults = input_trees[0];

	//	Switches for different plotting
	bool CV_Drift = false;
	bool Want_Physics_Params=false;
	bool Want_Extra_Physics_Param_Info=false;
	bool want_FC = true;
	bool high_res = true;	//	False method may only call plot fully once... still in development

	for( unsigned int i=0; i< other_params.size(); ++i )
	{
		if( other_params[i] == "-CV" )	   Want_Physics_Params = true;
		if( other_params[i] == "-RelCV" )   { Want_Physics_Params = true; CV_Drift = true; }
		if( other_params[i] == "-ExtraCV" ) { Want_Physics_Params = true; Want_Extra_Physics_Param_Info = true; }
		if( other_params[i] == "-NOFC" )    want_FC = false;
	}


	//	Strings that are relevent to this plot
	//	This is currently given by the user, but in time we can write an algorithm to determine this
	TString param1_gen = param1string + gen_suffix;
	TString param1_val = param1string + value_suffix;
	TString param1_err = param1string + error_suffix;
	TString param2_gen = param2string + gen_suffix;
	TString param2_val = param2string + value_suffix;
	TString param2_err = param2string + error_suffix;


	//	Get a list of ALL parameters within the file, this is given for free in the algorithms I wrote above

	//	Get a list of all branches in allresults
	vector<TString> all_parameters = get_free_parameters( allresults );
	//	Get a list of all branches in allresults with '_value' in their name
	vector<TString> all_parameter_values = postpend( all_parameters, value_suffix );
	vector<TString> all_parameter_errors = postpend( all_parameters, error_suffix );
	vector<TString> all_parameter_pulls = postpend( all_parameters, pull_suffix );


	//	Now to read in the Global Minima

	//	Store the Global fit corrdinate and value
	Float_t nll=0, Global_Best_NLL=0, x_point=0, y_point=0;

	//	Store the value of the 'nuisence' parameters as well
	Float_t* best_fit_values = new Float_t[all_parameter_values.size()];
	//	Construct an array to hold the best values for the rest of the parameters
	Float_t* best_fit_temp_values = new Float_t[all_parameter_values.size()];
	allresults->SetBranchAddress(NLL,&nll);

	//	Tell ROOT where to store the values
	for( unsigned short int i=0; i < all_parameter_values.size(); ++i )
	{
		allresults->SetBranchAddress(all_parameter_values[i],&best_fit_temp_values[i]);
	}


	//	Actually store the values in resident memory of the program
	cout << endl << "BEST FIT RESULTS AS DETERMINED FROM THE STANDARD FIT:" << endl;
	allresults->GetEntry(0);
	Global_Best_NLL = nll;
	for( unsigned short int i=0; i < all_parameter_values.size(); ++i )
	{
		best_fit_values[i] = best_fit_temp_values[i];
		if( string(all_parameter_values[i].Data()).compare(param1_val.Data()) == 0 ) x_point = best_fit_temp_values[i];
		if( string(all_parameter_values[i].Data()).compare(param2_val.Data()) == 0 ) y_point = best_fit_temp_values[i];
	}
	cout << endl;



	//	General Cuts to be applied for various plots

	//	Fit_Status == 3
	TString Fit_Cut = "(abs(" + Fit_Status + "-3.0)<"+Double_Tolerance+")";


	//	param1_gen == notgen	&&	param1_err == 0
	TString Param_1_Cut = "(abs(" + param1_gen + "-" + notgen +")<" + Double_Tolerance + ")&&(abs("+ param1_err + "-0.0)<"+ Double_Tolerance + ")";
	//	param2_gen == notgen	&&	param2_err == 0
	TString Param_2_Cut = "(abs(" + param2_gen + "-" + notgen + ")<"+ Double_Tolerance + ")&&(abs(" +param2_err +"-0.0)<" + Double_Tolerance + ")";

	//	Combine the individual Cuts
	TString Fit_Cut_String = Param_1_Cut + "&&" + Param_2_Cut + "&&" + Fit_Cut;

	//TString Fit_Cut_String = Fit_Cut;


	//	Toys have a defined generation Value, Fits to Data DO NOT

	//	param1_gen != notgen
	TString p1isatoy = "(abs(" + param1_gen + "-" + notgen + ")>" + Double_Tolerance + ")";
	//	param2_gen != notgen
	TString p2isatoy = "(abs(" + param2_gen + "-" + notgen + ")>" + Double_Tolerance + ")";

	//	Combine the individual Cuts
	TString Toy_Cut_String = p1isatoy + "&&" + p2isatoy;




	//	Check for Toys in the file and wether I should run the FC code

	allresults->Draw( NLL, Toy_Cut_String, "goff" );
	bool Has_Toys = allresults->GetSelectedRows() > 0;
	Has_Toys = false;

	cout << endl << "NUMBER OF TOYS IN FILE:\t" << allresults->GetSelectedRows() << endl;

	if( int(allresults->GetEntries() - allresults->GetSelectedRows()) == 0 )
	{
		cerr << "SERIOUS ERROR:\tSOMETHING HAS REALLY GOTTEN SCREWED UP!" << endl;
		exit(-3498);
	}


	//	Fit values for the global fit are now stored in:
	//
	//	Global_Best_NLL,	best_fit_values,	x_point,	y_point

	//	Tell the user
	cout << "GLOBAL DATA BEST FIT NLL:\t" << setprecision(10) << Global_Best_NLL << "\tAT:\tX:" << setprecision(10) << x_point << "\tY:\t" <<setprecision(10)<< y_point << endl;

	//	Check wether the minima as defined from the Global fit is the true minima within the phase-space

	//Check_Minima( allresults, Fit_Cut_String, &Global_Best_NLL, NLL, Double_Tolerance, param1_val, param2_val );


	TString NLL_Min;		//	Of course ROOT doesn't have USEFUL constructors!
	NLL_Min+=Global_Best_NLL;

	//	Plot the distribution of successfully fitted grid points for the PLL scan
	//	NB:	For FC this will likely saturate due to multiple layers of fits
	cout << endl << "FOUND UNIQUE GRID POINTS, PLOTTING" << endl;

	LL2D_Grid( allresults, Fit_Cut_String, param1string, param2string, rand_gen, "LL", outputdir );

	//	Construct a plot string for the NLL plot and plot it

	//	PLL part of Draw_String
	TString PLL = "(" + NLL + "-" + NLL_Min + ")";
	//	Gridding part of Draw_String
	TString Param1_Param2 = ":" + param1_val + ":" + param2_val;


	//	Draw String for PLL
	TString NLL_DrawString = PLL + Param1_Param2;

	cout << endl << "PLOTTING NLL VARIATION" << endl;
	//	PLL plot
	TH2* pllhist=NULL;
	TGraph2D* pllgraph=NULL;

	if( high_res )	{
		pllhist = Plot_From_Cut( allresults, NLL_DrawString, Fit_Cut_String, rand_gen, param1string, param2string );
	} else {
		pllgraph = Plot_From_Cut_lo( allresults, NLL_DrawString, Fit_Cut_String, rand_gen, param1string, param2string );
	}

	//	Array storing the addresses of all of the Physics Plots still in memory
	TH2** All_Physics_Plots = new TH2*[all_parameter_values.size()];
	TH2** All_Physics_Errors = new TH2*[all_parameter_values.size()];
	TH2** All_Physics_Pulls = new TH2*[all_parameter_values.size()];
	//	Making use of switch
	if( Want_Physics_Params )
	{
		//	Physics Plots
		Physics_Plots( all_parameter_values, best_fit_values, allresults, rand_gen, Param1_Param2, CV_Drift, All_Physics_Plots, Fit_Cut_String);
		cout << endl << "Finalising CV Plots" << endl << endl;
		Finalize_Physics_Plots( All_Physics_Plots, all_parameter_values, param1string, param2string, outputdir, CV_Drift );
		if( Want_Extra_Physics_Param_Info )
		{
			//      Physics Plots
			//      
			//      Passing false as Highly unlikely that I will 'ever' want to watch error/pull drift rather than the absolute value

			Float_t* null_pointer=NULL;
			Physics_Plots( all_parameter_errors, null_pointer, allresults, rand_gen, Param1_Param2, false, All_Physics_Errors, Fit_Cut_String);
			Physics_Plots( all_parameter_pulls, null_pointer, allresults, rand_gen, Param1_Param2, false, All_Physics_Pulls, Fit_Cut_String);
			cout << endl << "Finalising Extra CV Plots" << endl << endl;
			Finalize_Physics_Plots( All_Physics_Errors, all_parameter_errors, param1string, param2string, outputdir, false );
			Finalize_Physics_Plots( All_Physics_Pulls, all_parameter_pulls, param1string, param2string, outputdir, false );
		}
	}


	TFile* new_FC_Output = NULL;
	TTree* FC_Output = new TTree( "FC_Output", "FC_Output" );;
	TH2* FC_Plot = NULL;
	//	Making use of switch
	if( Has_Toys && want_FC )
	{
		//	FC Plot
		cout << "FOUND TOYS IN FILE, PLOTTING FC" <<endl;
		TString FC_Tuple_File( outputdir+ "/FC_Tuple.root" );
		new_FC_Output = new TFile( FC_Tuple_File, "RECREATE" );

		FC_Plot = FC_TOYS( allresults, Fit_Cut_String, param1string, param2string, NLL, Fit_Cut, Global_Best_NLL, FC_Output, Double_Tolerance, rand_gen );
		FC_Output->Write();

		//LL2D_Grid( FC_Output, Fit_Cut, param1_val, param2_val, rand_gen, "FC" );
	}

	//	Contours to be used in plotting
	int cont_num = 4;
	//int cont_num = 3;

	double* pllconts = new double[unsigned(cont_num)];
	//4 normal levels
	pllconts[0] = 1.15;
	pllconts[1] = 2.31;
	pllconts[2] = 3.;
	pllconts[3] = 4.61;
	//1, 2, 3 sigma levels
	//pllconts[0] = 1.15;
	//pllconts[1] = 3.09;
	//pllconts[2] = 5.915;

	double* fcconts = new double[unsigned(cont_num)];
	//4 normal levels
	fcconts[0] = 0.68;
	fcconts[1] = 0.9;
	fcconts[2] = 0.95;
	fcconts[3] = 0.99;
	//1, 2, 3 sigma levels
	//fcconts[0] = 0.6827;
	//fcconts[1] = 0.9545;
	//fcconts[2] = 0.9973;

	double* confs = new double[unsigned(cont_num)];
	//4 normal levels
	confs[0] = 68;
	confs[1] = 90;
	confs[2] = 95;
	confs[3] = 99;
	//1, 2, 3, sigma levels
	//confs[0] = 68.27;
	//confs[1] = 95.45;
	//confs[2] = 99.73;

	cout <<endl<< "SAVING GRAPHS" << endl;

	if( high_res )
	{
		Plot_Styled_Contour( pllhist, cont_num, pllconts, confs, outputdir, "Likelihood Profile" );
	} else {
		//Plot_Styled_Contour2( pllgraph, cont_num, pllconts, confs, outputdir, "Likelihood Profile" );
	}

	if( Has_Toys && want_FC )	//	If FC was generated
	{
		//	Thinking of implementing a PLL Plot for the FC as it's easier to plot that (and looks more impressive)
		//TH2* FC_PLL = Invert_PLL( FC_Plot );
		Plot_Styled_Contour( FC_Plot, cont_num, fcconts, confs, outputdir, "FeldmanCousins Profile" );
		//Plot_Styled_Contour( FC_PLL, cont_num, fcconts, confs, outputdir, "FeldmanCousins Profile PLL" );

		//TString FC_Tuple_File( outputdir+ "FC_Tuple.root" );
		//TFile* new_FC_Output = new TFile( FC_Tuple_File, "RECREATE" );
		//FC_Output->Write();

		Plot_Both( pllhist, FC_Plot, cont_num, fcconts, pllconts, confs, outputdir );
		FC_Stats( FC_Output, param1string, param2string, rand_gen, outputdir );
		if( new_FC_Output != NULL )	new_FC_Output->Close();
	}

	return 0;
}

