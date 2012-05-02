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
//	System Headers
#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace::std;

int main( int argc, char* argv[] )
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
	cout << "Usage:" << endl;
	cout << "\t" << argv[0] << "\tSomeFile.root\t" << "phys_param1\tphys_param2\toutput_directory" << endl;
	cout << endl << "\teg:\t" << argv[0] << "\tLLcontourData.root\tdeltaGamma\tPhi_s\toutputdir"<<endl;
	cout << endl;
	cout << "6th (Optional) Option:"<<endl;
	cout << "\t\t\tCV\t\t(Plot the CV of 'nuisence' Physics Parameters)"<<endl;
	cout << "\t\t\tRelCV\t\t(Plot the variation in the CV of the 'nuisence' Physics Parameters)"<<endl;
	cout << "\t\t\tExtraCV\t\t(Plot all of the CV, errors and pulls for each 'nuisence' Physics Parameter)"<<endl;
	cout << "\t\t\tNOFC\t\t(Don't process FC data in a file if any is there)"<<endl;
	cout << endl;
	if( (argc != 5) && (argc != 6) ) exit(-1);
	//	Use UUID based seed from ROOT, just used for unique identification of ROOT objects
	TRandom3* rand_gen = new TRandom3(0);

	//	Setup the Canvas and such
	EdStyle* RapidFit_Style = new EdStyle();
	RapidFit_Style->SetStyle();

	//	By Design
	TString outputdir = argv[4];
	TString param1string = argv[2];
	TString param2string = argv[3];

	//	Make OutputDir
	gSystem->mkdir( outputdir );

	//	Open the File
	TFile * input = TFile::Open( argv[1] );
	TNtuple* allresults;
	input->GetObject("RapidFitResult", allresults);
	if(!allresults){
		input->GetObject("RapidFitResult/RapidFitResult",allresults);
		if(!allresults){
			cout << "Couldn't find ntuple RapidFitResult in TFile" << endl;
			exit(1);
		}
	}

	//	Switches for different plotting
	bool CV_Drift = false;
	bool Want_Physics_Params=false;
	bool Want_Extra_Physics_Param_Info=false;
	bool want_FC = true;
	bool high_res = true;	//	False method may only call plot fully once... still in development

	if( argc == 6 )
	{
		if( string(argv[5]) == "CV" )	   Want_Physics_Params = true;
		if( string(argv[5]) == "RelCV" )   { Want_Physics_Params = true; CV_Drift = true; }
		if( string(argv[5]) == "ExtraCV" ) { Want_Physics_Params = true; Want_Extra_Physics_Param_Info = true; }
		if( string(argv[5]) == "NOFC" )    want_FC = false;
	}

	//	Strings that are universally defined
	TString Double_Tolerance="0.000001";
	TString Fit_Status = "Fit_Status";
	TString NLL = "NLL";
	TString notgen = "-9999.";
	TString error_suffix = "_error";
	TString value_suffix = "_value";
	TString pull_suffix = "_pull";
	TString gen_suffix = "_gen";
	TString Copy_Option = "fast";


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
	vector<TString> all_parameters = get_branch_names( allresults );
	//	Get a list of all branches in allresults with '_value' in their name
	vector<TString> all_parameter_values = filter_names( all_parameters, value_suffix );
	vector<TString> all_parameter_errors = filter_names( all_parameters, error_suffix );
	vector<TString> all_parameter_pulls = filter_names( all_parameters, pull_suffix );


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
		cout << all_parameter_values[i] << ":\t" << best_fit_temp_values[i] << endl;
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
	int cont_num = 3;

	double* pllconts = new double[unsigned(cont_num)];
	pllconts[0] = 1.15; pllconts[1] = 2.36;
	pllconts[2] = 3.0;//  pllconts[3] = 4.61;

	double* fcconts = new double[unsigned(cont_num)];
	fcconts[0] = 0.68; fcconts[1] = 0.9;
	fcconts[2] = 0.95;// fcconts[3] = 0.99;

	double* confs = new double[unsigned(cont_num)];
	confs[0] = 68.0; confs[1] = 90.0;
	confs[2] = 95.0;// confs[3] = 99.0;

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

