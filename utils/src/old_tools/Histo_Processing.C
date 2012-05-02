//	ROOT Headers
#include "TFile.h"
#include "TH1.h"
#include "TAxis.h"
#include "TString.h"
#include "TPolyMarker.h"
#include "TPolyMarker3D.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "TList.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TPaletteAxis.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
//	RapidFit Header
#include "EdStyle.h"
//	RapidFit Utils Headers
#include "Histo_Processing.h"
#include "TString_Processing.h"
#include "NTuple_Processing.h"
//	System Headers
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits.h>
#include <float.h>

using namespace::std;

bool double_equals_test(Double_t first, Double_t second)
{
	return fabs(first-second) < 0.0001;
}

//  This has been adapted from the original code in RapidFits Statistics code
//  It is intended to take a histogram and automatically rebin according to this function
//  As such it uses as much inbuilt functionality in root as possible
//Return the ideal number of bins for a histogram of a vector of doubles
////Uses D. Scott's method, published 1979
int OptimumBinNumber( TH1* input_hist, int axis )
{
	double wanted_bins = GetOptimalBins( input_hist, axis );
	double existing_bins = 0.;
	if( axis == 1 ) existing_bins = input_hist->GetNbinsX();
	else if( axis == 2 )  existing_bins = input_hist->GetNbinsY();
	else if( axis == 3 )  existing_bins = input_hist->GetNbinsZ();
	int rebin_factor = int( ( (double) existing_bins ) / (( double)wanted_bins ) );
	if( rebin_factor <= 0 ) rebin_factor = 1;
	input_hist->Rebin( rebin_factor );

	return (int) wanted_bins;
}

//  Return the optimal number of bins that should be used if this data is gaussianly distributed
unsigned int GetOptimalBins( TH1* input_hist, int axis )
{
	double varience = input_hist->GetRMS( axis );
	varience = varience*varience;
	double num_entries = input_hist->GetEntries();
	double width = 3.49 * sqrt( varience ) * pow( num_entries, -(1.0/3.0) );
	double min_range = 0.;
	double max_range = 0.;
	if( axis == 1 )
	{
		min_range = input_hist->GetXaxis()->GetXmin();
		max_range = input_hist->GetXaxis()->GetXmax();
	}
	else if( axis == 2 )
	{
		min_range = input_hist->GetYaxis()->GetXmin();
		max_range = input_hist->GetYaxis()->GetXmax();
	}
	else if( axis == 3 )
	{
		min_range = input_hist->GetZaxis()->GetXmin();
		max_range = input_hist->GetZaxis()->GetXmax();
	}

	double range = max_range - min_range;
	double wanted_bins = ceil( range / width );

	//  Catch SERIOUS rounding errors and caes where the information has already been lost due to underbinning
	if( ( axis == 1 ) && ( wanted_bins > input_hist->GetNbinsX() ) )  return 0;
	else if( ( axis == 2 ) && ( wanted_bins > input_hist->GetNbinsY() ) )  return 0;
	else if( ( axis == 3 ) && ( wanted_bins > input_hist->GetNbinsZ() ) )  return 0;
	else return (unsigned) wanted_bins;

	return 0;
}

//  This function simply rebins a histogram to an optimal number of bins for a gaussian like distribution
void OptimallyRebin( TH1* input_hist, int axis )
{
	int temp = OptimumBinNumber( input_hist, axis );
	(void) temp;	return;
}

//      This will find all of the unique coordinates stored in a 1D vector of doubles
vector<vector<Double_t> > Unique_Coords( vector<Double_t> input )
{
	vector<Double_t> output_temp = input;
	sort( output_temp.begin(), output_temp.end() );
	vector<Double_t>::iterator len;
	len = unique( output_temp.begin(), output_temp.end(), double_equals_test );
	output_temp.resize( len - output_temp.begin() );
	vector<vector<Double_t> > output( 1, output_temp );
	return output;
}
//	This will find all of the unique coordinates stored in a 2D vector of pairs of doubles
vector<vector<Double_t> > Unique_Coords( vector<pair<Double_t,Double_t> > input )
{
	//      Tolerance of the DOUBLE
	double DT = 1E-5;

	//      Get the number of coordinates in the TPolyMarker3D
	int number_of_points = (int)input.size();

	//      Will populate and return to the user
	vector<vector<Double_t> > Returnable_Coord_Data;

	//      Used within a for loop in logic, externally created/destroyed
	bool add_point = true;
	bool temp_decision_1 = true;
	bool temp_decision_2 = true;

	//      Run over all of the points contained within the TPolyMarker3D object
	for( int i=0; i< number_of_points; ++i )
	{
		//      Assume we haven't seen this point yet
		add_point = true;
		//      Check the anzats
		for( unsigned int j=0; j< Returnable_Coord_Data.size(); ++j )
		{
			temp_decision_1 = true;
			temp_decision_2 = true;
			if( fabs( input[i].first - Returnable_Coord_Data[j][0] ) < DT ) temp_decision_1 = false;
			if( fabs( input[i].second - Returnable_Coord_Data[j][1] ) < DT ) temp_decision_2 = false;

			//      If all of the information is the same do NOT add the point to the new vector of data
			if( ( temp_decision_1 == temp_decision_2 ) && ( temp_decision_1 == false ) )
				add_point = false;
		}
		//      If we haven't seen this point yet add it to the array of points
		if( add_point )
		{
			//cout << Coord_Data_pointer[i].first << "\t" << Coord_Data_pointer[i].second << endl;
			vector<Double_t> temp_vector;
			temp_vector.push_back( input[i].first );
			temp_vector.push_back( input[i].second );
			Returnable_Coord_Data.push_back( temp_vector );
		}
	}

	//      Return a 2D vector of unique corrdinates
	return Returnable_Coord_Data;
}

//      This will find all of the unique coordinates stored in a n-D vector of doubles
vector<vector<Double_t> > Unique_Coords( vector<vector<Double_t> > input )
{
	//      Tolerance of the DOUBLE
	double DT = 1E-5;

	//      Get the number of coordinates in the TPolyMarker3D
	int number_of_points = (int)input.size();

	//      Will populate and return to the user
	vector<vector<Double_t> > Returnable_Coord_Data;

	//      Used within a for loop in logic, externally created/destroyed
	bool add_point = true;

	//      Run over all of the points contained within the TPolyMarker3D object
	for( int i=0; i< number_of_points; ++i )
	{
		//      Assume we haven't seen this point yet
		add_point = true;
		//      Check the anzats
		for( vector<vector<Double_t> >::iterator ret_i = Returnable_Coord_Data.begin(); ret_i != Returnable_Coord_Data.end(); ++ret_i )
		{
			vector<bool> decisions;
			vector<Double_t>::iterator coord_i = input[i].begin();
			vector<Double_t>::iterator coord_j = ret_i->begin();
			for( ; coord_i != input[i].end(); ++coord_i, ++coord_j )
			{
				if( fabs( *coord_i - *coord_j ) < DT ) decisions.push_back( false );
				else decisions.push_back( true );
			}
			for( vector<bool>::iterator dec_i = decisions.begin(); dec_i != decisions.end(); ++dec_i )
			{
				add_point = add_point && *dec_i;
			}
		}
		//      If we haven't seen this point yet add it to the array of points
		if( add_point )
		{
			//cout << Coord_Data_pointer[i].first << "\t" << Coord_Data_pointer[i].second << endl;
			vector<Double_t> temp_vector;
			for( vector<Double_t>::iterator wanted_i = input[i].begin(); wanted_i != input[i].end(); ++wanted_i )
			{
				temp_vector.push_back( *wanted_i );
			}
			Returnable_Coord_Data.push_back( temp_vector );
		}
	}

	//      Return a 2D vector of unique corrdinates
	return Returnable_Coord_Data;
}

//      This will find all of the unique coordinates stored in a TPolyMarker3D object
//      NB this was written due to the fact that ROOT thows away the contents of {3/4}D TTree->Draw() objects... (God only knows why)
vector<vector<Double_t> > Unique_Coords( TPolyMarker3D *pm )
{
	//      Tolerance of the DOUBLE
	double DT = 1E-5;

	//      Get the number of coordinates in the TPolyMarker3D
	int number_of_points = pm->GetN();
	//      Get the data contained in the TPolyMarker3D
	Float_t* Coord_Data_pointer = pm->GetP();

	//      Will populate and return to the user
	vector<vector<Double_t> > Returnable_Coord_Data;

	//      Used within a for loop in logic, externally created/destroyed
	bool add_point = true;
	bool temp_decision_1 = true;
	bool temp_decision_2 = true;
	bool temp_decision_3 = true;

	//      Run over all of the points contained within the TPolyMarker3D object
	for( int i=0; i< number_of_points; ++i )
	{
		//      Assume we haven't seen this point yet
		add_point = true;
		//      Check the anzats
		for( unsigned int j=0; j< Returnable_Coord_Data.size(); ++j )
		{
			temp_decision_1 = true;
			temp_decision_2 = true;
			temp_decision_3 = true;
			if( fabs( Coord_Data_pointer[i*3] - Returnable_Coord_Data[j][0] ) < DT ) temp_decision_1 = false;
			if( fabs( Coord_Data_pointer[i*3+1] - Returnable_Coord_Data[j][1] ) < DT ) temp_decision_2 = false;
			if( fabs( Coord_Data_pointer[i*3+2] - Returnable_Coord_Data[j][2] ) < DT ) temp_decision_3 = false;

			//      If all of the information is the same do NOT add the point to the new vector of data
			if( ( ( temp_decision_1 == temp_decision_2 ) && ( temp_decision_2 == temp_decision_3 ) ) && ( temp_decision_1 == false ) )
				add_point = false;
		}
		//      If we haven't seen this point yet add it to the array of points
		if( add_point )
		{
			//cout << Coord_Data_pointer[i*3] << "\t" << Coord_Data_pointer[i*3+1] << "\t" << Coord_Data_pointer[i*3+2] << endl;
			vector<Double_t> temp_vector;
			temp_vector.push_back( Coord_Data_pointer[i*3] );
			temp_vector.push_back( Coord_Data_pointer[i*3+1] );
			temp_vector.push_back( Coord_Data_pointer[i*3+2] );
			Returnable_Coord_Data.push_back( temp_vector );
		}
	}

	//      Return a 2D vector of unique 3D corrdinates of size npoints*3
	return Returnable_Coord_Data;
}

vector<vector<Double_t> > Plotter_Data( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom3* random )
{
	//      Plot the graph using TTree->Draw()
	//      The resulting graph (TH3) contains empty axis and a TPolyMarker3D object
	input_tree->SetEstimate(input_tree->GetEntries());  // Fix the size of the array of doubles to be created (There will never be more than this
	TDirectory* temp_path = gDirectory;
	TString canv_name = "Plot_Data_"; canv_name += random->Rndm()*1000;
	TCanvas* temp_canv = new TCanvas( canv_name, canv_name, 1680, 1050 );
	input_tree->Draw( Draw_String, Cut_String );

	cout << Draw_String << "\t" << Cut_String << endl;

	int numberofinputColumns = input_tree->GetPlayer()->GetDimension();

	if( numberofinputColumns == 3 )
	{
		//      Get the Points that have been plotted (TPolyMarker3D object named "TPolyMarker3D", see ROOTtalk)
		TPolyMarker3D *pm = new TPolyMarker3D( *(TPolyMarker3D*)gPad->FindObject("TPolyMarker3D") );
		double temp = random->Rndm();
		//temp_canv->Close();
		//temp_path->cd();
		TString Name = "TPoly3_";
		Name+=temp;

		if( pm!= NULL )
		{
			pm->SetName(Name);

			//      Get a list of ONLY unique coordinates due to the short comings of the interpolation within TGraph2D
			vector<vector<Double_t> > returnable_data = Unique_Coords( pm );

			return returnable_data;
		}
		cout << "Didn't pick up any points in Cut" << endl;
		vector<vector<Double_t> > dummy_return;
		return dummy_return;
	}
	else if ( numberofinputColumns == 2 )
	{
		Double_t* first = input_tree->GetV1();
		Double_t* second = input_tree->GetV2();

		int rows = (int)input_tree->GetSelectedRows();

		vector<pair<Double_t,Double_t> > to_be_sorted;

		for( int i=0; i< rows; ++i )
		{
			to_be_sorted.push_back( make_pair( first[(unsigned)i], second[(unsigned)i] ) );
		}

		vector<vector<Double_t> > returnable_data = Unique_Coords( to_be_sorted );
		return returnable_data;
	}
	else
	{
		Double_t* first = input_tree->GetV1();

		int rows = (int)input_tree->GetSelectedRows();

		vector<Double_t> to_be_sorted;

		for( int i=0; i< rows; ++i )
		{
			to_be_sorted.push_back( first[(unsigned)i] );
		}

		vector<vector<Double_t> > returnable_data = Unique_Coords( to_be_sorted );
		return returnable_data;
	}
}

TGraph2D* Plotter( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom3* random )
{
	//      Get a list of ONLY unique coordinates due to the short comings of the interpolation within TGraph2D
	vector<vector<Double_t> > Coord_Data = Plotter_Data( input_tree, Draw_String, Cut_String, random );

	//      Make a new EMPTY TGraph2D with a unique name
	TGraph2D* new_plot = new TGraph2D();
	TString Plot="Plot_";
	double temp = random->Rndm();
	Plot+=temp;
	new_plot->SetName(Plot);
	new_plot->SetTitle(Plot);
	//      Set Binning of Plot
	new_plot->SetNpx( int( ceil( sqrt( double(Coord_Data.size()) ) ) ) );
	new_plot->SetNpy( int( ceil( sqrt( double(Coord_Data.size()) ) ) ) );

	//      Add the data to the TGraph2D object as points
	for( unsigned int i=0; i< Coord_Data.size(); ++i )
	{
		//cout << i << "\t" << Coord_Data[i][0] << "\t" << Coord_Data[i][1] << "\t" << Coord_Data[i][2] << endl;
		new_plot->SetPoint( int(i), Coord_Data[i][0], Coord_Data[i][1], Coord_Data[i][2] );
	}

	//      Return the TGraph2D object
	return new_plot;
}


//      Produce Plotting Histograms from an input TTree
TH2D* Plot_From_Cut( TTree* wanted_tree, TString Draw_String, TString Cut_String, TRandom3* random, TString param1, TString param2 )
{
	//      Create a canvas
	TString Canvas_Name("Canvas_");
	double rand = random->Rndm();
	Canvas_Name+=rand;
	TCanvas* temp_canvas = new TCanvas( Canvas_Name, Canvas_Name );

	//      Plot the Graph
	TGraph2D* new_graph = Plotter( wanted_tree, Draw_String, Cut_String, random );

	//      Update the Canvas to initialize anything that requires this step... this is very much a ROOT thing
	temp_canvas->Update();

	//      Return the Histogram from within this graph for plotting
	TH2D* Returnable_Hist = new_graph->GetHistogram();
	if( Returnable_Hist != NULL )
	{
		Returnable_Hist->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param2 ) );
		Returnable_Hist->GetYaxis()->SetTitle( EdStyle::GetParamRootName( param1 ) );
		return Returnable_Hist;
	}
	return NULL;
}

//      Produce Plotting Histograms from an input TTree
TGraph2D* Plot_From_Cut_lo( TTree* wanted_tree, TString Draw_String, TString Cut_String, TRandom3* random, TString param1, TString param2 )
{
	//      Create a canvas
	TString Canvas_Name("Canvas_");
	double rand = random->Rndm();
	Canvas_Name+=rand;
	TCanvas* temp_canvas = new TCanvas( Canvas_Name, Canvas_Name );

	//      Plot the Graph
	TGraph2D* new_graph = Plotter( wanted_tree, Draw_String, Cut_String, random );

	//      Update the Canvas to initialize anything that requires this step... this is very much a ROOT thing
	temp_canvas->Update();

	//      Return the Histogram from within this graph for plotting
	new_graph->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param2 ) );
	new_graph->GetYaxis()->SetTitle( EdStyle::GetParamRootName( param1 ) );

	return new_graph;
}

//	Return a TTree object composed from a vector of vector of data objects, I provide a stupid branch naming scheme
TTree* vecvec2TTree( vector<vector<Float_t> > input_vec )
{

	TTree* new_tree = new TTree( "tree2Draw", "tree2Drw" );

	Float_t* Float_data = new Float_t[ input_vec[0].size() ];
	for( unsigned int i=0; i<input_vec[0].size(); ++i )
	{
		TString br_name("Branch_");
		br_name+=i;
		TString br_title(br_name);
		br_title.Append("/F");
		new_tree->Branch( br_name, &Float_data[i], br_title );
	}

	for( unsigned int i=0; i<input_vec.size(); ++i)
	{
		for( unsigned int j=0; j<input_vec[i].size(); ++j )
		{
			Float_data[j] = input_vec[i][j];
			new_tree->Fill();
		}
	}

	return new_tree;
}

//	This will be removed in future versions
//pair<TH2*,TPolyMarker3D>* Plot_From_Cut_lo( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom3* random, TString param1, TString param2 )
//{
//	vector<vector<Float_t> > Coord_Data = Plotter_Data( input_tree, Draw_String, Cut_String, random );
//
//	TRandom3* rand = new TRandom(0);
//
//	TString name("namez");
//	name+=rand->Rndm();
//
//	TCanvas* mc = new TCanvas( name, name, 1680, 1050 );
//
//	TTree* plotting_tree = vecvec2TTree( Coord_Data );
//
//	plotting_tree->Draw("Branch_0:Branch_1:Branch_2");
//
//	mc->Update();
//	TPolyMarker3D *pm = (TPolyMarker3D*)gPad->FindObject("TPolyMarker3D");
//	temp = rand->Rndm();
//	TString Name = "TPoly3_";
//	Name+=temp;
//	pm->SetName(Name);
//
//	TH2* returnable_hist = (TH2*) plotting_tree->GetHistogram();
//
//	returnable_hist->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param2 ) );
//	returnable_hist->GetYaxis()->SetTitle( EdStyle::GetParamRootName( param1 ) );
//	pair<TH2*,TPolyMarker3D> returnable_pair;
//	returnable_pair.first = returnable_hist;
//	returnable_pair.second = pm;
//	return returnable_pair;
//}

//      Produce plots for all Physics Parameters stored within the given TTree
//      In theory could also extract the parameters from the input TTree, but as the user likely knows what they want, ask-em
void Physics_Plots( vector<TString> all_parameter_values, Float_t* best_fit_values, TTree* input_tree, TRandom3* rand_gen, TString Param1_Param2, bool CV_Drift, TH2** Physics_Param_Plots, TString Cut_String)
{
	//      Construct a plot string for the physics parameters and plot them
	cout << endl << "STARTING PLOTS SHOWING THE VARIATION OF PHYSICS PARAMETERS THROUGHOUT THE SCANS" << endl;

	//      Store all of the plotting strings
	TString* Physics_Param_DrawString = new TString[unsigned(all_parameter_values.size())];
	//      Loop over all of the wanted parameters that have been passed
	for( unsigned int i=0; i<all_parameter_values.size(); ++i )
	{
		cout << endl << "PLOTTING: " << all_parameter_values[i] << "\t" << i+1 << " of: " << all_parameter_values.size() << endl;

		//      Construct Plotting String based on input
		Physics_Param_DrawString[i] = "(" + all_parameter_values[i];
		if( CV_Drift )
		{
			TString Best_Value;
			Best_Value+=best_fit_values[i];
			Physics_Param_DrawString[i] += "-" + Best_Value;
		}
		Physics_Param_DrawString[i] += ")";
		Physics_Param_DrawString[i] += Param1_Param2;

		//      Actually perform the plots using the same tool as before
		cout << Physics_Param_DrawString[i] << endl;
		cout << Cut_String << endl;
		Physics_Param_Plots[i] = Plot_From_Cut( input_tree, Physics_Param_DrawString[i], Cut_String, rand_gen );

		//      Plot the graph
		TString Canvas_Name("Plot_Me_");
		double rand_canv = rand_gen->Rndm();
		Canvas_Name+=rand_canv;
		TCanvas* plot_me = new TCanvas( Canvas_Name, Canvas_Name, 1680, 1050 );
		plot_me->SetTitle( "" );
		plot_me->SetName( all_parameter_values[i] );
		Physics_Param_Plots[i]->Draw("colz");
		plot_me->Update();
	}
}

//	Actually plot the Physics Parameter plots properly
void Finalize_Physics_Plots( TH2* All_Physics_Plots[], vector<TString> all_parameter_values, TString param1string, TString param2string, TString outputdir, bool CV_Drift )
{
	for( unsigned short int i=0; i < all_parameter_values.size(); ++i )
	{
		TString Name("Physics_Param_");
		Name+=i;

		TCanvas* final_physics_canvas = new TCanvas( Name, Name, 1680, 1050 );

		//	Use the RapidFit naming scheme to translate the parameter name back into latex characters
		All_Physics_Plots[i]->SetTitle("");
		All_Physics_Plots[i]->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param2string ) );
		All_Physics_Plots[i]->GetYaxis()->SetTitle( EdStyle::GetParamRootName( param1string ) );

		if( CV_Drift ) 	all_parameter_values[i].Append("_var");

		vector<TString> Draw_Styles; Draw_Styles.push_back("colz"); Draw_Styles.push_back("cont1z");

		for( unsigned int j=0; j< Draw_Styles.size(); ++j )
		{
			if( j==1 ) All_Physics_Plots[i]->SetContour(80);
			All_Physics_Plots[i]->Draw(Draw_Styles[j]);
			addLHCbLabel( all_parameter_values[i] )->Draw();

			final_physics_canvas->Update();

			//	See the object for FC Stats as to why this is here
			TPaletteAxis *palette = (TPaletteAxis*)All_Physics_Plots[i]->GetListOfFunctions()->FindObject("palette");
			if( palette != NULL )
			{
				palette->SetX1NDC(0.957);
				palette->SetX2NDC(0.962);
				palette->SetLabelSize((Float_t)0.02);
				palette->GetAxis()->SetTickSize(0);
				final_physics_canvas->Update();
			}
			//	Output Filename for the Plot
			TString Output_File( outputdir );
			Output_File+= "/" + all_parameter_values[i]+"_"+Draw_Styles[j]+"_";

			final_physics_canvas->Print( Output_File + ".png" );
			final_physics_canvas->Print( Output_File + ".pdf" );
		}
	}
	return;
}

//	Plot a 2D grid of points that contain data on a TGraph2D
//	This allows for failing grid points and misssing grid jobs to be visualised and compared on a similar scale to explain any anomolies in other plots
void LL2D_Grid( TTree* input_tree, TString Cut_String, TString param1_val, TString param2_val, TRandom3* random, TString Suffix, TString outputdir )
{
	TString Name("Canvas");
	double rand = random->Rndm();
	Name+=rand;
	TCanvas* GRID = new TCanvas(Name,Name,1680,1050);
	input_tree->SetEstimate(input_tree->GetEntries());  // Fix the size of the array of doubles to be created (There will never be more than this
	TString Draw_Str = param1_val + "_value:" + param2_val + "_value";
	input_tree->Draw( Draw_Str, Cut_String );
	TGraph* GRID_Graph = new TGraph( (int)input_tree->GetSelectedRows(), input_tree->GetV2(), input_tree->GetV1() );
	rand = random->Rndm();
	TString GName= "Coord";GName+=rand;
	//	Use Black squares to indicate if data was here or not, nothing more intelligent, nothing less
	GRID_Graph->SetName( GName );
	GRID_Graph->SetMarkerStyle(21);
	GRID_Graph->SetMarkerSize(3);
	GRID_Graph->SetTitle("");
	Name+=random->Rndm();
	GRID_Graph->SetName( Name );
	GRID_Graph->Draw("P");
	GRID->Update();
	//	use the RapidFit nameing scheme
	GRID_Graph->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param2_val ) );
	GRID_Graph->GetYaxis()->SetTitle( EdStyle::GetParamRootName( param1_val ) );
	GRID->Print( outputdir + "/Coordinate_Grid"+Suffix+".png");
}


//	Take a contour and plot on a given set of contours and provide, 4 plots:
//	Tempterature plot,
//	Contour Plot with 40 levels
//	Contour Plot with the reuqested contours
//	Black and White version of the above
//	in .png & .pdf
void Plot_Styled_Contour( TH2* input_hist, int cont_num, double* input_conts, double* confs, TString outputdir, TString Name )
{
	vector<TString> Plot_Type;
	Plot_Type.push_back( "Temp" );
	Plot_Type.push_back( "Cont" );
	Plot_Type.push_back( "Conf" );
	Plot_Type.push_back( "Pub" );

	TString Draw_String;
	TString Canvas_Name("Styled_Canvas_");

	TString Base_Name = outputdir + "/" + Name + "_";

	for( unsigned int i = 0; i < Plot_Type.size(); ++i )
	{
		//	Define Plot Type for the DrawString
		if( i == 0 )	Draw_String = "colz";
		if( i == 1 )	Draw_String = "cont1z";

		if( i == 2 || i == 3 )
		{
			Draw_String = "cont LIST";
			input_hist->SetContour( cont_num, input_conts );
		}

		Canvas_Name+=i;
		TCanvas* Styled_Output_Canvas = new TCanvas( Canvas_Name, Canvas_Name, 1680, 1050 );

		//	Actually call plot function to create EVERYTHING on the convas
		//	THIS SHOULD BE 1 COMMAND BUT ISN'T IN ROOT...
		input_hist->Draw( Draw_String );
		Styled_Output_Canvas->Update();

		//	Temperature Plot
		if( i == 0 )
		{
			input_hist->Draw( Draw_String );
			Styled_Output_Canvas->Update();
			TPaletteAxis *palette = (TPaletteAxis*)input_hist->GetListOfFunctions()->FindObject("palette");
			palette->SetX1NDC(0.957);
			palette->SetX2NDC(0.962);
			palette->SetLabelSize((Float_t)0.02);
			palette->GetAxis()->SetTickSize(0);
			Styled_Output_Canvas->Update();
		}

		//	Lots of Contours
		if( i == 1 )
		{
			input_hist->SetContour(40);	
		}


		//	User requested Contours
		if( i == 2 || i == 3 )
		{
			TObjArray *contObjArr = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

			TList* contLevel = NULL;
			TGraph* curv = NULL;
			TGraph* gc = NULL;
			double cl=0;
			TString confname;

			int TotalConts = contObjArr->GetSize();

			TLegend *leg = new TLegend(0.80,0.89,0.95,0.7);
			leg->SetHeader("Conf. Levels");
			leg->SetBorderSize(0);
			leg->SetFillStyle(0);

			input_hist->Draw("AXIS");

			for(int j = 0; j < TotalConts; ++j )
			{
				confname = "";
				cl = confs[j];
				//confname +=cl;
				stringstream lim;
				lim << setprecision(4) << cl;
				confname.Append( lim.str() );
				confname += "% C.L.";
				contLevel = (TList*)contObjArr->At(j);
				for(int k =0; k < contLevel->GetSize(); ++k)
				{
					curv = (TGraph*)contLevel->At(k);
					gc = (TGraph*)curv->Clone();
					int offset1=1, offset2=2;
					if( j+offset1 >= 5 ) ++offset1;
					if( j+offset2 >= 5 ) ++offset2;
					if( Plot_Type[i] == "Pub"  )	gc->SetLineStyle( Style_t(j+offset1) );
					if( Plot_Type[i] == "Conf" )	gc->SetLineColor( Color_t(j+offset2) );
					gc->Draw("L");
				}
				leg->AddEntry( gc, confname, "L");
			}

			leg->Draw();
		}

		TString Output_Name = Base_Name + Plot_Type[i];

		//	Actually Output the graphs to disk
		addLHCbLabel( Name )->Draw();
		input_hist->SetTitle("");
		Styled_Output_Canvas->Update();

		double phi[1] = {-0.036};
		double phiE[1] = {0.002};
		double dg[1] = {0.087};
		double dgE[1] = {0.021};
		TGraphErrors* sm = new TGraphErrors( 1, phi, dg, phiE, dgE);
		sm->SetMarkerStyle(21);

		sm->Draw("PLSAME");
		Styled_Output_Canvas->Update();

		Styled_Output_Canvas->Print( Output_Name + ".png" );
		Styled_Output_Canvas->Print( Output_Name + ".pdf" );

	}

	//	Save the LL Contour that was user requested as a TGraph2D Object in a ROOT file

	//	Lets conserve the efforts of this plotting tool :D
	cout << endl << "Writing TH2D object for comparisons with other scans" << endl << endl;

	TString Hist_FileName = outputdir+"/"+Name+".root";

	TFile * output = new TFile( Hist_FileName, "RECREATE" );

	input_hist->Write();
	output->Close();

	// Return to default
	input_hist->SetContour(20);
	return;
}


//	Plot 2 canvased on top of each other
//	This allows for 2 sets of contours (of equal number) to be chosen for each dataset
//	
//	TODO:	Allow for the Legend to be customised, easy step just takes time to code up
//
//	This adopts the namescheme to plot a FC scan atop a 2DLL plot but is generically written enough to allow 2 group TH2 objects to be placed atop each other
//
void Plot_Both( TH2* pllhist, TH2* FC_Plot, int nconts, double* fcconts, double *llconts, double* confs, TString outputdir, TString Legend_Name_1, TString Legend_Name_2 )
{
	TCanvas* Temp_1 = new TCanvas("Temp","Temp",1680,1050);  

	TList* contLevel = NULL;
	TGraph* Line     = NULL;
	vector<TGraph*> Contour_Lines;

	//TString Legend_Name_1 = "NLL";
	//TString Legend_Name_2 = "FC";

	//	Construct the Legend
	//TLegend *leg = new TLegend(0.75,0.89,0.95,0.99);
	TLegend *leg = new TLegend(0.80,0.89,0.95,0.7);
	leg->SetHeader("Conf. Levels");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);

	//	Construct the Contours in the LL plot
	pllhist->SetContour( nconts, llconts );
	pllhist->Draw("cont LIST");
	Temp_1->Update();

	//	Get the Contours in the LL Plot
	TObjArray *LL_Contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

	//	Loop over all contours
	for( int i = 0; i < LL_Contours->GetSize(); ++i )
	{
		//	Name that contour...
		TString confname = "";
		stringstream lim;
		lim << setprecision(4) << confs[i];
		confname.Append( lim.str() );
		confname.Append( "% C.L. " + Legend_Name_1 );

		//	Get the List of lines making up this contour
		contLevel = (TList*) LL_Contours->At(i);

		//	Loop over all lines constructing this contour (max 4)
		for(int j =0; j < contLevel->GetSize(); ++j)
		{
			//	Current line
			Line = (TGraph*) contLevel->At(j);

			//	Change the contour line Style
			TGraph *gc = (TGraph*) Line->Clone();
			gc->SetLineStyle( Style_t(i+1) );

			//	Store this line
			Contour_Lines.push_back(gc);

			//	Add an entry in the Legend for it
			if( j==0 )  leg->AddEntry( gc, confname, "L");
		}
	}

	//	Construct the Legend
	FC_Plot->SetContour(nconts,fcconts);
	FC_Plot->Draw("cont LIST");
	Temp_1->Update();

	//	Get the Contours in the FC Plot
	TObjArray *FC_Contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	Temp_1->Update();

	//	Loop over all contours
	for(int i = 0; i < FC_Contours->GetSize(); ++i)
	{
		//	Name that contour...
		TString confname = "";
		confname += confs[i];
		confname.Append( "% C.L. " + Legend_Name_2 );

		//	Get the List of lines making up this contour
		contLevel = (TList*) FC_Contours->At(i);

		//	Loop over all lines constructing this contour (max 4)
		for(int j =0; j<contLevel->GetSize(); ++j)
		{
			//	Current line
			Line = (TGraph*) contLevel->At(j);

			//	Set the line Color
			TGraph *gc = (TGraph*) Line->Clone();
			gc->SetLineColor( Color_t(i+2) );

			//	Store this line
			Contour_Lines.push_back(gc);

			//	Add an entry in the Legend for it
			if(j==0)  leg->AddEntry( gc, confname, "L");
		}
	}

	TCanvas* Overlay_Output = new TCanvas( "Output_Overlay", "Output_Overlay", 1680, 1050);

	//	First construct the Axis
	pllhist->Draw("AXIS");

	//	Now plot all of the lines from all of the contours
	for( unsigned int i = 0; i < Contour_Lines.size(); ++i )
		Contour_Lines[i]->Draw("L SAME");

	//	Add the rest of the details to the plot
	TString Label_Str( Legend_Name_1 + " & " + Legend_Name_2 + " Overlay" );
	addLHCbLabel( Label_Str )->Draw();
	leg->Draw();
	Overlay_Output->Update();

	TString Overlay_FileName( outputdir + "/LL_FC_Overlay");

	Overlay_Output->Print( Overlay_FileName + ".png" );
	Overlay_Output->Print( Overlay_FileName + ".pdf" );

	//      Return the plots back to their input state
	pllhist->SetContour(20);
	FC_Plot->SetContour(20);
}

//	Actually Perform an FC analysis on a RapidFit dataset containing the data
//
//	The data is extraced from the input_tree using various cutstrings built up from data in the 2DLL data contained within the dataset
//
//	the FC_Output object has output data placed in as the analysis is performed which allows for the FC_Stat object to simply review the output
//	this stored things such as job efficiency/stability
TH2D* FC_TOYS( TTree* input_tree, TString Fit_Cut_String, TString param1, TString param2, TString NLL, TString Fit_Cut, double NLL_Global_Best, TTree* FC_Output, TString Double_Tolerance, TRandom3* random )
{

	TString param1_gen = param1 + "_gen";
	TString param2_gen = param2 + "_gen";
	TString param1_val = param1 + "_value";
	TString param2_val = param2 + "_value";
	TString param1_err = param1 + "_error";
	TString param2_err = param2 + "_error";

	TTree* Grid = Cut( input_tree, Fit_Cut_String, random );
	//	The TTree object from the PLL fit contains information on the corrdinates of the fit
	Float_t Param_1_Coord=0, Param_2_Coord=0, NLL_Local_Best=0;
	Grid->SetBranchAddress( param1_val, &Param_1_Coord );
	Grid->SetBranchAddress( param2_val, &Param_2_Coord );
	Grid->SetBranchAddress( NLL, &NLL_Local_Best );

	//	Easiest to store the output in a new TTree
	//	Want to store the CL calculated at every point and the efficiency of Toys at this coordinate
	Float_t CL=0.;
	Int_t Toy_Num=0, Successful_Toys=0, Processed_Toys=0, Fit_Status_Val=3;
	TString CL_Branch = "CL";
	TString CL_Branch_Type = "CL/F";
	TString Success_Branch = "Success";
	TString Success_Branch_Type = "Success/I";
	TString Total_Branch = "Total";
	TString Total_Branch_Type = "Total/I";
	TString Processed_Toys_Branch = "Processed_Toys";
	TString Processed_Toys_Branch_Type = "Processed_Toys/I";
	TString param1_val_type = param1_val + "/F";
	TString param2_val_type = param2_val + "/F";
	TString Fit_Status_Str = "Fit_Status";
	TString Fit_Status_Type = Fit_Status_Str + "/I";
	FC_Output->Branch( CL_Branch, &CL, CL_Branch_Type );
	FC_Output->Branch( Success_Branch, &Successful_Toys, Success_Branch_Type );
	FC_Output->Branch( Total_Branch, &Toy_Num, Total_Branch_Type );
	FC_Output->Branch( param1_val, &Param_1_Coord, param1_val_type );
	FC_Output->Branch( param2_val, &Param_2_Coord, param2_val_type );
	FC_Output->Branch( Processed_Toys_Branch, &Processed_Toys, Processed_Toys_Branch_Type );
	FC_Output->Branch( Fit_Status_Str, &Fit_Status_Val, Fit_Status_Type );

	vector<vector<Float_t> > Used_Coordinate;

	//	Move the object definitions outside of the loop
	bool Add_Point=true;
	bool decision_1=true;
	bool decision_2=true;
	int j=0;
	int Floated_Toy_Num=0;
	int Fixed_Toy_Num=0;

	//double rand=0;

	//Find the toys generated at the previously discovered gridpoints
	for( int i = 0; i < Grid->GetEntries(); ++i )
	{
		vector<Double_t> Floated_Toys_NLL;
		vector<Double_t> Fixed_Toys_NLL;

		//	Read in information about this Grid_Coordinate, coordinate and the NLL from the local best fit
		Grid->GetEntry( i );

		Add_Point = true;
		for( j=0; j < (int)Used_Coordinate.size(); ++j )
		{
			decision_1 = true;
			decision_2 = true;
			if( fabs( Param_1_Coord  - Used_Coordinate[(unsigned)j][0] ) < 1E-5 ) { decision_1 = false; }
			if( fabs( Param_2_Coord  - Used_Coordinate[(unsigned)j][1] ) < 1E-5 ) { decision_2 = false; }
			if( ( decision_1 == false ) && ( decision_2 == false ) ) { Add_Point = false; }
		}
		if( Add_Point ){
			vector<Float_t> temp_Coord;
			temp_Coord.push_back( Param_1_Coord );
			temp_Coord.push_back( Param_2_Coord );
			Used_Coordinate.push_back( temp_Coord );
		} else { continue; }

		//	Cut on only toys at this Grid coord
		TString Param_1_Coord_Str;
		Param_1_Coord_Str+=Param_1_Coord;
		TString Param_2_Coord_Str;
		Param_2_Coord_Str+=Param_2_Coord;
		TString Param_1_Grid = "(abs(" + param1_gen + "-" + Param_1_Coord_Str + ")<" + Double_Tolerance + ")";
		TString Param_2_Grid = "(abs(" + param2_gen + "-" + Param_2_Coord_Str + ")<" + Double_Tolerance + ")";

		TString Toys_At_Grid_Point = Param_1_Grid + "&&" + Param_2_Grid;

		Floated_Toys_NLL = *(Buffer_Branch( input_tree, Toys_At_Grid_Point, NLL ));
		Toy_Num = (int)Floated_Toys_NLL.size();

		Toys_At_Grid_Point += "&&" + Fit_Cut;

		Floated_Toys_NLL = *(Buffer_Branch( input_tree, Toys_At_Grid_Point, NLL ));
		Successful_Toys = (int)Floated_Toys_NLL.size();

		TString float_param_1 = "(abs(" + param1_err + ")>" + Double_Tolerance + ")";
		TString float_param_2 = "(abs(" + param2_err + ")>" + Double_Tolerance + ")";

		TString Floated_Toys_At_Grid_Point = Toys_At_Grid_Point + "&&" + float_param_1 + "&&" + float_param_2;

		TString fixed_param_1 = "(abs(" + param1_err + ")<" + Double_Tolerance + ")";
		TString fixed_param_2 = "(abs(" + param2_err + ")<" + Double_Tolerance + ")";

		TString Fixed_Toys_At_Grid_Point = Toys_At_Grid_Point + "&&" + fixed_param_1 + "&&" + fixed_param_2;

		//	Get all toys at this grid point that are Floated

		Floated_Toys_NLL = *(Buffer_Branch( input_tree, Floated_Toys_At_Grid_Point, NLL ));

		//	Get all toys at this grid point that are Fixed
		Fixed_Toys_NLL = *(Buffer_Branch( input_tree, Fixed_Toys_At_Grid_Point, NLL ));

		Floated_Toy_Num = int( Floated_Toys_NLL.size() );
		Fixed_Toy_Num = int( Fixed_Toys_NLL.size() );

		if( Floated_Toy_Num != Fixed_Toy_Num )
		{
			cerr << endl << "NUMBER OF FIXED AND FLOATED TOYS AT COORDINATE:\t" << Param_1_Coord << ":" << Param_2_Coord <<endl;
			cerr << "ARE DIFFERENT, USING THE SAME SAMBLE SIZE OF EACH WHICH REDUCES THE ACCURACY" << endl << endl;
			if( Floated_Toy_Num < Fixed_Toy_Num ) Fixed_Toy_Num = Floated_Toy_Num;
			else Floated_Toy_Num = Fixed_Toy_Num;
		}

		//	By definition here!
		Processed_Toys = Floated_Toy_Num;

		if( Fixed_Toy_Num != 0 ){

			//Loop over the toys, pulling out the NLL ratio
			UInt_t smaller_Toys = 0;

			Double_t Ratio = NLL_Local_Best - NLL_Global_Best;
			for( j = 0; j < Fixed_Toy_Num; ++j){

				//THE LINE BELOW IS THE FELDMAN-COUSINS ORDERING METHOD USED BY CDF/HEIDELBERG: 
				//if the toyratio is smaller than the data ratio at this point, increment:
				if ( (Fixed_Toys_NLL[(unsigned)j]-Floated_Toys_NLL[(unsigned)j]) < Ratio )
				{
					++smaller_Toys;
				}
			}

			//The C.L. is the percentage of toys that were smaller
			CL = Float_t(smaller_Toys)/Float_t(Fixed_Toy_Num);

		} else {
			//	THIS IS SERIOUS, tell the user about it!
			cerr << "\t" << Param_1_Coord << ":" << Param_2_Coord << "\tWARNING: NO TOYS FOUND HERE! " << endl;
			CL = +9999.;	//	This plots a spike here which shows on contours
		}
		cout << Processed_Toys << "\tTOYS PROCESSED AT:\t" << setprecision(4) << setw(10) << Param_1_Coord << "   :" << setw(10) << Param_2_Coord << endl;
		//	Store the relevent information for plotting in the FC_Output TTree
		FC_Output->Fill();
	}

	TString FCName="FC_Plot_";
	double temp = random->Rndm();
	FCName+=temp;

	TCanvas* FC_Plot = new TCanvas( FCName, FCName, 1680, 1050 );
	TString Param1_Param2 = ":" + param1_val + ":" + param2_val;

	FC_Output->Draw( "CL" + Param1_Param2 );
	FC_Plot->Update();
	TPolyMarker3D *pm = (TPolyMarker3D*)gPad->FindObject("TPolyMarker3D");
	temp = random->Rndm();
	TString Name = "TPoly3_";
	Name+=temp;
	pm->SetName(Name);

	vector<vector<Double_t> > Coord_Data = Unique_Coords( pm );

	TGraph2D* new_plot = new TGraph2D();
	TString Plot="FC_Plot_";
	temp = random->Rndm();
	Plot+=temp;
	new_plot->SetName(Plot);
	new_plot->SetTitle(Plot);
	//	Set Binning of Plot
	new_plot->SetNpx( int( ceil( sqrt( double(Coord_Data.size()) ) ) ) );
	new_plot->SetNpy( int( ceil( sqrt( double(Coord_Data.size()) ) ) ) );

	for( unsigned int i=0; i< Coord_Data.size(); ++i )
	{
		//cout << i << "\t" << Coord_Data[i][0] << "\t" << Coord_Data[i][1] << "\t" << Coord_Data[i][2] << endl;
		new_plot->SetPoint( int(i), Coord_Data[i][0], Coord_Data[i][1], Coord_Data[i][2] );
	}

	TH2D* Returnable_Hist = new_plot->GetHistogram();

	Returnable_Hist->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param2 ) );
	Returnable_Hist->GetYaxis()->SetTitle( EdStyle::GetParamRootName( param1 ) );

	return Returnable_Hist;
}

//	Plot some stats from the FC_Output ttree which stores some data from the full FC calculations
void FC_Stats( TTree* FC_Output, TString param1, TString param2, TRandom3* rand, TString outputdir )
{
	TString Name("FC_STAT_");
	Name+= rand->Rndm();
	vector<TString> unique;
	TString Name2("FC_STAT_");
	Name2+= rand->Rndm();
	TString Name3("FC_STAT_");
	Name3+= rand->Rndm();
	unique.push_back( Name2 );
	unique.push_back( Name3 );
	TCanvas* FC_Stat_Canvas = new TCanvas( Name, Name, 1680, 1050 );

	TString Total_Toy_Map = param1 + "_value:" + param2 + "_value:" + "Processed_Toys";
	TString Toy_Efficiency_Map = param1 + "_value:" + param2 + "_value:" + "(Processed_Toys/(Total/2))";

	vector<TString> all_plot_str;
	all_plot_str.push_back( Total_Toy_Map );
	all_plot_str.push_back( Toy_Efficiency_Map );
	vector<TString> Map_Names;
	Map_Names.push_back( "FC_Total_Toys.pdf" );
	Map_Names.push_back( "FC_Toy_Efficiency.pdf" );

	for( unsigned int i=0; i< all_plot_str.size(); ++i )
	{
		FC_Output->Draw( all_plot_str[i], "goff" );
		TGraph2D* FC_Stat_Graph = new TGraph2D( (int)FC_Output->GetSelectedRows(), FC_Output->GetV2(), FC_Output->GetV1(), FC_Output->GetV3() );
		FC_Stat_Graph->SetName( unique[i] );
		FC_Stat_Graph->SetTitle( unique[i] );
		FC_Stat_Graph->Draw();
		FC_Stat_Canvas->Update();
		FC_Stat_Graph->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param2 ) );
		FC_Stat_Graph->GetYaxis()->SetTitle( EdStyle::GetParamRootName( param1 ) );
		FC_Stat_Graph->Draw("P");
		FC_Stat_Canvas->Update();
		TH1* FC_Stat_Hist = FC_Stat_Graph->GetHistogram();
		unique[i]+=rand->Rndm();
		FC_Stat_Hist->SetName( unique[i] );
		FC_Stat_Hist->SetTitle( unique[i] );
		FC_Stat_Canvas->cd();
		FC_Stat_Hist->Draw("colz");
		FC_Stat_Canvas->Update();

		//	Due to some part of the construction the labels for the temp key (right of the plot) are plotted off canvas
		//	This is likely some net effect of doing something 'not quite ROOT'
		//
		//	Hence you can either resize the TFrame object containing the actual area where the plot is made,
		//	OR:
		//	Move and reisze the axis (which is something that is much easier and will likely always be less painful

		TPaletteAxis* palette = NULL;

		//	Get the PaletteAxis for temperature graphs which is created as part of the Update axis
		palette = (TPaletteAxis*)FC_Stat_Hist->GetListOfFunctions()->FindObject("palette");
		if( palette != NULL )
		{
			TString pal_name("pal_");
			pal_name+=rand->Rndm();
			palette->SetName(pal_name);
			//	Move it to a slightly better position
			palette->SetX1NDC(0.957);
			palette->SetX2NDC(0.962);
			//	Reduce the text size of the axis to allow the numbers to be shown on part of the plot
			palette->SetLabelSize((Float_t)0.02);
			//	Remove the pointless axis markings...
			palette->GetAxis()->SetTickSize(0);
			//	Store the output
		}
		FC_Stat_Canvas->Update();
		FC_Stat_Canvas->Print( outputdir + "/" + Map_Names[i] );

		//delete FC_Stat_Hist;
		//delete FC_Stat_Graph;
	}
	//delete FC_Throw;
	//delete FC_Stat_Canvas;
	//if( palette != NULL ) delete palette;
	//return;
}


pair<vector<double>,vector<double> > LL_Plot_Histo( TTree* input_TTree, TString Cut_String, double Global_Best_NLL, TString NLL, TString param )
{
	TString cut_val; cut_val+=Global_Best_NLL;
	TString Draw_String = "("+NLL+"-"+cut_val+"):"+param;

	input_TTree->SetEstimate(input_TTree->GetEntries() );
	input_TTree->Draw( Draw_String, Cut_String, "goff" );

	//      Return the Histogram from within this graph for plotting
	TH1* Returnable_Hist = (TH1*) input_TTree->GetHistogram();
	pair<vector<double>,vector<double> > output;
	if( Returnable_Hist != NULL )
	{
		int selected = (int)input_TTree->GetSelectedRows();
		double* param_plot = input_TTree->GetV1();
		double* NLL_plot = input_TTree->GetV2();
		vector<double> temp1, temp2;
		for( int i=0; i< selected; ++i )
		{
			temp1.push_back( param_plot[i] );
			temp2.push_back( NLL_plot[i] );
		}
		output.first = temp1;
		output.second = temp2;
	}
	return output;
}

bool Sort_first_Double( pair<double,double> one_pair, pair<double,double> two_pair )
{
	bool x_larger = one_pair.first < two_pair.first ;
	if( x_larger ) return true;
	else return false;
}

bool Sort_second_Double( pair<double,double> one_pair, pair<double,double> two_pair )
{
	bool y_larger = one_pair.second < two_pair.second ;
	if( y_larger ) return true;
	else return false;
}

bool Unique_2D_Double( pair<double,double> one_pair, pair<double,double> two_pair )
{
	bool x_same = fabs(one_pair.first - two_pair.first) < 1E-3 ;
	bool y_same = fabs(one_pair.second - two_pair.second) < 1E-3 ;
	return x_same && y_same;
}

TGraph* LL_Plot( TTree* input_TTree, TString Cut_String, double Global_Best_NLL, TString NLL, TString param, TRandom3* rand )
{
	TGraph* new_graph = NULL;
	pair<vector<double>,vector<double> > data = LL_Plot_Histo( input_TTree, Cut_String, Global_Best_NLL, NLL, param );

	vector<pair<double,double> > filter = reparam( data );
	sort( filter.begin(), filter.end(), Sort_first_Double );
	sort( filter.begin(), filter.end(), Sort_second_Double );

	vector<pair<double,double> >::iterator it;
	it = unique( filter.begin(), filter.end(), Unique_2D_Double );

	filter.resize( unsigned(it - filter.begin()) );
	pair<vector<double>,vector<double> > plottable = reparam( filter );

	unsigned int size = (unsigned) plottable.first.size();
	double* temp1 = new double[size]; double* temp2 = new double[size];
	for( unsigned int i=0; i< size; ++i )
	{
		temp1[i] = plottable.first[i];
		temp2[i] = plottable.second[i];
	}
	new_graph = new TGraph( (Int_t)size, temp2, temp1 );
	TString Name="Graph";
	Name+=rand->Rndm();
	new_graph->SetName( Name );

	return new_graph;
}


//  The gamma distribution coded up in root is the more general form of that found on wikipedia (there's a surprise)
//  
//  Using the root definition:                                  wiki:
//                              gamma = mean^2 / sigma^2                k     = mu^2 / sigma^2  
//                              beta  = sigma^2 / mean                  theta = sigma^2 / mu
//
//              For:            mu == 0                         The 2 conditions above are ONLY valid for this condition

TF1* gamma_func( int OutputLevel )
{
	//      the gamma function in ROOT has issues with verbosity, let's silence it and then return the verbosity back at the end

	streambuf *nullbuf=NULL, *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;
	ofstream filestr;
	filestr.open ("/dev/null");
	//      If the user wanted silence we point the Std Output Streams to /dev/null
	if( OutputLevel <= -1 )
	{
		cout_bak = cout.rdbuf();
		cerr_bak = cerr.rdbuf();
		clog_bak = clog.rdbuf();
		nullbuf = filestr.rdbuf();
		cout.rdbuf(nullbuf);
		cerr.rdbuf(nullbuf);
		clog.rdbuf(nullbuf);
	}

	TF1* output = new TF1( "gammaf", "[0]*TMath::GammaDist( x, ([1]*[1])/([2]*[2]), 0, ([2]*[2])/[1] )" );

	//  TF1 inherits from TFormula so we can use it's functions to rename the parameters to have consistancy with gaus / landau functions
	output->SetParName( 0, "Constant" );
	output->SetParName( 1, "Mean" );
	output->SetParName( 2, "Sigma" );

	//      Reset Std Output Streams
	if( OutputLevel <= -1 )
	{
		cout.rdbuf(cout_bak);
		cerr.rdbuf(cerr_bak);
		clog.rdbuf(clog_bak);
	}

	return output;
}

TF1* landau_func()
{
	//  Want plotting consistancy so have defined my own landau function
	TF1 *land = new TF1( "mylandau", "[0]*TMath::Landau( x, [1], [2] )" );
	land->SetParName( 0, "Constant" );
	land->SetParName( 1, "Mean" );
	land->SetParName( 2, "Sigma" );
	return land;
}

TString Best_Fit_Function( TH1* local_histogram, int OutputLevel )
{
	//	the gamma function in ROOT has issues with verbosity, let's silence it and then return the verbosity back at the end

	streambuf *nullbuf=NULL, *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;
	ofstream filestr;
	filestr.open ("/dev/null");
	//      If the user wanted silence we point the Std Output Streams to /dev/null
	if( OutputLevel <= -1 )
	{
		cout_bak = cout.rdbuf();
		cerr_bak = cerr.rdbuf();
		clog_bak = clog.rdbuf();
		nullbuf = filestr.rdbuf();
		cout.rdbuf(nullbuf);
		cerr.rdbuf(nullbuf);
		clog.rdbuf(nullbuf);
	}


	TString Fit_Options( "Q" );			// Reduce the verbosity slamming the user during this internal procedure
	TF1* my_landau = landau_func(); (void) my_landau;
	TF1* my_gamma = gamma_func(); (void) my_gamma; // this is used from the global state not here... this displeases me...

	TFitResult* result = NULL;

	//  Try Gaussian function
	result = new TFitResult( local_histogram->Fit ( "gaus", Fit_Options ) );
	Double_t chi_2_gaus ( local_histogram->GetFunction ( "gaus" )->GetChisquare() );
	if( result->Status() != 0 )	chi_2_gaus = DBL_MAX;

	//  The GammaDist function originally complained of invalid results (A LOT) and giving it sensible starting points was a way around this
	my_gamma->SetParameters( local_histogram->GetFunction( "gaus" )->GetParameter( 0 ), local_histogram->GetFunction( "gaus" )->GetParameter( 1 ), local_histogram->GetFunction( "gaus" )->GetParameter( 2 ) );

	//  Try Landau function
	result = new TFitResult( local_histogram->Fit ( "mylandau", Fit_Options ) );
	Double_t chi_2_landau ( local_histogram->GetFunction ( "mylandau" )->GetChisquare() );
	if( result->Status() != 0 )	chi_2_landau = DBL_MAX;

	result = new TFitResult( local_histogram->Fit ( "gammaf", Fit_Options ) );
	Double_t chi_2_gamma_f ( local_histogram->GetFunction( "gammaf" )->GetChisquare() );
	if( result->Status() != 0 )	chi_2_gamma_f = DBL_MAX;

	TString fit_type;

	//  Determine which function gives the lowest fit chi squared
	if( fabs(chi_2_landau) < fabs(chi_2_gaus) )  fit_type.Append( "mylandau" );
	else  if( fabs(chi_2_gamma_f) < fabs(chi_2_gaus) )  fit_type.Append( "gammaf" );
	else  fit_type.Append( "gaus" );

	//      Reset Std Output Streams
	if( OutputLevel <= -1 )
	{
		cout.rdbuf(cout_bak);
		cerr.rdbuf(cerr_bak);
		clog.rdbuf(clog_bak);
	}

	return fit_type;
}

void Silent_Fit( TH1* input_histo, TString fit_type, int OutputLevel )
{
	//      the gamma function in ROOT has issues with verbosity, let's silence it and then return the verbosity back at the end

	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;
	//      If the user wanted silence we point the Std Output Streams to /dev/null
	if( OutputLevel <= -1 )
	{
		cout_bak = cout.rdbuf();
		cerr_bak = cerr.rdbuf();
		clog_bak = clog.rdbuf();
		//  Redirect the errors to the empty void of nullness
		freopen("/dev/null", "w", stderr);
	}
	//  Redirect the errors to the empty void of nullness
	freopen("/dev/null", "w", stderr);
	input_histo->Fit( fit_type, "Q" );
	//      Reset Std Output Streams
	if( OutputLevel <= -1 )
	{
		cout.rdbuf(cout_bak);
		cerr.rdbuf(cerr_bak);
		clog.rdbuf(clog_bak);
	}
}

void Silent_Draw( TCanvas* c1, TH1* input_histo, TString options, int OutputLevel )
{
	//      the gamma function in ROOT has issues with verbosity, let's silence it and then return the verbosity back at the end

	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;
	//      If the user wanted silence we point the Std Output Streams to /dev/null
	if( OutputLevel <= -1 )
	{
		cout_bak = cout.rdbuf();
		cerr_bak = cerr.rdbuf();
		clog_bak = clog.rdbuf();
		//  Redirect the errors to the empty void of nullness
		freopen("/dev/null", "w", stderr);
	}
	input_histo->Draw(options);
	c1->Update();
	//      Reset Std Output Streams
	if( OutputLevel <= -1 )
	{
		cout.rdbuf(cout_bak);
		cerr.rdbuf(cerr_bak);
		clog.rdbuf(clog_bak);
	}
}

void Silent_Print( TCanvas* c1, TString Print_String, int OutputLevel )
{
	//      the gamma function in ROOT has issues with verbosity, let's silence it and then return the verbosity back at the end

	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;
	//      If the user wanted silence we point the Std Output Streams to /dev/null
	if( OutputLevel <= -1 )
	{
		cout_bak = cout.rdbuf();
		cerr_bak = cerr.rdbuf();
		clog_bak = clog.rdbuf();
		//  Redirect the errors to the empty void of nullness
		freopen("/dev/null", "w", stderr);
	}
	c1->Update();
	c1->Print( Print_String );
	//      Reset Std Output Streams
	if( OutputLevel <= -1 )
	{
		cout.rdbuf(cout_bak);
		cerr.rdbuf(cerr_bak);
		clog.rdbuf(clog_bak);
	}
}

TH1* Get_Histo( TTree* input_tree, TString draw_str, TString weight_str, TRandom3* rand )
{
	TString rand_str; rand_str+=rand->Rndm();
	input_tree->Draw( draw_str, weight_str, "goff" );
	TH1* output_histo = (TH1*) input_tree->GetHistogram();
	output_histo->SetName(output_histo->GetName()+rand_str);
	return output_histo;
}

TGraph* Get_Graph( TTree* input_tree, TString draw_str, TString weight_str, TRandom3* rand )
{
	TString rand_str; rand_str+=rand->Rndm();
	TGraph* new_graph = new TGraph( Get_Histo( input_tree, draw_str, weight_str, rand ) );
	new_graph->SetName("Graph_"+rand_str);
	return new_graph;
}

TH1* Get_TH1( vector<Double_t> input, TRandom3* rand, int bins )
{
	TString rand_num; rand_num += rand->Rndm();
	TH1* returnable_hist = new TH1D( "TH1_"+rand_num, "TH1_"+rand_num, bins, get_minimum( input ), get_maximum( input ) );
	vector<Double_t>::iterator index_i = input.begin();
	vector<Double_t>::iterator index_e = input.end();
	for( ; index_i != index_e; ++index_i )
	{
		returnable_hist->Fill( *index_i );
	}
	return returnable_hist;
}

TH2* Get_TH2( vector<vector<Double_t> > input, TRandom3* rand, int bins1, int bins2 )
{
	TString rand_num; rand_num += rand->Rndm();
	vector<Double_t> X_data = input[0];
	vector<Double_t> Y_data = input[1];
	TH2* returnable_hist = new TH2D( "TH2_"+rand_num, "TH2_"+rand_num, bins1, get_minimum(X_data), get_maximum(X_data), bins2, get_minimum(Y_data), get_maximum(Y_data) );
	vector<Double_t>::iterator index_x = X_data.begin();
	vector<Double_t>::iterator index_y = Y_data.begin();
	vector<Double_t>::iterator index_e = X_data.end();
	for( ; index_x != index_e; ++index_x, ++index_y )
	{
		returnable_hist->Fill( *index_x, *index_y );
	}
	return returnable_hist;
}

TH3* Get_TH3( vector<vector<Double_t> > input, TRandom3* rand, int bins1, int bins2, int bins3 )
{
	TString rand_num; rand_num += rand->Rndm();
	vector<Double_t> X_data = input[0];
	vector<Double_t> Y_data = input[1];
	vector<Double_t> Z_data = input[2];
	TH3* returnable_hist = new TH3D( "TH3_"+rand_num, "TH3_"+rand_num, bins1, get_minimum(X_data), get_maximum(X_data), bins2, get_minimum(Y_data), get_maximum(Y_data), bins3, get_minimum(Z_data), get_maximum(Z_data) );
	vector<Double_t>::iterator index_x = X_data.begin();
	vector<Double_t>::iterator index_y = Y_data.begin();
	vector<Double_t>::iterator index_z = Z_data.begin();
	vector<Double_t>::iterator index_e = X_data.end();
	for( ; index_x != index_e; ++index_x, ++index_y, ++index_z )
	{
		returnable_hist->Fill( *index_x, *index_y, *index_z );
	}
	return returnable_hist;
}

