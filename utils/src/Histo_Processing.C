//	ROOT Headers
#include "TFile.h"
#include "TH1.h"
#include "TAxis.h"
#include "TString.h"
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

using namespace::std;

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
	input_hist->Rebin( int( ( (double) existing_bins ) / (( double)wanted_bins ) ) );

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
	if( axis == 1 ){
		min_range = input_hist->GetXaxis()->GetXmin();
		max_range = input_hist->GetXaxis()->GetXmax();	}
	else if( axis == 2 ){
		min_range = input_hist->GetYaxis()->GetXmin();
		max_range = input_hist->GetYaxis()->GetXmax();	}
	else if( axis == 3 ){
		min_range = input_hist->GetZaxis()->GetXmin();
		max_range = input_hist->GetZaxis()->GetXmax();	}

		double range = max_range - min_range;
		double wanted_bins = ceil( range / width );

		//  Catch SERIOUS rounding errors and caes where the information has already been lost due to underbinning
		if( ( axis == 1 ) && ( wanted_bins > input_hist->GetNbinsX() ) )  return 0;
		else if( ( axis == 2 ) && ( wanted_bins > input_hist->GetNbinsY() ) )  return 0;
		else if( ( axis == 3 ) && ( wanted_bins > input_hist->GetNbinsZ() ) )  return 0;
		else return (int) wanted_bins;

		return 0;
}

//  This function simply rebins a histogram to an optimal number of bins for a gaussian like distribution
void OptimallyRebin( TH1* input_hist, int axis )
{
	int temp = OptimumBinNumber( input_hist, axis );
	(void) temp;	return;
}


//      This will find all of the unique corrdinates stored in a TPolyMarker3D object
//      NB this was written due to the fact that ROOT thows away the contents of {3/4}D TTree->Draw() objects... (God only knows why)
vector<vector<Float_t> > Unique_Coords( TPolyMarker3D *pm )
{
	//      Tolerance of the DOUBLE
	double DT = 1E-5;

	//      Get the number of coordinates in the TPolyMarker3D
	int number_of_points = pm->GetN();
	//      Get the data contained in the TPolyMarker3D
	Float_t* Coord_Data_pointer = pm->GetP();

	//      Will populate and return to the user
	vector<vector<Float_t> > Returnable_Coord_Data;

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
			vector<Float_t> temp_vector;
			temp_vector.push_back( Coord_Data_pointer[i*3] );
			temp_vector.push_back( Coord_Data_pointer[i*3+1] );
			temp_vector.push_back( Coord_Data_pointer[i*3+2] );
			Returnable_Coord_Data.push_back( temp_vector );
		}
	}

	//      Return a 2D vector of unique 3D corrdinates of size npoints*3
	return Returnable_Coord_Data;
}

vector<vector<Float_t> > Plotter_Data( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom3* random )
{
        //      Plot the graph using TTree->Draw()
	//      The resulting graph (TH3) contains empty axis and a TPolyMarker3D object
	input_tree->SetEstimate(input_tree->GetEntries());  // Fix the size of the array of doubles to be created (There will never be more than this
	input_tree->Draw( Draw_String, Cut_String );
                                      
	//      Get the Points that have been plotted (TPolyMarker3D object named "TPolyMarker3D", see ROOTtalk)
	TPolyMarker3D *pm = (TPolyMarker3D*)gPad->FindObject("TPolyMarker3D");
	double temp = random->Rndm();
	TString Name = "TPoly3_";
	Name+=temp;
	pm->SetName(Name);

	//      Get a list of ONLY unique coordinates due to the short comings of the interpolation within TGraph2D
	vector<vector<Float_t> > returnable_data = Unique_Coords( pm );

	return returnable_data;
}

TGraph2D* Plotter( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom3* random )
{
	//      Get a list of ONLY unique coordinates due to the short comings of the interpolation within TGraph2D
	vector<vector<Float_t> > Coord_Data = Plotter_Data( input_tree, Draw_String, Cut_String, random );

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
	Returnable_Hist->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param2 ) );
	Returnable_Hist->GetYaxis()->SetTitle( EdStyle::GetParamRootName( param1 ) );

	return Returnable_Hist;
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

void Finalize_Physics_Plots( TH2* All_Physics_Plots[], vector<TString> all_parameter_values, TString param1string, TString param2string, TString outputdir, bool CV_Drift )
{
	for( unsigned short int i=0; i < all_parameter_values.size(); ++i )
	{
		TString Name("Physics_Param_");
		Name+=i;

		TCanvas* final_physics_canvas = new TCanvas( Name, Name, 1680, 1050 );

		All_Physics_Plots[i]->SetTitle("");
		All_Physics_Plots[i]->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param2string ) );
		All_Physics_Plots[i]->GetYaxis()->SetTitle( EdStyle::GetParamRootName( param1string ) );

		if( CV_Drift ) 	all_parameter_values[i].Append("_var");

		All_Physics_Plots[i]->Draw("colz");
		addLHCbLabel( all_parameter_values[i] )->Draw();

		final_physics_canvas->Update();

		//	Output Filename for the Plot
		TString Output_File( outputdir );
		Output_File+= "/" + all_parameter_values[i] + ".png";

		final_physics_canvas->Print( Output_File );
	}
	return;
}

TH1* LL2D_Grid( TTree* input_tree, TString Cut_String, TString param1_val, TString param2_val, TRandom3* random, TString Suffix )
{
	TString Name("Canvas");
	double rand = random->Rndm();
	Name+=rand;
	TCanvas* GRID = new TCanvas(Name,Name,1680,1050);
	input_tree->SetEstimate(input_tree->GetEntries());  // Fix the size of the array of doubles to be created (There will never be more than this
	TString Draw_Str = param1_val + ":" + param2_val;
	input_tree->Draw( Draw_Str, Cut_String );
	//GRID->Print("Coordinate_Grid"+Suffix+".png");
	(void) GRID;
	(void) Suffix;
	return input_tree->GetHistogram();
}

void Plot_Styled_Contour2( TGraph2D* input_graph, int cont_num, double* input_conts, double* confs, TString outputdir, TString Name )
{
	(void) input_conts; (void) cont_num; 
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

		if( i == 0 )	Draw_String = "colz";
		if( i == 1 )	Draw_String = "cont1z";

		if( i == 2 || i == 3 )
		{
			Draw_String = "cont LIST";
			//input_graph->SetContour( cont_num, input_conts );
		}

		Canvas_Name+=i;
		TCanvas* Styled_Output_Canvas = new TCanvas( Canvas_Name, Canvas_Name, 1680, 1050 );

		input_graph->Draw( Draw_String );
		Styled_Output_Canvas->Update();

		if( i == 0 )
		{
			input_graph->Draw( Draw_String );
			Styled_Output_Canvas->Update();
		}

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

			input_graph->Draw("AXIS");

			for(int j = 0; j < TotalConts; ++j )
			{
				confname = "";
				cl = confs[j];
				confname +=cl;
				confname += "% C.L.";
				contLevel = input_graph->GetContourList(cl);
				//contLevel = (TList*)contObjArr->At(j);
				for(int k =0; k < contLevel->GetSize(); ++k)
				{
					curv = (TGraph*)contLevel->At(k);
					gc = (TGraph*)curv->Clone();
					if( Plot_Type[i] == "Pub"  )	gc->SetLineStyle( Style_t(j+1) );
					if( Plot_Type[i] == "Conf" )	gc->SetLineColor( Color_t(j+2) );
					gc->Draw("L");
				}
				leg->AddEntry( gc, confname, "L");
			}

			leg->Draw();
		}

		TString Output_Name = Base_Name + Plot_Type[i];

		addLHCbLabel( Name )->Draw();
		input_graph->SetTitle("");
		Styled_Output_Canvas->Update();

		Styled_Output_Canvas->Print( Output_Name + ".png" );
		Styled_Output_Canvas->Print( Output_Name + ".pdf" );

	}

	//	Lets conserve the efforts of this plotting tool :D
	cout << endl << "Writing TH2D object for comparisons with other scans" << endl << endl;

	TString Hist_FileName = outputdir+"/"+Name+".root";

	TFile * output = new TFile( Hist_FileName, "RECREATE" );

	input_graph->Write();
	output->Close();

	// Return to default
	//input_hist->SetContour(20);
	return;
}

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

		if( i == 0 )	Draw_String = "colz";
		if( i == 1 )	Draw_String = "cont1z";

		if( i == 2 || i == 3 )
		{
			Draw_String = "cont LIST";
			input_hist->SetContour( cont_num, input_conts );
		}

		Canvas_Name+=i;
		TCanvas* Styled_Output_Canvas = new TCanvas( Canvas_Name, Canvas_Name, 1680, 1050 );

		input_hist->Draw( Draw_String );
		Styled_Output_Canvas->Update();

		if( i == 0 )
		{
			input_hist->Draw( Draw_String );
			Styled_Output_Canvas->Update();
		}

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
				confname +=cl;
				confname += "% C.L.";
				contLevel = (TList*)contObjArr->At(j);
				for(int k =0; k < contLevel->GetSize(); ++k)
				{
					curv = (TGraph*)contLevel->At(k);
					gc = (TGraph*)curv->Clone();
					if( Plot_Type[i] == "Pub"  )	gc->SetLineStyle( Style_t(j+1) );
					if( Plot_Type[i] == "Conf" )	gc->SetLineColor( Color_t(j+2) );
					gc->Draw("L");
				}
				leg->AddEntry( gc, confname, "L");
			}

			leg->Draw();
		}

		TString Output_Name = Base_Name + Plot_Type[i];

		addLHCbLabel( Name )->Draw();
		input_hist->SetTitle("");
		Styled_Output_Canvas->Update();

		Styled_Output_Canvas->Print( Output_Name + ".png" );
		Styled_Output_Canvas->Print( Output_Name + ".pdf" );

	}

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


void Plot_Both( TH2* pllhist, TH2* FC_Plot, int nconts, double* fcconts, double *llconts, double* confs, TString outputdir )
{
	TCanvas* Temp_1 = new TCanvas("Temp","Temp",1680,1050);  

	TList* contLevel = NULL;
	TGraph* Line     = NULL;
	vector<TGraph*> Contour_Lines;

	//	Construct the Legend
	TLegend *leg = new TLegend(0.75,0.89,0.95,0.7);
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
		confname += confs[i];
		confname += "% C.L. (PLL)";

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
		confname += "% C.L. (FC)";

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
	addLHCbLabel( "LL & FC Overlay" )->Draw();
	leg->Draw();
	Overlay_Output->Update();

	TString Overlay_FileName( outputdir + "/LL_FC_Overlay");

	Overlay_Output->Print( Overlay_FileName + ".png" );
	Overlay_Output->Print( Overlay_FileName + ".pdf" );

	//      Return the plots back to their input state
	pllhist->SetContour(20);
	FC_Plot->SetContour(20);
}

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
	Float_t CL=0., Toy_Num=0., Successful_Toys=0., Processed_Toys=0., Fit_Status=3.;
	TString CL_Branch = "CL";
	TString Success_Branch = "Success";
	TString Total_Branch = "Total";
	TString Processed_Toys_Branch = "Processed_Toys";
	FC_Output->Branch( CL_Branch, &CL );
	FC_Output->Branch( Success_Branch, &Successful_Toys );
	FC_Output->Branch( Total_Branch, &Toy_Num );
	FC_Output->Branch( param1_val, &Param_1_Coord );
	FC_Output->Branch( param2_val, &Param_2_Coord );
	FC_Output->Branch( Processed_Toys_Branch, &Processed_Toys );
	FC_Output->Branch( "Fit_Status", &Fit_Status );

	vector<vector<Float_t> > Used_Coordinate;


	bool Add_Point=true;
	bool decision_1=true;
	bool decision_2=true;

	//double rand=0;

	//Find the toys generated at the previously discovered gridpoints
	for( int i = 0; i < Grid->GetEntries(); ++i )
	{
		vector<Double_t> Floated_Toys_NLL;
		vector<Double_t> Fixed_Toys_NLL;

		//	Read in information about this Grid_Coordinate, coordinate and the NLL from the local best fit
		Grid->GetEntry( i );

		Add_Point = true;
		for( unsigned int j=0; j < Used_Coordinate.size(); ++j )
		{
			decision_1 = true;
			decision_2 = true;
			if( fabs( Param_1_Coord  - Used_Coordinate[j][0] ) < 1E-5 ) { decision_1 = false; }
			if( fabs( Param_2_Coord  - Used_Coordinate[j][1] ) < 1E-5 ) { decision_2 = false; }
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

		Floated_Toys_NLL = Get_Data( input_tree, Toys_At_Grid_Point, NLL );
		Toy_Num = Float_t( Floated_Toys_NLL.size() );

		Toys_At_Grid_Point += "&&" + Fit_Cut;

		Floated_Toys_NLL = Get_Data( input_tree, Toys_At_Grid_Point, NLL );
		Successful_Toys = Float_t( Floated_Toys_NLL.size() );

		TString float_param_1 = "(abs(" + param1_err + ")>" + Double_Tolerance + ")";
		TString float_param_2 = "(abs(" + param2_err + ")>" + Double_Tolerance + ")";

		TString Floated_Toys_At_Grid_Point = Toys_At_Grid_Point + "&&" + float_param_1 + "&&" + float_param_2;


		TString fixed_param_1 = "(abs(" + param1_err + ")<" + Double_Tolerance + ")";
		TString fixed_param_2 = "(abs(" + param2_err + ")<" + Double_Tolerance + ")";

		TString Fixed_Toys_At_Grid_Point = Toys_At_Grid_Point + "&&" + fixed_param_1 + "&&" + fixed_param_2;

		//	Get all toys at this grid point that are Floated

		Floated_Toys_NLL = Get_Data( input_tree, Floated_Toys_At_Grid_Point, NLL );
		//TTree* Floated_Toys = input_tree->CopyTree( Floated_Toys_At_Grid_Point, "fast", input_tree->GetEntries() );
		//TString Floated_Name("Floated_");
		//rand = random->Rndm();
		//Floated_Name+=rand;
		//Floated_Toys->SetName(Floated_Name);

		//	Get all toys at this grid point that are Fixed
		Fixed_Toys_NLL = Get_Data( input_tree, Fixed_Toys_At_Grid_Point, NLL );
		//TTree* Fixed_Toys = input_tree->CopyTree( Fixed_Toys_At_Grid_Point, "fast", input_tree->GetEntries() );
		//TString Fixed_Name("Fixed_");
		//rand = random->Rndm();
		//Fixed_Name+=rand;
		//Fixed_Toys->SetName(Fixed_Name);

		int Floated_Toy_Num = int( Floated_Toys_NLL.size() );
		int Fixed_Toy_Num = int( Fixed_Toys_NLL.size() );

		if( Floated_Toy_Num != Fixed_Toy_Num )
		{
			cerr << endl << "NUMBER OF FIXED AND FLOATED TOYS AT COORDINATE:\t" << Param_1_Coord << ":" << Param_2_Coord <<endl;
			cerr << "ARE DIFFERENT, USING THE SAME SAMBLE SIZE OF EACH WHICH REDUCES THE ACCURACY" << endl << endl;
			if( Floated_Toy_Num < Fixed_Toy_Num ) Fixed_Toy_Num = Floated_Toy_Num;
			else Floated_Toy_Num = Fixed_Toy_Num;
		}

		//	By definition here!
		Processed_Toys = Float_t(Floated_Toy_Num);

		//	Remember coordinates stored in param1gridpoints and param2gridpoints
		//	With NLL values stored in NLL_Local_Best
		//	This again would be nicer&easier if associated data was stored in the one object


		//	Using GetEntry to look at each toy at each point so tell it where to store the value we get
		//Float_t Floated_NLL=0, Fixed_NLL=0;
		//Floated_Toys->SetBranchAddress( NLL, &Floated_NLL );
		//Fixed_Toys->SetBranchAddress( NLL, &Fixed_NLL );

		if( Fixed_Toy_Num != 0 ){

			//Loop over the toys, pulling out the NLL ratio
			UInt_t smaller_Toys = 0;

			Double_t Ratio = NLL_Local_Best - NLL_Global_Best;
			for(unsigned short int j = 0; j < Fixed_Toy_Num; ++j){

				//	VERY MEMORY INTENSIVE!!!
				//	Floated_Toys->GetEntry(j);
				//	Fixed_Toys->GetEntry(j);

				//THE LINE BELOW IS THE FELDMAN-COUSINS ORDERING METHOD USED BY CDF/HEIDELBERG: 
				//if the toyratio is smaller than the data ratio at this point, increment:
				if ( (Fixed_Toys_NLL[j]-Floated_Toys_NLL[j]) < Ratio )
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
		cout << Processed_Toys << "\tTOYS PROCESSED AT:\t" << setprecision(4) << Param_1_Coord << "\t:\t" << Param_2_Coord << endl;
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

	vector<vector<Float_t> > Coord_Data = Unique_Coords( pm );

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
