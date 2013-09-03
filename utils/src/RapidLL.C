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
#include "RapidLL.h"
#include "TTree_Processing.h"
#include "Histo_Processing.h"
#include "DoFCAnalysis.h"
#include "StringOperations.h"
#include "RapidFit_Output_File.h"
#include "Mathematics.h"
//	System Headers
#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace::std;


// 	Used to be 3 but whatever the EB wants the EB gets...

unsigned int RapidLL::GetFunctionLineWidth()
{
	return EdStyle::GetLHCbFunctionLineWidth();
}

unsigned int RapidLL::GetAxisWidth()
{
	return EdStyle::GetLHCbAxisLineWidth();
}

/*int RapidLL::PlotRapidFitLL( TTree* input_tree, TString controlled_parameter, TRandom3* rand, vector<string> other_params )
  {
  (void) other_params;

  vector<string> controlled_parameters( 1, controlled_parameter.Data() );

  vector<TString> free_parameters = 

  TString Draw_String = controlled_parameter+value_suffix + ":NLL";

  TString Cut_String = "Fit_Status==3";

  vector<vector<double> > plotting_data = TTree_Processing::Plotter_Data( input_tree, Draw_String, Cut_String, rand );

  return 0;
  }*/

void RapidLL::Help()
{
	cout << endl << "RapidLL which plots the LL scans accepts the following additional Runtime Arguments" << endl;
	cout << endl;
	cout << "--addLHCb" << "\t\t" << "This triggers the tool to add the LHCb logo to plots" << endl;
	cout << endl;
	cout << "--isFinal" << "\t\t" << "This replaces the LHCb Preliminary text box with just an LHCb Stamp on plots" << endl;
	cout << endl;
}

TGraph* RapidLL::PlotRapidLL( TString controlled_parameter, TTree* input_tree, TRandom3* rand_gen, vector<string> other_params )
{
	gStyle->SetTitleOffset((Float_t)0.8,"Y");
	gStyle->SetPadLeftMargin( (Float_t)0.14 );

	gROOT->UseCurrentStyle();
	gROOT->ForceStyle( true );

	bool isFinal=false;
	string isFinalString="--isFinal";
	isFinal = StringOperations::VectorContains( other_params, isFinalString ) != -1;

	bool addLHCb=false;
	string addLHCbString="--addLHCb";
	addLHCb = StringOperations::VectorContains( other_params, addLHCbString ) != -1;

	TString param_string = controlled_parameter;
	TString param_val = param_string+value_suffix;

	input_tree->Draw( param_val,"","goff",1,0 );
	double notgenvalue = input_tree->GetV1()[0];
	TString notgen; notgen+=notgenvalue;

	TString output_path("RapidFit_LLScan_");
	output_path.Append( controlled_parameter );
	output_path.Append( "_" + StringOperations::TimeString() );

	//      Make OutputDir
	TDirectory* output_dir=NULL;

	TString pwd = TString( gSystem->WorkingDirectory() );

	if( gSystem->mkdir( output_path ) < 0 )
	{
		int could_make=-1;
		for( unsigned int i=0; could_make < 0; ++i )
		{
			TString new_output_path = output_path;
			new_output_path.Append("_"); new_output_path+=i;
			could_make = gSystem->mkdir( output_path );
			if( could_make >= 0 )
			{
				output_path=new_output_path;
				break;
			}
		}
	}
	gSystem->cd( output_path );
	output_dir=gDirectory;


	//	Setup Objects for Plotting

	//	Quick(ish) way to establish the smallest NLL from a good fit (by definition this is the global result)
	TTree* temp_tree = input_tree->CopyTree( "NLL>=0","fast", input_tree->GetEntries(),0 );
	double true_min_NLL = temp_tree->GetMinimum("NLL");

	//	General Cuts to be applied for various plots

	//	Fit_Status == 3
	TString Fit_Cut = "(abs(" + Fit_Status + "-3.0)<"+double_tolerance+")";

	//	Combine the individual Cuts
	TString Fit_Cut_String = Fit_Cut;

	TString NLL_Min;		//	Of course ROOT doesn't have USEFUL constructors!
	NLL_Min+=true_min_NLL;

	gStyle->SetLineWidth( (Width_t)RapidLL::GetAxisWidth() );
	gROOT->UseCurrentStyle();
	gROOT->ForceStyle( true );

	//	Do the NLL plotting

	cout << "Plotting NLL variation\n\n";

	//	This function is amazing, it even returns TGraphs with unique names so as to not cause problems in ROOT namespace :D
	TGraph* drawn_histo = LL_Plot( input_tree, Fit_Cut_String, true_min_NLL, NLL, param_val, rand_gen );

	TGraph* color_graph = LL_Plot( input_tree, Fit_Cut_String, true_min_NLL, NLL, param_val, rand_gen );


	TString Name_Base="RapidLL_"+param_string;
	TString Name=Name_Base; Name.Append("_"); Name+=rand_gen->Rndm();

	gSystem->cd( output_path );

	TCanvas* new_canvas = EdStyle::RapidFitCanvas( Name, Name );

	drawn_histo->Draw( "AC*" );
	new_canvas->Update();
	TPaveText* text_1 = NULL;


	if( addLHCb )
	{
		text_1 = EdStyle::LHCbLabel();
		if( !isFinal )
		{
			text_1->AddText("");
			text_1->AddText( "Preliminary" );
		}

		text_1->SetFillStyle(0);
		text_1->Draw("SAME");
	}

	//      Now that the axis exist we can worry about labeling them
	drawn_histo->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param_string ) + " " + EdStyle::GetParamRootUnit( param_string ) );
	drawn_histo->GetYaxis()->SetTitle( EdStyle::GetParamRootName( TString("LLscan") ) );
	new_canvas->Update();

	Histogram_Processing::Silent_Print( new_canvas, Name_Base+"_pub.C");
	Histogram_Processing::Silent_Print( new_canvas, Name_Base+"_pub.png");
	Histogram_Processing::Silent_Print( new_canvas, Name_Base+"_pub.pdf");

	drawn_histo->Draw( "AC" );
	new_canvas->SetTitle("");
	if( text_1 != NULL ) text_1->Draw("SAME");
	new_canvas->Update();
	Histogram_Processing::Silent_Print( new_canvas, Name_Base+"_pub_nopoints.C");
	Histogram_Processing::Silent_Print( new_canvas, Name_Base+"_pub_nopoints.png");
	Histogram_Processing::Silent_Print( new_canvas, Name_Base+"_pub_nopoints.pdf");


	color_graph->Draw("AC*");
	new_canvas->Update();
	if( text_1 != NULL ) text_1->Draw("SAME");

	//      Now that the axis exist we can worry about labeling them
	color_graph->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param_string ) + " " + EdStyle::GetParamRootUnit( param_string ) );
	color_graph->GetYaxis()->SetTitle( EdStyle::GetParamRootName( TString("LLscan") ) );
	new_canvas->Update();

	Histogram_Processing::Silent_Print( new_canvas, Name_Base+"_conf.C");
	Histogram_Processing::Silent_Print( new_canvas, Name_Base+"_conf.png");
	Histogram_Processing::Silent_Print( new_canvas, Name_Base+"_conf.pdf");

	color_graph->Draw("AC");
	new_canvas->SetTitle("");
	if( text_1 != NULL ) text_1->Draw("SAME");
	new_canvas->Update();
	Histogram_Processing::Silent_Print( new_canvas, Name_Base+"_conf_nopoints.C");
	Histogram_Processing::Silent_Print( new_canvas, Name_Base+"_conf_nopoints.png");
	Histogram_Processing::Silent_Print( new_canvas, Name_Base+"_conf_nopoints.pdf");

	TString nuisance_dir = "nuisance_params_"; nuisance_dir.Append( StringOperations::TimeString() );

	gSystem->mkdir( nuisance_dir );
	gSystem->cd( nuisance_dir );

	vector<string> temp_vec( 1, controlled_parameter.Data() );

	vector<TString> nuisance_parameters = RapidFit_Output_File::get_free_non_scanned_parameters( input_tree, temp_vec );

	cout << "Plotting Variation in:" << endl;

	for( vector<TString>::iterator param_i = nuisance_parameters.begin(); param_i != nuisance_parameters.end(); ++param_i )
	{
		TCanvas* param_c = EdStyle::RapidFitCanvas( "param_canv_"+*param_i, "param_canv_"+*param_i );

		cout << *param_i << endl;

		TGraph* param_graph = LL_Plot( input_tree, Fit_Cut_String, 0., *param_i+value_suffix, param_val, rand_gen );
		TGraph* err_graph = LL_Plot( input_tree, Fit_Cut_String, 0., *param_i+error_suffix, param_val, rand_gen );

		param_graph->Draw("AC*");

		param_c->SetTitle("");

		param_c->Update();
		if( text_1 != NULL ) text_1->Draw("SAME");


		//      Now that the axis exist we can worry about labeling them
		param_graph->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param_string + value_suffix ) + " " + EdStyle::GetParamRootUnit( param_string ) );
		param_graph->GetYaxis()->SetTitle( EdStyle::GetParamRootName( *param_i + value_suffix ) + " " + EdStyle::GetParamRootUnit( *param_i ) );

		TString param_file = *param_i;

		Histogram_Processing::Silent_Print( param_c, param_file+".C" );
		Histogram_Processing::Silent_Print( param_c, param_file+".pdf" );
		Histogram_Processing::Silent_Print( param_c, param_file+".png" );


		err_graph->Draw("AC*");

		param_c->SetTitle("");

		param_c->Update();
		if( text_1 != NULL ) text_1->Draw("SAME");

		err_graph->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param_string + value_suffix ) + " " + EdStyle::GetParamRootUnit( param_string ) );
		err_graph->GetYaxis()->SetTitle( EdStyle::GetParamRootName( *param_i + error_suffix ) + " " + EdStyle::GetParamRootUnit( *param_i ) );

		Histogram_Processing::Silent_Print( param_c, param_file+"_"+error_suffix+".C" );
		Histogram_Processing::Silent_Print( param_c, param_file+"_"+error_suffix+".pdf" );
		Histogram_Processing::Silent_Print( param_c, param_file+"_"+error_suffix+".png" );
	}


	gSystem->cd( pwd );

	return color_graph;
}

pair<vector<double>,vector<double> > RapidLL::LL_Plot_Histo( TTree* input_TTree, TString Cut_String, double Global_Best_NLL, TString NLL, TString param )
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



TGraph* RapidLL::LL_Plot( TTree* input_TTree, TString Cut_String, double Global_Best_NLL, TString NLL, TString param, TRandom3* rand )
{
	TGraph* new_graph = NULL;
	pair<vector<double>,vector<double> > data = LL_Plot_Histo( input_TTree, Cut_String, Global_Best_NLL, NLL, param );

	vector<pair<double,double> > filter = reparam( data );
	sort( filter.begin(), filter.end(), Mathematics::Sort_first_Double );
	sort( filter.begin(), filter.end(), Mathematics::Sort_second_Double );

	vector<pair<double,double> >::iterator it;
	it = unique( filter.begin(), filter.end(), Mathematics::Unique_2D_Double );

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
	string rand_str_( Name.Data() );
	replace( rand_str_.begin(), rand_str_.end(), '.', '_' );
	new_graph->SetName( rand_str_.c_str() );
	new_graph->SetTitle("");

	new_graph->SetLineWidth( RapidLL::GetFunctionLineWidth() );

	return new_graph;
}

void RapidLL::OverlayMutliplePlots( TMultiGraph* GraphsToOverlay )
{
	gStyle->SetLineWidth( (Width_t)RapidLL::GetAxisWidth() );
	gROOT->UseCurrentStyle();
	gROOT->ForceStyle( true );

	TString OverlayName="OverlayPlots";
	gSystem->mkdir( OverlayName );
	gSystem->cd( OverlayName );

	string timeStamp = StringOperations::TimeString();

	TCanvas* newOverlay = EdStyle::RapidFitCanvas( "OverlayCanvas", "OverlayCanvas" );

	TLegend* thisLegend = EdStyle::LHCbLegend();
	GraphsToOverlay->Draw("A");
	newOverlay->Update();

	for( unsigned int i=0; i< GraphsToOverlay->GetListOfGraphs()->Capacity(); ++i )
	{
		TGraph* thisGraph = (TGraph*) GraphsToOverlay->GetListOfGraphs()->At(i);
		thisGraph->SetLineColor( (Color_t)(i+1) );
		thisGraph->SetMarkerColor( (Color_t)(i+1) );
		thisGraph->SetFillColor( kWhite );
	}
	GraphsToOverlay->Draw("PC");
	newOverlay->Update();
	for( unsigned int i=0; i< GraphsToOverlay->GetListOfGraphs()->Capacity(); ++i )
	{
		TGraph* thisGraph = (TGraph*) GraphsToOverlay->GetListOfGraphs()->At(i);
		TString Name="#Delta Log Likelihood Function ";Name+=(i+1);
		thisLegend->AddEntry( thisGraph, Name );
	}
	thisLegend->Draw();
	newOverlay->Update();

	GraphsToOverlay->GetXaxis()->SetTitle( ((TGraph*)GraphsToOverlay->GetListOfGraphs()->At(0) )->GetXaxis()->GetTitle() );
	GraphsToOverlay->GetYaxis()->SetTitle( ((TGraph*)GraphsToOverlay->GetListOfGraphs()->At(0) )->GetYaxis()->GetTitle() );
	newOverlay->Update();

	TString Name="Overlay_Graph";
	newOverlay->Print(Name+".png");
	newOverlay->Print(Name+".pdf");
	newOverlay->Print(Name+".C");
	Name.Append("_");Name.Append(timeStamp);
	newOverlay->Print(Name+".png");
	newOverlay->Print(Name+".pdf");
	newOverlay->Print(Name+".C");


	TCanvas* newOverlay2 = EdStyle::RapidFitCanvas( "OverlayCanvas2", "OverlayCanvas2" );

	GraphsToOverlay->Draw("A");
	newOverlay2->Update();

	GraphsToOverlay->Draw("C");
	thisLegend->Draw();
	newOverlay2->Update();

	Name="Overlay_Graph_noPoints";
	newOverlay2->Print(Name+".png");
	newOverlay2->Print(Name+".pdf");
	newOverlay2->Print(Name+".C");
	Name.Append("_");Name.Append(timeStamp);
	newOverlay2->Print(Name+".png");
	newOverlay2->Print(Name+".pdf");
	newOverlay2->Print(Name+".C");

	return;
}

