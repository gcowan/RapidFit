
#include "TLegend.h"
#include "TString.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TGraphErrors.h"

#include "Rapid2DLL.h"
#include "Histo_Processing.h"
#include "RapidFit_Output_File.h"
#include "TTree_Processing.h"
#include "StringOperations.h"
#include "EdStyle.h"

#include <vector>
#include <string>

using namespace::std;


//	The line widths were 3 and now are 2 and are always whatever the EB decides they will be tomorrrow...

unsigned int Rapid2DLL::GetFunctionLineWidth()
{
	return EdStyle::GetLHCbFunctionLineWidth();
}

unsigned int Rapid2DLL::GetAxisWidth()
{
	return EdStyle::GetLHCbAxisLineWidth();
}

//	This was 0.04 because it looked sensible but the EB font size is 0.060...

double Rapid2DLL::GetLegTextSize()
{
	return EdStyle::GetLHCbTextSize();
}

unsigned int Rapid2DLL::GetColors( unsigned int input )
{
	vector<unsigned int> colors;
	colors.push_back( 2 );
	colors.push_back( 3 );
	colors.push_back( 4 );
	colors.push_back( 6 );
	colors.push_back( 7 );
	colors.push_back( 8 );
	return colors[input];
}

unsigned int Rapid2DLL::GetStyle( unsigned int input )
{
	vector<unsigned int> style;
	style.push_back( 1 );
	style.push_back( 2 );
	style.push_back( 3 );
	style.push_back( 4 );
	style.push_back( 6 );
	style.push_back( 7 );
	style.push_back( 8 );
	return style[input];
}

/*
vector<pair<double,TString> > Rapid2DLL::GetContour( TString input )
{
	vector<pair<double,TString> > returnable;
	if( input == "2DLL" )
	{
		returnable.push_back( make_pair( 1.15, "68 % C.L." ) );
		returnable.push_back( make_pair( 2.36, "90 % C.L." ) );
		returnable.push_back( make_pair( 3., "95 % C.L." ) );
	}
	if( input == "2DLL-99" )
	{
		returnable.push_back( make_pair( 1.15, "68 % C.L." ) );
		returnable.push_back( make_pair( 2.36, "90 % C.L." ) );
		returnable.push_back( make_pair( 3., "95 % C.L." ) );
		returnable.push_back( make_pair( 4.61, "99 % C.L." ) );
	}
	if( input == "FC" )
	{
		returnable.push_back( make_pair( 0.6827, "68 % C.L." ) );
		returnable.push_back( make_pair( 0.90, "90 % C.L." ) );
		returnable.push_back( make_pair( 0.95, "95 % C.L." ) );
	}
	if( input == "FC-99" )
	{
		returnable.push_back( make_pair( 0.6827, "68 % C.L." ) );
		returnable.push_back( make_pair( 0.90, "90 % C.L." ) );
		returnable.push_back( make_pair( 0.95, "95 % C.L." ) );
		returnable.push_back( make_pair( 0.99, "99 % C.L." ) );
	}
	return returnable;
}
*/

vector<pair<double,TString> > Rapid2DLL::GetContour( TString input )
{
	vector<pair<double,TString> > returnable;
	if( input == "2DLL" )
	{
		returnable.push_back( make_pair( 1.15, "68 % CL" ) );
		returnable.push_back( make_pair( 2.36, "90 % CL" ) );
		returnable.push_back( make_pair( 3., "95 % CL" ) );
	}
	if( input == "2DLL-99" )
	{
		returnable.push_back( make_pair( 1.15, "68 % CL" ) );
		returnable.push_back( make_pair( 2.36, "90 % CL" ) );
		returnable.push_back( make_pair( 3., "95 % CL" ) );
		returnable.push_back( make_pair( 4.61, "99 % CL" ) );
	}
	if( input == "FC" )
	{
		returnable.push_back( make_pair( 0.6827, "68 % CL" ) );
		returnable.push_back( make_pair( 0.90, "90 % CL" ) );
		returnable.push_back( make_pair( 0.95, "95 % CL" ) );
	}
	if( input == "FC-99" )
	{
		returnable.push_back( make_pair( 0.6827, "68 % CL" ) );
		returnable.push_back( make_pair( 0.90, "90 % CL" ) );
		returnable.push_back( make_pair( 0.95, "95 % CL" ) );
		returnable.push_back( make_pair( 0.99, "99 % CL" ) );
	}
	return returnable;
}

void Rapid2DLL::get_Plotting_Data( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom* rand, vector<vector<Double_t> >& nll_data, vector<vector<Double_t> >& nll_data_rotated, vector<vector<Double_t> >& coords, vector<vector<Double_t> >& coords_rotated )
{
	vector<vector<Double_t> > nll_data_uniq = TTree_Processing::Plotter_Data( input_tree, Draw_String, Cut_String, rand );

	vector<double> temp( 3, 0. ), temp2( 2, 0. );

	for( vector<vector<double> >::iterator nll_i = nll_data_uniq.begin(); nll_i != nll_data_uniq.end(); ++nll_i )
	{
		temp[0] = (*nll_i)[2];
		temp[1] = (*nll_i)[1];
		temp[2] = (*nll_i)[0];
		nll_data.push_back( temp );
		temp2[0] = (*nll_i)[2];
		temp2[1] = (*nll_i)[1];
		coords.push_back( temp2 );
		temp2[0] = (*nll_i)[1];
		temp2[1] = (*nll_i)[2];
		coords_rotated.push_back( temp2 );
		temp[0] = (*nll_i)[1];
		temp[1] = (*nll_i)[2];
		temp[2] = (*nll_i)[0];
		nll_data_rotated.push_back( temp );
	}

	nll_data = rotate( nll_data );
	nll_data_rotated = rotate( nll_data_rotated );
	coords = rotate( coords );
	coords_rotated = rotate( coords_rotated );
}


int Rapid2DLL::PlotRapidFit2DLL( TString controlled_parameter1, TString controlled_parameter2, TTree* input_tree, TRandom3* rand, vector<string> other_params )
{
	//gStyle->SetPadLeftMargin( (Float_t)0.15 );
	//gStyle->SetTitleOffset((Float_t)0.85,"Y");

	gROOT->UseCurrentStyle();
	gROOT->ForceStyle( true );

	//double CV_val_X = RapidFit_Output_File::GetCV_val( input_tree, controlled_parameter1 );
	//double CV_val_Y = RapidFit_Output_File::GetCV_val( input_tree, controlled_parameter1 );
	double CV_val_NLL = RapidFit_Output_File::GetCV_val( input_tree, "NLL", false );

	TString NLL_val_str; NLL_val_str+=CV_val_NLL;
	TString Draw_String( controlled_parameter1+value_suffix+":"+controlled_parameter2+value_suffix+":(NLL-"+NLL_val_str+")");
	TString Coord_String(controlled_parameter1+value_suffix+":"+controlled_parameter2+value_suffix );
	TString Coord_String_Rotate(controlled_parameter2+value_suffix+":"+controlled_parameter1+value_suffix );

	TString Cut_String( "(Fit_Status == 3)" );

	cout << "Getting Data from File" << endl;


	//	Objects to hold the data for plotting
	//	nll_data         =    dG vs phi vs nll
	//	nll_data_rotated =    phi vs dG ns nll
	//	coords           dG vs phi and acts as a map of plots which have/haven't failed
	vector<vector<Double_t> > nll_data, nll_data_rotated, coords, coords_rotated;

	Rapid2DLL::get_Plotting_Data( input_tree, Draw_String, Cut_String, rand, nll_data, nll_data_rotated, coords, coords_rotated );

	TString timeStr( StringOperations::TimeString() );
	gSystem->mkdir( "RapidFit_2DLL_"+timeStr );
	gSystem->cd( "RapidFit_2DLL_"+timeStr );

	TGraph* coord_graph = Histogram_Processing::Get_TGraph( coords, (TRandom*)rand );
	TGraph* r_coord_graph = Histogram_Processing::Get_TGraph( coords_rotated, (TRandom*)rand );	

	TString TCanvas_Name_c("TCanvas_"); TCanvas_Name_c+=rand->Rndm();
	TCanvas_Name_c = StringOperations::Clean( TCanvas_Name_c );
	TCanvas* coord_c = EdStyle::RapidFitCanvas( TCanvas_Name_c, TCanvas_Name_c );
	coord_graph->Draw("AP");
	coord_c->Update();
	Histogram_Processing::Silent_Print( coord_c, "coords.pdf" );

	TString r_TCanvas_Name_c("TCanvas_"); r_TCanvas_Name_c+=rand->Rndm();
	r_TCanvas_Name_c = StringOperations::Clean( r_TCanvas_Name_c );
	TCanvas* r_coord_c = EdStyle::RapidFitCanvas( r_TCanvas_Name_c, r_TCanvas_Name_c );
	r_coord_graph->Draw("AP");
	r_coord_c->Update();
	Histogram_Processing::Silent_Print( r_coord_c, "coords_rotated.pdf" );

	TGraph2D* nll_graph = Histogram_Processing::Get_TGraph2D( nll_data, rand );
	TGraph2D* nll_graph_rotated = Histogram_Processing::Get_TGraph2D( nll_data_rotated, rand );

	TH2* test_histo = Histogram_Processing::Get_TH2( nll_data, rand );

	test_histo->GetXaxis()->SetTitle( EdStyle::GetParamRootName(controlled_parameter1) + " " + EdStyle::GetParamRootUnit(controlled_parameter1) );
	test_histo->GetYaxis()->SetTitle( EdStyle::GetParamRootName(controlled_parameter2) + " " + EdStyle::GetParamRootUnit(controlled_parameter2) );
	test_histo->GetZaxis()->SetTitle( "DLL" );

	TFile* test_file = new TFile( "test_file.root", "UPDATE" );
	test_histo->Write("", TObject::kOverwrite);
	test_file->Close();

	cout << "Interpolating Data" << endl;

	nll_graph->SetNpx( 100 );
	nll_graph->SetNpy( 100 );
	nll_graph_rotated->SetNpx( 100 );
	nll_graph_rotated->SetNpy( 100 );

	TH2* nll_hist = nll_graph->GetHistogram();
	TH2* nll_hist_rotated = nll_graph_rotated->GetHistogram();

	nll_hist->GetXaxis()->SetTitle( EdStyle::GetParamRootName(controlled_parameter1) + " " + EdStyle::GetParamRootUnit(controlled_parameter1) );
	nll_hist->GetYaxis()->SetTitle( EdStyle::GetParamRootName(controlled_parameter2) + " " + EdStyle::GetParamRootUnit(controlled_parameter2) );
	nll_hist->GetZaxis()->SetTitle( "DLL" );
	nll_hist_rotated->GetXaxis()->SetTitle( EdStyle::GetParamRootName(controlled_parameter2) + " " + EdStyle::GetParamRootUnit(controlled_parameter2) );
	nll_hist_rotated->GetYaxis()->SetTitle( EdStyle::GetParamRootName(controlled_parameter1) + " " + EdStyle::GetParamRootUnit(controlled_parameter1) );
	nll_hist_rotated->GetZaxis()->SetTitle( "DLL" );

	TString root_filename( controlled_parameter1+"_"+controlled_parameter2+".root");
	TString root_filename2( controlled_parameter2+"_"+controlled_parameter1+".root");

	TFile* root_file1 = new TFile( root_filename, "UPDATE" );
	nll_hist->Write("",TObject::kOverwrite);
	root_file1->Close();

	TFile* root_file2 = new TFile( root_filename2, "UPDATE" );
	nll_hist_rotated->Write("",TObject::kOverwrite);
	root_file2->Close();

	cout << "Constructing Contours" << endl;

	vector<double> cont_levels = return_first( Rapid2DLL::GetContour( "2DLL" ) );

	vector<TMultiGraph*> nll_contours = Histogram_Processing::GetContoursFromTH2( nll_hist, cont_levels, rand );
	vector<TMultiGraph*> nll_contours_rotated = Histogram_Processing::GetContoursFromTH2( nll_hist_rotated, cont_levels, rand );

	cout << "Plotting Contours" << endl;

	gStyle->SetLineWidth( (Width_t)Rapid2DLL::GetAxisWidth() );
	gROOT->UseCurrentStyle();
	gROOT->ForceStyle( true );

	unsigned int cont_num=0;

	vector<pair<TMultiGraph*,TString> > named_nll_contours;
	cont_num=0;
	for( vector<TMultiGraph*>::iterator cont_i = nll_contours.begin(); cont_i != nll_contours.end(); ++cont_i, ++cont_num )
	{
		named_nll_contours.push_back( make_pair( *cont_i, return_second( Rapid2DLL::GetContour( "2DLL" ) )[cont_num] ) );
	}
	TString filename( controlled_parameter1 + "_" + controlled_parameter2 + ".pdf" );
	Rapid2DLL::Plot_Contours( controlled_parameter1, controlled_parameter2, nll_hist, named_nll_contours, filename, rand, other_params );
	filename = TString( controlled_parameter1 + "_" + controlled_parameter2 + ".png" );
	Rapid2DLL::Plot_Contours( controlled_parameter1, controlled_parameter2, nll_hist, named_nll_contours, filename, rand, other_params );
	filename = TString( controlled_parameter1 + "_" + controlled_parameter2 + ".C" );
	Rapid2DLL::Plot_Contours( controlled_parameter1, controlled_parameter2, nll_hist, named_nll_contours, filename, rand, other_params );

	vector<pair<TMultiGraph*,TString> > named_nll_contours_rotated;
	cont_num=0;
	for( vector<TMultiGraph*>::iterator cont_i = nll_contours_rotated.begin(); cont_i != nll_contours_rotated.end(); ++cont_i, ++cont_num )
	{
		named_nll_contours_rotated.push_back( make_pair( *cont_i, return_second( Rapid2DLL::GetContour( "2DLL" ) )[cont_num] ) );
	}
	TString filename_rotated( controlled_parameter2 + "_" + controlled_parameter1 + ".pdf" );
	Rapid2DLL::Plot_Contours( controlled_parameter2, controlled_parameter1, nll_hist_rotated, named_nll_contours_rotated, filename_rotated, rand, other_params );
	filename_rotated = TString( controlled_parameter2 + "_" + controlled_parameter1 + ".png" );
	Rapid2DLL::Plot_Contours( controlled_parameter2, controlled_parameter1, nll_hist_rotated, named_nll_contours_rotated, filename_rotated, rand, other_params );
	filename_rotated = TString( controlled_parameter2 + "_" + controlled_parameter1 + ".C" );
	Rapid2DLL::Plot_Contours( controlled_parameter2, controlled_parameter1, nll_hist_rotated, named_nll_contours_rotated, filename_rotated, rand, other_params );

	cout << endl;


	gStyle->SetPadRightMargin( (Float_t)0.15 );
	gROOT->UseCurrentStyle();
	gROOT->ForceStyle( true );

	cout << "Plotting Variation in Nuisence Parameters!" << endl << endl;

	Rapid2DLL::Plot_Free_Parameters( input_tree, controlled_parameter1, controlled_parameter2, rand );

	return 0;
}

void Rapid2DLL::Help()
{
	cout << endl << "Rapid2DLL can Accept the following arguments for producing 2DLL plots" << endl;
	cout << "--isFinal" << "\t\t" << "As in producing 1DLL this decides whether we are adding 'LHCb' or 'LHCb Preliminary' to the Plots" << endl;
	cout << endl;
	cout << "--addPhis" << "\t\t" << "Add the Phi_s Theory" << endl;
	cout << endl;
}

void Rapid2DLL::Plot_Contours( TString controlled_parameter1, TString controlled_parameter2, TH1* nll_hist, vector<pair<TMultiGraph*,TString> > nll_contours,
		TString filename, TRandom* rand, vector<string> other_params )
{
	if( rand == NULL ) rand = gRandom;

	//TPaveText* label = NULL;

	string addLHCb="--isFinal";
	string addSMPhisString="--addPhis";
	bool addSMPhis = StringOperations::VectorContains( &other_params, &addSMPhisString ) != -1;

	/*
	   if( StringOperations::VectorContains( &other_params, &addLHCb ) != -1 )		label = Histogram_Processing::addLHCbLabel( "", true );
	   else										label = Histogram_Processing::addLHCbLabel( "", false );
	   */

	TPaveText* label = EdStyle::LHCbLabel();
	if( StringOperations::VectorContains( &other_params, &addLHCb ) == -1 )
	{
		label->AddText( "" );
		label->AddText( "Preliminary" );
	}

	//TLegend* leg = EdStyle::LHCbLegend();//0.65, 0.65, 0.9, 0.9 );

	TLegend* leg = new TLegend( 0.675, 0.725, 0.9, 0.925 );
	leg->SetFillColor( kWhite );
	leg->SetFillStyle( EdStyle::GetTransparentFillStyle() );
	leg->SetTextSize( 0.050 );//EdStyle::GetLHCbTextSize() );
	leg->SetTextFont( EdStyle::GetLHCbFont() );

	//TString TCanvas_Namea("TCanvas_");TCanvas_Namea+=rand->Rndm();
	//TCanvas* c1a = EdStyle::RapidFitCanvas( TCanvas_Namea, TCanvas_Namea );
	//nll_hist->Draw("AXIS");
	//c1a->Update();

	TString TCanvas_Name("TCanvas_"); TCanvas_Name+=rand->Rndm();
	TCanvas_Name = StringOperations::Clean( TCanvas_Name );
	TCanvas* c1 = EdStyle::RapidFitCanvas( TCanvas_Name, TCanvas_Name );

	nll_hist->Draw("AXIS");

	//nll_contours.back().first->Draw("AXIS");

	//c1->Update();

	//nll_contours.back().first->GetXaxis()->SetRangeUser( nll_hist->GetXaxis()->GetXmin(), nll_hist->GetXaxis()->GetXmax() );
	//nll_contours.back().first->GetYaxis()->SetRangeUser( nll_hist->GetYaxis()->GetXmin(), nll_hist->GetYaxis()->GetXmax() );

	c1->Update();

	unsigned int cont_num=0;
	for( vector<pair<TMultiGraph*,TString> >::iterator cont_i = nll_contours.begin(); cont_i != nll_contours.end(); ++cont_i, ++cont_num )
	{
		cont_i->first->Draw("C");
		TList* this_cont = cont_i->first->GetListOfGraphs();
		for( unsigned int i=0; i< (unsigned)this_cont->GetSize(); ++i )
		{
			TGraph* this_part = (TGraph*) this_cont->At( i );
			this_part->SetLineColor( (Color_t)Rapid2DLL::GetColors( cont_num ) );
			this_part->SetLineStyle( (Style_t)Rapid2DLL::GetStyle( cont_num ) );
			this_part->SetLineWidth( (Width_t)Rapid2DLL::GetFunctionLineWidth() );
		}
		leg->AddEntry( cont_i->first->GetListOfGraphs()->First() , nll_contours[cont_num].second, "L" );
	}
	nll_hist->GetXaxis()->SetTitle( EdStyle::GetParamRootName(controlled_parameter1) + " " + EdStyle::GetParamRootUnit(controlled_parameter1) );
	nll_hist->GetYaxis()->SetTitle( EdStyle::GetParamRootName(controlled_parameter2) + " " + EdStyle::GetParamRootUnit(controlled_parameter2) );
	c1->Update();
	leg->Draw();
	if( label )	label->Draw();
	if( addSMPhis )
	{
		double* px=NULL;
		double* xerr=NULL;
		double* py=NULL;
		double* yerr=NULL;
		if( controlled_parameter1 == "Phi_s" )
		{
			px = new double[1]; px[0] = -0.036;//0.087;
			xerr = new double[1]; xerr[0] = 0.002;//0.021;
			py = new double[1]; py[0] = 0.087;//0.036;
			yerr = new double[1]; yerr[0] = 0.021;//0.002;
		}
		else
		{
			px = new double[1]; px[0] = 0.087;
			xerr = new double[1]; xerr[0] = 0.021;
			py = new double[1]; py[0] = -0.036;
			yerr = new double[1]; yerr[0] = 0.002;
		}

		if( ( controlled_parameter1 == "Phi_s" || controlled_parameter1 == "deltaGamma" )&&( controlled_parameter2 == "Phi_s" || controlled_parameter2 == "deltaGamma" ) )
		{
			TGraphErrors* CV = new TGraphErrors( 1, px, py, xerr, yerr );
			CV->SetName( "Standard Model" );
			leg->AddEntry( CV, "Standard Model", "PE" );
			CV->Draw("LP SAME");
		}
	}
	c1->Update();
	Histogram_Processing::Silent_Print( c1, filename );
	return;
}

void Rapid2DLL::Plot_Free_Parameters( TTree* input_tree, TString controlled_parameter1, TString controlled_parameter2, TRandom3* rand )
{
	vector<string> controlled_parameter_name;
	controlled_parameter_name.push_back( controlled_parameter1.Data() );
	controlled_parameter_name.push_back( controlled_parameter2.Data() );

	vector<TString> free_param_names = RapidFit_Output_File::get_free_non_scanned_parameters( input_tree, controlled_parameter_name );

	for( vector<TString>::iterator param_i = free_param_names.begin(); param_i != free_param_names.end(); ++param_i )
	{
		cout << "Plotting: " << *param_i << endl;
		TString free_value(*param_i); free_value.Append( value_suffix );

		TString Draw_String( controlled_parameter1+value_suffix+":"+controlled_parameter2+value_suffix+":"+free_value );
		TString Cut_String("Fit_Status==3");

		vector<vector<double> > param_data_raw = TTree_Processing::Plotter_Data( input_tree, Draw_String, Cut_String, rand );
		vector<double> temp(3, 0. ); vector<vector<double> > param_data;

		for( vector<vector<double> >::iterator val_i = param_data_raw.begin(); val_i != param_data_raw.end(); ++val_i )
		{
			temp[0] = (*val_i)[2];
			temp[1] = (*val_i)[1];
			temp[2] = (*val_i)[0];
			param_data.push_back( temp );
		}

		param_data = rotate(param_data);

		TGraph2D* param_graph = Histogram_Processing::Get_TGraph2D( param_data, rand );
		TH1* param_hist = param_graph->GetHistogram();

		param_hist->SetContour( 40 );
		TString canvas_name("TCanvas_"); canvas_name+=rand->Rndm();
		canvas_name = StringOperations::Clean( canvas_name );
		TCanvas* c1 = EdStyle::RapidFitCanvas( canvas_name, canvas_name );

		c1->SetTitle( "Variation in " + *param_i );

		param_hist->Draw("cont1z");
		c1->Update();
		TString filename( *param_i+"_cont1z.pdf" );
		c1->Print( filename );
		filename = TString( *param_i+"_cont1z.png" );
		c1->Print( filename );
		filename = TString( *param_i+"_cont1z.C" );
		c1->Print( filename );
		param_hist->Draw("colz");
		c1->Update();
		c1->Print( *param_i+"_colz.pdf" );
		c1->Print( *param_i+"_colz.png" );
		c1->Print( *param_i+"_colz.C" );
		c1->Close();
	}

	return;
}



