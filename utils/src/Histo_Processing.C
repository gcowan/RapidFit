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
#include "TRandom.h"
#include "TROOT.h"
#include "TPaletteAxis.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
//	RapidFit Header
#include "EdStyle.h"
//	RapidFit Utils Headers
#include "Histo_Processing.h"
#include "StringOperations.h"
#include "TTree_Processing.h"
#include "Mathematics.h"
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

//  This has been adapted from the original code in RapidFits Statistics code
//  It is intended to take a histogram and automatically rebin according to this function
//  As such it uses as much inbuilt functionality in root as possible
//Return the ideal number of bins for a histogram of a vector of doubles
////Uses D. Scott's method, published 1979
int Histogram_Processing::OptimumBinNumber( TH1* input_hist, int axis )
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
unsigned int Histogram_Processing::GetOptimalBins( TH1* input_hist, int axis )
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
void Histogram_Processing::OptimallyRebin( TH1* input_hist, int axis )
{
	int temp = OptimumBinNumber( input_hist, axis );
	(void) temp;	return;
}


TString Histogram_Processing::Best_Fit_Function( TH1* local_histogram, int OutputLevel )
{
	//	the gamma function in ROOT has issues with verbosity, let's silence it and then return the verbosity back at the end

	streambuf *nullbuf=NULL, *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;
	ofstream filestr;
	filestr.open ("/dev/null");

	int thisErr = gErrorIgnoreLevel;

	//      If the user wanted silence we point the Std Output Streams to /dev/null
	if( OutputLevel <= -1 )
	{
		cout_bak = cout.rdbuf();
		cerr_bak = cerr.rdbuf();
		clog_bak = clog.rdbuf();
		nullbuf = filestr.rdbuf();
		//freopen("/dev/null", "w", stderr);
		cout.rdbuf(nullbuf);
		cerr.rdbuf(nullbuf);
		clog.rdbuf(nullbuf);
		gErrorIgnoreLevel = kFatal;
	}


	TString Fit_Options( "Q" );			// Reduce the verbosity slamming the user during this internal procedure
	TF1* my_landau = Mathematics::landau_func(); (void) my_landau;
	TF1* my_gamma = Mathematics::gamma_func(); (void) my_gamma; // this is used from the global state not here... this displeases me...

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

	gErrorIgnoreLevel = thisErr;

	return fit_type;
}

void Histogram_Processing::Silent_Fit( TH1* input_histo, TString fit_type, int OutputLevel )
{
	//      the gamma function in ROOT has issues with verbosity, let's silence it and then return the verbosity back at the end

	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL, *nullbuf=NULL;
	ofstream filestr;
	filestr.open ("/dev/null");

	int thisErr = gErrorIgnoreLevel;
	//      If the user wanted silence we point the Std Output Streams to /dev/null
	if( OutputLevel <= -1 )
	{
		cout_bak = cout.rdbuf();
		cerr_bak = cerr.rdbuf();
		clog_bak = clog.rdbuf();
		//  Redirect the errors to the empty void of nullness
		nullbuf = filestr.rdbuf();
		//freopen("/dev/null", "w", stderr);
		cout.rdbuf(nullbuf);
		cerr.rdbuf(nullbuf);
		clog.rdbuf(nullbuf);
		gErrorIgnoreLevel = kFatal;
	}
	//  Redirect the errors to the empty void of nullness
	//freopen("/dev/null", "w", stderr);
	input_histo->Fit( fit_type, "Q" );
	//      Reset Std Output Streams
	if( OutputLevel <= -1 )
	{
		cout.rdbuf(cout_bak);
		cerr.rdbuf(cerr_bak);
		clog.rdbuf(clog_bak);
		gErrorIgnoreLevel = thisErr;
	}
}

void Histogram_Processing::Silent_Draw( TCanvas* c1, TH1* input_histo, TString options, int OutputLevel )
{
	//      the gamma function in ROOT has issues with verbosity, let's silence it and then return the verbosity back at the end

	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;

	int thisErr = gErrorIgnoreLevel;
	//      If the user wanted silence we point the Std Output Streams to /dev/null
	if( OutputLevel <= -1 )
	{
		cout_bak = cout.rdbuf();
		cerr_bak = cerr.rdbuf();
		clog_bak = clog.rdbuf();
		//  Redirect the errors to the empty void of nullness
		//freopen("/dev/null", "w", stderr);
		cout.rdbuf(0);
		cerr.rdbuf(0);
		clog.rdbuf(0);
		gErrorIgnoreLevel = kFatal;
	}
	input_histo->Draw(options);
	c1->Update();
	//      Reset Std Output Streams
	if( OutputLevel <= -1 )
	{
		cout.rdbuf(cout_bak);
		cerr.rdbuf(cerr_bak);
		clog.rdbuf(clog_bak);
		thisErr = gErrorIgnoreLevel;
	}
}

void Histogram_Processing::Silent_Print( TCanvas* c1, TString Print_String, int OutputLevel )
{
	//      the gamma function in ROOT has issues with verbosity, let's silence it and then return the verbosity back at the end

	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;

	int thisErr = gErrorIgnoreLevel;
	//      If the user wanted silence we point the Std Output Streams to /dev/null
	if( OutputLevel <= -1 )
	{
		cout_bak = cout.rdbuf();
		cerr_bak = cerr.rdbuf();
		clog_bak = clog.rdbuf();
		//  Redirect the errors to the empty void of nullness
		//freopen("/dev/null", "w", stderr);
		gErrorIgnoreLevel = kFatal;
	}
	c1->Update();
	c1->Print( Print_String );
	//      Reset Std Output Streams
	if( OutputLevel <= -1 )
	{
		cout.rdbuf(cout_bak);
		cerr.rdbuf(cerr_bak);
		clog.rdbuf(clog_bak);
		gErrorIgnoreLevel = thisErr;
	}
}

TH1* Histogram_Processing::Get_Histo( TTree* input_tree, TString draw_str, TString weight_str, TRandom* rand )
{
	if( rand == NULL ) rand = gRandom;
	TString rand_str; rand_str+=rand->Rndm();
	string rand_num( rand_str.Data() );
	replace( rand_num.begin(), rand_num.end(), '.', '_' );
	rand_str = rand_num.c_str();
	input_tree->Draw( draw_str, weight_str, "goff" );
	TH1* output_histo = (TH1*)((TH1*)(input_tree->GetHistogram())->Clone(output_histo->GetName()+rand_str));
	return output_histo;
}

TGraph* Histogram_Processing::Get_Graph( TTree* input_tree, TString draw_str, TString weight_str, TRandom* rand )
{
	if( rand == NULL ) rand = gRandom;
	TString rand_str; rand_str+=rand->Rndm();
        string rand_num( rand_str.Data() );
        replace( rand_num.begin(), rand_num.end(), '.', '_' );
        rand_str = rand_num.c_str();
	TGraph* new_graph = new TGraph( Get_Histo( input_tree, draw_str, weight_str, rand ) );
	new_graph->SetName("Graph_"+rand_str);
	return new_graph;
}

TH1* Histogram_Processing::Get_TH1( const vector<Double_t>& input, TRandom* rand, int bins, double X_MIN, double X_MAX )
{
	if( rand == NULL ) rand = gRandom;
	if( input.empty() )
	{
		cout << "NO 1D INPUT DATA TO PLOT" << endl;
		return NULL;
	}
	double input_min = get_minimum( input );
	double input_max = get_maximum( input );
	if( X_MIN > -DBL_MAX ) input_min = X_MIN;
	if( X_MAX < DBL_MAX ) input_max = X_MAX;
	TString rand_num; rand_num += rand->Rndm();
        string rand_str( rand_num.Data() );
        replace( rand_str.begin(), rand_str.end(), '.', '_' );
        rand_num = rand_str.c_str();
	TH1* returnable_hist = new TH1D( "TH1_"+rand_num, "TH1_"+rand_num, bins, input_min, input_max );
	vector<Double_t>::const_iterator index_i = input.begin();
	vector<Double_t>::const_iterator index_e = input.end();
	for( ; index_i != index_e; ++index_i )
	{
		returnable_hist->Fill( *index_i );
	}
	return returnable_hist;
}

TH1* Histogram_Processing::Plot_1D( const vector<double>& input, TString Filename, TString Options, TRandom* rand, int bins, double min, double max )
{
	if( rand == NULL ) rand = gRandom;
	if( input.empty() ) return NULL;

	TH1* temp_histo = Histogram_Processing::Get_TH1( input, rand, bins, min, max );

	TString canv_name("canv_"); canv_name+=rand->Rndm();
        string rand_num( canv_name.Data() );
        replace( rand_num.begin(), rand_num.end(), '.', '_' );
        canv_name = rand_num.c_str();

	TCanvas* c1 = EdStyle::RapidFitCanvas( canv_name, canv_name );

	temp_histo->Draw( Options );

	c1->Update();
	c1->Print(Filename);
	c1->Close();
	delete c1;

	return temp_histo;
}

TGraph* Histogram_Processing::Get_TGraph( const vector<vector<Double_t> >& input, TRandom* rand )
{
	if( rand == NULL ) rand = gRandom;
	if( input.empty() ) return NULL;

	TGraph* returnable_graph = new TGraph( input[0].size(), &(input[0][0]), &(input[1][0]) );

	TString TGraph_Name("TGraph_"); TGraph_Name+=rand->Rndm();
        string rand_num( TGraph_Name.Data() );
        replace( rand_num.begin(), rand_num.end(), '.', '_' );
        TGraph_Name = rand_num.c_str();

	returnable_graph->SetTitle( TGraph_Name );
	returnable_graph->SetName( TGraph_Name );

	return returnable_graph;
}

TGraph2D* Histogram_Processing::Get_TGraph2D( const vector<vector<Double_t> >& input, TRandom* rand )
{
	if( rand == NULL ) rand = gRandom;
	if( input.empty() ) return NULL;
	if( input.size() != 3 ) return NULL;

	TGraph2D* returnable_TGraph2D = new TGraph2D( input[0].size(), const_cast<Double_t*>(&(input[0][0])), const_cast<Double_t*>(&(input[1][0])), const_cast<Double_t*>(&(input[2][0])) );

	TString TGraph2D_Name("TGraph2D_"); TGraph2D_Name+=rand->Rndm();
        string rand_num( TGraph2D_Name.Data() );
        replace( rand_num.begin(), rand_num.end(), '.', '_' );
        TGraph2D_Name = rand_num.c_str();

	returnable_TGraph2D->SetName( TGraph2D_Name );
	returnable_TGraph2D->SetTitle( TGraph2D_Name );

	return returnable_TGraph2D;

}

TH2* Histogram_Processing::Get_TH2( const vector<vector<Double_t> >& input, TRandom* rand, int bins1, int bins2 )
{
	if( rand == NULL ) rand = gRandom;
	if( input.empty() )
	{
		cout << "NO 2D INPUT DATA TO PLOT" << endl;
		return NULL;
	}
	TString rand_num; rand_num += rand->Rndm();
        string rand_str( rand_num.Data() );
        replace( rand_str.begin(), rand_str.end(), '.', '_' );
        rand_num = rand_str.c_str();

	vector<Double_t> X_data = input[0];
	vector<Double_t> Y_data = input[1];
	vector<Double_t> weight;
	if( input.size() == 3 )
	{
		weight = input[2];
	}
	else
	{
		weight = vector<double>( input[0].size(), 1. );
	}
	TH2* returnable_hist = new TH2D( "TH2_"+rand_num, "TH2_"+rand_num, bins1, get_minimum(X_data), get_maximum(X_data), bins2, get_minimum(Y_data), get_maximum(Y_data) );
	vector<Double_t>::const_iterator index_x = X_data.begin();
	vector<Double_t>::const_iterator index_y = Y_data.begin();
	vector<Double_t>::const_iterator index_e = X_data.end();
	vector<Double_t>::const_iterator index_w = weight.begin();
	for( ; index_x != index_e; ++index_x, ++index_y, ++index_w )
	{
		returnable_hist->Fill( *index_x, *index_y, *index_w );
	}
	return returnable_hist;
}

TH2* Histogram_Processing::Plot_2D( const vector<double>& X, const vector<double>& Y, TString Filename, TString Option, TRandom* rand )
{
	if( rand == NULL ) rand = gRandom;
	if( X.empty() || Y.empty() ) return NULL;
	if( X.size() != Y.size() ) return NULL;

	vector<vector<double> > data; data.push_back( X ); data.push_back( Y );

	TH2* temp_plot = Histogram_Processing::Get_TH2( data, rand );

	TString canv_name("Canv_"); canv_name += rand->Rndm();
        string rand_num( canv_name.Data() );
        replace( rand_num.begin(), rand_num.end(), '.', '_' );
        canv_name = rand_num.c_str();

	TCanvas* c1 = EdStyle::RapidFitCanvas( canv_name, canv_name );

	temp_plot->Draw( Option );

	c1->Update();
	c1->Print( Filename );

	c1->Close();

	return temp_plot;
}

TH3* Histogram_Processing::Get_TH3( const vector<vector<Double_t> >& input, TRandom* rand, int bins1, int bins2, int bins3 )
{
	if( rand == NULL ) rand = gRandom;
	if( input.empty() )
	{
		cout << "NO 3D INPUT DATA TO PLOT" << endl;
		return NULL;
	}
	TString rand_num; rand_num += rand->Rndm();
        string rand_str( rand_num.Data() );
        replace( rand_str.begin(), rand_str.end(), '.', '_' );
        rand_num = rand_str.c_str();

	vector<Double_t> X_data = input[0];
	vector<Double_t> Y_data = input[1];
	vector<Double_t> Z_data = input[2];
	TH3* returnable_hist = new TH3D( "TH3_"+rand_num, "TH3_"+rand_num, bins1, get_minimum(X_data), get_maximum(X_data), bins2, get_minimum(Y_data), get_maximum(Y_data), bins3, get_minimum(Z_data), get_maximum(Z_data) );
	vector<Double_t>::const_iterator index_x = X_data.begin();
	vector<Double_t>::const_iterator index_y = Y_data.begin();
	vector<Double_t>::const_iterator index_z = Z_data.begin();
	vector<Double_t>::const_iterator index_e = X_data.end();
	for( ; index_x != index_e; ++index_x, ++index_y, ++index_z )
	{
		returnable_hist->Fill( *index_x, *index_y, *index_z );
	}
	return returnable_hist;
}

TH1* Histogram_Processing::Get_Histo( const vector<vector<Double_t> >& input, TRandom* rand, int bins1, int bins2, int bins3 )
{
	if( rand == NULL ) rand = gRandom;
	if( input.empty() ) return NULL;
	if( input.size() == 1 ) return Histogram_Processing::Get_TH1( input[0], rand, bins1 );
	if( input.size() == 2 ) return (TH1*)Histogram_Processing::Get_TH2( input, rand, bins1, bins2 );
	if( input.size() == 3 ) return (TH1*)Histogram_Processing::Get_TH3( input, rand, bins1, bins2, bins3 );
	return NULL;
}

TH3* Histogram_Processing::Plot_3D( const vector<double>& X, const vector<double>& Y, const vector<double>& Z, TString Filename, TString Option, TRandom* rand )
{
	if( rand == NULL ) rand = gRandom;
	if( X.empty() || (Y.empty() || Z.empty() ) )	return NULL;
	if( (X.size() != Y.size()) || (Y.size() != Z.size()) || (X.size() != Z.size()) ) return NULL;

	vector<vector<double> > data; data.push_back( X ); data.push_back( Y ); data.push_back( Z );

	TH3* temp_plot = Histogram_Processing::Get_TH3( data, rand );

	TString canv_name("Canv_"); canv_name+=rand->Rndm();
        string rand_str( canv_name.Data() );
        replace( rand_str.begin(), rand_str.end(), '.', '_' );
        canv_name = rand_str.c_str();

	TCanvas* c1 = EdStyle::RapidFitCanvas( canv_name, canv_name );

	temp_plot->Draw( Option );

	c1->Update();
	c1->Print(Filename);
	c1->Close();

	delete c1;

	return temp_plot;
}

TH1* Histogram_Processing::Plot( const vector<vector<double> >& input, TString Filename, TString Option, TRandom* rand )
{
	if( rand = NULL ) rand = gRandom;
	if( input.empty() ) return NULL;
	if( input.size() == 1 ) return Histogram_Processing::Plot_1D( input[0], Filename, Option, rand );
	if( input.size() == 2 ) return (TH1*) Histogram_Processing::Plot_2D( input[0], input[1], Filename, Option, rand );
	if( input.size() == 3 ) return (TH1*) Histogram_Processing::Plot_3D( input[0], input[1], input[2], Filename, Option, rand );
	return NULL;
}

//	Originally Written by Conor Fitzpatrick
TPaveText* Histogram_Processing::addLHCbLabel(TString footer, bool final/*, bool DATA*/)
{
	//TPaveText * label = new TPaveText(0.18, 0.77, 0.18, 0.88,"BRNDC");
	TPaveText* label = new TPaveText( gStyle->GetPadLeftMargin() + 0.08, 
			0.87 - gStyle->GetPadTopMargin(),
			gStyle->GetPadLeftMargin() + 0.30,
			0.95 - gStyle->GetPadTopMargin(), "BRNDC");

			label->SetFillStyle(0);         //Transparent i.e. Opacity of 0 :D
			label->SetBorderSize(0);
			label->SetTextAlign(11);
			label->SetTextSize(Float_t(0.07));
			TString labelstring( "LHCb" );//"#splitline{LHCb}{#scale[1.0]{Preliminary}}" );
			//if( DATA ) labeltstring.Append( " Data" );
			//if( !DATA ) labeltstring.Append( " Simulation" );

			label->AddText( labelstring );
			if( !final )	label->AddText( "#scale[1.0]{Preliminary}" );

			//label->AddText("#sqrt{s} = 7TeV " + footer );
			return label;
}

vector<TMultiGraph*> Histogram_Processing::GetContoursFromTH2( TH2* input_th2, const vector<double>& contour_list, TRandom* rand )
{
	if( rand == NULL ) rand = gRandom;
	vector<TMultiGraph*> returnable_Contours;

	TString TCanvas_Name("TCanvas_");TCanvas_Name+=rand->Rndm();
        string rand_str_( TCanvas_Name.Data() );
        replace( rand_str_.begin(), rand_str_.end(), '.', '_' );
        TCanvas_Name = rand_str_.c_str();
	TCanvas* c1 = EdStyle::RapidFitCanvas( TCanvas_Name, TCanvas_Name );

	input_th2->SetContour( contour_list.size(), &(contour_list[0]) );

	input_th2->Draw("cont LIST");
	c1->Update();

	TObjArray* generated_contours = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");

	if( generated_contours == NULL ) return returnable_Contours;

	TString rand_str; rand_str+=rand->Rndm();
	string rand_str2( rand_str.Data() );
	replace( rand_str2.begin(), rand_str2.end(), '.', '_' );
	rand_str = rand_str2.c_str();

	for( unsigned int i=0; i< generated_contours->GetSize(); ++i )
	{
		TString Contour_Name( "Contour_"+rand_str+"_" ); Contour_Name+=i;

		TMultiGraph* this_contour = new TMultiGraph( Contour_Name, Contour_Name );

		TList* contour_parts = (TList*) generated_contours->At(i);

		for( unsigned int j=0; j< contour_parts->GetSize(); ++j )
		{
			TString Part_Name( "Contour_"+rand_str+"_" ); Part_Name+=i;
			Part_Name.Append("_"); Part_Name+=j;
			TGraph* this_part = (TGraph*) contour_parts->At(j);

			TGraph* copy_this_part = new TGraph( *this_part );

			copy_this_part->SetName(Part_Name);
			copy_this_part->SetTitle(Part_Name);

			this_contour->Add( copy_this_part );
		}

		returnable_Contours.push_back( this_contour );
	}

	c1->Close();

	return returnable_Contours;
}

