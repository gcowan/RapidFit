/* This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it
 * and/or modify it under the terms of the Do What The Fuck You Want
 * To Public License, Version 2, as published by Sam Hocevar. See
 * http://sam.zoy.org/wtfpl/COPYING for more details. */

//  Disclaimer: I have copied this from my personal code repo and have had to copy a LOT of routines from a 'library' I am working on,
//  		as such the version I include in RapidFit will likely never get maintained, however I have written this particular program to be
//  		file structure agnostic so that it should 'just work'(tm) with whatever toy study results you throw at it
//
//  		If you want a more upto date version of this particular script
//  		(and likely one which is slightly faster as I plan to opimize several file access routines soon)
//  		just drop me an email rcurrie@cern.ch and I will send you my upto date copy, which does require you have cmake installed and running

//  Root Headers
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
//	Easier with TFitResultPtr but that is 'too new' for some users...
#include "TFitResult.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TMath.h"
#include "TROOT.h"
#include "TColor.h"
#include "TStyle.h"
#include "TKey.h"
#include "TNtuple.h"
//	RapidFit Utils Headers
#include "TString_Processing.h"
#include "Histo_Processing.h"
#include "NTuple_Processing.h"
//  System Headers
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <cstdlib>

#define TAB_LENGTH 8

using namespace std;

short int mystyle_line_width(){  return 1;  }

//  Return the string of the filename
string True_Filename_str( string FilePath )
{
	size_t found;  //  Remove file extention
	found = FilePath.find_last_of(".");
	string temp_string = FilePath.substr(0,found);
	//  Remove relative path of the file
	found = temp_string.find_last_of("/\\");
	return temp_string.substr(found+1);
}

//  Return the filename of the file in the FilePath
TString True_Filename( string FilePath )
{
	return TString( True_Filename_str( FilePath ) );
}

int tab_formatter(TString input, int tab_number)
{
	string input_str(input.Data());
	//  Want tab_number of tabs after this string, i.e. how many to align words of long length
	return int( tab_number - floor(input_str.size()/TAB_LENGTH) );
}

//  Return False if the string new_name isn't in the vector of stings
//  found_names
bool check_found( string *new_name, vector<string> *found_names )
{
	for( unsigned int i=0; i < found_names->size(); i++)
	{
		if( new_name->compare((*found_names)[i]) == 0 ) return true;
	}
	return false;
}

bool check_found( string *new_name, vector<TString> *found_names )
{
	vector<string> temp_strings;
	for(unsigned int index=0; index < found_names->size(); index++ )
	{
		temp_strings.push_back((*found_names)[index].Data());
	}
	return check_found( new_name, &temp_strings);
}


//  From the root current path look for all keys (objects in root file) and loop over them
//  For each one that is actually an ntuple store it's name and the number of events it has
//  Slightly stupid and ineffient as it always has to decend into the file stuture
//  This shouldn't be a huge problem as I never anticipate finding files with hundereds of files
void get_number_events(int *number_o_events, TString current_path, TString top_dir, TString *relative_path, vector<string> *found_names)
{
	// goto current search path
	gDirectory->cd(current_path);
	TDirectory *current_path_directory = gDirectory;
	TString current_path_str(current_path);
	// get all keys present in current path
	TIter nextkey( current_path_directory->GetListOfKeys() );
	TKey *key, *oldkey=0;

	//cout << current_path << endl;
	while ( (key = (TKey*)nextkey()) )
	{
		//  Loop over all keys (objects) in the current path
		//  Leave on the last key
		if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
		//  Should already be here?
		//  Read in the object
		TObject *obj = key->ReadObj();

		// Get the name of the object in current directory
		string obj_name((char*) obj->GetName());
		TString temp(obj_name);

		//  If this is a Tuple object we want to store some info about it
		if ( obj->IsA()->InheritsFrom( TNtuple::Class()) || obj->IsA()->InheritsFrom( TTree::Class()) )
		{
			//  Store the full name of the object and it's full path to determine if we've already found it
			string full_name(relative_path->Data());
			full_name.append(obj_name);
			// If this is an ntuple we foud ignore it
			if( !check_found(&full_name,found_names) )
			{
				//  If this is a new ntuple save it's name and event_number
				TTree* found = (TTree*)gDirectory->Get( temp );
				*number_o_events = int(found->GetEntries());
				relative_path->Append(obj_name);
				found_names->push_back((*relative_path).Data());
				break;
			}
		} else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {

			//  If this is a directory we want to decend down into it
			current_path_str.Append(obj_name);
			current_path_str.Append("/");
			relative_path->Append(obj_name);
			relative_path->Append("/");
			get_number_events(number_o_events, current_path_str, top_dir, relative_path, found_names);

		}
	}
}

//  A Wrapper for the get_number_events script to be a slight tad more user friendly
//  Input	File name
//  Output	a vector containing the events in each ntuple
//		a vector containing the location of each ntuple within the root file
void get_tuple_infos(TString File_Name, vector<int>* number_o_events, vector<TString>* relative_path)
{

	TFile *first_source = new TFile( File_Name, "READ" );
	//  Should probably put an abort here if file does not exist,
	//  Will copy this from another chunk of code when I have time

	TString top_dir( gDirectory->GetPath() );
	TString current_directory( gDirectory->GetPath() );

	int zero=-1;
	number_o_events->push_back(zero);
	vector<string> found_names;
	//  Holder for a new name we expect to find
	found_names.push_back("");
	TString newstr;
	//  Holder for the path of this object within the file
	relative_path->push_back(newstr);
	//  Find the first ntuple and store it's number of events and address
	get_number_events( &(*number_o_events)[0], current_directory, top_dir, &(*relative_path)[0], &found_names);

	//  Keep searching for new TTree objects until the get_number_events returns -1 events
	//  indicating that it's finished
	for( unsigned int i=1; (*number_o_events)[i-1]!=-1; i++ )
	{
		int temp=-1;
		number_o_events->push_back(temp);	//Holder for number of events
		TString newstr;
		relative_path->push_back(newstr);	//Holder for path to new ntuple
		get_number_events( &(*number_o_events)[i], current_directory, top_dir, &(*relative_path)[i], &found_names);
		//cout << "called" << i << " " << number_o_events[i] << "events found " << endl;
	}

	//  If No ntuple found, remove the holders
	if( (*number_o_events)[number_o_events->size()-1] == -1 )
	{
		number_o_events->pop_back();
		relative_path->pop_back();
	}
	first_source->Close();
}

void get_tuple_infos(TString File_Name, vector<TString>* relative_path)
{
	vector<int> tuple_numbers;
	get_tuple_infos( File_Name, &tuple_numbers, relative_path);
}

vector<TString> get_tuple_infos(TString File_Name)
{
	vector<int> tuple_numbers;
	vector<TString> relative_path;
	get_tuple_infos( File_Name, &tuple_numbers, &relative_path);
	return relative_path;
}

//  Get the name of the TTree object from it's path within a file
void get_tuple_name( TString relative_path, TString *tuple_name )
{
	size_t found;
	string temp(relative_path.Data());
	found=temp.find_last_of("/\\");
	temp.substr(0,found);
	tuple_name->Append(temp.substr(found+1));
}

TString get_tuple_name( TString relative_path )
{
	size_t found;
	string temp(relative_path.Data());
	found=temp.find_last_of("/\\");
	temp.substr(0,found);
	return TString(temp.substr(found+1));
}

//  Return the absolute path
void get_absolute_path( TString relative_path, TString file_name, TString *absolute_path)
{
	TFile* temp = new TFile( file_name, "READ" );
	absolute_path->Append( gDirectory->GetPath() );
	absolute_path->Append( relative_path );
	temp->Close();
}

TH1F* get_branch( TFile *Input_File, TString Relative_Tuple_Name, TString Branch_Name )
{
	Input_File->cd();
	TTree* temp_tree = (TTree*) gDirectory->Get( Relative_Tuple_Name );
	temp_tree->Draw( Branch_Name );
	return (TH1F*) temp_tree->GetHistogram()->Clone();
}

void tabs_TStr( TString input, int tab_number, TString* output)
{
	for(int tabs=0; tabs < tab_formatter(input, tab_number); tabs++)
	{
		output->Append("\t");
	}
}

void set_style()
{
	TStyle* mystyle = new TStyle("newstyle","newstyle");
	mystyle->SetPalette(1,0);
	mystyle->SetCanvasColor(10);
	mystyle->SetTitleColor(4);
	mystyle->SetStatColor(10);
	mystyle->SetTitleColor(1);
	mystyle->SetPadColor(10);
	mystyle->SetLineColor(1);
	mystyle->SetLineWidth( mystyle_line_width() );
	mystyle->SetFuncWidth( mystyle_line_width() );
	mystyle->SetStatW(float(0.2));
	mystyle->SetStatH(float(0.2));
	mystyle->SetLabelColor(1,"xy");
	//mystyle->SetOptStat(11);
	mystyle->SetOptStat(0);
	//  gStyle->SetOptStat("neMR");
	mystyle->SetOptFit(111);
	mystyle->SetFillColor(0);
	mystyle->SetPadColor(10);
	mystyle->SetTitleOffset(float(1.3),"xyz");
	mystyle->SetEndErrorSize(8);
	mystyle->SetFrameFillColor(10);
	mystyle->SetTitleFillColor(10);
	mystyle->SetTitleFontSize(float(0.07));
	mystyle->SetTitleBorderSize(1);
	mystyle->SetTitleX(0.2f);
	mystyle->cd();
	mystyle->SetTextSize(5);
	//mystyle->SetStatFontSize(22);
	mystyle->SetStatX(float(0.9));
	mystyle->SetStatY(float(0.9));
	gROOT->ForceStyle();

}

//  Some external constants to use in code

// 5 colors x 5 styles should be enough for most complex plots, if not your doing it wrong

vector<short int> get_line_colors( ){
	vector<short int> my_line_colors;
	my_line_colors.push_back(4);
	my_line_colors.push_back(2);
	my_line_colors.push_back(6);
	my_line_colors.push_back(7);
	my_line_colors.push_back(8);
	return my_line_colors;
}

vector<short int> get_line_styles( ){
	vector<short int> my_line_styles;
	my_line_styles.push_back(1);
	my_line_styles.push_back(2);
	my_line_styles.push_back(3);
	my_line_styles.push_back(9);
	my_line_styles.push_back(10);
	return my_line_styles;
}

int main ( int argc, char* argv[] )
{
	//  Redirect the errors to the empty void of nullness
	freopen("/dev/null", "w", stderr);

	cout << "\n\nTINTER Is Not a TuplE Reader:\n\n";

	cout<<"\t\t'########:'####:'##::: ##:'########:'########:'########::"<<endl;
	cout<<"\t\t... ##..::. ##:: ###:: ##:... ##..:: ##.....:: ##.... ##:"<<endl;
	cout<<"\t\t::: ##::::: ##:: ####: ##:::: ##:::: ##::::::: ##:::: ##:"<<endl;
	cout<<"\t\t::: ##::::: ##:: ## ## ##:::: ##:::: ######::: ########::"<<endl;
	cout<<"\t\t::: ##::::: ##:: ##. ####:::: ##:::: ##...:::: ##.. ##:::"<<endl;
	cout<<"\t\t::: ##::::: ##:: ##:. ###:::: ##:::: ##::::::: ##::. ##::"<<endl;
	cout<<"\t\t::: ##::::'####: ##::. ##:::: ##:::: ########: ##:::. ##:"<<endl;
	cout<<"\t\t:::..:::::....::..::::..:::::..:::::........::..:::::..::"<<endl;

	cout << "\n\n\nUsage: " << argv[0] << "\tfilename.root\t"<< "tuple_name(optional)\t" << "branch_name(requires a tuple_name)\t" << "fit_function(requires a branch)"<< endl << endl;

	cout << "eg:   tinter pullPlots.root gamma value gaus \t\t\t gives you a plot for gamma values fitted with a Gaussian\n\n";
	//cout << "\tThese plot(s) has been automatically rebinned, fitted correctly and nicely presented in a new directory: ./pullPlots_output :D " <<endl;
	//cout << "\n\n\tThis has been 'optimised' to be able to generate the toy study plots in pdf format for Note2 in Moriond." << endl;

	cout << "\n\n\t\t\tEnjoy :D\tI will now give you some verbose output about the file: "<<  argv[1] << "\n\n" << endl;

	//  Set the output style before I forget
	set_style();

	//  Take the first argument open that file
	TString Input_Filename ( argv[1] ); 
	TFile *input_pull = new TFile ( Input_Filename,"READ" );

	TString outputdir = True_Filename( string( argv[1] ) );
	outputdir.Append("_output");
	gSystem->mkdir( outputdir );

	vector<short int> line_colors = get_line_colors();

	gStyle->SetFuncColor(line_colors[1]);
	gROOT->ForceStyle();

	//  Lookup the Existing ntuples for reference
	vector<TString> Param_tree_path = get_tuple_infos ( Input_Filename );

	//  Inspect the user input for any wanted ntuples

	if ( ( argc >= 3 )  && ( string ( argv[2] ).compare ( "*" ) != 0 ) )
	{
		string pass(argv[2]);
		if ( check_found ( &pass, &Param_tree_path ) )
		{
			while ( !Param_tree_path.empty() )
			{
				Param_tree_path.pop_back();
			}
			Param_tree_path.push_back ( argv[2] );
		}
		else
		{
			cout << "Tuple with name: " << argv[2] << " not found in file " << argv[1] << " exiting." << endl;
			return -1;
		}
	}

	//  Holder for the absolute path of objects in root
	vector<TString> Absolute_Param_Paths ( Param_tree_path.size() );

	//  2D holder for all of the plots should we wish to keep them in scope for use later
	vector<vector<TH1F*> > Param_FitResult_Plots;
	//  2D holder for the names of each of the plots based on what is being plotted
	vector<vector<TString> > Param_FitResult_Names;

	vector<vector<bool> > left_handed_plot;

	//  Loop over all of the trees
	for ( unsigned int param_num=0; param_num < Param_tree_path.size(); param_num++ )
	{
		//  Increment the holder for the hisrogram pointers
		vector<TH1F*> temp;
		Param_FitResult_Plots.push_back ( temp );

		//  Get the name of the tree in this path and open it
		TString tree_name = get_tuple_name ( Param_tree_path[param_num] );

		//  Lookup the existing branches for reference
		get_absolute_path ( Param_tree_path[param_num], Input_Filename, &Absolute_Param_Paths[param_num] );
		vector<TString> output_branch_names = get_branch_names ( Absolute_Param_Paths[param_num] );

		//  Inspect the user input for any requested branch names

		if ( ( argc >= 4 ) && ( string ( argv[3] ).compare ( "*" ) != 0 ) )
		{
			string pass(argv[3]);
			//  We only want a the output_branch based on the user input
			if ( check_found ( &pass, &output_branch_names ) )
			{
				while ( !output_branch_names.empty() )
				{
					output_branch_names.pop_back();
				}
				output_branch_names.push_back ( argv[3] );
			}
			else
			{	//  Not found exit
				cout << "Branch with name: " << argv[3] << " not found in ttree object " << tree_name << " continuing to next tuple." << endl << endl;
				while ( !output_branch_names.empty() )
				{
					output_branch_names.pop_back();
				}
			}
		}

		vector<TString> Branch_Names;
		vector<bool> left_handed_vector;

		//  Loop over all of the branches within this TTree
		for ( unsigned int output_param_num=0; output_param_num < output_branch_names.size(); output_param_num++ )
		{
			bool left_handedness = false;

			//  Get the Branch of interest in the form of a Histogram
			TH1F* local_histogram = get_branch ( input_pull, Param_tree_path[param_num], output_branch_names[output_param_num] );

			//  Set an individual tree name
			TString Name_String ( tree_name );
			Name_String.Append ( "_" );
			Name_String.Append ( output_branch_names[output_param_num] );
			local_histogram->SetName ( Name_String );
			local_histogram->SetTitle ( Name_String );
			Branch_Names.push_back( Name_String );

			//  Read some simple data from the histogram due to the default binning
			double bin_width = local_histogram->GetBinWidth ( 0 );
			int most_probable_bin = local_histogram->GetMaximumBin();
			double most_probable_value = local_histogram->GetXaxis()->GetBinCenter ( most_probable_bin );

			//double param_rms = local_histogram->GetRMS ( 1 );

			//  If the user requested to fir with the off function, don't fit to the data
			TString fit_type;
			if ( ( argc>4 ) && ( string ( argv[4] ).compare ( "off" ) !=0 ) )  fit_type.Append ( argv[4] );

			//double bin_scale = local_histogram->GetBinContent ( local_histogram->GetMaximumBin() );


			//  Ignore unphysical RMS  it's the telltale sign of an empty TTree RapidFit spits out from fixed parameters in fit

			if( true )//GetOptimalBins( local_histogram ) != 0 )
			{
				//  
				TString format_1;
				tabs_TStr ( tree_name, 3, &format_1 );
				TString format_2;
				tabs_TStr ( output_branch_names[output_param_num], 2, &format_2 );

				cout << "Tree:\t" << tree_name << format_1 << "Branch_Name:\t" <<  output_branch_names[output_param_num] << format_2 << "most_probable_value:\t\t" << most_probable_value << " \\pm " << bin_width/2.0 << endl;

				//	Rebin
				int newnum = OptimumBinNumber( local_histogram );
				(void) newnum;

				TString Fit_Options ( "Q" );


				//  The gamma distribution coded up in root is the more general form of that found on wikipedia (there's a surprise)
				//  
				//  Using the root definition:					wiki:
				//  				gamma = mean^2 / sigma^2		k     = mu^2 / sigma^2	
				//  				beta  = sigma^2 / mean			theta = sigma^2 / mu
				//
				//  		For:		mu == 0				The 2 conditions above are ONLY valid for this condition


				//  More correct to fit to the gamma distribution than the landau function for those that require it
				TF1 *func = new TF1( "gammaf", "[0]*TMath::GammaDist( x, ([1]*[1])/([2]*[2]), 0, ([2]*[2])/[1] )" );

				//  TF1 inherits from TFormula so we can use it's functions to rename the parameters to have consistancy with gaus / landau functions
				func->SetParName( 0, "Constant" );
				func->SetParName( 1, "Mean" );
				func->SetParName( 2, "Sigma" );

				//  I have constructed the gamma function as gammaf in such a way that it is compatible with the gaus and landau built in functions within root
				//  This allows me to simply compare the dfferent fit functions to establish which provides the best fit solution based on chi_squared
				//  I will only fit and plot the best fitting function and the rest of the code is unaware of the changes as the mean and sigma are in the same place for each

				double mean=0,sigma=0,mean_err=0,sigma_err=0,error=0;

				//  If the user has not selected a special function just fit a gaussian
				//		if( ( output_branch_names[output_param_num] != "error" ) )  {
				if ( is_empty ( fit_type ) ){
					//		else {
					//  Not all errors are expected to be gaussian like so we have to test to find the best fit function

					//  Want plotting consistancy so have defined my own landau function
					TF1 *land = new TF1( "mylandau", "[0]*TMath::Landau( x, [1], [2] )" );
					land->SetParName( 0, "Constant" );
					land->SetParName( 1, "Mean" );
					land->SetParName( 2, "Sigma" );

					//  Try Landau function
					local_histogram->Fit ( "mylandau", Fit_Options );
					Double_t chi_2_landau ( local_histogram->GetFunction ( "mylandau" )->GetChisquare() );

					//  Try Gaussian function
					local_histogram->Fit ( "gaus", Fit_Options );
					Double_t chi_2_gaus ( local_histogram->GetFunction ( "gaus" )->GetChisquare() );

					//  The GammaDist function originally complained of invalid results (A LOT) and giving it sensible starting points was a way around this
					func->SetParameters( 1, local_histogram->GetFunction ( "gaus" )->GetParameter ( 1 ), local_histogram->GetFunction ( "gaus" )->GetParameter ( 2 ), 0 );

					local_histogram->Fit ( "gammaf", Fit_Options );
					Double_t chi_2_gamma_f ( local_histogram->GetFunction( "gammaf" )->GetChisquare() );

					//			cout << endl << chi_2_landau << "\t" << chi_2_gaus << "\t" << chi_2_gamma_f << endl << endl;

					//  Determine which function gives the lowest fit chi squared
					if( chi_2_landau < chi_2_gaus )  fit_type.Append( "mylandau" );
					else  if( chi_2_gamma_f < chi_2_gaus )  fit_type.Append( "gammaf" );
					else  fit_type.Append( "gaus" );
				}

				//  Perform fit again (I know it's slightly taxing but what the hey)
				TFitResult* result;

				//		if ( fit_type == "gaus" )  {  result = local_histogram->Fit ( fit_type, Fit_Options );  }
				//		else  {  result = local_histogram->Fit ( fit_type, Fit_Options );  }

				result = new TFitResult( local_histogram->Fit ( fit_type, Fit_Options ) );

				Int_t fitresult = result->Status();	//  Fit Status at end of Fit	i.e.	did it converge?

				if ( string ( argv[4] ).compare ( "off" ) != 0 )
				{
					//  Return 0 for a fit which doesn't converge correctly as this is a problem
					if ( ( fitresult != 0 ) )
					{
						mean = 0;
						error = 0;
					}
					else
					{
						//  Luckily the definition of the parameters for landau and gaus functions is the same for mean and sigma
						mean = local_histogram->GetFunction ( fit_type )->GetParameter ( 1 );
						mean_err = local_histogram->GetFunction ( fit_type )->GetParError ( 1 );
						sigma = local_histogram->GetFunction ( fit_type )->GetParameter ( 2 );
						sigma_err = local_histogram->GetFunction ( fit_type )->GetParError ( 2 );
						error = sigma / sqrt ( local_histogram->GetEntries() );
					}

					cout << "Tree:\t" << tree_name << format_1 << "Branch_Name:\t" <<  output_branch_names[output_param_num] << format_2 << "Fit_most_probable_value:\t" << mean << " \\pm " << error << "\t\tBest Fit_Func:\t" << fit_type << endl;
					cout << "Fit Results:" << endl;
					cout << "Mean:\t" << mean << "\\pm" << mean_err << "\tSigma:\t" << sigma << "\\pm" << sigma_err  << "\t\t\\sqrt{N}:\t" << sqrt( local_histogram->GetEntries() ) << endl << endl;

					cout << "Numerical Average in TBrowser:\t" << local_histogram->GetMean ( 1 ) << endl;
					cout << endl;
					Double_t x_min = local_histogram->GetXaxis()->GetXmin();
					Double_t x_max = local_histogram->GetXaxis()->GetXmax();
					if( ( mean-0.*sigma ) > ( ( x_max + x_min )/2.0 ) )  left_handedness = true;
				}
				} else if( output_branch_names[output_param_num] == "value" ){
					//  It only makes sense to talk about Fixed paramater 'value's as the error is now strictly a histogram binning error
					TString format_1;
					tabs_TStr ( tree_name, 3, &format_1 );
					TString format_2;
					tabs_TStr ( output_branch_names[output_param_num], 2, &format_2 );
					cout << "Tree:\t" << tree_name << format_1 << "was Fixed to: " << most_probable_value << "\\pm" << bin_width/2.0 << endl << endl;
				}
				//  Append the container of all of the histograms
				TH1F* temp2;
				Param_FitResult_Plots[param_num].push_back ( temp2 );

				left_handed_vector.push_back( left_handedness );

				//  This allows me to keep the code cleaner & understandable
				Param_FitResult_Plots[param_num][output_param_num] = ( TH1F* ) local_histogram->Clone();
			}
			left_handed_plot.push_back( left_handed_vector );
			Param_FitResult_Names.push_back( Branch_Names );
			}


			printf("\n\n Plotting Graphs Now:\n\n");

			//  Apologies for the poor naming convntion of pointers here as a plotting function was a very last minute design feature
			//  This will plot the branch data for each ttree object as a new file that was Free in the fit as well as a horizontal plot of each of the branches side by side
			//  1680x1050 was chosen as this is the max res of the latest mbp machines. i.e. full-screen for me, I provide the .C should anyone wish to reproduce the output again

			vector<TString> Extensions;
			Extensions.push_back( TString("png") );
			Extensions.push_back( TString("C") );
			Extensions.push_back( TString("pdf") );
			TCanvas* temp_multi;
			for( unsigned int count1=0; count1<Param_FitResult_Plots.size(); count1++ )
			{
				TString Param_Name = get_tuple_name ( Param_tree_path[count1] );
				int sub_plots = (int) Param_FitResult_Plots[count1].size();
				temp_multi = new TCanvas( Param_Name, Param_Name, 1680*3, 1050 );
				temp_multi->Divide( sub_plots, int(1),0,0,0 );
				vector<TCanvas*> sub_canvas;

				for( unsigned int count2=0; count2<Param_FitResult_Plots[count1].size(); count2++ )
				{
					TCanvas* local_canvas = new TCanvas( Param_FitResult_Names[count1][count2], Param_FitResult_Names[count1][count2], 1680, 1050);
					//double Hist_RMS = Param_FitResult_Plots[count1][count2]->GetRMS ( 1 );
					//double bin_scale = Param_FitResult_Plots[count1][count2]->GetBinContent ( Param_FitResult_Plots[count1][count2]->GetMaximumBin() );
					Param_FitResult_Plots[count1][count2]->SetLineColor( line_colors[0] );
					Param_FitResult_Plots[count1][count2]->SetLineWidth( mystyle_line_width() );
					Param_FitResult_Plots[count1][count2]->GetYaxis()->SetTitle("Entries");
					if( ( (Param_FitResult_Names[count1][count2] == "gamma_value") || (Param_FitResult_Names[count1][count2] == "gamma_error") || (Param_FitResult_Names[count1][count2] == "deltaGamma_value") || (Param_FitResult_Names[count1][count2] == "deltaGamma_error") ) )
					{
						Param_FitResult_Plots[count1][count2]->GetXaxis()->SetTitle( TString( Param_FitResult_Plots[count1][count2]->GetXaxis()->GetTitle()) +" ps^{-1}") ;
					}
					Param_FitResult_Plots[count1][count2]->SetTitle("");
					if( true )//GetOptimalBins( Param_FitResult_Plots[count1][count2] ) != 0 )
					{
						for( unsigned int print_count=0; print_count < Extensions.size(); print_count++ )
						{
							if( left_handed_plot[count1][count2] )
							{  gStyle->SetStatX(float(0.46));  gStyle->SetStatY(float(0.9));  gStyle->SetTitleX(0.8f);  }
							gROOT->ForceStyle();
							Param_FitResult_Plots[count1][count2]->Draw();
							local_canvas->Print(outputdir+"/"+Param_FitResult_Names[count1][count2]+"."+Extensions[ print_count ]);
							if( left_handed_plot[count1][count2] )
							{  gStyle->SetStatX(float(0.9));  gStyle->SetStatY(float(0.9));  gStyle->SetTitleX(0.2f);  }
							gROOT->ForceStyle();
						}
					}
					sub_canvas.push_back( local_canvas );
					temp_multi->cd( count2+1 );
					Param_FitResult_Plots[count1][count2]->Draw();
					temp_multi->Update();
				}
				if( GetOptimalBins( Param_FitResult_Plots[count1][0] ) != 0 )
				{
					for( unsigned int print_count=0; print_count < Extensions.size(); print_count++ )
					{
						temp_multi->Print(outputdir+"/"+Param_Name+"_all"+"."+Extensions[ print_count ]);
					}
				}
				delete temp_multi;
			}

			printf("Finished Successfully, Output Plots stored in: %s\n\n",outputdir.Data());
		}
