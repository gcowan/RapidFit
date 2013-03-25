
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
#include <sstream>
#include <iomanip>
#include <fstream>

using namespace::std;

void ToyStudyAnalysis::LatexDocHeader( stringstream& latex )
{
	latex << "\\documentclass[a4paper,10pt]{article}" << endl;
	latex << "\\usepackage[utf8]{inputenc}" << endl;
	latex << "\\title{}" << endl;
	latex << "\\author{}" << endl;
	latex << "\\date{}" << endl;
	latex << endl;
	latex << "\\pdfinfo{" << endl;
	latex << "/Title    ()" << endl;
	latex << "/Author   ()" << endl;
	latex << "/Creator  ()" << endl;
	latex << "/Producer ()" << endl;
	latex << "/Subject  ()" << endl;
	latex << "/Keywords ()" << endl;
	latex << "}" << endl;
	latex << endl;
	latex << "\\usepackage{amsmath}" << endl;
	latex << "\\usepackage{caption}" << endl;
	latex << "\\usepackage{graphicx}" << endl;
	latex << "\\usepackage{float}" << endl;
	latex << "\\usepackage{subfig}" << endl;
	latex << endl;
	latex << "\\begin{document}" << endl;
	latex << endl;
}

void ToyStudyAnalysis::LatexDocFooter( stringstream& latex )
{
	latex << endl;
	latex << "\\end{document}" << endl;
	latex << endl;
}

void ToyStudyAnalysis::TableHeader( stringstream& latex )
{
	latex << endl;
	latex << "\\let\\oldpm\\pm" << endl;
	latex << "\\renewcommand{\\pm}{\\ensuremath{\\oldpm}}" << endl;
	latex << "\\begin{table}[h]" << endl;
	latex << "\\begin{center}" << endl;
	latex << "\\begin{tabular}{@{}|l|r|r|@{}}" << endl;
	latex << "\\hline" << endl;
	latex << "Paramater & Sensitivity & Pull Bias \\\\" << endl;
	latex << "\\hline" << endl;
}


void ToyStudyAnalysis::TableFooter( stringstream& latex )
{
	latex << "\\hline" << endl;
	latex << "\\end{tabular}" << endl;
	latex << "\\caption{Some Caption}" << endl;
	latex << "\\label{thisTable}" << endl;
	latex << "\\end{center}" << endl;
	latex << "\\end{table}" << endl;
	latex << "\\renewcommand{\\pm}{\\oldpm}" << endl;
	latex << endl;
}

void ToyStudyAnalysis::ImageSplit( stringstream& latex )
{
	latex << "\\caption{Some Caption for The Image}" << endl;
	latex << "\\label{fig:Some_Label}" << endl;
	latex << "\\end{figure}" << endl;
	latex << "\\clearpage" << endl;
	latex << endl;
	latex << "\\begin{figure}" << endl;
	latex << "\\ContinuedFloat" << endl;
}

void ToyStudyAnalysis::ProcessImageContent( stringstream& latex, vector<TString> all_parameter_plots )
{
	//      3 x each parameter, value, error, pull
	int num_params = (int)all_parameter_plots.size()/3;

	latex << "\\begin{figure}" << endl;

	for( int i=0; i< num_params; ++i )
	{
		if( (i%5 == 0) && (i != 0) )
		{
			ImageSplit( latex );
		}
		string param_name = EdStyle::Remove_Suffix( all_parameter_plots[i] ).Data();
		string latex_name = string(EdStyle::GetParamLatexName( param_name ) );
		latex << "\\subfloat[" << latex_name << " error, value and pull ]{" << endl;
		latex << "\\includegraphics[width=0.3\\textwidth]{" << param_name << "_error_c_thru.png}" << endl;
		latex << "\\includegraphics[width=0.3\\textwidth]{" << param_name << "_value_c_thru.png}" << endl;
		latex << "\\includegraphics[width=0.3\\textwidth]{" << param_name << "_pull_c_thru.png}" << endl;
		latex << "\\label{fig:" << param_name << "}}" << endl;
		latex << "\\newline" << endl;
	}

	latex << "\\caption{Some Caption for The Image}" << endl;
	latex << "\\label{fig:Some_Label}" << endl;
	latex << "\\end{figure}" << endl;
	latex << endl;
}

void ToyStudyAnalysis::ProcessTableContent( stringstream& latex, vector<TString> all_parameter_plots, vector<pair<pair<double,double>,pair<double,double> > > table_content )
{
	//	3 x each parameter, value, error, pull
	int num_params = (int)all_parameter_plots.size()/3;

	latex << setprecision(3);

	for( int i=0; i<num_params; ++i  )
	{
		string param_name = EdStyle::Remove_Suffix( all_parameter_plots[i] ).Data();

		latex << string(EdStyle::GetParamLatexName( param_name ) ) << " & ";
		latex << table_content[i+num_params].first.first << " & ";
		latex << table_content[i+num_params+num_params].first.first << " \\pm ";
		latex << table_content[i+num_params+num_params].first.second << "\\\\" << endl;
	}
}

void ToyStudyAnalysis::Help()
{
	cout << endl << "Toy Studies currently accept the extra run time arguments:" << endl << endl;
	cout << "--allData" << "\t\t" << "This option stops the tool automatically filtering toys which have a parameter with a pull > 5 sigma from being used" << endl;
	cout << "\t\t\t\t" << "Currently this useful to keep turned on as the toys are post-processed in a binned fasion rather than an un-binned fit" << endl;
	cout << endl;
	cout << "--runWithAnyToys" << "\t" << "This option is used to force the tool to run when you have < 100 toys in the input file" << endl;
	cout << "\t\t\t\t" << "This is known to likely cause problems in poorly behaved studies due to low stats" << endl;
	cout << endl;
}

int ToyStudyAnalysis::Toy_Study( TTree* input_tree, TRandom3* rand_gen, vector<string> OtherOptions )
{
	bool noPullCuts=false;

	for( unsigned int i=0; i< OtherOptions.size(); ++i )
	{
		string thisOption = OtherOptions[i];

		if( thisOption == "--allData" )
		{
			noPullCuts = true;
		}
	}

	string runAlways="--runWithAnyToys";

	if( input_tree->GetEntries() < 100 )
	{
		cout << endl << endl;
		cout << "You will likely end up with errors due to the low number of entries in this study" << endl;
		cout << input_tree->GetEntries() << endl;
		cout << "If you are serious about continuing then re-run with --runWithAnyToys" << endl;
		cout << "Try using hadd to merge the contents of smaller studies which are part of a larger run" << endl << endl;
		cout << "i.e. hadd myStudy.root file1.root file2.root ..." << endl;
		cout << endl;
		if( StringOperations::VectorContains( OtherOptions, runAlways ) == -1 )
		{
			cout << "Silently continuing without running Toy Study analysis" << endl;
			return -1;
		}
		else
		{
			cout << "Performing Toy Study despite having less than 100 toys in the study... good luck" << endl;
		}
	}

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(111);

	TString Output_Dir;

	Output_Dir = "ToyStudy_Output_"; Output_Dir.Append( StringOperations::TimeString() );

	gSystem->mkdir( Output_Dir );

        vector<TString> free_param = RapidFit_Output_File::get_free_parameters( input_tree );

	vector<TString> all_parameters, all_parameter_values, all_parameter_errors, all_parameter_pulls;

	for( unsigned int i=0; i< free_param.size(); ++i )
	{
		//cout << free_param[i] << " " << i << endl;
		string this_param = free_param[i].Data();
		all_parameters.push_back( this_param );
		all_parameter_values.push_back( this_param+string(value_suffix.Data()) );
		all_parameter_errors.push_back( this_param+string(error_suffix.Data()) );
		all_parameter_pulls.push_back( this_param+string(pull_suffix.Data()) );
	}

	/*
	 *	OLD WAY THAT PLOTTED TOO MUCH USELESS FREE IGNORED PARAMETERS
	 *
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
	}*/

	vector<vector<TString> > all; all.push_back( all_parameter_values ); all.push_back( all_parameter_errors ); all.push_back( all_parameter_pulls );

	vector<TString> all_parameter_plots = concatonnate( all );

	cout << all_parameter_plots.size() << " Plots to Draw" << endl;

	TString all_pulls_lt5;
	if( !noPullCuts )
	{
		for( unsigned int j=0; j< all_parameter_plots.size(); ++j )
		{
			string param_name = EdStyle::Remove_Suffix( all_parameter_plots[j] ).Data();
			//      Results at 5 sigma are heavily biased
			all_pulls_lt5.Append( "&&(abs("+TString(param_name.c_str())+"_pull)<=5.)" );
		}
	}

	stringstream latex;

	LatexDocHeader( latex );

	stringstream latexTable;

	TableHeader( latexTable );

	vector<pair<pair<double,double>,pair<double,double> > > table_content;

	//	Read in the number of 'good' rows in the file
	TString cut_str("(Fit_Status==3)");
	input_tree->Draw( all_parameter_plots[0], cut_str, "goff" );
	int total_rows = (int)input_tree->GetSelectedRows();

	//	Read in the number of rows biased by <= 5 sigma
	TString cut_str2 = cut_str; cut_str2.Append( all_pulls_lt5 );
	input_tree->Draw( all_parameter_plots[0], cut_str2, "goff" );
	int usable_rows = (int)input_tree->GetSelectedRows();

	//	if more than 90% of the data is lost then the data is strongly biased and cannot construct a sensible cut on the toys to clean things up.
	if( ((double)usable_rows)<=(0.1*(double)total_rows)  )
	{
		cerr << endl << "All of your data is biased by > 5 sigma, please check you understand what is going on!" << endl;
		cerr << "Removing cut that removes all fit results with a pull > 5" << endl << endl;
		noPullCuts = true;
	}		//	50% 'Turn Off Cut' test
	else if( ((double)usable_rows)<=(0.5*(double)total_rows)  )
	{
		cerr << endl << "More than half of your data is biased by > 5 sigma, check you understand what is going on!" << endl << endl;
	}

	//	90% Cut test			   &&		50% 'Turn Off Cut' test
	if( ((double)usable_rows)<=(0.9*(double)total_rows)&&((double)usable_rows)>(0.5*(double)total_rows)  )
	{
		cerr << "More than 10\% of data has been removed by the cut 'All pulls <=5'. To replot with all data use the option --allData" << endl << endl;
		cerr << "Total Fit Results: " << total_rows << "\tNumber of Fit Results Cut by demanding 'All pulls <=5': " << total_rows-usable_rows << endl;
	}

	if( usable_rows!=total_rows )
	{
		cerr << "This Tool has thrown Away " << total_rows - usable_rows << " by insisting that all toys are <= 5 sigma." << endl;
		cerr << endl << "To use ALL DATA regardless of outliers run with the additional option: --allData" << endl << endl;
	}

	//      Results at >5 sigma are heavily biased, fit is by definition badly behaved
	if( !noPullCuts ) cut_str.Append( all_pulls_lt5 );

	cout << endl;
	for( unsigned int j=0; j< all_parameter_plots.size(); ++j )
	{

		cout << "Plotting: " << all_parameter_plots[j] << setw(20) << " " <<  "\r" << flush;

		input_tree->Draw( all_parameter_plots[j], cut_str, "goff" );

		TH1* input_histo = (TH1*) input_tree->GetHistogram();
		TString Histo_Name( all_parameter_plots[j] );
		Histo_Name+=rand_gen->Rndm();
		input_histo->SetName( Histo_Name );

		Histogram_Processing::OptimallyRebin( input_histo );

		TString fit_type = Histogram_Processing::Best_Fit_Function( input_histo );

		Histogram_Processing::Silent_Fit( input_histo, fit_type );

		TString Canvas_Name("Canvas");	Canvas_Name+=rand_gen->Rndm();
		TCanvas* c1 = EdStyle::RapidFitCanvas( Canvas_Name, Canvas_Name );
		Histogram_Processing::Silent_Draw( c1, input_histo );
		TAxis* x_axis = input_histo->GetXaxis();
		string units;
		string suffix = EdStyle::Get_Suffix( string(all_parameter_plots[j].Data()) );
		if( suffix != EdStyle::Get_Suffix(string(pull_suffix.Data())) )
		{
			units = " " + EdStyle::GetParamRootUnit( EdStyle::Remove_Suffix(all_parameter_plots[j]) );
		}
		//cout << endl << suffix << "  " << pull_suffix << endl;
		x_axis->SetTitle( EdStyle::GetParamRootName( EdStyle::Remove_Suffix(all_parameter_plots[j]) ) + " " 
				+ EdStyle::Get_Suffix( all_parameter_plots[j] )
				+ units.c_str() );
		input_histo->SetTitle( "Toy Study " + EdStyle::GetParamRootName( all_parameter_plots[j] ) + " distribution" );

		c1->Update();

		Histogram_Processing::Silent_Print( c1, Output_Dir+"/"+all_parameter_plots[j]+".png" );
		Histogram_Processing::Silent_Print( c1, Output_Dir+"/"+all_parameter_plots[j]+".pdf" );
		Histogram_Processing::Silent_Print( c1, Output_Dir+"/"+all_parameter_plots[j]+".C" );

		TPaveStats* thisStats = (TPaveStats*)input_histo->GetListOfFunctions()->FindObject("stats");
		thisStats->SetFillStyle(3023);
		thisStats->SetTextColor(1);

		thisStats->Draw();

		c1->Update();

		TF1* fitted = input_histo->GetFunction( fit_type );

		table_content.push_back( make_pair(
					make_pair( fitted->GetParameter(1), fitted->GetParError(1) ),
					make_pair( fitted->GetParameter(2), fitted->GetParError(2) )
					) );

		Histogram_Processing::Silent_Print( c1, Output_Dir+"/"+all_parameter_plots[j]+"_c_thru.png" );
		Histogram_Processing::Silent_Print( c1, Output_Dir+"/"+all_parameter_plots[j]+"_c_thru.pdf" );
		Histogram_Processing::Silent_Print( c1, Output_Dir+"/"+all_parameter_plots[j]+"_c_thru.C" );
	}
	cout << "Finished" << setw(40) << " " << endl;
	cout << endl;
	ProcessTableContent( latexTable, all_parameter_plots, table_content );

	TableFooter( latexTable );

	ofstream outFileTable;
	outFileTable.open( Output_Dir+"/ToyStudy_Table.tex");
	outFileTable << latexTable.str();
	outFileTable.close();

	latex << "\\input{ToyStudy_Table}" << endl;
	latex << endl;
	latex << "\\clearpage" << endl;

	stringstream latexImage;


	ProcessImageContent( latexImage, all_parameter_plots );

	ofstream outFileImage;
	outFileImage.open( Output_Dir+"/ToyStudy_Image.tex");
	outFileImage << latexImage.str();
	outFileImage.close();

	latex << "\\input{ToyStudy_Image}" << endl;
	latex << endl;

	LatexDocFooter( latex );

	ofstream outFile;
	outFile.open( Output_Dir+"/ToyStudy_LaTeX.tex" );
	outFile << latex.str();
	outFile.close();

	cout << "Compiling Latex Summary..." << endl;

	gSystem->cd( Output_Dir );

	TString command="pdflatex -interaction=batchmode ToyStudy_LaTeX.tex > /dev/null";
	//cout << command << endl;
	system( command );

	gSystem->cd("..");

	return 0;
}

