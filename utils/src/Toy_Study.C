
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

void LatexDocHeader( stringstream& latex )
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
	latex << "\\usepackage{subcaption}" << endl;
	latex << "\\usepackage{graphicx}" << endl;
	latex << endl;
	latex << "\\begin{document}" << endl;
	latex << endl;
}

void LatexDocFooter( stringstream& latex )
{
	latex << endl;
	latex << "\\end{document}" << endl;
	latex << endl;
}

void TableHeader( stringstream& latex )
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


void TableFooter( stringstream& latex )
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

void ImageSplit( stringstream& latex )
{
	latex << "\\caption{Some Caption for The Image}" << endl;
	latex << "\\label{fig:Some_Label}" << endl;
	latex << "\\end{figure}" << endl;
	latex << "\\clearpage" << endl;
	latex << endl;
	latex << "\\begin{figure}" << endl;
	latex << "\\ContinuedFloat" << endl;
}

void ProcessImageContent( stringstream& latex, vector<TString> all_parameter_plots )
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
		latex << "\\begin{subfigure}[" << latex_name << "]{\\textwidth}" << endl;
		latex << "{\\includegraphics[width=0.3\\textwidth]{" << param_name << "_error_c_thru.png} }" << endl;
		latex << "{\\includegraphics[width=0.3\\textwidth]{" << param_name << "_value_c_thru.png} }" << endl;
		latex << "{\\includegraphics[width=0.3\\textwidth]{" << param_name << "_pull_c_thru.png} }" << endl;
		latex << "\\caption{some sub-caption " << latex_name << "}" << endl;
		latex << "\\label{fig:" << param_name << "}" << endl;
		latex << "\\end{subfigure}" << endl;
		latex << "\\newline" << endl;
	}

	latex << "\\caption{Some Caption for The Image}" << endl;
	latex << "\\label{fig:Some_Label}" << endl;
	latex << "\\end{figure}" << endl;
	latex << endl;
}

void ProcessTableContent( stringstream& latex, vector<TString> all_parameter_plots, vector<pair<pair<double,double>,pair<double,double> > > table_content )
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

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(111);

	TString Output_Dir;

	Output_Dir = "Output" ;

	gSystem->mkdir( Output_Dir );

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
	}

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
		TCanvas* c1 = new TCanvas( Canvas_Name, Canvas_Name, 1680, 1050 );
		Histogram_Processing::Silent_Draw( c1, input_histo );
		TAxis* x_axis = input_histo->GetXaxis();
		x_axis->SetTitle( EdStyle::GetParamRootName( EdStyle::Remove_Suffix(all_parameter_plots[j]) ) + " " 
				+ EdStyle::Get_Suffix( all_parameter_plots[j] ) + " "
				+ EdStyle::GetParamRootUnit( EdStyle::Remove_Suffix(all_parameter_plots[j]) ) );
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

	return 0;
}

