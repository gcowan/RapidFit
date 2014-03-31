
#include "EdStyle.h"
#include "TTree_Processing.h"
#include "StringOperations.h"
#include "RapidFit_Output_File.h"
#include "ROOT_File_Processing.h"

#include "TTree.h"
#include "TString.h"

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <fstream>

using namespace::std;

string buildLatexTable( vector<TString> params, vector<double> CV1, vector<double> Err1, vector<double> CV2, vector<double> Err2 );

int main( int argc, char* argv[] )
{

	if( argc != 3 )
	{
		cout << "USAGE:" << endl;
		cout << "\t\t" << argv[0] << "\t" << "RapidFitOutput_1.root" << "\t" << "RapidFitOutput_2.root" << endl;
		cout << endl;
		exit(-1);
	}

	TTree* RapidFitResult_1 = ROOT_File_Processing::GetFirstTree( argv[1], "READ" );
	TTree* RapidFitResult_2 = ROOT_File_Processing::GetFirstTree( argv[2], "READ" );

	vector<TString> FreeParams_Result1 = RapidFit_Output_File::get_free_parameters( RapidFitResult_1 );
	vector<TString> FreeParams_Result2 = RapidFit_Output_File::get_free_parameters( RapidFitResult_2 );

	if( FreeParams_Result1.size() != FreeParams_Result2.size() )
	{
		cout << endl;
		cout << "Free Number of Params different between 2 fits" << endl;
		cout << endl;
	}

	vector<TString> common;
	for( unsigned int i=0; i< FreeParams_Result1.size(); ++i )
	{
		if( StringOperations::VectorContains( FreeParams_Result2, FreeParams_Result1[i] ) != -1 )
		{
			common.push_back( FreeParams_Result1[i] );
		}
	}

	vector<Double_t> CV_1;
	vector<Double_t> CV_2;
	vector<Double_t> Err_1;
	vector<Double_t> Err_2;
	vector<Double_t> CV_diffs;
	vector<Double_t> Err_diffs;

	for( unsigned int i=0; i< common.size(); ++i )
	{
		TString valueStr( common[i] );
		valueStr.Append( "_value" );
		TString errorStr( common[i] );
		errorStr.Append( "_error" );


		vector<Double_t>* thisValue1 = TTree_Processing::Buffer_Branch( RapidFitResult_1, valueStr );
		vector<Double_t>* thisError1 = TTree_Processing::Buffer_Branch( RapidFitResult_1, errorStr );

		vector<Double_t>* thisValue2 = TTree_Processing::Buffer_Branch( RapidFitResult_2, valueStr );
		vector<Double_t>* thisError2 = TTree_Processing::Buffer_Branch( RapidFitResult_2, errorStr );

		double diff_CV = thisValue1->at( 0 ) - thisValue2->at( 0 );
		double diff_Err = thisError1->at( 0 ) - thisError2->at( 0 );

		CV_1.push_back( thisValue1->at( 0 ) );
		CV_2.push_back( thisValue2->at( 0 ) );
		Err_1.push_back( thisError1->at( 0 ) );
		Err_2.push_back( thisError2->at( 0 ) );

		CV_diffs.push_back( diff_CV );
		Err_diffs.push_back( diff_Err );

		delete thisValue1;
		delete thisError1;
		delete thisValue2;
		delete thisError2;
	}

	cout << endl;

	cout << setprecision( 8 ) << fixed << endl;
	cout << setw( 20 ) << "Parameter     " << "   " << setw( 29 ) << " Result1      " << setw( 29 ) << " Result2      ";
	cout << setw( 21 ) << " Diff in CV " << setw( 21 ) << " Error Ratio " << "  " << setw( 23 ) << " sqrt(|Err1^2-Err2^2|) ";
	cout << setw( 25 ) << " sqrt(|Err1^2-Err2^2|)/Err1 " << "  " << setw( 25 ) << " sqrt(|Err1^2-Err2^2|)/Err2 " << "  ";
	cout << setw( 25 ) << "sqrt( Diff^2+DiffErr^2 )" << endl;

	cout << endl;
	for( unsigned int i=0; i< common.size(); ++i )
	{
		double ErrDiff = sqrt( fabs(Err_1[i]*Err_1[i] - Err_2[i]*Err_2[i]) );
		cout << setw( 20 ) << common[i] << " : " << setw( 13 ) << CV_1[i] << " ± " << setw( 13 ) << Err_1[i] << "\t" << setw( 13 ) << CV_2[i] << " ± " << setw( 13 ) << Err_2[i] << "   ";
		cout << setw( 15 ) << CV_diffs[i] << "   " <<  setw( 15 ) << Err_1[i] / Err_2[i] << "  " << setw( 23 ) << ErrDiff;
		cout << setw( 25 ) << ErrDiff/Err_1[i] << "  " << setw( 25 ) << ErrDiff/Err_2[i];
		cout << "  " << setw( 23 ) << sqrt( CV_diffs[i]*CV_diffs[i] + ErrDiff*ErrDiff );
		cout << endl;
	}
	cout << endl;

	string thisTable = buildLatexTable( common, CV_1, Err_1, CV_2, Err_2 );

	cout << thisTable;

	cout << endl << endl;

	return 0;

}

string buildLatexTable( vector<TString> params, vector<double> CV1, vector<double> Err1, vector<double> CV2, vector<double> Err2 )
{

	stringstream thisTable;

	thisTable << setprecision( 8 ) << fixed << endl;

	thisTable << "\\documentclass[a4paper,10pt]{article}" << endl;
	thisTable << "\\usepackage[utf8]{inputenc}" << endl;
	thisTable << "\\usepackage[landscape]{geometry}" << endl;
	thisTable << "\\usepackage{amsmath}" << endl;
	thisTable << endl;
	thisTable << "\\begin{document}" << endl;
	thisTable << endl;
	thisTable << "\\renewcommand{\\arraystretch}{1.25}" << endl;
	thisTable << "\\begin{table}" << endl;
	thisTable << "\\hspace*{-3.cm}\\begin{tabular}{@{}l|c|c|c|c|c|c|c|c@{}}" << endl;
	thisTable << setw( 30 ) << "Parameter   & " << "   " << setw( 29 ) << " Result1    & " << setw( 29 ) << " Result2    & ";
	thisTable << setw( 15 ) << " Diff in CV & " << setw( 10 ) << " Error Ratio & " << setw( 20 ) << " $\\sqrt{ \\left|E1^2-E2^2\\right| }$";
	thisTable << setw( 25 ) << " & $\\dfrac{\\sqrt{ \\left|E1^2-E2^2\\right| } }{ E1 }$ " << " & " << setw( 25 ) << " $\\dfrac{\\sqrt{ \\left|E1^2-E2^2\\right| } }{ E2 }$ ";
	thisTable << setw( 23 ) << " & $\\sqrt{ \\left(\\Delta{CV}\\right)^2+\\Delta{E^2} }$ " << "   \\\\ \\hline" << endl;

	for( unsigned int i=0; i< params.size(); ++i )
	{

		double ErrDiff = sqrt( fabs(Err1[i]*Err1[i] - Err2[i]*Err2[i]) );
		double CVDiff = CV1[i]-CV2[i];
		thisTable << setw( 30 ) << EdStyle::GetParamLatexName( params[i] ) << " & " << setw( 13 ) << CV1[i] << " $\\pm$ " << setw( 13 ) << Err1[i] << " & ";
		thisTable << setw( 13 ) << CV2[i] << " $\\pm$ " << setw( 13 ) << Err2[i] << " & ";
		thisTable << setw( 13 ) << CVDiff << " & " << setw( 13 ) << Err1[i] / Err2[i] << " & " << setw( 13 ) << ErrDiff;
		if( fabs( ErrDiff / Err1[i] ) > 0.1 || fabs( ErrDiff / Err2[i] ) > 0.1 ) thisTable << " & \\textbf{" << setw( 13 ) << ErrDiff / Err1[i] << "} & \\textbf{" << setw( 13 ) << ErrDiff / Err2[i] << "}";
		else thisTable << " & " << setw( 13 ) << ErrDiff / Err1[i] << " & " << setw( 13 ) << ErrDiff / Err2[i];
		thisTable << " & " << setw( 13 ) << sqrt( CVDiff*CVDiff + ErrDiff*ErrDiff );
		thisTable << "   \\\\";
		thisTable << endl;

	}

	thisTable << "\\end{tabular}" << endl;
	thisTable << "\\end{table}" << endl;
	thisTable << endl;
	thisTable << "\\end{document}" << endl;
	thisTable << endl;

	return thisTable.str();
}


