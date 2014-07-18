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
#include "ROOT_File_Processing.h"


#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
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

	TString Output_Dir;
	Output_Dir = "ToyStudy_Output_"; Output_Dir.Append( StringOperations::TimeString() );
	gSystem->mkdir( Output_Dir );

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


	for( unsigned int i=0; i< common.size(); ++i )
	{
		TString valueStr( common[i] );
		valueStr.Append( "_value" );

		TString errorStr( common[i] );
		errorStr.Append( "_error" );

		TString pullStr( common[i] );
		pullStr.Append( "_pull" );

		vector<Double_t>* thisValue1 = TTree_Processing::Buffer_Branch( RapidFitResult_1, valueStr );
		vector<Double_t>* thisError1 = TTree_Processing::Buffer_Branch( RapidFitResult_1, errorStr );
		vector<Double_t>* thisPull1 = TTree_Processing::Buffer_Branch( RapidFitResult_1, pullStr );

		vector<Double_t>* thisValue2 = TTree_Processing::Buffer_Branch( RapidFitResult_2, valueStr );
		vector<Double_t>* thisError2 = TTree_Processing::Buffer_Branch( RapidFitResult_2, errorStr );
		vector<Double_t>* thisPull2 = TTree_Processing::Buffer_Branch( RapidFitResult_2, pullStr );

		TString ValueName(valueStr);

		TH1D* deltaValue = new TH1D(valueStr,valueStr,1000,-1.0,1.0);
		deltaValue->SetName( ValueName );


		cout << "Number of entries in vectors: " << thisValue1->size() << endl;
		for(unsigned int j=0; j< thisValue1->size(); j++){

			deltaValue->Fill(thisValue1->at( j ) - thisValue2->at( j ));
			//deltaError->Fill(thisError1->at( j ) - thisError2->at( j ));
			//deltaPull->Fill(thisPull1->at( j ) - thisPull2->at( j ));



		}

		delete thisValue1;
		delete thisError1;
		delete thisPull1;
		delete thisValue2;
		delete thisError2;
		delete thisPull2;


		Histogram_Processing::OptimallyRebin( deltaValue );

		TString ValueFitType = Histogram_Processing::Best_Fit_Function( deltaValue );
		Histogram_Processing::Silent_Fit( deltaValue, ValueFitType );
		TString Canvas_Name(valueStr+"_Canvas");  
		TCanvas* c1 = EdStyle::RapidFitCanvas( Canvas_Name, Canvas_Name );
		Histogram_Processing::Silent_Draw( c1, deltaValue );
		TAxis* x_axis = deltaValue->GetXaxis();
		string units;
		units = " " + EdStyle::GetParamRootUnit( EdStyle::Remove_Suffix(valueStr) );
		x_axis->SetTitle( EdStyle::GetParamRootName( EdStyle::Remove_Suffix(valueStr )) + " " + EdStyle::Get_Suffix(valueStr) + units.c_str() );
		deltaValue->SetTitle( "Toy Study " + EdStyle::GetParamRootName( valueStr ) + " difference" );

				c1->Update();

				Histogram_Processing::Silent_Print( c1, Output_Dir+"/"+valueStr+".png" );
				Histogram_Processing::Silent_Print( c1, Output_Dir+"/"+valueStr+".pdf" );
				Histogram_Processing::Silent_Print( c1, Output_Dir+"/"+valueStr+".C" );




				}

				}


