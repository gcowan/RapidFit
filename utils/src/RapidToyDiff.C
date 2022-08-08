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

void MakePlot(TH1* hist, TString Name, TString dir){

	Histogram_Processing::OptimallyRebin( hist );

	TString FitType = Histogram_Processing::Best_Fit_Function( hist );
	Histogram_Processing::Silent_Fit( hist, FitType );
	TString Canvas_Name(Name+"_Canvas");  
	TCanvas* c1 = EdStyle::RapidFitCanvas( Canvas_Name, Canvas_Name );
	Histogram_Processing::Silent_Draw( c1, hist );
	TAxis* x_axis = hist->GetXaxis();
	string units;
	units = " " + EdStyle::GetParamRootUnit( EdStyle::Remove_Suffix(Name) );
	if(!EdStyle::GetParamRootName(Name).Contains("ull")){
	x_axis->SetTitle( EdStyle::GetParamRootName( EdStyle::Remove_Suffix(Name )) + " " + EdStyle::Get_Suffix(Name) + units.c_str() );
	}else{
	x_axis->SetTitle( EdStyle::GetParamRootName( EdStyle::Remove_Suffix(Name )) + " " + EdStyle::Get_Suffix(Name) );
	}
	hist->SetTitle( "Toy Study " + EdStyle::GetParamRootName( Name ) + " difference" );

	c1->Update();

	Histogram_Processing::Silent_Print( c1, dir+"/"+Name+".png" );
	Histogram_Processing::Silent_Print( c1, dir+"/"+Name+".pdf" );
	Histogram_Processing::Silent_Print( c1, dir+"/"+Name+".C" );

	TPaveStats* thisStats = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
	thisStats->SetFillStyle(3023);
	thisStats->SetTextColor(1);

	thisStats->Draw();

	c1->Update();

	TF1* fitted = hist->GetFunction( FitType );
	(void) fitted;

	Histogram_Processing::Silent_Print( c1, dir+"/"+Name+"_c_thru.png" );
	Histogram_Processing::Silent_Print( c1, dir+"/"+Name+"_c_thru.pdf" );
	Histogram_Processing::Silent_Print( c1, dir+"/"+Name+"_c_thru.C" );





}

int main( int argc, char* argv[] )
{

	if( argc != 3 )
	{
		cout << "USAGE:" << endl;
		cout << "\t\t" << argv[0] << "\t" << "RapidFitOutput_1.root" << "\t" << "RapidFitOutput_2.root" << endl;
		cout << endl;
		exit(-1);
	}

	        //      Setup the Canvas and such
		              EdStyle* RapidFit_Style = new EdStyle();
		                      RapidFit_Style->SetStyle();


	     gStyle->SetOptStat(0);
	             gStyle->SetOptFit(111);

	TString Output_Dir;
	Output_Dir = "DiffToyStudy_Output_"; Output_Dir.Append( StringOperations::TimeString() );
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

		TString ptpullStr(common[i]);
		ptpullStr.Append("_perToyPull");

		vector<Double_t>* thisValue1 = TTree_Processing::Buffer_Branch( RapidFitResult_1, valueStr );
		vector<Double_t>* thisError1 = TTree_Processing::Buffer_Branch( RapidFitResult_1, errorStr );
		vector<Double_t>* thisPull1 = TTree_Processing::Buffer_Branch( RapidFitResult_1, pullStr );

		vector<Double_t>* thisValue2 = TTree_Processing::Buffer_Branch( RapidFitResult_2, valueStr );
		vector<Double_t>* thisError2 = TTree_Processing::Buffer_Branch( RapidFitResult_2, errorStr );
		vector<Double_t>* thisPull2 = TTree_Processing::Buffer_Branch( RapidFitResult_2, pullStr );



		TH1D* deltaValue = new TH1D(valueStr,valueStr,1000,1.0,-1.0);
		deltaValue->SetName( valueStr );
		
		TH1D* deltaError = new TH1D(errorStr,errorStr,1000,1.0,-1.0);
		deltaError->SetName( errorStr );

		TH1D* deltaPull = new TH1D(pullStr,pullStr,1000,1.0,-1.0);
		deltaPull->SetName( pullStr );

		TH1D* perToyPull = new TH1D(ptpullStr,ptpullStr,1000,1.0,-1.0);
		perToyPull->SetName( ptpullStr );
	

		cout << "Number of entries in vectors: " << thisValue1->size() << endl;
		for(unsigned int j=0; j< thisValue1->size(); j++){

			deltaValue->Fill(thisValue1->at( j ) - thisValue2->at( j ));
			deltaError->Fill(thisError1->at( j ) - thisError2->at( j ));

			deltaPull->Fill(thisPull1->at( j ) - thisPull2->at( j ));
			perToyPull->Fill((thisValue1->at( j ) - thisValue2->at( j ))/(thisError1->at( j )));



		}


		MakePlot(deltaValue, valueStr, Output_Dir);
		MakePlot(deltaError, errorStr, Output_Dir);
		MakePlot(deltaPull, pullStr, Output_Dir);
		MakePlot(perToyPull, ptpullStr, Output_Dir);

		delete thisValue1;
		delete thisError1;
		delete thisPull1;
		delete thisValue2;
		delete thisError2;
		delete thisPull2;




	}

}



