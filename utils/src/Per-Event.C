
#include "TMath.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "Histo_Processing.h"
#include "ROOT_File_Processing.h"
#include "TTree_Processing.h"
#include "Template_Functions.h"
#include "RapidFit_Output_File.h"
#include "StringProcessing.h"

#include <string>
#include <map>
#include <vector>

using namespace::std;


int main( int argc, char* argv[] )
{
	if( argc != 3 )
	{
		cout << "Usage: " << argv[0] << " nTuple1.root  nTuple2.root" << endl;
		exit(0);
	}

	string nTuple1Name=argv[1];//"Scale_woSWave.root";
	string nTuple2Name=argv[2];//"Scale_wSWave.root";

	string FitStatusName="Fit_Status";

	string OutputFileName = "OffSetStudy_Output.root";

	TFile* outputFile = new TFile( OutputFileName.c_str(), "RECREATE" );

	TTree* outputTree = new TTree( "TwoNtuplePerEventDiff", "TwoNtuplePerEventDiff" );

	gSystem->mkdir( "PerEventOutput" );

	TTree* nTupleOne = ROOT_File_Processing::GetFirstTree( nTuple1Name );
	TTree* nTupleTwo = ROOT_File_Processing::GetFirstTree( nTuple2Name );

	vector<TString> free_param = RapidFit_Output_File::get_free_parameters( nTupleOne );
	vector<TString> free_param2 = RapidFit_Output_File::get_free_parameters( nTupleTwo );

	vector<string> freeParam; for( unsigned int i=0; i< free_param.size(); ++i ) freeParam.push_back( free_param[i].Data() );
	vector<string> freeParam2; for( unsigned int i=0; i< free_param2.size(); ++i ) freeParam2.push_back( free_param2[i].Data() );

	vector<string> allNames = StringProcessing::CombineUniques( freeParam, freeParam2 );

	for( unsigned int name_i=0; name_i< allNames.size(); ++name_i )
	{
		string ValName = allNames[name_i]; ValName.append( "_value" );
		string GenName = allNames[name_i]; GenName.append( "_gen" );
		string ErrName = allNames[name_i]; ErrName.append( "_error" );

		vector<double>* nTupleOne_Val = TTree_Processing::Buffer_Branch( nTupleOne, ValName );
		vector<double>* nTupleOne_Gen = TTree_Processing::Buffer_Branch( nTupleOne, GenName );
		vector<double>* nTupleOne_Err = TTree_Processing::Buffer_Branch( nTupleOne, ErrName );
		vector<double>* nTupleOne_FitStat = TTree_Processing::Buffer_Branch( nTupleOne, FitStatusName );

		vector<double>* nTupleTwo_Val = TTree_Processing::Buffer_Branch( nTupleTwo, ValName );
		vector<double>* nTupleTwo_Gen = TTree_Processing::Buffer_Branch( nTupleTwo, GenName );
		vector<double>* nTupleTwo_Err = TTree_Processing::Buffer_Branch( nTupleTwo, ErrName );
		vector<double>* nTupleTwo_FitStat = TTree_Processing::Buffer_Branch( nTupleTwo, FitStatusName );

		vector<double> OneDiff;
		vector<double> TwoDiff;
		vector<double> AbsDiff;

		vector<double> OnePull;
		vector<double> TwoPull;

		vector<double> Pull;
		vector<double> PropPull;

		vector<double> PullDiff;

		for( unsigned int i=0; i< nTupleOne_Val->size(); ++i )
		{
			if( fabs((*nTupleOne_FitStat)[i]-3.)<1E-5 && fabs((*nTupleTwo_FitStat)[i]-3.)<1E-5 )
			{
				double OneDiffVal = (*nTupleOne_Gen)[i] - (*nTupleOne_Val)[i];
				double TwoDiffVal = (*nTupleTwo_Gen)[i] - (*nTupleTwo_Val)[i];

				double DiffVal = (*nTupleTwo_Val)[i] - (*nTupleOne_Val)[i];

				double OnePullVal = OneDiffVal / (*nTupleOne_Err)[i];
				double TwoPullVal = TwoDiffVal / (*nTupleTwo_Err)[i];

				double PropErr = sqrt( (*nTupleOne_Err)[i]*(*nTupleOne_Err)[i]
						+ (*nTupleTwo_Err)[i]*(*nTupleTwo_Err)[i] );

				double AvErr = 0.5*( (*nTupleOne_Err)[i] + (*nTupleTwo_Err)[i] );

				double PullVal = DiffVal / AvErr;

				double PullPropVal = DiffVal / PropErr;

				double PullDiffVal = OnePullVal - TwoPullVal;

				OneDiff.push_back( OneDiffVal );
				TwoDiff.push_back( TwoDiffVal );
				AbsDiff.push_back( DiffVal );
				OnePull.push_back( OnePullVal );
				TwoPull.push_back( TwoPullVal );
				Pull.push_back( PullVal );
				PropPull.push_back( PullPropVal );
				PullDiff.push_back( PullDiffVal );
			}
			else
			{
				OneDiff.push_back( -9999. );
				TwoDiff.push_back( -9999. );
				AbsDiff.push_back( -9999. );
				OnePull.push_back( -9999. );
				TwoPull.push_back( -9999. );
				Pull.push_back( -9999. );
				PropPull.push_back( -9999. );
				PullDiff.push_back( -9999. );
			}
		}

		string Param_Diff1 = allNames[name_i]; Param_Diff1.append( "_value_gen_diff_nTuple1" );
		string Param_Diff2 = allNames[name_i]; Param_Diff2.append( "_value_gen_diff_nTuple2" );
		string DiffBranch = allNames[name_i]; DiffBranch.append( "_value_value_diff" );
		string Param_Pull1 = allNames[name_i]; Param_Pull1.append( "_calc_pull" );
		string Param_Pull2 = allNames[name_i]; Param_Pull2.append( "_calc_pull" );
		string Pull_Name = allNames[name_i]; Pull_Name.append( "_DiffPull" );
		string PropPull_Name = allNames[name_i]; PropPull_Name.append( "_DiffPull_wPropErr" );
		string PullDiff_Name = allNames[name_i]; PullDiff_Name.append( "_PullDiff" );

		TTree_Processing::AddBranch( outputTree, Param_Diff1, OneDiff );
		TTree_Processing::AddBranch( outputTree, Param_Diff2, TwoDiff );
		TTree_Processing::AddBranch( outputTree, DiffBranch, AbsDiff );
		TTree_Processing::AddBranch( outputTree, Param_Pull1, OnePull );
		TTree_Processing::AddBranch( outputTree, Param_Pull2, TwoPull );
		TTree_Processing::AddBranch( outputTree, Pull_Name, Pull );
		TTree_Processing::AddBranch( outputTree, PropPull_Name, PropPull );
		TTree_Processing::AddBranch( outputTree, PullDiff_Name, PullDiff );

		outputFile->cd();
		outputTree->Write("",TObject::kOverwrite);

		string cut = PullDiff_Name; cut.append(">-9999.");
		vector<double> DiffPlotData;
		for( unsigned int i=0; i< PullDiff.size(); ++i ) if( PullDiff[i] > -9999. ) DiffPlotData.push_back( PullDiff[i] );

		vector<double> AbsDiffData;
		for( unsigned int i=0; i< AbsDiff.size(); ++i ) if( AbsDiff[i] > -9999. ) AbsDiffData.push_back( AbsDiff[i] );

		string DiffPlotFileName = "PerEventOutput/"; DiffPlotFileName.append( PullDiff_Name ); DiffPlotFileName.append(".pdf");
		TString ThisCanv("ThisCanv");ThisCanv+=name_i;
		TCanvas* c0 = new TCanvas( ThisCanv, ThisCanv );
		TH1* thisTH1 = Histogram_Processing::Plot_1D( DiffPlotData, DiffPlotFileName, "PE9" );
		c0->cd();
		thisTH1->Draw( "PE9" );
		c0->Update();
		c0->Print( DiffPlotFileName.c_str() );

		TString FuncName="newfunc_"; FuncName+=name_i;
		//TF1* newFunc = new TF1( FuncName, "[0]*exp(-((x-[1])*(x-[1])/[2])) + [3]*exp(-((x-[4])*(x-[4])/[5]))" );
		TF1* newFunc = new TF1( FuncName, "[0]*exp(-((x-[1])*(x-[1])/[2])) + [3]*exp(-((x-[1])*(x-[1])/[4]))" );
		newFunc->SetParameter( 0, 100);
		//newFunc->SetParameter( 1, 0.05);
		newFunc->SetParameter( 1, thisTH1->GetBinCenter( thisTH1->GetMaximumBin()) );
		newFunc->SetParameter( 2, 1.);
		newFunc->SetParameter( 3, 500);
		//newFunc->SetParameter( 4, -0.1);
		//newFunc->SetParameter( 5, 50.);
		newFunc->SetParameter( 4, 50. );

		TString CanvName = "Canvas_"; CanvName+=name_i;
		TCanvas* thisCanv = new TCanvas( CanvName, CanvName );

		thisTH1->Fit( newFunc, "M" );

		thisCanv->Update();

		TString X_Title = allNames[name_i]; X_Title.Append( " Difference in measured Pull" );
		thisTH1->GetYaxis()->SetTitle( "# of Results" );
		thisTH1->GetXaxis()->SetTitle( X_Title );

		thisCanv->Update();

		TString FitFileName = "PerEventOutput/"; FitFileName.Append( PullDiff_Name ); FitFileName.Append( "_Fit.pdf" );
		thisCanv->Print( FitFileName );

		string DiffFileName = "PerEventOutput/Abs_"; DiffFileName.append( PullDiff_Name ); DiffFileName.append(".pdf");
		TString ThisCanv2("ThisCanv"); ThisCanv2+=name_i; ThisCanv2.Append("_2");
		TCanvas* c1 = new TCanvas( ThisCanv2, ThisCanv2 );
		c1->cd();

		TH1* thisTH1_2 = Histogram_Processing::Plot_1D( AbsDiffData, DiffFileName, "PE9" );
		c1->cd();
		thisTH1_2->Draw( "PE9" );
		c1->Update();
		thisTH1_2->GetYaxis()->SetTitle( "# of Results" );
		X_Title = allNames[name_i]; X_Title.Append( " Absolute Difference" );
		thisTH1_2->GetXaxis()->SetTitle( X_Title );
		c1->Update();

		TString FuncName2="newfunc2_"; FuncName2+=name_i;
		TF1* newFunc2 = new TF1( FuncName2, "[0]*exp(-((x-[1])*(x-[1])/[2])) + [3]*exp(-((x-[1])*(x-[1])/[4]))" );
                newFunc2->SetParameter( 0, 200.);
                //newFunc2->SetParameter( 1, 0.05);
                newFunc2->SetParameter( 1, thisTH1_2->GetBinCenter( thisTH1_2->GetMaximumBin()) );
                newFunc2->SetParameter( 2, 0.3);
                newFunc2->SetParameter( 3, 1000.);
                //newFunc2->SetParameter( 4, -0.1);
                //newFunc2->SetParameter( 5, 50.);
                newFunc2->SetParameter( 4, 0.01 );

		thisTH1_2->Fit( newFunc2, "M" );
		c1->Update();

		c1->Print( DiffFileName.c_str() );
	}

	/*
	TH1* PhisOneGraph = Histogram_Processing::Get_TH1( PhisOneDiff, NULL, 100, -0.5, 0.5 );
	TH1* PhisTwoGraph = Histogram_Processing::Get_TH1( PhisTwoDiff, NULL, 100, -0.5, 0.5 );

	TH1* PhisDiffGraph = Histogram_Processing::Get_TH1( PhisDiff, NULL, 100, -0.5, 0.5 );

	TH1* PhisOnePullGraph = Histogram_Processing::Get_TH1( PhisOnePull, NULL, 100, -5, 5 );
	TH1* PhisTwoPullGraph = Histogram_Processing::Get_TH1( PhisTwoPull, NULL, 100, -5, 5 );

	TH1* PhisPullGraph = Histogram_Processing::Get_TH1( PhisPull, NULL, 100, -0.75, 0.75 );
	TH1* PhisPropPullGraph = Histogram_Processing::Get_TH1( PhisPropPull, NULL, 100, -0.75, 0.75 );

	TH1* PhisPullDiffGraph = Histogram_Processing::Get_TH1( PhisPullDiff, NULL, 100, -5, 5 );

	TH1* PhisPullDiffDisk = Histogram_Processing::Get_TH1( PhisPullDiff, NULL, 10000, -5, 5 );

	TCanvas* c1 = new TCanvas( "DiffOneGraph", "DiffOneGraph" );
	PhisOneGraph->Draw("EP9L");
	c1->Update(); c1->Print("DiffPhisOne.pdf");
	TCanvas* c2 = new TCanvas( "DiffTwoGraph", "DiffTwoGraph" );
	PhisTwoGraph->Draw("EP9L");
	c2->Update(); c2->Print("DiffPhisTwo.pdf");
	
	TCanvas* c3 = new TCanvas( "DiffGraph", "DiffGraph" );
	PhisDiffGraph->Draw("EP9L");
	c3->Update(); c3->Print("DiffGraph.pdf");

	TCanvas* c4 = new TCanvas( "PullOneGraph", "PullOneGraph" );
	PhisOnePullGraph->Draw("EP9L");
	c4->Update(); c4->Print("PullPhisOne.pdf");
	TCanvas* c5 = new TCanvas( "PullTwoGraph", "PullTwoGraph" );
	PhisTwoPullGraph->Draw("EP9L");
	c5->Update(); c5->Print("PullPhisTwo.pdf");

	TCanvas* c6 = new TCanvas( "PhisPullGraph", "PhisPullGraph" );
	PhisPullGraph->Draw("EP9L");
	c6->Update(); c6->Print("PhisPullGraph.pdf");
	TCanvas* c7 = new TCanvas( "PhisPropPullGraph", "PhisPropPullGraph" );
	PhisPropPullGraph->Draw("EP9L");
	c7->Update(); c7->Print("PhisPropPullGraph.pdf");

	TCanvas* c8 = new TCanvas( "PhisPullDiff", "PhisPullDiff" );
	PhisPullDiffGraph->Draw("EP9L");
	c8->Update(); c8->Print( "PhisPullDiff.pdf" );

	TFile* output = new TFile( "OFFSetStudy.root", "RECREATE" );
	PhisPullDiffDisk->Write("",TObject::kOverwrite);
	output->Write( "",TObject::kOverwrite);
	*/

	return 0;
}


