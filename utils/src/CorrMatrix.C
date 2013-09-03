
#include "TROOT.h"
#include "TTree.h"
#include "TH1.h"
#include "TMatrixDSym.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPaletteAxis.h"

#include "StringOperations.h"
#include "Histo_Processing.h"
#include "CorrMatrix.h"

#include "EdStyle.h"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using namespace::std;

void CorrMatrix::Help()
{
	cout << endl;
	cout << "Correlation Matricies are post-processed from a RapidFit file ONLY when requested." << endl;
	cout << "To request for Correlation Matricies to be plotted you should pass the runtime option of:" << endl;
	cout << endl << "--CorrMatrix" << "\t\t" << "This allows the Correlation Matrix code to be run over the file and the output to be post-processed and plotted" << endl;
	cout << endl;
}

void CorrMatrix::Analyse( const vector<TTree*> corr_trees, const vector<string> other_params )
{
	cout << "Processing Correlation Matricies" << endl;
	if( gDirectory->GetDirectory( "MatrixPlots" ) == 0 ) gSystem->mkdir( "MatrixPlots" );
	gSystem->cd( "MatrixPlots" );

	gStyle->SetPalette(1);
	if( StringOperations::VectorContains( other_params, "--noLegend" ) == -1 )
	{
		gStyle->SetPadRightMargin( (Float_t)0.2375 );
	}
	else
	{
		gStyle->SetPadRightMargin( (Float_t)0.15 );
	}
	gStyle->SetOptStat(0);
	Double_t lhcbTSize = 0.06;
	//gStyle->SetTitleSize( (Float_t)0.8*(Float_t)lhcbTSize, "xyz" );

	gROOT->UseCurrentStyle();
	gROOT->ForceStyle( true );

	vector<double>* thisMatrix = new vector<double>();
	vector<string>* thisNames = new vector<string>();

	vector<vector<vector<double> > > allData;

	cout << "Files with Correlations: " << corr_trees.size() << endl;
	for( unsigned int i=0; i< corr_trees.size(); ++i )
	{
		if( corr_trees[i] != NULL )
		{
			//cout << "File: " << i+1 << endl;
			corr_trees[i]->GetEntries();
			corr_trees[i]->SetBranchAddress( "MartrixElements", &thisMatrix );
			corr_trees[i]->SetBranchAddress( "MartrixNames", &thisNames );
			for( unsigned int j=0; j< corr_trees[i]->GetEntries(); ++j )
			{
				corr_trees[i]->GetEntry(j);
				//cout << "Matrix Elements: " << thisMatrix->size() << endl;
				//for( unsigned int k=0; k< thisNames->size(); ++k )
				//{
				//	cout << (*thisNames)[k] << endl;
				//}
				int dim = (int)sqrt(thisMatrix->size());
				//cout << dim << " x " << dim << endl;
				TMatrixDSym* newMatrix = new TMatrixDSym( dim );
				vector<vector<double> > thisMatrixData( dim, vector<double>(dim, 0.) );
				for( unsigned int k1=0; k1 < dim; ++k1 )
				{
					double drow = (*thisMatrix)[ k1*dim + k1 ];
					for( unsigned int k2=0; k2< dim; ++k2 )
					{
						double dcol = (*thisMatrix)[ k2*dim + k2 ];
						double covariance = (*thisMatrix)[ k1*dim + k2 ];
						//cout << "k1: " << k1 << " k2: " << k2 << "\t" << covariance / sqrt(fabs(drow * dcol)) << endl;
						double thisVal = covariance / sqrt(fabs(drow * dcol));
						(*newMatrix)( k1, k2 ) = thisVal;
						thisMatrixData[k1][k2] = thisVal;
					}
				}
				allData.push_back( thisMatrixData );
				TString Name="Matrix_"; Name+=j;
				TCanvas* c1 = EdStyle::RapidFitCanvas( Name, Name );
				newMatrix->Draw("colz");
				c1->Update();
				TH1* matrixPlot = (TH1*) c1->GetPrimitive("TMatrixDBase")->Clone("newMatrixPlot");
				matrixPlot->SetTitle("Correlation Matrix");
				c1->SetTitle("Correlation Matrix");

				if( StringOperations::VectorContains( other_params, "--noLegend" ) != -1 )
				{
					for( unsigned int k=0; k< thisNames->size(); ++k )
					{
						matrixPlot->GetXaxis()->SetBinLabel( k+1, EdStyle::GetParamRootName( TString((*thisNames)[k].c_str()) ) );
						matrixPlot->GetYaxis()->SetBinLabel( k+1, EdStyle::GetParamRootName( TString((*thisNames)[k].c_str()) ) );
					}
				}
				else
				{
					TLegend* thisLegend = new TLegend( 0.875, 0., 1., 1. );
					thisLegend->SetFillColor( kWhite );
					thisLegend->SetFillStyle( 3001 );
					for( unsigned int k=0; k< thisNames->size(); ++k )
					{
						TString num; num+=k;
						matrixPlot->GetXaxis()->SetBinLabel( k+1, num );
						matrixPlot->GetYaxis()->SetBinLabel( k+1, num );

						TString thisEntry = num;
						thisEntry.Append(" "); thisEntry.Append( EdStyle::GetParamRootName( TString((*thisNames)[k].c_str()) ) );
						thisLegend->AddEntry((TObject*)0, thisEntry, "");
					}
					thisLegend->Draw("SAME");
				}

				if( StringOperations::VectorContains( other_params, "--noLegend" ) != -1 )
				{
					matrixPlot->GetXaxis()->SetLabelSize( (Float_t)0.6*(Float_t)lhcbTSize );
					matrixPlot->GetYaxis()->SetLabelSize( (Float_t)0.6*(Float_t)lhcbTSize );
					matrixPlot->GetZaxis()->SetLabelSize( (Float_t)0.6*lhcbTSize );
				}
				else
				{
					matrixPlot->GetXaxis()->SetLabelSize( (Float_t)0.8*(Float_t)lhcbTSize );
					matrixPlot->GetYaxis()->SetLabelSize( (Float_t)0.8*(Float_t)lhcbTSize );
					matrixPlot->GetZaxis()->SetLabelSize( (Float_t)0.8*lhcbTSize );
				}
				matrixPlot->LabelsOption("v");
				matrixPlot->GetZaxis()->SetRangeUser( -1., 1. );
				matrixPlot->Draw("colz");
				c1->Update();
				// Because this is obvious!...
				TPaletteAxis *palette = (TPaletteAxis*)matrixPlot->GetListOfFunctions()->FindObject("palette");
				palette->SetLabelSize( matrixPlot->GetZaxis()->GetLabelSize() );
				c1->Update();
				Histogram_Processing::Silent_Print( c1 , Name+".pdf" );
				Histogram_Processing::Silent_Print( c1 , Name+".C" );
			}
		}
		cout << endl;
	}

	unsigned int dim = allData[0].size();
	vector<vector<TH1*> > allHistos;
	if( corr_trees[0]->GetEntries() > 5 )
	{
		for( unsigned int i=0; i< dim; ++i )
		{
			for( unsigned int j=0; i< dim; ++j )
			{
				vector<double> thisElement;
				for( unsigned int k=0; k< allData.size(); ++k )
				{
					thisElement.push_back( allData[k][i][j] );
				}
				TH1* thisHisto = Histogram_Processing::Get_TH1( thisElement );
			}
		}
	}
	gStyle->SetPadRightMargin( (Float_t)0.08 );
	gSystem->cd("..");

	EdStyle::SetStyle();
}

