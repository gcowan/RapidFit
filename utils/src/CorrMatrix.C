
#include "TTree.h"
#include "TH1.h"
#include "TMatrixDSym.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"

#include "Histo_Processing.h"
#include "CorrMatrix.h"

#include "EdStyle.h"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using namespace::std;

void CorrMatrix::Analyse( const vector<TTree*> corr_trees )
{
	cout << "Processing Correlation Matricies" << endl;
	gSystem->mkdir( "MatrixPlots" );
	gSystem->cd( "MatrixPlots" );
	gStyle->SetPalette(1);
	gStyle->SetPadRightMargin( (Float_t)0.15 );
	gStyle->SetOptStat(0);
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
				TCanvas* c1 = new TCanvas( Name, Name, 1680, 1050 );
				newMatrix->Draw("colz");
				c1->Update();
				TH1* matrixPlot = (TH1*) c1->GetPrimitive("TMatrixDBase");
				matrixPlot->SetTitle("Correlation Matrix");
				c1->SetTitle("Correlation Matrix");
				for( unsigned int k=0; k< thisNames->size(); ++k )
				{
					matrixPlot->GetXaxis()->SetBinLabel( k+1, (*thisNames)[k].c_str() );
					matrixPlot->GetYaxis()->SetBinLabel( k+1, (*thisNames)[k].c_str() );
				}
				matrixPlot->GetZaxis()->SetRangeUser( 0., 1. );
				c1->Update();
				Histogram_Processing::Silent_Print( c1 , Name+".pdf" );
			}
		}
		cout << endl;
	}
	unsigned int dim = allData[0].size();
	vector<vector<TH1*> > allHistos;
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
	gStyle->SetPadRightMargin( (Float_t)0.08 );
	gSystem->cd("..");

	EdStyle::SetStyle();
}

