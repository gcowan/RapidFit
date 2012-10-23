
#include "TTree.h"
#include "TH1.h"
#include "TMatrixDSym.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "CorrMatrix.h"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using namespace::std;

void CorrMatrix::Analyse( const vector<TTree*> corr_trees )
{
	gStyle->SetPalette(1);
	gStyle->SetPadRightMargin( (Float_t)0.15 );
	vector<double>* thisMatrix = new vector<double>();
	vector<string>* thisNames = new vector<string>();
	cout << "Files with Correlations: " << corr_trees.size() << endl;
	for( unsigned int i=0; i< corr_trees.size(); ++i )
	{      
		if( corr_trees[i] != NULL )
		{
			cout << "File: " << i+1 << endl;
			corr_trees[i]->GetEntries();
			corr_trees[i]->SetBranchAddress( "MartrixElements", &thisMatrix );
			corr_trees[i]->SetBranchAddress( "MartrixNames", &thisNames );
			for( unsigned int j=0; j< corr_trees[i]->GetEntries(); ++j )
			{
				corr_trees[i]->GetEntry(j);
				cout << "Matrix Elements: " << thisMatrix->size() << endl;
				for( unsigned int k=0; k< thisNames->size(); ++k )
				{
					cout << (*thisNames)[k] << endl;
				}
				TMatrixDSym* newMatrix = new TMatrixDSym( (int)sqrt(thisMatrix->size()), &((*thisMatrix)[0]) );
				TString Name="Matrix_"; Name+=j;
				TCanvas* c1 = new TCanvas( Name, Name, 1680, 1050 );
				newMatrix->Draw("colz");
				c1->Update();
				TH1* matrixPlot = (TH1*) c1->GetPrimitive("TMatrixDBase");
				for( unsigned int k=0; k< thisNames->size(); ++k )
				{
					matrixPlot->GetXaxis()->SetBinLabel( k+1, (*thisNames)[k].c_str() );
					matrixPlot->GetYaxis()->SetBinLabel( k+1, (*thisNames)[k].c_str() );
				}
				c1->Update();
				c1->Print(Name+".pdf");
			}
		}
		cout << endl;
	}      
	gStyle->SetPadRightMargin( (Float_t)0.08 );
}

