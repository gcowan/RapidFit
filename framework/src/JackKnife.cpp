/**
  @namespace JackKnife 

  Namespace allowing jackknifing procedure to be carried out 

  @author Greig A Cowan greig.cowan@cern.ch
  @date 2012-05-25
 */

//	ROOT Headers
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TFile.h"
//	RapidFit Headers
#include "JackKnife.h"
#include "MemoryDataSet.h"
#include "ObservableContinuousConstraint.h"
#include "RapidFitIntegrator.h"
#include "XMLConfigReader.h"
#include "PDFWithData.h"
#include "IPDF.h"
#include "IDataSet.h"
#include "PhaseSpaceBoundary.h"
#include "EdStyle.h"
#include "DataSetConfiguration.h"
#include "FitResult.h"
#include "FitAssembler.h"
#include "ResultFormatter.h"
//	System Headers
#include <math.h>
#include <iostream>
#include <set>
#include <algorithm>
#include <time.h>

using namespace::std;

namespace JackKnife
{
	void jackknife( XMLConfigReader * xmlFile, MinimiserConfiguration * theMinimiser, 
		FitFunctionConfiguration * theFunction, ParameterSet* argumentParameterSet, vector<string> CommandLineParam, int start, int stop )
	{
		cout << "Starting JackKnife" << endl;
		PDFWithData  * pdfAndData  = xmlFile->GetPDFsAndData()[0];
		ParameterSet * parSet      = xmlFile->GetFitParameters( CommandLineParam );
		PhaseSpaceBoundary * phase = xmlFile->GetPhaseSpaceBoundaries()[0];
		pdfAndData->SetPhysicsParameters( parSet );
			
		MemoryDataSet * dataset = (MemoryDataSet*) pdfAndData->GetDataSet();
		vector<IDataSet*> data;
		data.push_back(dataset);
	
		IPDF * pdf = pdfAndData->GetPDF();
		vector<IPDF*> vectorPDFs;
		vectorPDFs.push_back(pdf);
	
		// Repeat nominal fit to get the nominal value of the physics parameter	
		FitResult * nominal = FitAssembler::DoFit( theMinimiser, theFunction, argumentParameterSet, vectorPDFs, data, xmlFile->GetConstraints() );
		double nominal_value = nominal->GetResultParameterSet()->GetResultParameter("tau")->GetValue();
		cout << nominal_value << endl;
		double jackknifed_value = 0.;

		MemoryDataSet * subset = new MemoryDataSet( phase );
		int nData = dataset->Yield();	
                string fileName = ResultFormatter::GetOutputFolder();
                fileName.append("/jackknife.root");         
		TFile * outputFile = new TFile("jackknife.root", "RECREATE");	
                TNtuple * jack = new TNtuple("jack", "jackknifed - nominal", "diff_jackknifed_nominal:reco_time:true_time");

		double reco_time = 0.;
		double true_time = 0.;

		for ( int i = start; i < stop; i++ )
		{
			for ( int j = 0; j < nData; j++ )
			{
				if ( j == i ) {
					reco_time = dataset->GetDataPoint( j )->GetObservable( "time" )->GetValue();
					true_time = dataset->GetDataPoint( j )->GetObservable( "truetime" )->GetValue();
					continue;
				}
				subset->AddDataPoint( dataset->GetDataPoint( j ) );
			}
			cout << "Creating subset " << i << " containing " << subset->Yield() << " candidates" << endl;
			vector<IDataSet*> vectorData;
			vectorData.push_back(subset);
			FitResult * result = FitAssembler::DoFit( theMinimiser, theFunction, argumentParameterSet, vectorPDFs, vectorData, xmlFile->GetConstraints() );
			if ( result->GetFitStatus() == 3 ) {
				jackknifed_value = result->GetResultParameterSet()->GetResultParameter("tau")->GetValue();
				jack->Fill( (Float_t)(jackknifed_value - nominal_value), (Float_t)reco_time, (Float_t)true_time);
			}
			subset->Clear();
			delete result;
		}
		outputFile->Write();
		outputFile->Close();
	}

	void plotUstatistic( IPDF * pdf, IDataSet * data, PhaseSpaceBoundary * phase, string plot )
	{
		(void) pdf; (void) phase;
		unsigned int nD = (unsigned) data->GetDataNumber();
		int bins = 30;
		double level = (double)nD/(double)bins;
		TH1D * distances = new TH1D( "distances", "U", bins, 0., 1.);
		distances->Sumw2();

		char buff[10];
		sprintf( buff, "%f", level );
		TF1 * line = new TF1("line", buff, 0, 1);
		line->SetLineColor(kBlue);

		TCanvas* ca = new TCanvas("ca","canvas", 3000, 3000);
		gStyle->SetOptStat(0);
		distances->SetTitle("");
		distances->GetXaxis()->SetNdivisions(508);
		distances->GetXaxis()->SetTitle("U");
		distances->Draw();
		line->Draw("same");
		ca->Update();
		ca->SaveAs(plot.c_str());
		//delete distances;
		//delete ca;
	}

}
