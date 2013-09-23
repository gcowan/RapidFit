/**
  @namespace GoodnessOfFit

  Namespace holding some common GOF tools

  @author Greig A Cowan greig.cowan@cern.ch
  @date 2011-04-06
 */

//	ROOT Headers
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1D.h"
#include "TFile.h"
#include "TVectorD.h"
//	RapidFit Headers
#include "GoodnessOfFit.h"
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
#include <iomanip>

using namespace::std;

namespace GoodnessOfFit
{

	double gofLoop( XMLConfigReader * xmlFile, MinimiserConfiguration * theMinimiser, FitFunctionConfiguration * theFunction, ParameterSet* argumentParameterSet, vector<string> CommandLineParam, int nData )
	{
		cout << "Starting GOF" << endl;

		double pvalue = 0.;

		for ( int i = 0; i < 1; i ++ )
		{
			cout << "Iteration " << i << endl;
			// Set up the PDF parameters
			ParameterSet* parSet = xmlFile->GetFitParameters( CommandLineParam );

			// First, generate our toy data, fit to it, use PDF to generate large MC and then calculate corresponding p-value from permutation
			vector<double> dataPvalue;
			GoodnessOfFit::generateFitAndCalculatePvalue( xmlFile, parSet, theMinimiser, theFunction, argumentParameterSet, nData, 1, &dataPvalue );

			// Now need to account for any potential fit bias so repeat above procedure N times to get distribution of p-values
			vector<double> distOfPvalues;
			GoodnessOfFit::generateFitAndCalculatePvalue( xmlFile, parSet, theMinimiser, theFunction, argumentParameterSet, nData, 50, &distOfPvalues );

			// Finally, compare dataPvalue with distribution of pvalues to get the actual p-value of the fit, which correctly accounts for the bias
			pvalue = GoodnessOfFit::getPvalue( dataPvalue[0], distOfPvalues );
		}
		return pvalue;
	}

	double getPvalue( double datavalue, vector<double> distribution )
	{
		int count = 0;
		vector<double>::iterator iter;
		for ( iter = distribution.begin(); iter != distribution.end(); ++iter ){
			//cout << datavalue << " " << *iter << endl;
			if ( datavalue < *iter ) count++;
		}
		double pvalue = count/double(distribution.size());
		return pvalue;
	}

	double fitDataCalculatePvalue( XMLConfigReader * xmlFile, MinimiserConfiguration * theMinimiser, FitFunctionConfiguration * theFunction, ParameterSet* argumentParameterSet, FitResult * result)
	{
		//	Unused parameters, keep gcc happy
		(void) theMinimiser; (void) theFunction; (void) argumentParameterSet;
		cout << "Starting GOF for data" << endl;
		vector<PDFWithData *> pdfAndData = xmlFile->GetPDFsAndData();
		vector<PDFWithData *>::iterator iter;
		vector<IDataSet*> data;

		// Take the fitted parameters
		ParameterSet* parSetFromFit;
		parSetFromFit = result->GetResultParameterSet()->GetDummyParameterSet();

		IPDF * pdf = NULL;
		PhaseSpaceBoundary * phase = NULL;
		DataSetConfiguration * dataConfig = NULL;
		MemoryDataSet * mc = NULL;
		unsigned int ulTime = 0;
		unsigned int nData = 0;
		vector<IDataSet*> mcData;

		// Generate large sample of MC using the fitted PDF parameters
		for ( iter = pdfAndData.begin(); iter != pdfAndData.end(); ++iter )
		{
			data.push_back( (*iter)->GetDataSet() );
			nData = (unsigned)data.back()->GetDataNumber();
			phase = data.back()->GetBoundary();

			(*iter)->SetPhysicsParameters( parSetFromFit );
			ulTime = static_cast<unsigned int>( time( NULL ));
			pdf = (*iter)->GetPDF();
			pdf->SetRandomFunction( (int)ulTime );
			pdf->SetMCCacheStatus( false );

			dataConfig = (*iter)->GetDataSetConfig();
			dataConfig->SetSource( "Foam" );

			mc = (MemoryDataSet*)dataConfig->MakeDataSet( phase, pdf, 5*(int)nData );
			mcData.push_back( mc );
		}

		// Finally calculate the p-value relative to the large MC sample
		double pvalue = GoodnessOfFit::pValueFromPoint2PointDissimilarity( data[0], mcData[0] );
		cout << "p-value " << pvalue << endl;
		return pvalue;
	}

	void generateFitAndCalculatePvalue( XMLConfigReader * xmlFile, ParameterSet* parSet, MinimiserConfiguration * theMinimiser, FitFunctionConfiguration * theFunction, ParameterSet* argumentParameterSet, int nData, int repeats, vector<double> * pvalues)
	{
		PDFWithData * pdfAndData = xmlFile->GetPDFsAndData()[0];
		pdfAndData->SetPhysicsParameters( parSet );
		IPDF * pdf = pdfAndData->GetPDF();
		vector<IPDF*> vectorPDFs;
		vectorPDFs.push_back(pdf);
		PhaseSpaceBoundary * phase = xmlFile->GetPhaseSpaceBoundaries()[0];
		EdStyle * greigFormat = new EdStyle();
		greigFormat->SetStyle();

		ParameterSet* parSetFromFit = NULL;

		// Set up to be able to generate some MC data
		DataSetConfiguration * dataConfig = pdfAndData->GetDataSetConfig();
		dataConfig->SetSource( "Foam" );
		bool model = false;
		double pvalue = 0.;
		unsigned int ulTime = static_cast<unsigned int>( time( NULL ));
		for ( int i = 0; i < repeats; )
		{
			cout << "Ensemble " << i << endl;

			// First, generate some data from this PDF
			pdfAndData->SetPhysicsParameters( parSet );
			pdf->SetRandomFunction( (int)ulTime );
			pdf->SetMCCacheStatus( false );
			MemoryDataSet * subset = (MemoryDataSet*)dataConfig->MakeDataSet( phase, pdf, nData );
			vector<IDataSet*> vectorData;
			vectorData.push_back(subset);

			if ( model ) {
				// Use the same PDF for the large MC dataset (the Model approach)
				pdfAndData->SetPhysicsParameters( parSet );
			}
			else {
				// Fit the generated data (the Fit I approach)
				FitResult * gofResult = FitAssembler::DoFit( theMinimiser, theFunction, argumentParameterSet, vectorPDFs, vectorData, xmlFile->GetConstraints() );
				if ( gofResult->GetFitStatus() == 3 ) {
					delete parSetFromFit;
					parSetFromFit = gofResult->GetResultParameterSet()->GetDummyParameterSet();
					pdfAndData->SetPhysicsParameters( parSetFromFit );
					i++; // let us go to the next iteration of the loop so that we always get "repeats" good results
				}
			}
			// Generate large sample of MC
			ulTime = static_cast<unsigned int>( time( NULL ));
			pdf->SetRandomFunction( (int)ulTime );
			pdf->SetMCCacheStatus( false );
			MemoryDataSet * mcData = (MemoryDataSet*)dataConfig->MakeDataSet( phase, pdf, 10*nData );

			// Finally calculate the p-value relative to the large MC sample
			pvalue = GoodnessOfFit::pValueFromPoint2PointDissimilarity( subset, mcData );
			//pvalue = GoodnessOfFit::pValueFromPoint2PointDissimilarity( data, mcData );
			pvalues->push_back(pvalue);
			//delete subset;
			//delete mcData;
		}
		parSet = parSetFromFit;  // Need this so that we have the correct parameters during the second call of this function
	}

	void plotUstatistic( IPDF * pdf, IDataSet * data, PhaseSpaceBoundary * phase, string plot )
	{
		unsigned int nD = (unsigned) data->GetDataNumber();
		int bins = 30;
		double level = (double)nD/(double)bins;
		TH1D * distances = new TH1D( "distances", "U", bins, 0., 1.);
		distances->Sumw2();
		calculateUstatisticNum(pdf, data, phase, distances);

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

	void calculateUstatistic( IPDF * pdf, IDataSet * data, PhaseSpaceBoundary * phase, TH1D * distances)
	{
		double pdfValue = 0.;
		double distance = 0.;
		double smallest_distance = 0., sd = 0.;
		double U = 0;
		DataPoint * event_i = 0;
		DataPoint * event_j = 0;
		DataPoint * closest = 0;
		bool firstEvent = true;
		size_t dimension = (pdf->GetPrototypeDataPoint()).size();
		if ( dimension == 7 ) dimension = 4;
		dimension = 2;
		cout << "number of dimensions: " << dimension << endl;
		unsigned int nD = (unsigned)data->GetDataNumber();
		int count = 0;
		for (unsigned int i = 0; i < nD; i++) {
			if ( !(i % 100) ) cout << i << endl;
			event_i = data->GetDataPoint( (int)i );
			pdfValue = pdf->Evaluate( event_i )/pdf->Integral( event_i, phase );
			firstEvent = true;
			for (unsigned int j = 0; j < nD; j++) {
				if ( j == i ) continue;
				event_j = data->GetDataPoint( (int)j );
				distance = getDistance( event_i, event_j );
				if ( firstEvent ) {
					smallest_distance = distance;
					firstEvent = false;
				}
				if (distance < smallest_distance) {
					smallest_distance = distance;
					closest = event_j;
				}
			}
			sd = smallest_distance;
			vector<double> vectorOfDistances = getDistances( event_i, closest );
			PhaseSpaceBoundary * tempPhase = new PhaseSpaceBoundary( phase->GetAllNames() );
			copyPhaseSpaceBoundary( tempPhase, phase );
			updatePhaseSpaceBoundary( event_i, tempPhase, phase, vectorOfDistances );
			if ( vectorOfDistances[0] > vectorOfDistances[1] ) count++;

			//double distToCosThetaBoundary = 1. - fabs(event_i->GetObservable( "cosTheta" )->GetValue()) ;
			//double distToCosPsiBoundary = 1. - fabs(event_i->GetObservable( "cosPsi" )->GetValue()) ;
			//double distToTimeBoundary = event_i->GetObservable( "time" )->GetValue() - 0.3;
			//double distToMassBoundary = min(fabs(event_i->GetObservable( "mass" )->GetValue() - 5200.),fabs(event_i->GetObservable( "mass" )->GetValue() - 5550.)) ;
			//(void) distToTimeBoundary; (void) distToMassBoundary;	//Unused
			double f1 = 1.;
			double f2 = 1.;
			//if( sd > distToCosThetaBoundary ) f1 = 1. - acos(distToCosThetaBoundary/sd)/TMath::Pi() ;
			//if( sd > distToCosPsiBoundary ) f2 = 1. - acos(distToCosPsiBoundary/sd)/TMath::Pi() ;
			//if( sd > distToTimeBoundary ) f1 = 1. - acos(distToTimeBoundary/sd)/TMath::Pi() ;
			//if( sd > distToMassBoundary ) f2 = 1. - acos(distToMassBoundary/sd)/TMath::Pi() ;

			if ( dimension == 1 ) U = exp(-1.*nD*2.                            *     sd     * pdfValue * f1 );
			if ( dimension == 2 ) U = exp(-1.*nD         *     TMath::Pi()     * pow(sd, 2) * pdfValue * f1 * f2 );
			if ( dimension == 3 ) U = exp(-1.*nD*4./3.   *     TMath::Pi()     * pow(sd, 3) * pdfValue );
			if ( dimension == 4 ) U = exp(-1.*nD*1./2.   *pow( TMath::Pi(), 2) * pow(sd, 4) * pdfValue * f1 * f2 );
			if ( dimension == 5 ) U = exp(-1.*nD*8./15.  *pow( TMath::Pi(), 2) * pow(sd, 5) * pdfValue );
			if ( dimension == 6 ) U = exp(-1.*nD*1./6.   *pow( TMath::Pi(), 3) * pow(sd, 6) * pdfValue );
			if ( dimension == 7 ) U = exp(-1.*nD*16./105.*pow( TMath::Pi(), 3) * pow(sd, 7) * pdfValue );
			/*if ( U < 1.e-02 ) {
			  cout << "pdfValue " << pdfValue << " smallest " << sd <<  " U " << U << " time " << event_i->GetObservable( "time" )->GetValue()
			  << " mass " << event_i->GetObservable( "mass" )->GetValue()
			  << " cosTheta " << event_i->GetObservable( "cosTheta" )->GetValue()
			  << " phi " << event_i->GetObservable( "phi" )->GetValue()
			  << " cosPsi " << event_i->GetObservable( "cosPsi" )->GetValue()  << endl;
			  }*/
			distances->Fill( U );
		}
		cout << "count " << count << endl;
	}

	void calculateUstatisticNum( IPDF * pdf, IDataSet * data, PhaseSpaceBoundary * phase, TH1D * distances)
	{
                double pdfValue = 0.;
		double distance = 0.;
		double smallest_distance = 0., sd = 0.;
		double U = 0;
		DataPoint * event_i = 0;
		DataPoint * event_j = 0;
		DataPoint * closest = 0;
		bool firstEvent = true;
		size_t dimension = (pdf->GetPrototypeDataPoint()).size();
		if ( dimension == 7 ) dimension = 5;
		cout << "number of dimensions: " << dimension << endl;
		unsigned int nD = (unsigned) data->GetDataNumber();
		RapidFitIntegrator * integrator = new RapidFitIntegrator( pdf, false, false );
		vector<string> doNotIntegrate;
		doNotIntegrate.push_back("fsig_sw");
		double volume = integrator->NumericallyIntegratePhaseSpace( phase, doNotIntegrate );
		delete integrator;
		for (unsigned int i = 0; i < nD; i++) {
			if ( !(i % 100) ) cout << i << endl;
			event_i = data->GetDataPoint( (int)i );
                        pdfValue = pdf->Evaluate( event_i );
			firstEvent = true;
			for (unsigned int j = 0; j < nD; j++) {
				if ( j == i ) continue;
				event_j = data->GetDataPoint( (int)j );
				distance = getDistance( event_i, event_j );
				if ( firstEvent ) {
					smallest_distance = distance;
					firstEvent = false;
				}
				if (distance < smallest_distance)
				{
					smallest_distance = distance;
					closest = event_j;
				}
			}
			sd = smallest_distance;
			if ( dimension == 1 ) U = exp(-1.*nD*2.                            *     sd     * pdfValue/volume );
                        if ( dimension == 2 ) U = exp(-1.*nD         *     TMath::Pi()     * pow(sd, 2) * pdfValue/volume );
                        if ( dimension == 3 ) U = exp(-1.*nD*4./3.   *     TMath::Pi()     * pow(sd, 3) * pdfValue/volume );
                        if ( dimension == 4 ) U = exp(-1.*nD*1./2.   *pow( TMath::Pi(), 2) * pow(sd, 4) * pdfValue/volume );
                        if ( dimension == 5 ) U = exp(-1.*nD*8./15.  *pow( TMath::Pi(), 2) * pow(sd, 5) * pdfValue/volume );
                        if ( dimension == 6 ) U = exp(-1.*nD*1./6.   *pow( TMath::Pi(), 3) * pow(sd, 6) * pdfValue/volume );
                        if ( dimension == 7 ) U = exp(-1.*nD*16./105.*pow( TMath::Pi(), 3) * pow(sd, 7) * pdfValue/volume );
			//cout << sd << " " << pdfValue << " " << volume << endl;
			distances->Fill( U );
		}
	}

	void copyPhaseSpaceBoundary( PhaseSpaceBoundary * newBoundary, PhaseSpaceBoundary * oldBoundary )
	{
		vector<string> names = oldBoundary->GetAllNames();
		vector<string>::iterator nameIter = names.begin();
		vector<double> values;
		double min, max;
		string unit;
		for ( ; nameIter != names.end(); ++nameIter )
		{
			IConstraint * constraint = oldBoundary->GetConstraint( *nameIter );
			unit = constraint->GetUnit();
			if ( constraint->IsDiscrete() )
			{
				values = constraint->GetValues();
				newBoundary->SetConstraint( *nameIter, values, unit );
			}
			else {
				min = constraint->GetMinimum();
				max = constraint->GetMaximum();
				newBoundary->SetConstraint( *nameIter, min, max, unit );
			}
		}
	}

	vector<double> getDistances( DataPoint * x, DataPoint * y )
	{
		vector<double> distances;
		vector<string> xListOfNames = x->GetAllNames();
		vector<string> yListOfNames = y->GetAllNames();
		vector<string>::iterator xIter = xListOfNames.begin();
		vector<string>::iterator yIter = yListOfNames.begin();
		Observable * xVar = 0;
		Observable * yVar = 0;
		while ( (xVar = x->GetObservable( *xIter ) ) && (yVar = y->GetObservable( *yIter ) ) ) {
			if ( (*xIter == "time" ) || (*xIter == "mass" )) {
				//cout << fabs(xVar->GetValue() - yVar->GetValue() ) << endl;
				distances.push_back( fabs(xVar->GetValue() - yVar->GetValue() ) );
			}
			else if ( (*xIter == "cosTheta") ) {
				distances.push_back( fabs(xVar->GetValue() - yVar->GetValue()) );
			}
			/*
			   else if(  (*xIter == "cosPsi") ) {
			   distances.push_back( fabs( xVar->GetValue() - yVar->GetValue() ) );
			   }
			   else if ( (*xIter == "phi") ) {
			   double diff = diffmod2pi( xVar->GetValue() - yVar->GetValue() );
			   distances.push_back( diff );
			   }
			 */
			++xIter, ++yIter;
			if ( xIter == xListOfNames.end() || yIter == yListOfNames.end() ) break;
		}
		return distances;
	}

	void updatePhaseSpaceBoundary( DataPoint * x, PhaseSpaceBoundary * newPhase, PhaseSpaceBoundary * oldPhase, vector<double> distances )
	{
		vector<string> xListOfNames = x->GetAllNames();
		vector<string>::iterator xIter = xListOfNames.begin();
		Observable * xVar = 0;
		double min = 0., max = 0.;
		int i = 0;
		char buffer[100];	(void) buffer;//	Unused

		while ( ( xVar = x->GetObservable( *xIter ) ) ) {
			if ( *xIter == "time" || *xIter == "cosTheta" || *xIter == "mass"){//|| *xIter == "cosPsi" || *xIter == "phi") {
				min = xVar->GetValue() - distances[ (unsigned)i ];
				max = xVar->GetValue() + distances[ (unsigned)i ];
				if ( min > oldPhase->GetConstraint( *xIter )->GetMinimum() ) ((ObservableContinuousConstraint *)newPhase->GetConstraint( *xIter ))->SetMinimum( min );
				if ( max < oldPhase->GetConstraint( *xIter )->GetMaximum() ) ((ObservableContinuousConstraint *)newPhase->GetConstraint( *xIter ))->SetMaximum( max );
				//cout << *xIter << " " << buffer << endl;
			}
			++xIter;
			if ( xIter == xListOfNames.end() ) break;
			i++;
			}
		}

		double pValueFromPoint2PointDissimilarity(IDataSet * data, IDataSet * mcData)
		{
            std::cout << "GC: calculating T stat for data" << std::endl;
            double T = calculateTstatistic( data, mcData );
			char buffer[20];
			sprintf( buffer, "Tdata = %f", T );
			cout << buffer << endl;

			int nPerm = 25;
			vector<double> Tvalues = permutation( data, mcData, nPerm );

            int count = 0;
			vector<double>::iterator Titer;
			for ( Titer = Tvalues.begin(); Titer != Tvalues.end(); ++Titer ){
				if ( T < *Titer ) count++;
			}

			double pvalue = count/double(nPerm);

	    string fileName = ResultFormatter::GetOutputFolder();	
	    fileName.append("/tvalues.root");
	    	
            TFile * outputFile = new TFile(fileName.c_str(), "RECREATE");
            TNtuple * ntuple = new TNtuple("tvalues", "tvalues", "T:Tdata:pvalue");
            for ( int i = 0; i < nPerm; i++ ) ntuple->Fill(Tvalues[i], T, pvalue);
            ntuple->Write();
            outputFile->Close();
            delete outputFile;

       		return pvalue;
        }

		double calculateTstatistic( IDataSet * data, IDataSet * mcData )
		{
			return sumEvents( data ) /*+ sumEvents( mcData )*/ - sumDataMCEvents( data, mcData );
		}

		vector<double> permutation( IDataSet * data, IDataSet * mc, int nPerm )
		{
            cout << "Performing boostrapping " << nPerm << " times" << endl;
			vector<double> bootstrappedTvalues;
			for ( int i = 0; i < nPerm; i++ ) {
                cout << "Boostrap #" << i << endl;
				double T = permutationCore( data, mc, i );
				bootstrappedTvalues.push_back(T);
				char buffer[20];
				sprintf( buffer, "Tperm%i = %f", i, T );
				cout << buffer << endl;
			}
			return bootstrappedTvalues;
		}

		double permutationCore( IDataSet * data, IDataSet * mc, int iteration )
		{
			unsigned int nD = (unsigned) data->GetDataNumber();
			unsigned int nMC = (unsigned) mc->GetDataNumber();
			// Append the two datasets
			MemoryDataSet * dataClone = new MemoryDataSet( data->GetBoundary() );
			for ( unsigned int i = 0; i < nD;  i++ ) dataClone->AddDataPoint( data->GetDataPoint((int)i) );
			for ( unsigned int j = 0; j < nMC; j++ ) dataClone->AddDataPoint( mc  ->GetDataPoint((int)j) );

			TRandom3 * rand = new TRandom3((unsigned)iteration);

			// Need to get some empty version of the dataset
			MemoryDataSet * tempMC   = new MemoryDataSet( mc->GetBoundary() );
			MemoryDataSet * tempData = new MemoryDataSet( data->GetBoundary() );

			set<unsigned int> eventNotAccepted;
			for ( unsigned int i = 0; i < nD + nMC; i++ ) eventNotAccepted.insert(i);

			unsigned int count = 0;
			set<unsigned int>::iterator iter;
			while ( count < nD ) {
				unsigned int random = (unsigned int)rand->Uniform(0., nD + nMC);
				iter = find(eventNotAccepted.begin(), eventNotAccepted.end(), random);
				if ( iter != eventNotAccepted.end() ) {
					tempData->AddDataPoint( dataClone->GetDataPoint((int)random) );
					eventNotAccepted.erase(random);
					count++;
				}
			}

			for ( iter = eventNotAccepted.begin(); iter != eventNotAccepted.end(); ++iter) {
				tempMC->AddDataPoint( dataClone->GetDataPoint((int)*iter) );
			}
			double T = calculateTstatistic( tempData, tempMC );
			//delete rand;
			//delete tempMC;
			//delete tempData;
			//delete dataClone;
			return T;
		}

		/*
		   vector< vector< double > > getMeanAndSigma( IDataSet * data ) {
		   vector< vector< double > > meanAndSigma;
		   int n = data->GetDataNumber();
		   for ( int i = 0; i < n; i++ ){

		   }
		   return meanAndSigma;
		   }
		 */

		double sumEvents( IDataSet * data )
		{
            int n = data->GetDataNumber();
			double distance = 0.;
            double weight_i = 1.;
            double weight_j = 1.;
            double sum_weights_i = 0.;
            double T = 0.;
			for ( int i = 0; i < n; i++ ){
				if (i % 100 == 0) cout << "sumEvents - Event # " << i << "\t\t" << setprecision(4) << 100.*(double)i/(double)n << "\% Complete\b\b\b\b\b\b\b\r\r\r\r\r\r\r\r\r\r\r";
                DataPoint * event_i = data->GetDataPoint(i);
                weight_i = getWeight( event_i );
                sum_weights_i += weight_i;
				for ( int j = i+1; j < n; j++){
					DataPoint * event_j = data->GetDataPoint(j);
					distance = getDistance( event_i, event_j );
                    weight_j = getWeight( event_j );
					//cout << "sumEvents " << distance << " " << event_j->GetObservable("mass")->GetValue() << " " << event_i->GetObservable("mass")->GetValue() << endl;
					//if ( distance == 0. ) cout << "sumEvents " << distance << " " << i << " " << j << endl;
					//T += 1./distance;  // This is an alternative function which could be used
					//T += exp( -distance * distance / (2.*0.01) ); // sigma = 0.1 is a nuisance parameter. 0.1 is chosen arbitrarily.
					//T += distance * distance;
					T += -log( distance + 1./n ) * (weight_i * weight_j);
				}
			}
            cout << endl;
			T *= 1. / (sum_weights_i * sum_weights_i);
			return T;
		}

		double sumDataMCEvents( IDataSet * data, IDataSet * mcData )
		{
			int nD = data->GetDataNumber();
			int nMC = mcData->GetDataNumber();
			double distance = 0.;
            double weight_i = 1.;
            double weight_j = 1.; // MC weight should always be 1 for an sFit since signal only
            double sum_weights_i = 0.;
			double T = 0.;
			for ( int i = 0; i < nD; i++){
				if (i % 100 == 0) cout << "sumDataMCEvents - DataEvent # " << i << "\t\t" << setprecision(4) << 100.*(double)i/(double)nD << "\% Complete\b\b\b\b\b\b\b\r\r\r\r\r\r\r\r\r\r\r";
				DataPoint* event_i = data->GetDataPoint(i);
                weight_i = getWeight( event_i );
                sum_weights_i += weight_i;
				for ( int j = 0; j < nMC; j++){
					DataPoint * event_j = mcData->GetDataPoint(j);
					distance = getDistance( event_i, event_j );
					//cout << "sumDataMCEvents " << distance << " " << event_j->GetObservable("mass")->GetValue() << " " << event_i->GetObservable("mass")->GetValue() << endl;
					//if ( distance == 0. ) cout << "sumDataMCEvents " << distance << " " << i << " " << j << endl;
					//T += 1./distance;
					//T += exp( -distance * distance / (2.*0.01) );
					//T += distance * distance;
					T += -log( distance + 1./nD ) * (weight_i * weight_j);
				}
			}
            cout << endl;
			T *= 1. / (sum_weights_i * nMC);
			return T;
		}

        double getWeight( DataPoint * x )
        {
       		vector<string> xListOfNames = x->GetAllNames();
			vector<string>::iterator xIter = xListOfNames.begin();
			Observable * xVar = 0;
			double xVal = 1.;
			while ( ( xVar = x->GetObservable( *xIter ) ) ) {
                if ( (*xIter == "fsig_sw") ) {
					xVal = xVar->GetValue();
				}
				++xIter;
				if ( xIter == xListOfNames.end() ) break;
            }
            return xVal;
        }

		double getDistance(DataPoint * x, DataPoint * y)
		{
			double distance = 0.;
			double diff = 0.;	(void) diff;	//	Unused
			vector<string> xListOfNames = x->GetAllNames();
			vector<string> yListOfNames = y->GetAllNames();
			vector<string>::iterator xIter = xListOfNames.begin();
			vector<string>::iterator yIter = yListOfNames.begin();
			Observable * xVar = 0;
			Observable * yVar = 0;
			double xVal = 0.;
			double yVal = 0.;
			char buffer[100];	(void) buffer;	//	Unused
			while ( (xVar = x->GetObservable( *xIter ) ) && (yVar = y->GetObservable( *yIter ) ) ) {
				//cout << *xIter << endl;
				//if ( (*xIter == "time") ) {
				//	xVal = ( xVar->GetValue() - 1.6 )/1.2;
			    //	yVal = ( yVar->GetValue() - 1.6 )/1.2;
				//}
				if ( (*xIter == "mass") ) {
					xVal = ( xVar->GetValue() - 5200. )/350.;
					yVal = ( yVar->GetValue() - 5200. )/350.;
				}
				else if ( (*xIter == "cosTheta1") ) {
					xVal = ( xVar->GetValue() );
					yVal = ( yVar->GetValue() );
				}
				else if(  (*xIter == "cosTheta2") ) {
					xVal = ( xVar->GetValue() );
					yVal = ( yVar->GetValue() );
				}
				else if ( (*xIter == "phi") ) {
					xVal = ( xVar->GetValue() )/3.14159;
					yVal = ( yVar->GetValue() )/3.14159;
				}
				else {
					xVal = 0.;
					yVal = 0.;
				}
				distance += ( (xVal - yVal)*(xVal - yVal) );
				//sprintf(buffer, " %f %f %f", xVar->GetValue(), yVar->GetValue(), sqrt(distance));
				//cout << (*xIter) << buffer << endl;
				++xIter, ++yIter;
				if ( xIter == xListOfNames.end() || yIter == yListOfNames.end() ) break;
			}
			return sqrt(distance);
		}

		double diffmod2pi( double input ) {
			double pi2 = 2.*TMath::Pi();
			double pi = TMath::Pi();
			double diff = fabs(input);
			if( diff >= pi ) diff = pi2 - diff;
			return diff;
		}


}
