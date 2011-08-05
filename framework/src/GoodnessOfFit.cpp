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
//	System Headers
#include <math.h>
#include <iostream>
#include <set>
#include <algorithm>
#include <time.h>

using namespace std;

namespace GoodnessOfFit
{

	double gofLoop( XMLConfigReader * xmlFile, MinimiserConfiguration * theMinimiser, FitFunctionConfiguration * theFunction, vector<ParameterSet*> argumentParameterSet, vector<string> CommandLineParam, int nData )
	{
		cout << "Starting GOF" << endl;

		double pvalue = 0.;
		for ( int i = 0; i < 1; i ++ )
		{
			cout << "Iteration " << i << endl;
			// Set up the PDF parameters
			vector< ParameterSet * > parSet = xmlFile->GetFitParameters( CommandLineParam );

			// First, generate our toy data, fit to it, use PDF to generate large MC and then calculate corresponding p-value from permutation
			vector<double> dataPvalue;
			GoodnessOfFit::generateFitAndCalculatePvalue( xmlFile, &parSet, theMinimiser, theFunction, argumentParameterSet, nData, 1, &dataPvalue );

			// Now need to account for any potential fit bias so repeat above procedure N times to get distribution of p-values
			vector<double> distOfPvalues;
			GoodnessOfFit::generateFitAndCalculatePvalue( xmlFile, &parSet, theMinimiser, theFunction, argumentParameterSet, nData, 100, &distOfPvalues );

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

	void generateFitAndCalculatePvalue( XMLConfigReader * xmlFile, vector<ParameterSet*> * parSet, MinimiserConfiguration * theMinimiser, FitFunctionConfiguration * theFunction, vector<ParameterSet*> argumentParameterSet, int nData, int repeats, vector<double> * pvalues)
	{
		PDFWithData * pdfAndData = xmlFile->GetPDFsAndData()[0];
		IPDF * pdf = pdfAndData->GetPDF();
		vector<IPDF*> vectorPDFs;
		vectorPDFs.push_back(pdf);
		PhaseSpaceBoundary * phase = xmlFile->GetPhaseSpaceBoundaries()[0];
		EdStyle * greigFormat = new EdStyle();
		greigFormat->SetStyle();

		vector< ParameterSet * > parSetFromFit;

		// Set up to be able to generate some MC data
		DataSetConfiguration * dataConfig = pdfAndData->GetDataSetConfig();
		dataConfig->SetSource( "Foam" );
		bool model = false;
		double pvalue = 0.;
		unsigned int ulTime = static_cast<unsigned int>( time( NULL ));
		for ( int i = 0; i < repeats; i++ ) {
			cout << "Ensemble " << i << endl;

			// First, generate some data from this PDF
			pdfAndData->SetPhysicsParameters( *parSet );
			pdf->SetRandomFunction( ulTime );
			pdf->SetMCCacheStatus( false );
			MemoryDataSet * subset = (MemoryDataSet*)dataConfig->MakeDataSet( phase, pdf, nData );
			vector<IDataSet*> vectorData;
			vectorData.push_back(subset);

			if ( model ) {
				// Use the same PDF for the large MC dataset (the Model approach)
				pdfAndData->SetPhysicsParameters( *parSet );
			}
			else {
				// Fit the generated data (the Fit I approach)
				FitResult * gofResult = FitAssembler::DoFit( theMinimiser, theFunction, argumentParameterSet, vectorPDFs, vectorData, xmlFile->GetConstraints() );
				parSetFromFit.clear();
				parSetFromFit.push_back(gofResult->GetResultParameterSet()->GetDummyParameterSet());
				pdfAndData->SetPhysicsParameters( parSetFromFit );
			}
			// Generate large sample of MC
			ulTime = static_cast<unsigned int>( time( NULL ));
			pdf->SetRandomFunction( ulTime );
			pdf->SetMCCacheStatus( false );
			MemoryDataSet * mcData = (MemoryDataSet*)dataConfig->MakeDataSet( phase, pdf, 10*nData );

			// Finally calculate the p-value relative to the large MC sample
			pvalue = GoodnessOfFit::pValueFromPoint2PointDissimilarity( subset, mcData );
			//pvalue = GoodnessOfFit::pValueFromPoint2PointDissimilarity( data, mcData );
			pvalues->push_back(pvalue);
			delete subset;
			delete mcData;
		}
		parSet = &parSetFromFit;  // Need this so that we have the correct parameters during the second call of this function
	}

	void plotUstatistic( IPDF * pdf, IDataSet * data, PhaseSpaceBoundary * phase, string plot )
	{
		unsigned int nD = data->GetDataNumber();
		int bins = 30;
		double level = nD/bins;
		TH1D * distances = new TH1D( "distances", "U", bins, 0., 1.);
		distances->Sumw2();
		calculateUstatistic(pdf, data, phase, distances);

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
		delete distances;
		delete ca;
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
		unsigned int nD = data->GetDataNumber();
		int count = 0;
		for (unsigned int i = 0; i < nD; i++) {
			if ( !(i % 100) ) cout << i << endl;
			event_i = data->GetDataPoint( i );
			pdfValue = pdf->Evaluate( event_i )/pdf->Integral( event_i, phase );
			firstEvent = true;
			for (unsigned int j = 0; j < nD; j++) {
				if ( j == i ) continue;
				event_j = data->GetDataPoint( j );
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
			double distToTimeBoundary = event_i->GetObservable( "time" )->GetValue() - 0.3;
			double distToMassBoundary = min(fabs(event_i->GetObservable( "mass" )->GetValue() - 5200.),fabs(event_i->GetObservable( "mass" )->GetValue() - 5550.)) ;
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
		double volume = 0.;
		double distance = 0.;
		double smallest_distance = 0.;
		double U = 0;
		DataPoint * event_i = 0;
		DataPoint * event_j = 0;
		DataPoint * closest = 0;
		bool firstEvent = true;
		size_t dimension = (pdf->GetPrototypeDataPoint()).size();
		if ( dimension == 7 ) dimension = 5;
		cout << "number of dimensions: " << dimension << endl;
		unsigned int nD = data->GetDataNumber();
		int count = 0;
		for (unsigned int i = 0; i < nD; i++) {
			if ( !(i % 100) ) cout << "Event " << i << endl;
			event_i = data->GetDataPoint( i );
			firstEvent = true;
			for (unsigned int j = 0; j < nD; j++) {
				if ( j == i ) continue;
				event_j = data->GetDataPoint( j );
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
			vector<double> vectorOfDistances = getDistances( event_i, closest );
			if ( vectorOfDistances[0] > vectorOfDistances[1] ) count++;
			PhaseSpaceBoundary * tempPhase = new PhaseSpaceBoundary( phase->GetAllNames() );
			copyPhaseSpaceBoundary( tempPhase, phase );
			updatePhaseSpaceBoundary( event_i, tempPhase, phase, vectorOfDistances ); 
			RapidFitIntegrator * integrator = new RapidFitIntegrator( pdf, true );
			volume = integrator->Integral( event_i, tempPhase, false );
			delete tempPhase;
			delete integrator;
			U = exp( -1. * nD * volume );
			double anaVolume = TMath::Pi() * pow(smallest_distance, 2);
			double anaU = exp(-1.*nD         *     TMath::Pi()     * pow(smallest_distance, 2) *   pdf->Evaluate( event_i ) );
			//double anaVolume = 2.*smallest_distance;
			//double anaU = exp(-1.*nD*2.                            *     smallest_distance     * pdf->Evaluate( event_i ) );
			double pdfVal = pdf->Evaluate( event_i );
			cout << "sd " << smallest_distance << " numVol: " << volume << " anaVol: " << anaVolume << " numU: " << U << " anaU: " << anaU << " pdfVal: " << pdfVal << endl;
			distances->Fill( U );
		}
		cout << "count " << count << endl;
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
		char buffer[100];

		while ( ( xVar = x->GetObservable( *xIter ) ) ) {
			if ( *xIter == "time" || *xIter == "cosTheta" || *xIter == "mass"){//|| *xIter == "cosPsi" || *xIter == "phi") {
				min = xVar->GetValue() - distances[ i ];
				max = xVar->GetValue() + distances[ i ];
				if ( min > oldPhase->GetConstraint( *xIter )->GetMinimum() ) ((ObservableContinuousConstraint *)newPhase->GetConstraint( *xIter ))->SetMinimum( min );
				if ( max < oldPhase->GetConstraint( *xIter )->GetMaximum() ) ((ObservableContinuousConstraint *)newPhase->GetConstraint( *xIter ))->SetMaximum( max );
				//sprintf(buffer, "%f %f %f %f %f %f", xVar->GetValue(), distances[i], min, max, newPhase->GetConstraint( *xIter )->GetMinimum(), newPhase->GetConstraint( *xIter )->GetMaximum());
				//cout << *xIter << " " << buffer << endl;
			}
			++xIter;
			if ( xIter == xListOfNames.end() ) break;
			i++;
			}
		}

		double pValueFromPoint2PointDissimilarity(IDataSet * data, IDataSet * mcData)
		{
			double T = calculateTstatistic( data, mcData );
			//char buffer[20];
			//sprintf( buffer, "T = %f", T );
			//cout << buffer << endl;

			int nPerm = 25;
			vector<double> Tvalues = permutation( data, mcData, nPerm );

			/*	
				TH1D * tHist = new TH1D("tvalues", "tvalues", 50, 1.0, 1.3);
				for ( int i = 0; i < nPerm; i++ ) tHist->Fill(Tvalues[i]);
				TFile * outputFile = new TFile("tvalues.root", "UPDATE");
				tHist->Write();
				outputFile->Close();
				delete outputFile;
			 */

			int count = 0;
			vector<double>::iterator Titer;
			for ( Titer = Tvalues.begin(); Titer != Tvalues.end(); ++Titer ){
				if ( T < *Titer ) count++;
			}

			double pvalue = count/double(nPerm);
			return pvalue;
		}

		double calculateTstatistic( IDataSet * data, IDataSet * mcData )
		{
			return sumEvents( data ) /*+ sumEvents( mcData )*/ - sumDataMCEvents( data, mcData );
		}

		vector<double> permutation( IDataSet * data, IDataSet * mc, int nPerm )
		{
			vector<double> bootstrappedTvalues;
			for ( int i = 0; i < nPerm; i++ ) {
				double T = permutationCore( data, mc, i );
				bootstrappedTvalues.push_back(T);
				//char buffer[20];
				//sprintf( buffer, "Tperm%i = %f", i, T );
				//cout << buffer << endl;
			}
			return bootstrappedTvalues;
		}

		double permutationCore( IDataSet * data, IDataSet * mc, int iteration )
		{
			unsigned int nD = data->GetDataNumber();
			unsigned int nMC = mc->GetDataNumber();
			// Append the two datasets
			MemoryDataSet * dataClone = new MemoryDataSet( data->GetBoundary() );
			for ( unsigned int i = 0; i < nD;  i++ ) dataClone->AddDataPoint( data->GetDataPoint(i) );
			for ( unsigned int j = 0; j < nMC; j++ ) dataClone->AddDataPoint( mc  ->GetDataPoint(j) );

			TRandom3 * rand = new TRandom3(iteration);

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
					tempData->AddDataPoint( dataClone->GetDataPoint(random) );
					eventNotAccepted.erase(random);
					count++;
				}
			}

			for ( iter = eventNotAccepted.begin(); iter != eventNotAccepted.end(); ++iter) {
				tempMC->AddDataPoint( dataClone->GetDataPoint(*iter) );
			}
			double T = calculateTstatistic( tempData, tempMC );
			delete rand;
			delete tempMC;
			delete tempData;
			delete dataClone;
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
			double T = 0.;
			for ( int i = 0; i < n; i++ ){
				DataPoint * event_i = data->GetDataPoint(i);
				for ( int j = i+1; j < n; j++){
					DataPoint * event_j = data->GetDataPoint(j);
					distance = getDistance( event_i, event_j );
					//cout << "sumEvents " << distance << " " << event_j->GetObservable("mass")->GetValue() << " " << event_i->GetObservable("mass")->GetValue() << endl;
					//if ( distance == 0. ) cout << "sumEvents " << distance << " " << i << " " << j << endl;
					//T += 1./distance;  // This is an alternative function which could be used
					//T += exp( -distance * distance / (2.*0.01) ); // sigma = 0.1 is a nuisance parameter. 0.1 is chosen arbitrarily.
					//T += distance * distance;
					T += -log( distance + 1./n );
				}
			}
			T *= 1./(n * n); 
			return T;
		}

		double sumDataMCEvents( IDataSet * data, IDataSet * mcData )
		{
			int nD = data->GetDataNumber();
			int nMC = mcData->GetDataNumber();
			double distance = 0.;
			double T = 0.;
			for ( int i = 0; i < nD; i++){
				DataPoint* event_i = data->GetDataPoint(i);
				for ( int j = 0; j < nMC; j++){
					DataPoint * event_j = mcData->GetDataPoint(j);
					distance = getDistance( event_i, event_j );
					//cout << "sumDataMCEvents " << distance << " " << event_j->GetObservable("mass")->GetValue() << " " << event_i->GetObservable("mass")->GetValue() << endl;
					//if ( distance == 0. ) cout << "sumDataMCEvents " << distance << " " << i << " " << j << endl;
					//T += 1./distance;
					//T += exp( -distance * distance / (2.*0.01) );
					//T += distance * distance;
					T += -log( distance + 1./nD );
				}
			}
			T *= 1./(nD * nMC);
			return T;
		}

		double getDistance(DataPoint * x, DataPoint * y)
		{
			double distance = 0.;
			double diff = 0.;
			vector<string> xListOfNames = x->GetAllNames();
			vector<string> yListOfNames = y->GetAllNames();
			vector<string>::iterator xIter = xListOfNames.begin();
			vector<string>::iterator yIter = yListOfNames.begin();
			Observable * xVar = 0;
			Observable * yVar = 0;
			double xVal = 0.;
			double yVal = 0.;
			char buffer[100];
			while ( (xVar = x->GetObservable( *xIter ) ) && (yVar = y->GetObservable( *yIter ) ) ) {
				if ( (*xIter == "time") ) {
					xVal = ( xVar->GetValue() - 1.6 )/1.2;
					yVal = ( yVar->GetValue() - 1.6 )/1.2;
				}
				else if ( (*xIter == "mass") ) {
					xVal = ( xVar->GetValue() - 5366. )/47.;
					yVal = ( yVar->GetValue() - 5366. )/47.;
				}
				else if ( (*xIter == "cosTheta") ) {
					xVal = ( xVar->GetValue() - 0. )/0.5;
					yVal = ( yVar->GetValue() - 0. )/0.5;
				}
				else if(  (*xIter == "cosPsi") ) {
					xVal = ( xVar->GetValue() - 0. )/0.5;
					yVal = ( yVar->GetValue() - 0. )/0.5;
				}
				else if ( (*xIter == "phi") ) {
					xVal = ( xVar->GetValue() - 0. )/(3.14159/2.);
					yVal = ( yVar->GetValue() - 0. )/(3.14159/2.);
				}
				distance += ( (xVal - yVal)*(xVal - yVal) );
				//sprintf(buffer, "%f %f %f", xVar->GetValue(), yVar->GetValue(), sqrt(distance));
				//cout << buffer << endl;
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
