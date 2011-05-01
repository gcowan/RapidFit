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
//	RapidFit Headers
#include "GoodnessOfFit.h"
#include "MemoryDataSet.h"
//	System Headers
#include <math.h>
#include <iostream>
#include <set>

using namespace std;

namespace GoodnessOfFit
{
	void plotUstatistic( IPDF * pdf, IDataSet * data, PhaseSpaceBoundary * phase, string plot )
	{
                unsigned int nD = data->GetDataNumber();
		int bins = 100;
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
		bool firstEvent = true;
		size_t dimension = (pdf->GetPrototypeDataPoint()).size();
		if ( dimension == 7 ) dimension = 2;
		cout << "number of dimensions: " << dimension << endl;
		unsigned int nD = data->GetDataNumber();
		for (unsigned int i = 0; i < nD; i++) {
			if ( !(i % 100) ) cout << i << endl;
			event_i = data->GetDataPoint( i );
			pdfValue = pdf->Evaluate( event_i )/pdf->Integral( event_i, phase );
			//cout << pdf->Integral( event_i, phase ) << endl;
			firstEvent = true;
			for (unsigned int j = 0; j < nD; j++) {
				if ( j == i ) continue;
				event_j = data->GetDataPoint( j );
				distance = getDistance( event_i, event_j );
				if ( firstEvent ) {
					smallest_distance = distance;
					firstEvent = false;
				}
				if (distance < smallest_distance) smallest_distance = distance;
			}
			sd = smallest_distance;
			if ( dimension == 1 ) U = exp(-1.*nD*2.                            *     sd     * pdfValue );
			if ( dimension == 2 ) U = exp(-1.*nD         *     TMath::Pi()     * pow(sd, 2) * pdfValue );
			if ( dimension == 3 ) U = exp(-1.*nD*4./3.   *     TMath::Pi()     * pow(sd, 3) * pdfValue );
			if ( dimension == 4 ) U = exp(-1.*nD*1./2.   *pow( TMath::Pi(), 2) * pow(sd, 4) * pdfValue );
			if ( dimension == 5 ) U = exp(-1.*nD*8./15.  *pow( TMath::Pi(), 2) * pow(sd, 5) * pdfValue );
			if ( dimension == 6 ) U = exp(-1.*nD*1./6.   *pow( TMath::Pi(), 3) * pow(sd, 6) * pdfValue );
			if ( dimension == 7 ) U = exp(-1.*nD*16./105.*pow( TMath::Pi(), 3) * pow(sd, 7) * pdfValue );
			if ( U < 1.e-02 ) {
				cout << "pdfValue " << pdfValue << " smallest " << sd <<  " U " << U << " time " << event_i->GetObservable( "time" )->GetValue()
				<< " mass " << event_i->GetObservable( "mass" )->GetValue()
				<< " cosTheta " << event_i->GetObservable( "cosTheta" )->GetValue()
				<< " phi " << event_i->GetObservable( "phi" )->GetValue() 
				<< " cosPsi " << event_i->GetObservable( "cosPsi" )->GetValue()  << endl;
			}
			distances->Fill( U );
		}
	}

	double pValueFromPoint2PointDissimilarity(IDataSet * data, IDataSet * mcData)
	{
		double T = calculateTstatistic( data, mcData );
		char buffer[20];
		sprintf( buffer, "T = %f", T );
		cout << buffer << endl;

		int nPerm = 20;
		vector<double> Tvalues = permutation( data, mcData, nPerm );

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
			char buffer[20];
			sprintf( buffer, "Tperm%i = %f", i, T );
			cout << buffer << endl;
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
		MemoryDataSet * tempData = new MemoryDataSet( data->GetBoundary() );
		MemoryDataSet * tempMC   = new MemoryDataSet( mc->GetBoundary() );

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

		for ( iter = eventNotAccepted.begin(); iter != eventNotAccepted.end(); ++iter) tempMC->AddDataPoint( dataClone->GetDataPoint(*iter) );

		double T = calculateTstatistic( tempData, tempMC );
		delete rand;
		delete dataClone;
		delete tempData;
		delete tempMC;
		return T;
	}

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
				//if ( distance == 0. ) cout << "sumEvents " << distance << " " << i << " " << j << endl;
				//T += 1./distance;  // This is an alternative function which could be used
				T += exp( -distance * distance / (2.*0.0001) ); // sigma = 0.1 is a nuisance parameter. 0.1 is chosen arbitrarily.
			}
		}
		//T *= 1./(n * (n - 1.));
		T *= 1./(n * n); // This is better behaved for low statistics
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
				//if ( distance == 0. ) cout << "sumDataMCEvents " << distance << " " << i << " " << j << endl;
				//T += 1./distance;
				T += exp( -distance * distance / (2.*0.0001) );
			}
		}
		T *= 1./(nD * nMC);
		return T;
	}

        double getDistance(DataPoint * x, DataPoint * y)
        {
                double distance = 0.;
                vector<string> xListOfNames = x->GetAllNames();
                vector<string> yListOfNames = y->GetAllNames();
                vector<string>::iterator xIter = xListOfNames.begin();
                vector<string>::iterator yIter = yListOfNames.begin();
                Observable * xVar = 0;
                Observable * yVar = 0;
                //char buffer[100];
                while ( (xVar = x->GetObservable( *xIter ) ) && (yVar = y->GetObservable( *yIter ) ) ) {
                        distance += (xVar->GetValue() - yVar->GetValue())*(xVar->GetValue() - yVar->GetValue());
                        //sprintf(buffer, "%f %f %f", xVar->GetValue(), yVar->GetValue(), sqrt(distance));
			//if ( distance == 0. ) cout << buffer << endl;
                        ++xIter, ++yIter;
                        if ( xIter == xListOfNames.end() || yIter == yListOfNames.end() ) break;
                }
                return sqrt(distance);
        }
}
