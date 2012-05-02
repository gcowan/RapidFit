// $Id: Bs2JpsiPhiMassBkgLL.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $

/** @class Bs2JpsiPhiMassBkgLL Bs2JpsiPhiMassBkgLL.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi mass background
 *
 *  @author Pete Clarke - copy to give two different PDFs for LL and prompt
 *  @date 2010-01-24
 */

#include "Bs2DsPi_lowmassbkg.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "Mathematics.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TH1.h"

#define DOUBLE_TOLERANCE 1E-6

PDF_CREATOR( Bs2DsPi_lowmassbkg );

//.....................................................
//Constructor
Bs2DsPi_lowmassbkg::Bs2DsPi_lowmassbkg(PDFConfigurator* configurator ) :

	// Physics parameters
	//alphaM_llName	( "alphaM_ll" )
	// Observables
	//, recoMassName  ( "mass" )

	massName  ( configurator->getName("mass") )
	,constraint_massName()
, histo(), xaxis(), nxbins(), xmin(), xmax(), deltax(), total_num_entries()
{

	cout << "Constructing PDF:Bs2DsPi_lowmassbkg" << endl ;


	{
		MakePrototypes();
	}

	//Find name of histogram needed to define 1-D mass distribution
	string fileName = configurator->getConfigurationValue( "MassDistributionHistogram" ) ;

	//Initialise depending upon whether configuration parameter was found
	if( fileName == "" )
	{
		cout << "   No MassDistributionHistogram found" << endl;
		exit(1);
	}
	else
	{
		cout << "   MassDistributionHistogram found: " << fileName << endl ;

		//Read in histo
		TFile* f =  TFile::Open(fileName.c_str());
		if( f == NULL ) exit(1) ;
		histo = new TH1D(  *(    (TH1D*)f ->Get("hist")      )     ); //(fileName.c_str()))); 

		xaxis = histo->GetXaxis();
		xmin = xaxis->GetXmin();
		xmax = xaxis->GetXmax();
		nxbins = histo->GetNbinsX();
		deltax = (xmax-xmin)/nxbins;

		//method for Checking whether histogram is sensible

		total_num_entries = histo->GetEntries();
		int total_num_bins = nxbins;
		int sum = 0;

		vector<int> zero_bins;

		//loop over each bin in histogram and print out how many zero bins there are
		for (int i=1; i < nxbins+1; i++)
		{
			double bin_content = histo->GetBinContent(i);
			//cout << "Bin content: " << bin_content << endl;
			if(bin_content<=0)
			{
				zero_bins.push_back(1);
			}
			//cout << " Zero bins " << zero_bins.size() << endl;}
			else if (bin_content>0)
			{
				sum += (int) bin_content;
			}
		}

		int average_bin_content = sum / total_num_bins;

		cout << "\n\n\t\t" << "****" << "For total number of bins " << total_num_bins << " there are " << zero_bins.size() << " bins containing zero events " << "****" << endl;
		cout <<  "\t\t\t" << "***" << "Average number of entries of non-zero bins: " << average_bin_content << "***" << endl;
		cout << endl;
		cout << endl;

		if ((xmax-xmin) < 2.)
		{
			cout << "In Bs2DsPi_lowmassbkg::Bs2DsPi_lowmassbkg: The full mass range is not used in this histogram" << endl;
			exit(1);
		}

		cout << "Finishing processing histo" << endl;
	}
}



//Make the data point and parameter set
void Bs2DsPi_lowmassbkg::MakePrototypes()
{
	allObservables.push_back( massName );
	constraint_massName = massName;
}

//Destructor
Bs2DsPi_lowmassbkg::~Bs2DsPi_lowmassbkg()
{
}


//..............................................................
//Main method to build the PDF return value
double Bs2DsPi_lowmassbkg::Evaluate(DataPoint * measurement)
{
	// Observable
	double mass = measurement->GetObservable( massName )->GetValue();


	int globalbin = histo->FindBin(mass);
	double num_entries_bin = histo->GetBinContent(globalbin);

	double value = num_entries_bin/total_num_entries;

	return value;

}


double Bs2DsPi_lowmassbkg::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	(void)measurement;
	double mhigh, mlow, sum ;

	IConstraint * massBound = boundary->GetConstraint( constraint_massName );
	if ( massBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on mass not provided in Bs2DsPi_lowmassbkg" << endl;
		return 1.0 ;
	}
	else
	{
		mlow = massBound->GetMinimum();
		mhigh = massBound->GetMaximum();
	}


	double histo_max=histo->GetXaxis()->GetBinUpEdge(nxbins);
	double histo_min=histo->GetXaxis()->GetBinLowEdge(1);

	double histo_range = histo_max - histo_min;
	double bin_width = histo_range / nxbins;

	int nbin_low =1;
	int nbin_high = nxbins;

	bool c0 = (mlow>mhigh) ;
	bool c1 = (mhigh>histo_max) || (mlow<histo_min) ;	
	if ( c0 || c1  )
	{
		cerr << "Mass boundaries aren't withing the background histogram range in Bs2DsPi_lowmassbkg" << endl;
		exit(1) ;
	}

	for( nbin_low=1; nbin_low <= nxbins; ++nbin_low ) {		
		if( histo->GetXaxis()->GetBinLowEdge(nbin_low)  > mlow ) break ;
	}
	nbin_low--;
	for( nbin_high=nxbins; nbin_high >= 1; --nbin_high ) {		
		if( histo->GetXaxis()->GetBinUpEdge(nbin_high)  < mhigh ) break ;
	}
	nbin_high++;

	//cout << "mlow occurs in bin" << nbin_low << endl;
	//cout << "mhigh occurs in bin" << nbin_high << endl;

	sum = 0;

	for(int i = nbin_low; i <= nbin_high; ++i){
		double bin_content = histo->GetBinContent(i);
		sum += (int) bin_content;
		//cout << "loop " << i << ", ";
	}



	//cout << "number of bins " << nxbins << endl;
	//cout << "sum of bins looped " << sum << endl;
	//cout << "histo range " << histo_range << endl;
	//cout << "bin_width " << bin_width << endl;	
	//cout << "total number enteries in histo " << total_num_entries << endl;	


	//cout << "end normalisation " << endl;


	double value = (sum * bin_width) / total_num_entries;


	//cout << "return " << value << endl;




	return value;
}

