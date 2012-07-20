// $Id: PerEventErrorHistogram.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class PerEventErrorHistogram PerEventErrorHistogram.cpp
 *
 *  PDF for background due to wrong PV association 
 *
 *  @author Greig Cowan
 *  @date 2012-05-29
 */

#include "PerEventErrorHistogram.h"
#include "Mathematics.h"
#include <iostream>
#include <fstream>
#include "math.h"
#include "TMath.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TH3D.h"
#include "TAxis.h"
#include "TH1.h"

PDF_CREATOR( PerEventErrorHistogram );

//Constructor
PerEventErrorHistogram::PerEventErrorHistogram( PDFConfigurator* config ) :

	// Observables
	 eventResolutionName   ( config->getName("eventResolution") )
	, xaxis()
	, nxbins()
	, xmin(), xmax()
	, deltax()
	, total_num_entries()
{

	cout << "Constructing PDF: PerEventErrorHistogram  " << endl ;

	//Make prototypes
	MakePrototypes();
	
	//Find name of histogram needed to define 3-D angular distribution
	string fileName = config->getConfigurationValue( "ErrorHistogram" ) ;
	string histName = config->getConfigurationValue( "HistogramName" ) ;

	//Initialise depending upon whether configuration parameter was found
	if( fileName == "" )
	{
		cout << "   No background histogram found: exiting " << endl ;
		exit(-1);
	}
	else
	{
		cout << "   Per-event error histogram requested: " << fileName << endl ;

		//File location
		ifstream input_file;
		input_file.open( fileName.c_str(), ifstream::in );
		input_file.close();

		bool local_fail = input_file.fail();

		if( !getenv("RAPIDFITROOT") && local_fail )
		{
			cerr << "\n" << endl;
			//cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cerr << "$RAPIDFITROOT NOT DEFINED, PLEASE DEFINE IT SO I CAN USE ACCEPTANCE DATA" << endl;
			//cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cerr << "\n" << endl;
			exit(-987);
		}

		string fullFileName;

		if( getenv("RAPIDFITROOT") )
		{
			string path( getenv("RAPIDFITROOT") ) ;

			cout << "RAPIDFITROOT defined as: " << path << endl;

			fullFileName = path+"/pdfs/configdata/"+fileName ;

			input_file.open( fullFileName.c_str(), ifstream::in );
			input_file.close();
		}
		bool elsewhere_fail = input_file.fail();

		if( elsewhere_fail && local_fail )
		{
			cerr << "\n\tFileName:\t" << fullFileName << "\t NOT FOUND PLEASE CHECK YOUR RAPIDFITROOT" << endl;
			cerr << "\t\tAlternativley make sure your XML points to the correct file, or that the file is in the current working directory\n" << endl;
			exit(-89);
		}

		if( fullFileName.empty() || !local_fail )
		{
			fullFileName = fileName;
		}

		//Read in histo
		TFile* f =  TFile::Open(fullFileName.c_str());
		histo = (TH1D*) f->Get(histName.c_str());

		// time
		xaxis = histo->GetXaxis();
		cout << "X axis Name: " << xaxis->GetName() << "\tTitle: " << xaxis->GetTitle() << endl;
		xmin = xaxis->GetXmin();
		xmax = xaxis->GetXmax();
		nxbins = xaxis->GetNbins();
		cout << "X axis Min: " << xmin << "\tMax: " << xmax << "\tBins: " << nxbins << endl;
		deltax = (xmax-xmin)/nxbins;

		//method for Checking whether histogram is sensible

		for ( int i = 0; i < histo->GetEntries(); i++ ) {
			total_num_entries += histo->GetBinContent(i);
		}
		int total_num_bins = nxbins ;
		int sum = 0;

		vector<int> zero_bins;

		//loop over each bin in histogram and print out how many zero bins there are
		for (int i=1; i < nxbins+1; ++i)
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

		cout << "\n\n\t\t" << "****" << "For total number of bins " << total_num_bins << " there are " << zero_bins.size() << " bins containing zero events and " << total_num_entries << "total num entries ****" << endl;
		cout <<  "\t\t\t" << "***" << "Average number of entries of non-zero bins: " << average_bin_content << "***" << endl;
		cout << endl;
		cout << endl;

		// Check.  This order works for both bases since phi is always the third one.
		if ((xmax-xmin) < 0.1 ) 
		{
			cout << "In PerEventErrorHistogram::PerEventErrorHistogram: The full angular range is not used in this histogram - the PDF does not support this case" << endl;
			exit(1);
		}

		cout << "Finishing processing histo" << endl;
	}
}
/*
// Copy
PerEventErrorHistogram::PerEventErrorHistogram( const PerEventErrorHistogram& input ) : BasePDF( (BasePDF) input ),
histo(input.histo), xaxis(input.xaxis), 
nxbins(input.nxbins), 
xmin(input.xmin), xmax(input.xmax), 
deltax(input.deltax), 
total_num_entries(input.total_num_entries), 
eventResolutionName( input.eventResolutionName ), 
eventResolution(input.eventResolution)
{


}
*/
//Destructor
PerEventErrorHistogram::~PerEventErrorHistogram()
{
}

//Make the data point and parameter set
void PerEventErrorHistogram::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( eventResolutionName );

	//Make the parameter set
	vector<string> parameterNames;
	allParameters = ParameterSet(parameterNames);
}

//Return a list of observables not to be integrated
vector<string> PerEventErrorHistogram::GetDoNotIntegrateList()
{
        vector<string> list;
        return list;
}


bool PerEventErrorHistogram::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	return isOK;
}

//Main method to build the PDF return value
double PerEventErrorHistogram::Evaluate(DataPoint * measurement)
{
	// Observable
	eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();
	return this->timeMassFactor();
}


// Normlisation
double PerEventErrorHistogram::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void)boundary;
	return buildPDFdenominator();
}

double PerEventErrorHistogram::buildPDFdenominator()
{
	return 1.;
}


//Angular distribution function
double PerEventErrorHistogram::timeMassFactor( )
{
	double returnValue = 0.;

	int globalbin = -1;
	int xbin = -1;
	double num_entries_bin = -1.;

		//Find global bin number for values of angles, find number of entries per bin, divide by volume per bin and normalise with total number of entries in the histogram
		xbin = xaxis->FindFixBin( eventResolution ); if( xbin > nxbins ) xbin = nxbins;
		
		globalbin = histo->GetBin( xbin );
		num_entries_bin = histo->GetBinContent(globalbin);

		if ( num_entries_bin < 0. ) num_entries_bin = 0.;

		//Angular factor normalized with phase space of histogram and total number of entries in the histogram
		returnValue = num_entries_bin / (deltax) / total_num_entries ;
	
	return returnValue;
}
