// $Id: LongLivedBkg_3Dangular.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class LongLivedBkg_3Dangular LongLivedBkg_3Dangular.cpp
 *
 *  PDF for  long lived background with 3D histogram angular description. A 3D root Histogram required provided by user and input in the xml file. The transversity angles in the order (x,y,z) --> (cosPsi, cosTheta, phi)
 *
 *  @author Ailsa Sparkes
 *  @date 2011-05-30
 */

#include "LongLivedBkg_3Dangular.h"
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

PDF_CREATOR( LongLivedBkg_3Dangular );

//.....................................................
//Constructor
LongLivedBkg_3Dangular::LongLivedBkg_3Dangular( PDFConfigurator* config ) :

	// Physics parameters
	f_LL1Name		( config->getName("f_LL1")  )
	, tauLL1Name		( config->getName("tau_LL1") )
	, tauLL2Name		( config->getName("tau_LL2") )
	//Detector parameters
	, timeResLL1FracName	( config->getName("timeResLL1Frac") )
	, sigmaLL1Name		( config->getName("sigma_LL1") )
	, sigmaLL2Name		( config->getName("sigma_LL2") )
	// Observables
	, timeName		( config->getName("time") )
	, cosThetaName		( config->getName("cosTheta") )
	, phiName		( config->getName("phi") )
	, cosPsiName		( config->getName("cosPsi") )
	//Other things to be initialised
	, timeconstName		( config->getName("time") )

, tauLL1(), tauLL2(), f_LL1(), sigmaLL(), sigmaLL1(), sigmaLL2(), timeResLL1Frac(), tlow(), thigh(), time(), cosTheta(),
	phi(), cosPsi(), histo(), xaxis(), yaxis(), zaxis(), nxbins(), nybins(), nzbins(), xmin(), xmax(), ymin(),
	ymax(), zmin(), zmax(), deltax(), deltay(), deltaz(), total_num_entries(), useFlatAngularDistribution(true)
{

	cout << "LongLivedBkg_3Dangular::  " ;

	MakePrototypes();

	//Find name of histogram needed to define 3-D angular distribution
	string fileName = config->getConfigurationValue( "AngularDistributionHistogram" ) ;

	//Initialise depending upon whether configuration parameter was found
	if( fileName == "" )
	{
		cout << "   No AngularDistributionHistogram found: using flat background " << endl ;
		useFlatAngularDistribution = true ;
	}
	else
	{
		cout << "   AngularDistributionHistogram requested: " << fileName << endl ;
		useFlatAngularDistribution = false ;

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
		histo = (TH3D*) f->Get("histo"); //(fileName.c_str())));

		xaxis = histo->GetXaxis();
		cout << "X axis Name: " << xaxis->GetName() << "\tTitle: " << xaxis->GetTitle() << endl;
		xmin = xaxis->GetXmin();
		xmax = xaxis->GetXmax();
		nxbins = histo->GetNbinsX();
		cout << "X axis Min: " << xmin << "\tMax: " << xmax << "\tBins: " << nxbins << endl;
		deltax = (xmax-xmin)/nxbins;

		yaxis = histo->GetYaxis();
		cout << "Y axis Name: " << yaxis->GetName() << "\tTitle: " << yaxis->GetTitle() << endl;
		ymin = yaxis->GetXmin();
		ymax = yaxis->GetXmax();
		nybins = histo->GetNbinsY();
		cout << "Y axis Min: " << ymin << "\tMax: " << ymax << "\tBins: " << nybins << endl;
		deltay = (ymax-ymin)/nybins;

		zaxis = histo->GetZaxis();
		cout << "Z axis Name: " << zaxis->GetName() << "\tTitle: " << zaxis->GetTitle() << endl;
		zmin = zaxis->GetXmin();
		zmax = zaxis->GetXmax();
		nzbins = histo->GetNbinsZ();
		cout << "Z axis Min: " << zmin << "\tMax: " << zmax << "\tBins: " << nzbins << endl;
		deltaz = (zmax-zmin)/nzbins;

		//method for Checking whether histogram is sensible

		total_num_entries = histo->GetEntries();
		int total_num_bins = nxbins * nybins * nzbins;
		int sum = 0;

		vector<int> zero_bins;

		//loop over each bin in histogram and print out how many zero bins there are
		for (int i=1; i < nxbins+1; ++i)
		{
			for (int j=1; j < nybins+1; ++j)
			{
				for (int k=1; k < nzbins+1; ++k)
				{

					double bin_content = histo->GetBinContent(i,j,k);
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
			}
		}

		int average_bin_content = sum / total_num_bins;

		cout << "\n\n\t\t" << "****" << "For total number of bins " << total_num_bins << " there are " << zero_bins.size() << " bins containing zero events " << "****" << endl;
		cout <<  "\t\t\t" << "***" << "Average number of entries of non-zero bins: " << average_bin_content << "***" << endl;
		cout << endl;
		cout << endl;

		if ((xmax-xmin) < 2. || (ymax-ymin) < 2. || (zmax-zmin) < 2.*TMath::Pi() )
		{
			cout << "In LongLivedBkg_3Dangular::LongLivedBkg_3Dangular: The full angular range is not used in this histogram - the PDF does not support this case" << endl;
			exit(1);
		}

		cout << "Finishing processing histo" << endl;
	}
}


//..................................................................
//Make the data point and parameter set
void LongLivedBkg_3Dangular::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( cosThetaName );
	allObservables.push_back( phiName );
	allObservables.push_back( cosPsiName );


	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( f_LL1Name );
	parameterNames.push_back( tauLL1Name );
	parameterNames.push_back( tauLL2Name );
	parameterNames.push_back( timeResLL1FracName );
	parameterNames.push_back( sigmaLL1Name );
	parameterNames.push_back( sigmaLL2Name );

	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}


LongLivedBkg_3Dangular::LongLivedBkg_3Dangular( const LongLivedBkg_3Dangular& input ) : BasePDF( (BasePDF) input ),
	f_LL1Name(input.f_LL1Name), tauLL1Name(input.tauLL1Name), tauLL2Name(input.tauLL2Name), timeResLL1FracName(input.timeResLL1FracName), sigmaLL1Name(input.sigmaLL1Name),
	sigmaLL2Name(input.sigmaLL2Name), timeName(input.timeName), timeconstName(input.timeconstName), tauLL1(input.tauLL1), tauLL2(input.tauLL2), f_LL1(input.f_LL1),
	sigmaLL(input.sigmaLL), sigmaLL1(input.sigmaLL1), sigmaLL2(input.sigmaLL2), timeResLL1Frac(input.timeResLL1Frac), tlow(input.tlow), thigh(input.thigh), time(input.time),
	histo(input.histo), xaxis(input.xaxis), yaxis(input.yaxis), zaxis(input.zaxis), nxbins(input.nxbins), nybins(input.nybins), nzbins(input.nzbins), xmin(input.xmin),
	xmax(input.xmax), ymin(input.ymin), ymax(input.ymax), zmin(input.zmin), zmax(input.zmax), deltax(input.deltax), deltay(input.deltay), deltaz(input.deltaz),
	total_num_entries(input.total_num_entries), useFlatAngularDistribution(input.useFlatAngularDistribution),
	cosThetaName(input.cosThetaName), phiName(input.phiName), cosPsiName(input.cosPsiName), cosTheta(input.cosTheta), phi(input.phi), cosPsi(input.cosPsi)
{
}

//................................................................
//Destructor
LongLivedBkg_3Dangular::~LongLivedBkg_3Dangular()
{
}

bool LongLivedBkg_3Dangular::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	f_LL1       = allParameters.GetPhysicsParameter( f_LL1Name )->GetValue();
	tauLL1      = allParameters.GetPhysicsParameter( tauLL1Name )->GetValue();
	tauLL2      = allParameters.GetPhysicsParameter( tauLL2Name )->GetValue();
	timeResLL1Frac = allParameters.GetPhysicsParameter( timeResLL1FracName )->GetValue();
	sigmaLL1    = allParameters.GetPhysicsParameter( sigmaLL1Name )->GetValue();
	sigmaLL2    = allParameters.GetPhysicsParameter( sigmaLL2Name )->GetValue();

	return isOK;
}

//..............................................................
//Main method to build the PDF return value
double LongLivedBkg_3Dangular::Evaluate(DataPoint * measurement)
{
	// Observable
	time = measurement->GetObservable( timeName )->GetValue();
	cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
	phi      = measurement->GetObservable( phiName )->GetValue();
	cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();

	double returnValue = 0;
	double val1=-1., val2=-1.;
	//Deal with propertime resolution
	if( timeResLL1Frac >= 0.9999 )
	{
		// Set the member variable for time resolution to the first value and calculate
		sigmaLL = sigmaLL1;
		returnValue =  this->buildPDFnumerator() ;
	}
	else
	{
		// Set the member variable for time resolution to the first value and calculate
		sigmaLL = sigmaLL1;
		val1 = this->buildPDFnumerator();
		// Set the member variable for time resolution to the second value and calculate
		sigmaLL = sigmaLL2;
		val2 = this->buildPDFnumerator();
		returnValue = (timeResLL1Frac*val1 + (1. - timeResLL1Frac)*val2) ;
	}

	if (returnValue <= 0)
	{
		cout << "PDF returns zero!" << endl;
	}

	return returnValue;
}


//.............................................................
// Core calculation of PDF value
double LongLivedBkg_3Dangular::buildPDFnumerator()
{
	// Sum of two exponentials, using the time resolution functions

	double returnValue = 0.;

	double val1=-1., val2=-1.;
	if( f_LL1 >= 0.9999 ) {
		if( tauLL1 <= 0 ) {
			cout << " In LongLivedBkg_3Dangular() you gave a negative or zero lifetime for tauLL1 " << endl ;
			exit(1) ;
		}
		returnValue = Mathematics::Exp(time, 1./tauLL1, sigmaLL);
	}
	else {
		if( (tauLL1 <= 0) ||  (tauLL2 <= 0) ) {
			cout << " In LongLivedBkg_3Dangular() you gave a negative or zero lifetime for tauLL1/2 " << endl ;
			exit(1) ;
		}
		val1 = Mathematics::Exp(time, 1./tauLL1, sigmaLL);
		val2 = Mathematics::Exp(time, 1./tauLL2, sigmaLL);
		returnValue = f_LL1 * val1 + (1. - f_LL1) * val2;
	}

	returnValue *= angularFactor();

	return returnValue;
}


//..............................................................
// Normlisation
double LongLivedBkg_3Dangular::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	//	Stupid gcc
	(void)measurement;

	IConstraint * timeBound = boundary->GetConstraint( timeconstName );
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		return -1.;
	}
	else
	{
		tlow = timeBound->GetMinimum();
		thigh = timeBound->GetMaximum();
	}

	double returnValue = 0;

	if( timeResLL1Frac >= 0.9999 )
	{
		// Set the member variable for time resolution to the first value and calculate
		sigmaLL = sigmaLL1;
		returnValue = buildPDFdenominator();
	}
	else
	{
		// Set the member variable for time resolution to the first value and calculate
		sigmaLL = sigmaLL1;
		double val1 = buildPDFdenominator();
		// Set the member variable for time resolution to the second value and calculate
		sigmaLL = sigmaLL2;
		double val2 = buildPDFdenominator();
		returnValue =  timeResLL1Frac*val1 + (1. - timeResLL1Frac)*val2;
	}

	return returnValue ;
}

//.............................................................
//
double LongLivedBkg_3Dangular::buildPDFdenominator()
{
	// Sum of two exponentials, using the time resolution functions

	double returnValue = 0;

	if( f_LL1 >= 0.9999 ) {
		if( tauLL1 <= 0 ) {
			cout << " In LongLivedBkg_3Dangular() you gave a negative or zero lifetime for tauLL1 " << endl ;
			exit(1) ;
		}
		returnValue = Mathematics::ExpInt(tlow, thigh, 1./tauLL1, sigmaLL);

	}
	else {
		if( (tauLL1 <= 0) ||  (tauLL2 <= 0) ) {
			cout << " In LongLivedBkg_3Dangular() you gave a negative or zero lifetime for tauLL1/2 " << endl ;
			exit(1) ;
		}
		double val1 = Mathematics::ExpInt(tlow, thigh, 1./tauLL1, sigmaLL);
		double val2 = Mathematics::ExpInt(tlow, thigh, 1./tauLL2, sigmaLL);
		returnValue = f_LL1 * val1 + (1. - f_LL1) * val2;
	}

	//This PDF only works for full angular phase space= 8pi factor included in the factors in the Evaluate() method - so no angular normalisation term.
	return returnValue ;

}


//................................................................
//Angular distribution function
double LongLivedBkg_3Dangular::angularFactor( )
{
	double returnValue=0.;

	int globalbin=-1;
	int xbin=-1, ybin=-1, zbin=-1;
	double num_entries_bin=-1.;

	if( useFlatAngularDistribution ) {
		returnValue = 1.0 / 8.0 / TMath::Pi() ;
	}
	else {

		//Find global bin number for values of angles, find number of entries per bin, divide by volume per bin and normalise with total number of entries in the histogram
		xbin = xaxis->FindFixBin( cosPsi ); if( xbin > nxbins ) xbin = nxbins;
		ybin = yaxis->FindFixBin( cosTheta ); if( ybin > nybins ) ybin = nybins;
		zbin = zaxis->FindFixBin( phi ); if( zbin > nzbins ) zbin = nzbins;

		globalbin = histo->GetBin( xbin, ybin, zbin );
		num_entries_bin = histo->GetBinContent(globalbin);

		//Angular factor normalized with phase space of histogram and total number of entries in the histogram
		returnValue = num_entries_bin / (deltax * deltay * deltaz) / total_num_entries ;
	}

	return returnValue;
}

