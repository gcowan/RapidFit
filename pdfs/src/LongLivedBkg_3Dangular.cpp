// $Id: LongLivedBkg_3Dangular.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class LongLivedBkg_3Dangular LongLivedBkg_3Dangular.cpp
 *
 *  PDF for  long lived background with 3D histogram angular description
 *
 *  @author Ailsa Sparkes
 *  @date 2011-05-30
 */

#include "LongLivedBkg_3Dangular.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TH3D.h"
#include "TAxis.h"
#include "TH1.h"


//.....................................................
//Constructor
LongLivedBkg_3Dangular::LongLivedBkg_3Dangular(PDFConfigurator config ) :

	// Physics parameters
	  f_LL1Name     ( make_pair(config.getName("f_LL1"),-1)  )
	, tauLL1Name	( make_pair(config.getName("tau_LL1"),-1) )
	, tauLL2Name	( make_pair(config.getName("tau_LL2"),-1) )
    //Detector parameters
	, timeResLL1FracName( make_pair(config.getName("timeResLL1Frac"),-1) )
	, sigmaLL1Name  ( make_pair(config.getName("sigma_LL1"),-1) )
	, sigmaLL2Name  ( make_pair(config.getName("sigma_LL2"),-1) )
	// Observables
	, timeName      ( make_pair(config.getName("time"),-1) )
	, cosThetaName	( make_pair(config.getName("cosTheta"),-1) )
	, phiName	    ( make_pair(config.getName("phi"),-1) )
	, cosPsiName	( make_pair(config.getName("cosPsi"),-1) )
	//Other things to be initialised


	, tauLL1(), tauLL2(), f_LL1(), sigmaLL(), sigmaLL1(), sigmaLL2(), timeResLL1Frac(), tlow(), thigh(), time(), cosTheta(),
	phi(), cosPsi(), histo(), xaxis(), yaxis(), zaxis(), nxbins(), nybins(), nzbins(), xmin(), xmax(), ymin(),
	ymax(), zmin(), zmax(), deltax(), deltay(), deltaz(), total_num_entries(), useFlatAngularDistribution(true)
{

	cout << "Constructing PDF:LongLivedBkg_3Dangular" << endl ;

	MakePrototypes();

	//Find name of histogram needed to define 3-D angular distribution
	string fileName = config.getConfigurationValue( "AngularDistributionHistogram" ) ;

	//Initialise depending upon whether configuration parameter was found
	if( fileName == "" ) {
		cout << "   No AngularDistributionHistogram found: using flat background " << endl ;
		useFlatAngularDistribution = true ;
	}
	else {
		cout << "   AngularDistributionHistogram found: " << fileName << endl ;
		useFlatAngularDistribution = false ;

		//Read in histo
		TFile* f =  TFile::Open(fileName.c_str());
		histo = new TH3D(  *(    (TH3D*)f ->Get("histo")      )     ); //(fileName.c_str())));

	xaxis = histo->GetXaxis();
        xmin = xaxis->GetXmin();
        xmax = xaxis->GetXmax();
        nxbins = histo->GetNbinsX();
        deltax = (xmax-xmin)/nxbins;

        yaxis = histo->GetYaxis();
        ymin = yaxis->GetXmin();
        ymax = yaxis->GetXmax();
        nybins = histo->GetNbinsY();
        deltay = (ymax-ymin)/nybins;

        zaxis = histo->GetZaxis();
        zmin = zaxis->GetXmin();
        zmax = zaxis->GetXmax();
        nzbins = histo->GetNbinsZ();
        deltaz = (zmax-zmin)/nzbins;

		total_num_entries = histo->GetEntries();


		if ((xmax-xmin) < 2. || (ymax-ymin) < 2. || (zmax-zmin) < 2.*TMath::Pi() ){
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
	allObservables.push_back( timeName.first );
	allObservables.push_back( cosThetaName.first );
	allObservables.push_back( phiName.first );
	allObservables.push_back( cosPsiName.first );


	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( f_LL1Name.first );
	parameterNames.push_back( tauLL1Name.first );
	parameterNames.push_back( tauLL2Name.first );
	parameterNames.push_back( timeResLL1FracName.first );
	parameterNames.push_back( sigmaLL1Name.first );
	parameterNames.push_back( sigmaLL2Name.first );

	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}


//................................................................
//Destructor
LongLivedBkg_3Dangular::~LongLivedBkg_3Dangular()
{
}

bool LongLivedBkg_3Dangular::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
        bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
		f_LL1       = allParameters.GetPhysicsParameter( &f_LL1Name )->GetValue();
        tauLL1      = allParameters.GetPhysicsParameter( &tauLL1Name )->GetValue();
        tauLL2      = allParameters.GetPhysicsParameter( &tauLL2Name )->GetValue();
		timeResLL1Frac = allParameters.GetPhysicsParameter( &timeResLL1FracName )->GetValue();
        sigmaLL1    = allParameters.GetPhysicsParameter( &sigmaLL1Name )->GetValue();
        sigmaLL2    = allParameters.GetPhysicsParameter( &sigmaLL2Name )->GetValue();

	return isOK;
}

//..............................................................
//Main method to build the PDF return value
double LongLivedBkg_3Dangular::Evaluate(DataPoint * measurement)
{
	// Observable
	time = measurement->GetObservable( &timeName )->GetValue();
	cosTheta = measurement->GetObservable( &cosThetaName )->GetValue();
	phi      = measurement->GetObservable( &phiName )->GetValue();
	cosPsi   = measurement->GetObservable( &cosPsiName )->GetValue();

	//Deal with propertime resolution
	if( timeResLL1Frac >= 0.9999 )
	{
		// Set the member variable for time resolution to the first value and calculate
		sigmaLL = sigmaLL1;
		return buildPDFnumerator() ;
	}
	else
	{
		// Set the member variable for time resolution to the first value and calculate
		sigmaLL = sigmaLL1;
		double val1 = buildPDFnumerator();
		// Set the member variable for time resolution to the second value and calculate
		sigmaLL = sigmaLL2;
		double val2 = buildPDFnumerator();
		return (timeResLL1Frac*val1 + (1. - timeResLL1Frac)*val2) ;
	}
}


//.............................................................
// Core calculation of PDF value
double LongLivedBkg_3Dangular::buildPDFnumerator()
{
	// Sum of two exponentials, using the time resolution functions

	double returnValue = 0;

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
		double val1 = Mathematics::Exp(time, 1./tauLL1, sigmaLL);
		double val2 = Mathematics::Exp(time, 1./tauLL2, sigmaLL);
		returnValue = f_LL1 * val1 + (1. - f_LL1) * val2;
	}

	return returnValue * angularFactor();
}


//..............................................................
// Normlisation
double LongLivedBkg_3Dangular::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	//	Stupid gcc
	(void)measurement;

	IConstraint * timeBound = boundary->GetConstraint( &timeName );
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

	if( useFlatAngularDistribution ) {
		return 1.0 / 8.0 / TMath::Pi() ;
	}
	else {

		//Find global bin number for values of angles, find number of entries per bin, divide by volume per bin and normalise with total number of entries in the histogram
		int globalbin = histo->FindBin(cosTheta, cosPsi, phi);
		double num_entries_bin = histo->GetBinContent(globalbin);

		//Angular factor normalized with phase space of histogram and total number of entries in the histogram
		double factor = num_entries_bin / (deltax * deltay * deltaz) / total_num_entries ;

		return factor;
	}

}
