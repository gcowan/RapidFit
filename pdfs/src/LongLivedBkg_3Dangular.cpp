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

//Constructor
LongLivedBkg_3Dangular::LongLivedBkg_3Dangular() :

	// Physics parameters
	  tauLL1Name	( "tau_LL1" )
	, tauLL2Name	( "tau_LL2" )
	, f_LL1Name      ( "f_LL1" )
        , sigmaLL1Name  ( "sigma_LL1" )
        , sigmaLL2Name  ( "sigma_LL2" )
        , timeResLL1FracName( "timeResLL1Frac" )

// Observables
        , timeName      ( "time" )
        , cosThetaName  ( "cosTheta" )
        , phiName       ( "phi" )
        , cosPsiName    ( "cosPsi" )

	, f(TFile::Open("output.root"))

{
	MakePrototypes();

	//Read in histo
	histo = new TH3D(*((TH3D*)f->Get("histo")));

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

	if ((xmax-xmin) < 2 || (ymax-ymin) < 2 || (zmax-zmin) < 2*TMath::Pi() ){
	cout << "The full angular range is not used in this histogram - the PDF LongLovedBkg_3Dangular does not support this case" << endl;
	exit(1);
	}

}

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

//Destructor
LongLivedBkg_3Dangular::~LongLivedBkg_3Dangular()
{
}

bool LongLivedBkg_3Dangular::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
        bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
        tauLL1      = allParameters.GetPhysicsParameter( tauLL1Name )->GetValue();
        tauLL2      = allParameters.GetPhysicsParameter( tauLL2Name )->GetValue();
        f_LL1       = allParameters.GetPhysicsParameter( f_LL1Name )->GetValue();
        sigmaLL1    = allParameters.GetPhysicsParameter( sigmaLL1Name )->GetValue();
        sigmaLL2    = allParameters.GetPhysicsParameter( sigmaLL2Name )->GetValue();
        timeResLL1Frac = allParameters.GetPhysicsParameter( timeResLL1FracName )->GetValue();

	return isOK;
}

//Calculate the function value
double LongLivedBkg_3Dangular::Evaluate(DataPoint * measurement)
{
	// Observable
	time = measurement->GetObservable( timeName )->GetValue();
	cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
	phi      = measurement->GetObservable( phiName )->GetValue();
	cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();


	if( timeResLL1Frac >= 0.9999 )
        {
                // Set the member variable for time resolution to the first value and calculate
                sigmaLL = sigmaLL1;
                return buildPDFnumerator() * angularFactor(cosTheta, cosPsi, phi);
        }
        else
        {
                // Set the member variable for time resolution to the first value and calculate
                sigmaLL = sigmaLL1;
                double val1 = buildPDFnumerator();
                // Set the member variable for time resolution to the second value and calculate
                sigmaLL = sigmaLL2;
                double val2 = buildPDFnumerator();
                return (timeResLL1Frac*val1 + (1. - timeResLL1Frac)*val2) * angularFactor(cosTheta, cosPsi, phi);
        }
}

double LongLivedBkg_3Dangular::buildPDFnumerator()
{
	// Sum of two exponentials, using the time resolution functions

	if( f_LL1 >= 0.9999 ) {
		if( tauLL1 <= 0 ) {
			cout << " In LongLivedBkg_3Dangular() you gave a negative or zero lifetime for tauLL1 " << endl ;
			exit(1) ;
		}
		double val = Mathematics::Exp(time, 1./tauLL1, sigmaLL);
		return val;
	}
	else {
		if( (tauLL1 <= 0) ||  (tauLL2 <= 0) ) {
			cout << " In LongLivedBkg_3Dangular() you gave a negative or zero lifetime for tauLL1/2 " << endl ;
			exit(1) ;
		}
		double val1 = Mathematics::Exp(time, 1./tauLL1, sigmaLL);
		double val2 = Mathematics::Exp(time, 1./tauLL2, sigmaLL);
		double val = f_LL1 * val1 + (1. - f_LL1) * val2;
		return val;
	}
}

double LongLivedBkg_3Dangular::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	//	Stupid gcc
	(void)measurement;
        IConstraint * timeBound = boundary->GetConstraint( timeName );
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

	//This PDF only works for full angular phase space= 8pi factor included in the Evaluate() method
	return returnValue ;
}

double LongLivedBkg_3Dangular::buildPDFdenominator()
{
	// Sum of two exponentials, using the time resolution functions

	if( f_LL1 >= 0.9999 ) {
		if( tauLL1 <= 0 ) {
			cout << " In LongLivedBkg_3Dangular() you gave a negative or zero lifetime for tauLL1 " << endl ;
			exit(1) ;
		}
		double val = Mathematics::ExpInt(tlow, thigh, 1./tauLL1, sigmaLL);
		return val;
	}
	else {
		if( (tauLL1 <= 0) ||  (tauLL2 <= 0) ) {
			cout << " In LongLivedBkg_3Dangular() you gave a negative or zero lifetime for tauLL1/2 " << endl ;
			exit(1) ;
		}
		double val1 = Mathematics::ExpInt(tlow, thigh, 1./tauLL1, sigmaLL);
		double val2 = Mathematics::ExpInt(tlow, thigh, 1./tauLL2, sigmaLL);
		double val = f_LL1 * val1 + (1. - f_LL1) * val2;
		return val;
	}

}


double LongLivedBkg_3Dangular::angularFactor(double cosTheta, double cosPsi, double phi )
{

	//Find global bin number for values of angles, find number of entries per bin, divide by volume per bin and normalise with total number of entries in the histogram
	int globalbin = histo->FindBin(cosTheta, cosPsi, phi);
	double num_entries_bin = histo->GetBinContent(globalbin);

	//Angular factor normalized with phase space of histogram and total number of entries in the histogram
	double factor = num_entries_bin / (deltax * deltay * deltaz) / total_num_entries ;

	return factor;

}
