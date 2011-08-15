// $Id: LongLivedBkg.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class LongLivedBkg LongLivedBkg.cpp
 *
 *  PDF for  long lived background with 3D histogram angular description
 *
 *  @author Ailsa Sparkes
 *  @date 2011-05-30
 */

#include "LongLivedBkg.h"
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
LongLivedBkg::LongLivedBkg(PDFConfigurator config ) :

	// Physics parameters
	  f_LL1Name		( config.getName("f_LL1")  )
	, tauLL1Name		( config.getName("tau_LL1") )
	, tauLL2Name		( config.getName("tau_LL2") )
    //Detector parameters
	, timeResLL1FracName	( config.getName("timeResLL1Frac") )
	, sigmaLL1Name		( config.getName("sigma_LL1") )
	, sigmaLL2Name		( config.getName("sigma_LL2") )
	// Observables
	, timeName		( config.getName("time") )
	//, cosThetaName		( config.getName("cosTheta") )
	//, phiName		( config.getName("phi") )
	//, cosPsiName		( config.getName("cosPsi") )
	//Other things to be initialised
	, timeconstName		( config.getName("time") )

	, tauLL1(), tauLL2(), f_LL1(), sigmaLL(), sigmaLL1(), sigmaLL2(), timeResLL1Frac(), tlow(), thigh(), time(),
	histo(), xaxis(), yaxis(), zaxis(), nxbins(), nybins(), nzbins(), xmin(), xmax(), ymin(),
	ymax(), zmin(), zmax(), deltax(), deltay(), deltaz(), total_num_entries(), useFlatAngularDistribution(true)
{

	cout << "LongLivedBkg::  " ;

	MakePrototypes();

	//Find name of histogram needed to define 3-D angular distribution
	string fileName = config.getConfigurationValue( "AngularDistributionHistogram" ) ;

	//Initialise depending upon whether configuration parameter was found
	useFlatAngularDistribution = true ;

}


//..................................................................
//Make the data point and parameter set
void LongLivedBkg::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
//	allObservables.push_back( cosThetaName );
//	allObservables.push_back( phiName );
//	allObservables.push_back( cosPsiName );


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


//................................................................
//Destructor
LongLivedBkg::~LongLivedBkg()
{
}

bool LongLivedBkg::SetPhysicsParameters( ParameterSet * NewParameterSet )
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
double LongLivedBkg::Evaluate(DataPoint * measurement)
{
	// Observable
	time = measurement->GetObservable( timeName )->GetValue();
//	cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
//	phi      = measurement->GetObservable( phiName )->GetValue();
//	cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();

	double returnValue = 0;
	//Deal with propertime resolution
	if( timeResLL1Frac >= 0.9999 )
	{
		// Set the member variable for time resolution to the first value and calculate
		sigmaLL = sigmaLL1;
		returnValue =  buildPDFnumerator() ;
	}
	else
	{
		// Set the member variable for time resolution to the first value and calculate
		sigmaLL = sigmaLL1;
		double val1 = buildPDFnumerator();
		// Set the member variable for time resolution to the second value and calculate
		sigmaLL = sigmaLL2;
		double val2 = buildPDFnumerator();
		returnValue = (timeResLL1Frac*val1 + (1. - timeResLL1Frac)*val2) ;
	}

	if (returnValue <= 0) cout << "PDF returns zero!" << endl;

	return returnValue;

}


//.............................................................
// Core calculation of PDF value
double LongLivedBkg::buildPDFnumerator()
{
	// Sum of two exponentials, using the time resolution functions

	double returnValue = 0;

	if( f_LL1 >= 0.9999 ) {
		if( tauLL1 <= 0 ) {
			cout << " In LongLivedBkg() you gave a negative or zero lifetime for tauLL1 " << endl ;
			exit(1) ;
		}
		returnValue = Mathematics::Exp(time, 1./tauLL1, sigmaLL);
	}
	else {
		if( (tauLL1 <= 0) ||  (tauLL2 <= 0) ) {
			cout << " In LongLivedBkg() you gave a negative or zero lifetime for tauLL1/2 " << endl ;
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
double LongLivedBkg::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
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
double LongLivedBkg::buildPDFdenominator()
{
	// Sum of two exponentials, using the time resolution functions

	double returnValue = 0;

	if( f_LL1 >= 0.9999 ) {
		if( tauLL1 <= 0 ) {
			cout << " In LongLivedBkg() you gave a negative or zero lifetime for tauLL1 " << endl ;
			exit(1) ;
		}
		returnValue = Mathematics::ExpInt(tlow, thigh, 1./tauLL1, sigmaLL);

	}
	else {
		if( (tauLL1 <= 0) ||  (tauLL2 <= 0) ) {
			cout << " In LongLivedBkg() you gave a negative or zero lifetime for tauLL1/2 " << endl ;
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
double LongLivedBkg::angularFactor( )
{

		return 1.0 / 8.0 / TMath::Pi() ;
}
