// $Id: LongLivedBkg.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class LongLivedBkg LongLivedBkg.cpp
 *
 *  PDF for  long lived background with 3D histogram angular description
 *
 *  @author Ailsa Sparkes
 *  @date 2011-05-30
 */

#include "TMath.h"
#include <cmath>

#include "LongLivedBkg.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH3D.h"
#include "TAxis.h"
#include "TH1.h"

PDF_CREATOR( LongLivedBkg );

//.....................................................
//Constructor
LongLivedBkg::LongLivedBkg(PDFConfigurator* configurator ) :
	// Physics parameters
	  f_LL1Name				( configurator->getName("f_LL1")  )
	, tauLL1Name			( configurator->getName("tau_LL1") )
	, tauLL2Name			( configurator->getName("tau_LL2") )
    //Detector parameters
	, timeResLL1FracName	( configurator->getName("timeResLL1Frac") )
	, sigmaLL1Name			( configurator->getName("sigma_LL1") )
	, sigmaLL2Name			( configurator->getName("sigma_LL2") )
	// Observables
	, timeName				( configurator->getName("time") )
	//Other things to be initialised
	, _useTimeAcceptance(false)
	, tauLL1(), tauLL2(), f_LL1() 
	, sigmaLL(), sigmaLL1(), sigmaLL2(), timeResLL1Frac()
	, tlow(), thigh(), time()
	, timeAcc(NULL)
{
	cout << "Constructing LongLivedBkg::  " << endl ;

	//...........................................
        // Configure to use time acceptance machinery 
        _useTimeAcceptance = configurator->isTrue( "UseTimeAcceptance" ) ;
        if( _useTimeAcceptance ) {
                        timeAcc = new SlicedAcceptance( "File" , configurator->getConfigurationValue( "TimeAcceptanceFile" ) ) ;
                        cout << "LongLivedBkg:: Constructing timeAcc: using file: " << configurator->getConfigurationValue( "TimeAcceptanceFile" ) << endl ;
        }

	//...............
	MakePrototypes();
	
}


//..................................................................
//Make the data point and parameter set
void LongLivedBkg::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( f_LL1Name );
	parameterNames.push_back( tauLL1Name );
	parameterNames.push_back( tauLL2Name );
	parameterNames.push_back( timeResLL1FracName );
	parameterNames.push_back( sigmaLL1Name );
	parameterNames.push_back( sigmaLL2Name );

	allParameters = ParameterSet(parameterNames);
}

//LongLivedBkg::LongLivedBkg( const LongLivedBkg& input ) : BasePDF( (BasePDF) input ),
//	f_LL1Name(input.f_LL1Name), tauLL1Name(input.tauLL1Name), tauLL2Name(input.tauLL2Name), timeResLL1FracName(input.timeResLL1FracName), sigmaLL1Name(input.sigmaLL1Name),
//	sigmaLL2Name(input.sigmaLL2Name), timeName(input.timeName), timeconstName(input.timeconstName), tauLL1(input.tauLL1), tauLL2(input.tauLL2), f_LL1(input.f_LL1),
//	sigmaLL(input.sigmaLL), sigmaLL1(input.sigmaLL1), sigmaLL2(input.sigmaLL2), timeResLL1Frac(input.timeResLL1Frac), tlow(input.tlow), thigh(input.thigh), time(input.time),
//	histo(input.histo), xaxis(input.xaxis), yaxis(input.yaxis), zaxis(input.zaxis), nxbins(input.nxbins), nybins(input.nybins), nzbins(input.nzbins), xmin(input.xmin),
//	xmax(input.xmax), ymin(input.ymin), ymax(input.ymax), zmin(input.zmin), zmax(input.zmax), deltax(input.deltax), deltay(input.deltay), deltaz(input.deltaz),
//	total_num_entries(input.total_num_entries), useFlatAngularDistribution(input.useFlatAngularDistribution), _useTimeAcceptance(input._useTimeAcceptance), timeAcc(NULL)
//{
//	timeAcc = new SlicedAcceptance( *(input.timeAcc) );
//}

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
	if(std::isnan(returnValue)) {
		cout << "PDF returns nan!  " << returnValue << endl;
		cout << "   f_LL1    "  << f_LL1 << endl;
		cout << "   tauLL1   "  << tauLL1 << endl;
		cout << "   tauLL2   "  << tauLL2 << endl;
		cout << "   sigmaLL  " << sigmaLL << endl;
		cout << "   time     "  << time << endl;
		exit(1);
	}

	
	if( _useTimeAcceptance ) returnValue = returnValue * timeAcc->getValue(time);
	
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
		returnValue = f_LL1 * val1/tauLL1 + (1. - f_LL1) * val2/tauLL2;
	}

	return returnValue;
}


//..............................................................
// Normlisation
double LongLivedBkg::Norm(PhaseSpaceBoundary * boundary)
{
	(void)boundary;

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

	return returnValue;
}

double LongLivedBkg::Normalisation(PhaseSpaceBoundary * boundary)
{
	// Use this if you want to ignore the time acceptance calculation
	//return Norm( measurement, boundary );

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

	double tlo_boundary = tlow;
	double thi_boundary = thigh;
	
	if( _useTimeAcceptance ) {
		//This loops over each time slice, does the normalisation between the limits, and accumulates
		for( unsigned int islice = 0; islice < timeAcc->numberOfSlices(); ++islice )
		{
			tlow = tlo_boundary > timeAcc->getSlice(islice)->tlow() ? tlo_boundary : timeAcc->getSlice(islice)->tlow() ;
			thigh = thi_boundary < timeAcc->getSlice(islice)->thigh() ? thi_boundary : timeAcc->getSlice(islice)->thigh() ;
			if( thigh > tlow ) returnValue += this->Norm( boundary ) * timeAcc->getSlice(islice)->height() ;
		}
	}
	else {
		returnValue = this->Norm( boundary );
	}

	tlow  = tlo_boundary;
	thigh = thi_boundary;

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
		returnValue = f_LL1 * val1/tauLL1 + (1. - f_LL1) * val2/tauLL2;
	}

	//This PDF only works for full angular phase space= 8pi factor included in the factors in the Evaluate() method - so no angular normalisation term.
	return returnValue ;
}


