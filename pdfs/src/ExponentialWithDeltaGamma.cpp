// $Id: ExponentialWithDeltaGamma.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class ExponentialWithDeltaGamma ExponentialWithDeltaGamma.cpp
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-11-13
 */

#include "ExponentialWithDeltaGamma.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"

PDF_CREATOR( ExponentialWithDeltaGamma );

//Constructor
ExponentialWithDeltaGamma::ExponentialWithDeltaGamma( PDFConfigurator* configurator) :

	// Physics parameters
	tauName	( configurator->getName("tau") )
	, deltaGammaName( configurator->getName("deltaGamma") )
	, sin2betaName	( configurator->getName("sin2beta") )
	, resScale1Name          ( configurator->getName("timeResolutionScale1") )
	, resScale2Name          ( configurator->getName("timeResolutionScale2") )
	, resScale3Name          ( configurator->getName("timeResolutionScale3") )
	, eventResolutionName   ( configurator->getName("eventResolution") )
	, sigma1Name	( configurator->getName("timeResolution1") )
	, sigma2Name	( configurator->getName("timeResolution2") )
	, timeRes1FracName( configurator->getName("timeResolution1Fraction") )
	, timeRes2FracName( configurator->getName("timeResolution2Fraction") )
	// Observables
	, timeName      ( configurator->getName("time") )
	//objects used in XML
	, tau(), deltaGamma(), sin2beta(), taueff()
	, sigma(), sigma1(), sigma2(), timeRes1Frac(), timeRes2Frac()
	, resolutionScale1(), resolutionScale2(), resolutionScale3()
	, tlow(), thigh(), time()
	, _useEventResolution(false)
	, _useTimeAcceptance(false)
, _numericIntegralForce(false)
{
	_useEventResolution = configurator->isTrue( "UseEventResolution" );
	_useTimeAcceptance  = configurator->isTrue( "UseTimeAcceptance" );
	_numericIntegralForce = configurator->isTrue( "UseNumericalIntegration" );
	_useTimeAcceptance = configurator->isTrue( "UseTimeAcceptance" );
	if( useTimeAcceptance() ) {
		if( configurator->hasConfigurationValue( "TimeAcceptanceType", "Upper" ) ) {
			timeAcc = new SlicedAcceptance( 0., 14.0, 0.011 );
			cout << "ExponentialWithDeltaGamma:: Constructing timeAcc: Upper time acceptance beta=0.015 [0 < t < 14] " << endl;
		}
		else if( configurator->getConfigurationValue( "TimeAcceptanceFile" ) != "" ) {
			timeAcc = new SlicedAcceptance( "File" , configurator->getConfigurationValue( "TimeAcceptanceFile" ) );
			cout << "ExponentialWithDeltaGamma:: Constructing timeAcc: using file: " << configurator->getConfigurationValue( "TimeAcceptanceFile" ) << endl;
		}
	}
	else {
		timeAcc = new SlicedAcceptance( -1.5, 14. );
		cout << "ExponentialWithDeltaGamma:: Constructing timeAcc: DEFAULT FLAT [-1.5 < t < 14]  " << endl;
	}
	MakePrototypes();
}

//Make the data point and parameter set
void ExponentialWithDeltaGamma::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	if( useEventResolution() )
	{
		allObservables.push_back( eventResolutionName );
		this->TurnCachingOff();
	}

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( tauName );
	parameterNames.push_back( deltaGammaName );
	parameterNames.push_back( sin2betaName );
	parameterNames.push_back( resScale1Name );
	parameterNames.push_back( resScale2Name );
	parameterNames.push_back( resScale3Name );
	parameterNames.push_back( timeRes1FracName );
	parameterNames.push_back( timeRes2FracName );
	if( ! useEventResolution() ) {
		parameterNames.push_back( sigma1Name );
		parameterNames.push_back( sigma2Name );
	}
	allParameters = ParameterSet(parameterNames);
}

//Return a list of observables not to be integrated
vector<string> ExponentialWithDeltaGamma::GetDoNotIntegrateList()
{
	vector<string> list;
	if( useEventResolution() ) list.push_back(eventResolutionName);
	return list;
}

//Destructor
ExponentialWithDeltaGamma::~ExponentialWithDeltaGamma()
{
}

bool ExponentialWithDeltaGamma::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	resolutionScale1 = allParameters.GetPhysicsParameter( resScale1Name )->GetValue();
	resolutionScale2 = allParameters.GetPhysicsParameter( resScale2Name )->GetValue();
	resolutionScale3 = allParameters.GetPhysicsParameter( resScale3Name )->GetValue();
	timeRes1Frac = allParameters.GetPhysicsParameter( timeRes1FracName )->GetValue();
	timeRes2Frac = allParameters.GetPhysicsParameter( timeRes2FracName )->GetValue();
	tau = allParameters.GetPhysicsParameter( tauName )->GetValue();
	deltaGamma = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	sin2beta = allParameters.GetPhysicsParameter( sin2betaName )->GetValue();
	// Currently fitting for DeltaGamma_d/Gamma_d with this parameterisation
	taueff = tau/(1. - sqrt(1.-sin2beta*sin2beta) * deltaGamma/2.);// * tau );
	if( ! useEventResolution() ) {
		sigma1    = allParameters.GetPhysicsParameter( sigma1Name )->GetValue();
		sigma2    = allParameters.GetPhysicsParameter( sigma2Name )->GetValue();
	}

	return isOK;
}

//Calculate the function value
double ExponentialWithDeltaGamma::Evaluate(DataPoint * measurement)
{
	// Observable
	time = measurement->GetObservable( timeName )->GetValue();
	if( useEventResolution() ) eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();

	double num = 0.;

	if( resolutionScale1 <= 0. ) {
		//This is the "code" to run with resolution=0
		sigma = 0.;
		num = buildPDFnumerator();
	}
	else if( useEventResolution() ) {
		// Event-by-event resolution has been selected
		if( timeRes1Frac >= 0.9999 )
		{
			// Set the member variable for time resolution to the first value and calculate
			sigma = eventResolution * resolutionScale1;
			num = buildPDFnumerator();
		}
		else
		{
			// Set the member variable for time resolution to the first value and calculate
			sigma = eventResolution * resolutionScale1;
			double val1 = buildPDFnumerator();
			// Set the member variable for time resolution to the second value and calculate
			sigma = eventResolution * resolutionScale2;
			double val2 = buildPDFnumerator();
			// Set the member variable for time resolution to the second value and calculate
			sigma = eventResolution * resolutionScale3;
			double val3 = buildPDFnumerator();
			num = timeRes1Frac*val1 + timeRes2Frac*val2 + (1. - timeRes1Frac - timeRes2Frac)*val3;
		}
	}
	else {
		if( timeRes1Frac >= 0.9999 )
		{
			// Set the member variable for time resolution to the first value and calculate
			sigma = sigma1;
			num = buildPDFnumerator();
		}
		else
		{
			// Set the member variable for time resolution to the first value and calculate
			sigma = sigma1;
			double val1 = buildPDFnumerator();
			// Set the member variable for time resolution to the second value and calculate
			sigma = sigma2;
			double val2 = buildPDFnumerator();
			num = timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
		}
	}
	if( useTimeAcceptance() ) num = num * timeAcc->getValue(time);
	return num;
}

double ExponentialWithDeltaGamma::buildPDFnumerator()
{
	// Sum of two exponentials, using the time resolution functions

	if( tau <= 0 ) {
		cout << " In ExponentialWithDeltaGamma() you gave a negative or zero lifetime for tau " << endl ;
		throw(10) ;
	}
	double val = Mathematics::Exp(time, 1./taueff, sigma);
	return val;
}

double ExponentialWithDeltaGamma::Normalisation( PhaseSpaceBoundary* boundary )
{
	if( _numericIntegralForce ) return -1.;

	double norm = 0.;
	if( useEventResolution() )
	{
		norm = -1.;
	}
	else
	{
		IConstraint * timeBound = boundary->GetConstraint( timeName );
		if ( timeBound->GetUnit() == "NameNotFoundError" )
		{
			cerr << "Bound on time not provided" << endl;
			norm = -1.;
		}
		else    
		{
			tlow = timeBound->GetMinimum();
			thigh = timeBound->GetMaximum();
		}
		if( timeRes1Frac >= 0.9999 )
		{
			// Set the member variable for time resolution to the first value and calculate
			sigma = sigma1;
			norm = buildPDFdenominator();
		}
		else
		{
			// Set the member variable for time resolution to the first value and calculate
			sigma = sigma1;
			double val1 = buildPDFdenominator();
			// Set the member variable for time resolution to the second value and calculate
			sigma = sigma2;
			double val2 = buildPDFdenominator();
			norm = timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
		}
	}
	return norm;
}

double ExponentialWithDeltaGamma::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	if( _numericIntegralForce ) return -1.;

	double norm = 0.;

	IConstraint * timeBound = boundary->GetConstraint( timeName );
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		norm = -1.;
	}
	else
	{
		tlow = timeBound->GetMinimum();
		thigh = timeBound->GetMaximum();
	}

	if( useEventResolution() )  {
		eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();
		if( timeRes1Frac >= 0.9999 )
		{
			// Set the member variable for time resolution to the first value and calculate
			sigma = eventResolution * resolutionScale1;
			norm = buildPDFdenominator();
		}
		else
		{
			// Set the member variable for time resolution to the first value and calculate
			sigma = eventResolution * resolutionScale1;
			double val1 = buildPDFdenominator();
			// Set the member variable for time resolution to the second value and calculate
			sigma = eventResolution * resolutionScale2;
			double val2 = buildPDFdenominator();
			// Set the member variable for time resolution to the second value and calculate
			sigma = eventResolution * resolutionScale3;
			double val3 = buildPDFdenominator();
			norm = timeRes1Frac*val1 + timeRes2Frac*val2 + (1. - timeRes1Frac - timeRes2Frac)*val3;
		}
	}
	else {
		norm = -1.;
	}
	return norm;
}

double ExponentialWithDeltaGamma::buildPDFdenominator()
{
	// Sum of two exponentials, using the time resolution functions

	if( taueff <= 0 ) {
		cout << " In ExponentialWithDeltaGamma() you gave a negative or zero lifetime for taueff " << endl ;
		throw(10) ;
	}

	double tlo_boundary = tlow;
	double thi_boundary = thigh;
	double val = 0;

	for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
	{
		tlow  = tlo_boundary > timeAcc->getSlice(islice)->tlow() ? tlo_boundary : timeAcc->getSlice(islice)->tlow();
		thigh = thi_boundary < timeAcc->getSlice(islice)->thigh() ? thi_boundary : timeAcc->getSlice(islice)->thigh();
		if( thigh > tlow ) val += Mathematics::ExpInt(tlow, thigh, 1./taueff, sigma) * timeAcc->getSlice(islice)->height();
	}

	tlow  = tlo_boundary;
	thigh = thi_boundary;
	return val;
	//return Mathematics::ExpInt(tlow, thigh, 1./tau, sigma);
}

