// $Id: DoubleExponential.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class DoubleExponential DoubleExponential.cpp
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-11-13
 */

#include "DoubleExponential.h"
#include "Mathematics.h"
#include <iostream>
#include <cmath>

using namespace::std;

PDF_CREATOR( DoubleExponential );

//Constructor
DoubleExponential::DoubleExponential( PDFConfigurator* configurator) :
	// Physics parameters
	tau1Name	( configurator->getName("tauLL1") )
	, tau2Name	( configurator->getName("tauLL2") )
	, fraction1Name	( configurator->getName("fractionLL1") )
	, eventResolutionName   ( configurator->getName("eventResolution") )
	, resScale1Name          ( configurator->getName("timeResolutionScale1") )
	, resScale2Name          ( configurator->getName("timeResolutionScale2") )
	, resScale3Name          ( configurator->getName("timeResolutionScale3") )
	, sigma1Name	( configurator->getName("timeResolution1") )
	, sigma2Name	( configurator->getName("timeResolution2") )
	, sigma3Name	( configurator->getName("timeResolution3") )
	, timeRes2FracName( configurator->getName("timeResolution2Fraction") )
	, timeRes3FracName( configurator->getName("timeResolution3Fraction") )
	, timeOffsetName( configurator->getName("timeOffset") )
	// Observables
	, timeName      ( configurator->getName("time") )
	//objects used in XML
	, tau1(), tau2(), fraction1()
	, sigma(), sigma1(), sigma2(), sigma3(), timeRes2Frac(), timeRes3Frac()
	, resolutionScale1(), resolutionScale2(), resolutionScale3()
	, timeOffset()
	, tlow(), thigh(), time()
	, _useEventResolution(false)
	, _useTimeAcceptance(false)
, _numericIntegralForce(false)
, _usePunziSigmat(false)
{
	_useEventResolution = configurator->isTrue( "UseEventResolution" );
	_useTimeAcceptance  = configurator->isTrue( "UseTimeAcceptance" );
	_numericIntegralForce = configurator->isTrue( "UseNumericalIntegration" );
        _usePunziSigmat = configurator->isTrue( "UsePunziSigmat" );
	if( useTimeAcceptance() ) {
		if( configurator->hasConfigurationValue( "TimeAcceptanceType", "Upper" ) ) {
			timeAcc = new SlicedAcceptance( 0., 14.0, 0.0033 );
			cout << "DoubleExponential:: Constructing timeAcc: Upper time acceptance beta=0.0033 [0 < t < 14] " << endl;
		}
		else if( configurator->getConfigurationValue( "TimeAcceptanceFile" ) != "" ) {
			timeAcc = new SlicedAcceptance( "File" , configurator->getConfigurationValue( "TimeAcceptanceFile" ) );
			cout << "DoubleExponential:: Constructing timeAcc: using file: " << configurator->getConfigurationValue( "TimeAcceptanceFile" ) << endl;
		}
	}
	else {
		timeAcc = new SlicedAcceptance( -25., 25. );
		cout << "DoubleExponential:: Constructing timeAcc: DEFAULT FLAT [-25 < t < 25]  " << endl;
	}
	MakePrototypes();
}


DoubleExponential::DoubleExponential( const DoubleExponential &copy ) :
        BasePDF( (BasePDF)copy )
        , tau1Name ( copy.tau1Name )
        , tau2Name ( copy.tau2Name )
        , fraction1Name ( copy.fraction1Name )
	, eventResolutionName   ( copy.eventResolutionName )
	, resScale1Name          ( copy.resScale1Name )
	, resScale2Name          ( copy.resScale2Name )
	, resScale3Name          ( copy.resScale3Name )
	, sigma1Name	( copy.sigma1Name )
	, sigma2Name	( copy.sigma2Name )
	, sigma3Name	( copy.sigma3Name )
	, timeRes2FracName( copy.timeRes2FracName )
	, timeRes3FracName( copy.timeRes3FracName )
	, timeOffsetName( copy.timeOffsetName )
	, timeName      ( copy.timeName )
	, tau1( copy.tau1 )
	, tau2( copy.tau2 )
	, fraction1( copy.fraction1 )
	, sigma( copy.sigma ), sigma1( copy.sigma1 ), sigma2( copy.sigma2 ), sigma3( copy.sigma3 )
	, timeRes2Frac( copy.timeRes2Frac), timeRes3Frac( copy.timeRes3Frac )
        , resolutionScale1( copy.resolutionScale1 ), resolutionScale2( copy.resolutionScale2 ), resolutionScale3( copy.resolutionScale3 )
        , timeOffset( copy.timeOffset )
	, tlow( copy.tlow ), thigh( copy.thigh ), time( copy.time )
	, _useEventResolution( copy._useEventResolution )
        , _useTimeAcceptance( copy._useTimeAcceptance )
	, _numericIntegralForce( copy._numericIntegralForce )
	, _usePunziSigmat( copy._usePunziSigmat )
{
        timeAcc = new SlicedAcceptance( *(copy.timeAcc) );
}


//Make the data point and parameter set
void DoubleExponential::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	if(useEventResolution())
	{
		allObservables.push_back( eventResolutionName );
		this->TurnCachingOff();
	}

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( tau1Name );
	parameterNames.push_back( tau2Name );
	parameterNames.push_back( fraction1Name );
	parameterNames.push_back( timeRes2FracName );
	parameterNames.push_back( timeRes3FracName );
	if( ! useEventResolution() ) {
		parameterNames.push_back( sigma1Name );
		parameterNames.push_back( sigma2Name );
		parameterNames.push_back( sigma3Name );
	}
	parameterNames.push_back( resScale1Name );
	parameterNames.push_back( resScale2Name );
	parameterNames.push_back( resScale3Name );
	parameterNames.push_back( timeOffsetName );
	allParameters = ParameterSet(parameterNames);
}

//Return a list of observables not to be integrated
vector<string> DoubleExponential::GetDoNotIntegrateList()
{
	vector<string> list;
	if( useEventResolution() && !_usePunziSigmat ) list.push_back(eventResolutionName);
	return list;
}

//Destructor
DoubleExponential::~DoubleExponential()
{
}

bool DoubleExponential::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	timeRes2Frac = allParameters.GetPhysicsParameter( timeRes2FracName )->GetValue();
	timeRes3Frac = allParameters.GetPhysicsParameter( timeRes3FracName )->GetValue();
	tau1 = allParameters.GetPhysicsParameter( tau1Name )->GetValue();
	tau2 = allParameters.GetPhysicsParameter( tau2Name )->GetValue();
	fraction1 = allParameters.GetPhysicsParameter( fraction1Name )->GetValue();
	if( ! useEventResolution() ) {
		sigma1    = allParameters.GetPhysicsParameter( sigma1Name )->GetValue();
		sigma2    = allParameters.GetPhysicsParameter( sigma2Name )->GetValue();
		sigma3    = allParameters.GetPhysicsParameter( sigma3Name )->GetValue();
	} 
	resolutionScale1 = allParameters.GetPhysicsParameter( resScale1Name )->GetValue();
	resolutionScale2 = allParameters.GetPhysicsParameter( resScale2Name )->GetValue();
	resolutionScale3 = allParameters.GetPhysicsParameter( resScale3Name )->GetValue();
	timeOffset       = allParameters.GetPhysicsParameter( timeOffsetName )->GetValue();

	return isOK;
}

//Calculate the function value
double DoubleExponential::Evaluate(DataPoint * measurement)
{
	// Observable
	time = measurement->GetObservable( timeName )->GetValue() - timeOffset;
	if( useEventResolution() ) eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();

	double num = 0.;

	if( resolutionScale1 <= 0. ) {
		//This is the "code" to run with resolution=0
		sigma = 0.;
		num = buildPDFnumerator();
	}
	else if( useEventResolution() ) {
		// Event-by-event resolution has been selected
		if( (1. - timeRes2Frac - timeRes3Frac ) >= 0.9999 )
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
			num = (1. - timeRes2Frac - timeRes3Frac)*val1 + timeRes2Frac*val2 + timeRes3Frac*val3;
		}
	}
	else {
		if(  (1. - timeRes2Frac - timeRes3Frac ) >= 0.9999 )
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
			// Set the member variable for time resolution to the second value and calculate
			sigma = sigma3;
			double val3 = buildPDFnumerator();
			num = (1. - timeRes2Frac - timeRes3Frac)*val1 + timeRes2Frac*val2 + timeRes3Frac*val3;
		}
	}
	Observable * timeObs = measurement->GetObservable( timeName );
	if( useTimeAcceptance() ) num = num * timeAcc->getValue(timeObs, timeOffset);
	//cout << eventResolution << endl;
	return num;
}

double DoubleExponential::buildPDFnumerator()
{
	// Sum of two exponentials, using the time resolution functions

	if( tau1 <= 0 || tau2 <= 0 ) {
		cout << " In DoubleExponential() you gave a negative or zero lifetime for tau " << endl ;
		throw(10) ;
	}
	double val = fraction1*Mathematics::Exp(time, 1./tau1, sigma);
	val += (1.-fraction1)*Mathematics::Exp(time, 1./tau2, sigma);
	return val;
}

double DoubleExponential::Normalisation( PhaseSpaceBoundary* boundary )
{
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
		if(  (1. - timeRes2Frac - timeRes3Frac ) >= 0.9999 )
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
                        sigma = sigma3;
                        double val3 = buildPDFdenominator();
                        norm = (1. - timeRes2Frac - timeRes3Frac)*val1 + timeRes2Frac*val2 + timeRes3Frac*val3;
		}
	}
	return norm;
}

double DoubleExponential::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
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
		if(  (1. - timeRes2Frac - timeRes3Frac ) >= 0.9999 )
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
			norm =  (1. - timeRes2Frac - timeRes3Frac )*val1 + timeRes2Frac*val2 + timeRes3Frac*val3;
		}
	}
	else {
		if(  (1. - timeRes2Frac - timeRes3Frac ) >= 0.9999 )
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
                        sigma = sigma3;
                        double val3 = buildPDFdenominator();
                        norm = (1. - timeRes2Frac - timeRes3Frac)*val1 + timeRes2Frac*val2 + timeRes3Frac*val3;
		}
	}
	return norm;
}

double DoubleExponential::buildPDFdenominator()
{
	// Sum of two exponentials, using the time resolution functions

	if( tau1 <= 0 || tau2 <= 0 ) {
		cout << " In DoubleExponential() you gave a negative or zero lifetime for tau " << endl ;
		throw(10) ;
	}

	double tlo_boundary = tlow;
	double thi_boundary = thigh;
	double val = 0;

	for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
	{
		tlow  = tlo_boundary > timeAcc->getSlice(islice)->tlow() ? tlo_boundary : timeAcc->getSlice(islice)->tlow();
		thigh = thi_boundary < timeAcc->getSlice(islice)->thigh() ? thi_boundary : timeAcc->getSlice(islice)->thigh();
		if( thigh > tlow ) {
			val += fraction1*Mathematics::ExpInt(tlow, thigh, 1./tau1, sigma) * timeAcc->getSlice(islice)->height();
			val += (1.-fraction1)*Mathematics::ExpInt(tlow, thigh, 1./tau2, sigma) * timeAcc->getSlice(islice)->height();
		}
	}

	tlow  = tlo_boundary;
	thigh = thi_boundary;
	return val;
	//return Mathematics::ExpInt(tlow, thigh, 1./tau, sigma);
}

