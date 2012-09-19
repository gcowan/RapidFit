// $Id: Exponential.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Exponential Exponential.cpp
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-11-13
 */

#include "Exponential.h"
#include "Mathematics.h"
#include "SlicedAcceptance.h"
#include <iostream>
#include <cmath>

using namespace::std;

PDF_CREATOR( Exponential );

//Constructor
Exponential::Exponential( PDFConfigurator* configurator) :
	// Physics parameters
	tauName	( configurator->getName("tau") )
	, eventResolutionName   ( configurator->getName("eventResolution") )
	, resScale1Name          ( configurator->getName("timeResolutionScale1") )
	, resScale2Name          ( configurator->getName("timeResolutionScale2") )
	, resScale3Name          ( configurator->getName("timeResolutionScale3") )
	, sigma1Name	( configurator->getName("timeResolution1") )
	, sigma2Name	( configurator->getName("timeResolution2") )
	, sigma3Name	( configurator->getName("timeResolution3") )
	, timeRes2FracName( configurator->getName("timeResolution2Fraction") )
	, timeRes3FracName( configurator->getName("timeResolution3Fraction") )
	// Observables
	, timeName      ( configurator->getName("time") )
	//objects used in XML
	, tau(), gamma(), sigmaNum(0), sigma(), sigma1(), sigma2(), sigma3(), timeRes2Frac(), timeRes3Frac()
	, resolutionScale1(), resolutionScale2(), resolutionScale3(), _dataPoint(NULL)
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
			cout << "Exponential:: Constructing timeAcc: Upper time acceptance beta=0.0033 [0 < t < 14] " << endl;
		}
		else if( configurator->getConfigurationValue( "TimeAcceptanceFile" ) != "" ) {
			timeAcc = new SlicedAcceptance( "File" , configurator->getConfigurationValue( "TimeAcceptanceFile" ) );
			cout << "Exponential:: Constructing timeAcc: using file: " << configurator->getConfigurationValue( "TimeAcceptanceFile" ) << endl;
		}
	}
	else {
		timeAcc = new SlicedAcceptance( -25., 25. );
		cout << "Exponential:: Constructing timeAcc: DEFAULT FLAT [-25 < t < 25]  " << endl;
	}
	this->MakePrototypes();
	this->prepareTimeInt();
}


Exponential::Exponential( const Exponential &copy ) :
	BasePDF( (BasePDF)copy )
	, tauName ( copy.tauName )
	, eventResolutionName   ( copy.eventResolutionName )
	, resScale1Name          ( copy.resScale1Name )
	, resScale2Name          ( copy.resScale2Name )
	, resScale3Name          ( copy.resScale3Name )
	, sigma1Name	( copy.sigma1Name )
	, sigma2Name	( copy.sigma2Name )
	, sigma3Name	( copy.sigma3Name )
	, timeRes2FracName( copy.timeRes2FracName )
	, timeRes3FracName( copy.timeRes3FracName )
	, timeName      ( copy.timeName )
	, tau( copy.tau ), gamma( copy.gamma ), sigmaNum( copy.sigmaNum ), _intexpIntObs_vec(), _dataPoint( copy._dataPoint )
	, sigma( copy.sigma ), sigma1( copy.sigma1 ), sigma2( copy.sigma2 ), sigma3( copy.sigma3 )
	, timeRes2Frac( copy.timeRes2Frac), timeRes3Frac( copy.timeRes3Frac )
	, resolutionScale1( copy.resolutionScale1 ), resolutionScale2( copy.resolutionScale2 ), resolutionScale3( copy.resolutionScale3 )
	, tlow( copy.tlow ), thigh( copy.thigh ), time( copy.time )
	, _useEventResolution( copy._useEventResolution )
	, _useTimeAcceptance( copy._useTimeAcceptance )
	, _numericIntegralForce( copy._numericIntegralForce )
	, _usePunziSigmat( copy._usePunziSigmat )
{
	timeAcc = new SlicedAcceptance( *(copy.timeAcc) );
	for( unsigned int i=0; i< copy._intexpIntObs_vec.size(); ++i )
	{
		_intexpIntObs_vec.push_back( PseudoObservable(copy._intexpIntObs_vec[i]) );
	}
	//cout << "making copy " << tau << endl;
}


//Make the data point and parameter set
void Exponential::MakePrototypes()
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
	parameterNames.push_back( tauName );
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
	allParameters = ParameterSet(parameterNames);
}

//Return a list of observables not to be integrated
vector<string> Exponential::GetDoNotIntegrateList()
{
	vector<string> list;
	if( useEventResolution() && !_usePunziSigmat ) list.push_back(eventResolutionName);
	return list;
}

//Destructor
Exponential::~Exponential()
{
}

bool Exponential::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	timeRes2Frac = allParameters.GetPhysicsParameter( timeRes2FracName )->GetValue();
	timeRes3Frac = allParameters.GetPhysicsParameter( timeRes3FracName )->GetValue();
	tau = allParameters.GetPhysicsParameter( tauName )->GetValue();
	gamma = 1./tau;
	if( ! useEventResolution() ) {
		sigma1    = allParameters.GetPhysicsParameter( sigma1Name )->GetValue();
		sigma2    = allParameters.GetPhysicsParameter( sigma2Name )->GetValue();
		sigma3    = allParameters.GetPhysicsParameter( sigma3Name )->GetValue();
	}
	resolutionScale1 = allParameters.GetPhysicsParameter( resScale1Name )->GetValue();
	resolutionScale2 = allParameters.GetPhysicsParameter( resScale2Name )->GetValue();
	resolutionScale3 = allParameters.GetPhysicsParameter( resScale3Name )->GetValue();

	return isOK;
}

//Calculate the function value
double Exponential::Evaluate(DataPoint * measurement)
{
	_dataPoint=measurement;
	// Observable
	Observable* timeObs = measurement->GetObservable( timeName );
	time = timeObs->GetValue();
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
	if( useTimeAcceptance() ) num = num * timeAcc->getValue(timeObs);
	//cout << eventResolution << endl;
	return num;
}

double Exponential::buildPDFnumerator()
{
	// Sum of two exponentials, using the time resolution functions

	if( tau <= 0 ) {
		PDF_THREAD_LOCK
		cout << " In Exponential() you gave a negative or zero lifetime for tau " << endl ;
		PDF_THREAD_UNLOCK
		throw(10) ;
	}
	//double val = Mathematics::Exp(time, 1./tau, sigma);
	double val = Mathematics::Exp(time, gamma, sigma);
	return val;
}


double Exponential::Normalisation( PhaseSpaceBoundary* boundary )
{
	(void) boundary;
	return -1.;
	//if( _dataPoint == NULL ) return -1.;
	/*
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
			PDF_THREAD_LOCK
			cerr << "Bound on time not provided" << endl;
			PDF_THREAD_UNLOCK
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
	return norm;*/
}

double Exponential::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	_dataPoint = measurement;
	if( _numericIntegralForce ) return -1.;

	double norm = 0.;

	IConstraint * timeBound = boundary->GetConstraint( timeName );
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		PDF_THREAD_LOCK
		cerr << "Bound on time not provided" << endl;
		PDF_THREAD_UNLOCK
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
			sigmaNum=0;
			norm = buildPDFdenominator();
		}
		else
		{
			// Set the member variable for time resolution to the first value and calculate
			sigma = eventResolution * resolutionScale1;
			sigmaNum=0;
			double val1 = buildPDFdenominator();
			// Set the member variable for time resolution to the second value and calculate
			sigma = eventResolution * resolutionScale2;
			sigmaNum=1;
			double val2 = buildPDFdenominator();
			// Set the member variable for time resolution to the second value and calculate
			sigma = eventResolution * resolutionScale3;
			sigmaNum=2;
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
			sigmaNum=0;
			double val1 = buildPDFdenominator();
			// Set the member variable for time resolution to the second value and calculate
			sigma = sigma2;
			sigmaNum=1;
			double val2 = buildPDFdenominator();
			sigma = sigma3;
			sigmaNum=2;
			double val3 = buildPDFdenominator();
			norm = (1. - timeRes2Frac - timeRes3Frac)*val1 + timeRes2Frac*val2 + timeRes3Frac*val3;
		}
	}
	return norm;
}

double Exponential::buildPDFdenominator()
{
	// Sum of two exponentials, using the time resolution functions

	if( tau <= 0 )
	{
		PDF_THREAD_LOCK
		cout << " In Exponential() you gave a negative or zero lifetime for tau " << endl ;
		PDF_THREAD_UNLOCK
		throw(10);
	}

	double tlo_boundary = tlow;
	double thi_boundary = thigh;
	double val = 0;

	for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
	{
		AcceptanceSlice* thisSlice = timeAcc->getSlice(islice);
		tlow  = tlo_boundary > thisSlice->tlow() ? tlo_boundary : thisSlice->tlow();
		thigh = thi_boundary < thisSlice->thigh() ? thi_boundary : thisSlice->thigh();
		//if( thigh > tlow ) val += Mathematics::ExpInt(tlow, thigh, 1./tau, sigma) * timeAcc->getSlice(islice)->height();
		vector<double> input(4, 0.); input[0]=tlow; input[1]=thigh; input[2]=gamma; input[3]=sigma;
		if( thigh > tlow )
		{
			//val += Mathematics::ExpInt(tlow, thigh, gamma, sigma) * timeAcc->getSlice(islice)->height();
			unsigned int index = islice + timeAcc->numberOfSlices()*(unsigned)sigmaNum;
			val += _dataPoint->GetPseudoObservable( _intexpIntObs_vec[ index ], input ) * thisSlice->height();
		}
	}

	tlow  = tlo_boundary;
	thigh = thi_boundary;
	return val;
	//return Mathematics::ExpInt(tlow, thigh, 1./tau, sigma);
}

void Exponential::prepareTimeInt()
{
	for( unsigned int j=0; j< 3; ++j )
	{
		for( unsigned int i=0; i< timeAcc->numberOfSlices(); ++i )
		{
			TString thisPseudoName="time_expInt_";
			thisPseudoName+=i+j*timeAcc->numberOfSlices();
			_intexpIntObs_vec.push_back( PseudoObservable( thisPseudoName.Data() ) );
			_intexpIntObs_vec.back().AddFunction( Mathematics::ExpInt_Wrapper );
		}
	}
}

