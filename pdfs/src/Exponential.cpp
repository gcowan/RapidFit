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
        , timeOffsetName                ( configurator->getName("timeOffset") )
	// Observables
	, timeName      ( configurator->getName("time") )
	//objects used in XML
	, tau(), gamma(), sigmaNum(0), sigma(), sigma1(), sigma2(), sigma3(), timeRes2Frac(), timeRes3Frac()
	, resolutionScale1(), resolutionScale2(), resolutionScale3(), _dataPoint(NULL)
	, tlow(), thigh(), time()
	, _useEventResolution(false), _useTimeAcceptance(false), _numericIntegralForce(false), _usePunziSigmat(false), _useSteppedProjection(false),
	_intexpIntObs_vec(), timeAcc( NULL ), eventResolution(), timeOffset()
{
	_useEventResolution = configurator->isTrue( "UseEventResolution" );
	_useTimeAcceptance  = configurator->isTrue( "UseTimeAcceptance" );
	_numericIntegralForce = configurator->isTrue( "UseNumericalIntegration" );
	_usePunziSigmat = configurator->isTrue( "UsePunziSigmat" );
	_useSteppedProjection = configurator->isTrue( "UseSteppedProjection" );
	if( useTimeAcceptance() )
	{
		if( configurator->hasConfigurationValue( "TimeAcceptanceType", "Upper" ) )
		{
			timeAcc = new SlicedAcceptance( 0., 14.0, 0.0033 );
			cout << "Exponential:: Constructing timeAcc: Upper time acceptance beta=0.0033 [0 < t < 14] " << endl;
		}
		else if( configurator->getConfigurationValue( "TimeAcceptanceFile" ) != "" )
		{
			timeAcc = new SlicedAcceptance( "File" , configurator->getConfigurationValue( "TimeAcceptanceFile" ) );
			cout << "Exponential:: Constructing timeAcc: using file: " << configurator->getConfigurationValue( "TimeAcceptanceFile" ) << endl;
		}
	}
	else
	{
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
	, timeOffsetName         ( copy.timeOffsetName )
	, sigma1Name	( copy.sigma1Name )
	, sigma2Name	( copy.sigma2Name )
	, sigma3Name	( copy.sigma3Name )
	, timeRes2FracName( copy.timeRes2FracName )
	, timeRes3FracName( copy.timeRes3FracName )
	, timeName      ( copy.timeName )
	, tau( copy.tau ), gamma( copy.gamma ), sigmaNum( copy.sigmaNum ), _intexpIntObs_vec(), _dataPoint( NULL )
	, sigma( copy.sigma ), sigma1( copy.sigma1 ), sigma2( copy.sigma2 ), sigma3( copy.sigma3 )
	, timeRes2Frac( copy.timeRes2Frac), timeRes3Frac( copy.timeRes3Frac )
	, resolutionScale1( copy.resolutionScale1 ), resolutionScale2( copy.resolutionScale2 ), resolutionScale3( copy.resolutionScale3 )
	, timeOffset( copy.timeOffset )
	, tlow( copy.tlow ), thigh( copy.thigh ), time( copy.time ), _useEventResolution( copy._useEventResolution )
	, _useTimeAcceptance( copy._useTimeAcceptance ), _numericIntegralForce( copy._numericIntegralForce ), _usePunziSigmat( copy._usePunziSigmat )
	, _useSteppedProjection( copy._useSteppedProjection ), eventResolution( copy.eventResolution ), timeAcc( NULL )
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
	if( useEventResolution() )
	{
		allObservables.push_back( eventResolutionName );
		this->TurnCachingOff();
	}

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( tauName );

	if( ! useEventResolution() )
	{
		parameterNames.push_back( timeRes2FracName );
		parameterNames.push_back( timeRes3FracName );
		parameterNames.push_back( sigma1Name );
		parameterNames.push_back( sigma2Name );
		parameterNames.push_back( sigma3Name );
		parameterNames.push_back( resScale2Name );
		parameterNames.push_back( resScale3Name );
		parameterNames.push_back( timeOffsetName );
	}
	parameterNames.push_back( resScale1Name );
	parameterNames.push_back( timeOffsetName );
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
	if( timeAcc != NULL ) delete timeAcc;
}

bool Exponential::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	timeOffset       = allParameters.GetPhysicsParameter( timeOffsetName )->GetValue();
	tau = allParameters.GetPhysicsParameter( tauName )->GetValue() - timeOffset;
	gamma = 1./tau;
	if( ! useEventResolution() )
	{
		timeRes2Frac = allParameters.GetPhysicsParameter( timeRes2FracName )->GetValue();
		timeRes3Frac = allParameters.GetPhysicsParameter( timeRes3FracName )->GetValue();
		sigma1    = allParameters.GetPhysicsParameter( sigma1Name )->GetValue();
		sigma2    = allParameters.GetPhysicsParameter( sigma2Name )->GetValue();
		sigma3    = allParameters.GetPhysicsParameter( sigma3Name )->GetValue();
		resolutionScale1 = allParameters.GetPhysicsParameter( resScale1Name )->GetValue();
		resolutionScale2 = allParameters.GetPhysicsParameter( resScale2Name )->GetValue();
		resolutionScale3 = allParameters.GetPhysicsParameter( resScale3Name )->GetValue();
	}
	else
	{
		resolutionScale1 = allParameters.GetPhysicsParameter( resScale1Name )->GetValue();
	}     

	return isOK;
}

double Exponential::EvaluateComponent( DataPoint* input, ComponentRef* thisRef )
{
	(void) thisRef;		//	This PDF has no concept of a component, it is JUST an exponential
	Observable* timeObs = input->GetObservable( timeName );
	if( _useSteppedProjection )
	{
		if( timeAcc->GetIsSorted() )
		{
			unsigned int binNum = timeAcc->findSliceNum( timeObs, timeOffset );
			double t_low = timeAcc->getSlice( binNum )->tlow();
			double t_high = timeAcc->getSlice( binNum )->thigh();
			double frac=0.5;
			//frac = ( exp(-gamma*t_low) - exp(-gamma*(t_low+0.5*(t_high-t_low))) ) / ( exp(-gamma*t_low) - exp(-gamma*t_high) );
			frac = ( exp(-gamma*(t_low+0.5*(t_high-t_low))) - exp(-gamma*t_low) ) / ( exp(-gamma*t_high) - exp(-gamma*t_low) );
			double mean = t_low+(1.-frac)*( t_high-t_low );
			DataPoint* thisPoint = new DataPoint( *input );
			Observable* newTimeObs = new Observable( timeName, mean, " " );
			thisPoint->SetObservable( timeName, newTimeObs );
			double returnable = this->Evaluate( thisPoint );
			delete newTimeObs;
			delete thisPoint;
			return returnable;
		}
		else
		{
			cerr << "Acceptance Histogram has to be in an ordered format (histogram-like) in order to do this!" << endl;
			cerr << "Disabling the Stepped Projections" << endl;
			_useSteppedProjection = false;
			return this->Evaluate( input );
		}
	}
	else
	{
		return this->Evaluate( input );
	}
	return 0.;
}

//Calculate the function value
double Exponential::Evaluate(DataPoint * measurement)
{
	_dataPoint=measurement;
	// Observable
	Observable* timeObs = measurement->GetObservable( timeName );
	time = timeObs->GetValue();

	thigh = measurement->GetPhaseSpaceBoundary()->GetConstraint( timeName )->GetMaximum();
	tlow = measurement->GetPhaseSpaceBoundary()->GetConstraint( timeName )->GetMinimum();

	double num = 0.;

	if( useEventResolution() )
	{
		eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();
		sigma = eventResolution * resolutionScale1;
		double numer = buildPDFnumerator();
		sigmaNum=0;
		double denom = buildPDFdenominator();

		num = numer / denom;
	}
	else if( resolutionScale1 <= 0. )
	{
		//This is the "code" to run with resolution=0
		sigma = 0.;
		double numer = buildPDFnumerator();
		sigmaNum=0;
		double denom = buildPDFdenominator();

		num = numer / denom;
	}
	else
	{
		if(  (1. - timeRes2Frac - timeRes3Frac ) >= 0.9999 )
		{
			// Set the member variable for time resolution to the first value and calculate
			sigma = sigma1;
			double numer = buildPDFnumerator();
			sigmaNum = 0;
			double denom = buildPDFdenominator();

			num = numer / denom;
		}
		else
		{
			// Set the member variable for time resolution to the first value and calculate
			sigma = sigma1;
			double val1 = buildPDFnumerator();
			sigmaNum=0;
			double denom1 = buildPDFdenominator();
			// Set the member variable for time resolution to the second value and calculate
			sigma = sigma2;
			double val2 = buildPDFnumerator();
			sigmaNum=1;
			double denom2 = buildPDFdenominator();
			// Set the member variable for time resolution to the second value and calculate
			sigma = sigma3;
			double val3 = buildPDFnumerator();
			sigmaNum=2;	
			double denom3 = buildPDFdenominator();

			num = (1. - timeRes2Frac - timeRes3Frac)*val1/denom1 + timeRes2Frac*val2/denom2 + timeRes3Frac*val3/denom3;
		}
	}

	//if( num >= 1E3 )
	//{
	//	cout << num << "\t\t" << timeAcc->getValue(timeObs) << endl;
	//	exit(0);
	//}

	if( useTimeAcceptance() ) num *= timeAcc->getValue( timeObs, timeOffset );

	//if( useTimeAcceptance() ) num = num * timeAcc->getValue(timeObs);
	//cout << eventResolution << endl;
	return num;
}

double Exponential::buildPDFnumerator()
{
	// Sum of two exponentials, using the time resolution functions
	if( tau <= 0 )
	{
		PDF_THREAD_LOCK
		cout << " In Exponential() you gave a negative or zero lifetime for tau " << endl ;
		PDF_THREAD_UNLOCK
		throw(10) ;
	}

	double val = Mathematics::Exp(time, gamma, sigma);
	return val;
}

double Exponential::Normalisation( DataPoint * measurement, PhaseSpaceBoundary * boundary )
{
	(void) measurement; (void) boundary;
	return 1.;
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

	vector<double> input(4, 0.); input[0]=tlow; input[1]=thigh; input[2]=gamma; input[3]=sigma;

	for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
	{
		AcceptanceSlice* thisSlice = timeAcc->getSlice(islice);
		double this_tlow  = tlo_boundary > thisSlice->tlow() ? tlo_boundary : thisSlice->tlow();
		double this_thigh = thi_boundary < thisSlice->thigh() ? thi_boundary : thisSlice->thigh();
		input[0]=this_tlow; input[1]=this_thigh;
		if( this_thigh > this_tlow )
		{
			//val += Mathematics::ExpInt(this_tlow, this_thigh, 1./tau, sigma) * timeAcc->getSlice(islice)->height();
			unsigned int index = islice + timeAcc->numberOfSlices()*(unsigned)sigmaNum;
			val += _dataPoint->GetPseudoObservable( _intexpIntObs_vec[ index ], input ) * thisSlice->height();
		}
	}

	return val;
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

