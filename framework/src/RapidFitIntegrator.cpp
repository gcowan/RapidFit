/**
  @class RapidFitIntegrator

  A numerical integrator to be used for projecting the PDF, and as a fallback if the PDF does not provide its own normalisation
  This class uses two integrator classes provided by root: AdaptiveIntegratorMultiDim and GaussLegendreIntegrator.
  Both of these assume the function to be integrated is well behaved.

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-8
  */

//	ROOT Headers
#include "RVersion.h"
//	RapidFit Headers
#include "RapidFitIntegrator.h"
#include "StringProcessing.h"
#include "StatisticsFunctions.h"
#include "ClassLookUp.h"
#include "NormalisedSumPDF.h"
#include "DebugClass.h"
#include "Threading.h"
#include "MultiThreadedFunctions.h"
#include "MemoryDataSet.h"
//	System Headers
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <float.h>
#include <pthread.h>
#include <exception>
#include <pthread.h>

#ifdef __RAPIDFIT_USE_GSL
//	GSL for MC integration
#include <gsl/gsl_qrng.h>
//	GSLMCIntegrator from ROOT
#include "Math/GSLMCIntegrator.h"
#endif

pthread_mutex_t GSL_DATAPOINT_MANAGEMENT_THREADLOCK;
pthread_mutex_t GSL_DATAPOINT_GET_THREADLOCK;

pthread_mutex_t check_settings_lock;

pthread_mutex_t gsl_mutex;

pthread_mutex_t one_dim_lock;
pthread_mutex_t multi_dim_lock;

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

using namespace::std;

//1% tolerance on integral value
//const double INTEGRAL_PRECISION_THRESHOLD = 0.01;

//Constructor with correct argument
RapidFitIntegrator::RapidFitIntegrator( IPDF * InputFunction, bool ForceNumerical, bool UsePseudoRandomIntegration ) :
	ratioOfIntegrals(-1.), fastIntegrator(NULL), functionToWrap(InputFunction), multiDimensionIntegrator(NULL), oneDimensionIntegrator(NULL),
	functionCanIntegrate(false), haveTestedIntegral(false), num_threads(4),
	RapidFitIntegratorNumerical( ForceNumerical ), obs_check(false), checked_list(), debug(new DebugClass(false) ),
	pseudoRandomIntegration( UsePseudoRandomIntegration ), GSLFixedPoints( __DEFAULT_RAPIDFIT_FIXEDINTEGRATIONPOINTS )
{
	multiDimensionIntegrator = new AdaptiveIntegratorMultiDim();
#if ROOT_VERSION_CODE > ROOT_VERSION(5,28,0)
	multiDimensionIntegrator->SetAbsTolerance( __DEFAULT_RAPIDFIT_INTABSTOL );      //      Absolute error for things such as plots
	multiDimensionIntegrator->SetRelTolerance( __DEFAULT_RAPIDFIT_INTRELTOL );
	multiDimensionIntegrator->SetMaxPts( __DEFAULT_RAPIDFIT_MAXINTEGRALSTEPS );		//	These are the defaults, and it's unlikely you will be able to realistically push the integral without using "double double"'s
#endif

	ROOT::Math::IntegrationOneDim::Type type = ROOT::Math::IntegrationOneDim::kGAUSS;
	oneDimensionIntegrator = new IntegratorOneDim(type);
}

RapidFitIntegrator::RapidFitIntegrator( const RapidFitIntegrator& input ) : ratioOfIntegrals( input.ratioOfIntegrals ),
	fastIntegrator( NULL ), functionToWrap( input.functionToWrap ), multiDimensionIntegrator( NULL ), oneDimensionIntegrator( NULL ),
	pseudoRandomIntegration(input.pseudoRandomIntegration), functionCanIntegrate( input.functionCanIntegrate ), haveTestedIntegral( true ),
	RapidFitIntegratorNumerical( input.RapidFitIntegratorNumerical ), obs_check( input.obs_check ), checked_list( input.checked_list ),
	debug((input.debug==NULL)?NULL:new DebugClass(*input.debug)), num_threads(input.num_threads), GSLFixedPoints( input.GSLFixedPoints )
{
	//	We don't own the PDF so no need to duplicate it as we have to be told which one to use
	//if( input.functionToWrap != NULL ) functionToWrap = ClassLookUp::CopyPDF( input.functionToWrap );

	//	Only Construct the Integrators if they are required
	if( input.fastIntegrator != NULL ) fastIntegrator = new FoamIntegrator( *(input.fastIntegrator) );
	if( input.multiDimensionIntegrator != NULL )
	{
		multiDimensionIntegrator = new AdaptiveIntegratorMultiDim();// *(input.multiDimensionIntegrator) );//new AdaptiveIntegratorMultiDim();
		//      These functions only exist with ROOT > 5.27 I think, at least they exist in 5.28/29
#if ROOT_VERSION_CODE > ROOT_VERSION(5,28,0)
		multiDimensionIntegrator->SetAbsTolerance( __DEFAULT_RAPIDFIT_INTABSTOL );      //      Absolute error for things such as plots
		multiDimensionIntegrator->SetRelTolerance( __DEFAULT_RAPIDFIT_INTRELTOL );
		multiDimensionIntegrator->SetMaxPts( __DEFAULT_RAPIDFIT_MAXINTEGRALSTEPS );
#endif
	}
	if( input.oneDimensionIntegrator != NULL )
	{
		//ROOT::Math::IntegrationOneDim::Type type = ROOT::Math::IntegrationOneDim::kGAUSS;
		//(void) type;	//	For the global state of ROOT
		//oneDimensionIntegrator = new IntegratorOneDim( type );//*(input.oneDimensionIntegrator) );
		oneDimensionIntegrator = new IntegratorOneDim( *(input.oneDimensionIntegrator) );
	}
	else
	{
		ROOT::Math::IntegrationOneDim::Type type = ROOT::Math::IntegrationOneDim::kGAUSS;
		oneDimensionIntegrator = new IntegratorOneDim(type);
	}
}

void RapidFitIntegrator::SetUpIntegrator( const RapidFitIntegratorConfig* config )
{
	this->SetNumThreads( config->numThreads );
	this->SetUseGSLIntegrator( config->useGSLIntegrator );
	this->SetFixedIntegralPoints( config->FixedIntegrationPoints );
	this->SetMaxIntegrationSteps( config->MaxIntegrationSteps );
	this->SetIntegrationAbsTolerance( config->IntegrationRelTolerance );
	this->SetIntegrationRelTolerance( config->IntegrationAbsTolerance );
}

void RapidFitIntegrator::SetMaxIntegrationSteps( const unsigned int input )
{
#if ROOT_VERSION_CODE > ROOT_VERSION(5,28,0)
	multiDimensionIntegrator->SetMaxPts( input );
#endif
}

void RapidFitIntegrator::SetIntegrationAbsTolerance( const double input )
{
#if ROOT_VERSION_CODE > ROOT_VERSION(5,28,0)
	multiDimensionIntegrator->SetAbsTolerance( input );
#endif
}

void RapidFitIntegrator::SetIntegrationRelTolerance( const double input )
{
#if ROOT_VERSION_CODE > ROOT_VERSION(5,28,0)
	multiDimensionIntegrator->SetRelTolerance( input );
#endif
}

void RapidFitIntegrator::SetNumThreads( const unsigned int input )
{
	num_threads = input;
}

bool RapidFitIntegrator::GetUseGSLIntegrator() const
{
	return pseudoRandomIntegration;
}

void RapidFitIntegrator::SetUseGSLIntegrator( const bool input )
{
	pseudoRandomIntegration = input;
	if( debug != NULL )
	{
		if( debug->DebugThisClass( "RapidFitIntegrator" ) )
		{
			if( input ) cout << "Requesting GSL." << endl;
		}
	}
}

void RapidFitIntegrator::SetFixedIntegralPoints( const unsigned int input )
{
	GSLFixedPoints = input;

	if( debug != NULL )
	{
		if( debug->DebugThisClass( "RapidFitIntegrator" ) )
		{
			cout << "Setting the Fixed number of Integration Points" << endl;
		}
	}
}

//	Don't want the projections to be insanely accurate
void RapidFitIntegrator::ProjectionSettings()
{
	//	These functions only exist with ROOT > 5.27 I think, at least they exist in 5.28/29
#if ROOT_VERSION_CODE > ROOT_VERSION(5,28,0)
	multiDimensionIntegrator->SetAbsTolerance( 1E-4 );	//	Absolute error for things such as plots
	multiDimensionIntegrator->SetRelTolerance( 1E-4 );
	multiDimensionIntegrator->SetMaxPts( 10000 );
#endif
}

RapidFitIntegrator::~RapidFitIntegrator()
{
	if( multiDimensionIntegrator != NULL ) delete multiDimensionIntegrator;
	if( oneDimensionIntegrator != NULL ) delete oneDimensionIntegrator;
	if( fastIntegrator != NULL ) delete fastIntegrator;
	if( debug != NULL ) delete debug;
}

//Return the integral over all observables
double RapidFitIntegrator::Integral( DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary )
{
	NewDataPoint->SetPhaseSpaceBoundary( NewBoundary );
	if( functionToWrap == NULL )
	{
		cerr << "WHAT ARE YOU DOING TO ME, YOUR TRYING TO INTEGRATE OVER A NULL OBJECT!!!" << endl;
		exit(-69);
	}

	double PDF_test_result = 0.;

	//cout << "R:here1" << endl;

	try
	{
		//	cout << "R:here1b" << endl;
		//	cout << functionToWrap << endl;
		PDF_test_result = functionToWrap->Integral( NewDataPoint, NewBoundary );
		//	cout << "R:here1b2" << endl;
	}
	catch(...)
	{
		cerr << "RapidFitIntegrator: PDF " << functionToWrap->GetLabel() << " has failed the basic test of not crashing, check your PhaseSpace and/or return -1 for Numerical Integration" << endl;
		throw(1280);
	}

	//cout << "PDF: " << PDF_test_result << endl;

	functionCanIntegrate = !functionToWrap->GetNumericalNormalisation();

	if( functionCanIntegrate == false )
	{
		//	PDF doesn't know how to Normalise no need to check PDF!
		haveTestedIntegral = true;
	}

	//cout << "R:here" << endl;

	bool cacheEnabled = functionToWrap->GetCachingEnabled();
	//cout << "integral " << TestIntegral( NewDataPoint, NewBoundary ) << endl;

	//Test the integration method the function has provided
	if( haveTestedIntegral )
	{
		//If the function has been tested already, use the result
		if( functionCanIntegrate && !RapidFitIntegratorNumerical )
		{
			return PDF_test_result;
		}
		else
		{
			double return_value = -5.;
			if( functionToWrap->CacheValid( NewDataPoint, NewBoundary ) )
			{
				return_value = PDF_test_result;
			}
			else
			{
				//Make a list of observables not to integrate
				vector<string> dontIntegrate = this->DontNumericallyIntegrateList( NewDataPoint );

				return_value = this->NumericallyIntegrateDataPoint( NewDataPoint, NewBoundary, dontIntegrate );

				if( cacheEnabled )
				{
					functionToWrap->SetCache( return_value, NewDataPoint, NewBoundary );
				}
			}
			return return_value;
		}
	}
	else
	{
		//cout << "R:here2" << endl;
		double returnVal = TestIntegral( NewDataPoint, NewBoundary );
		return returnVal;
	}
	return -1;
}

vector<string> RapidFitIntegrator::DontNumericallyIntegrateList( const DataPoint* NewDataPoint, vector<string> input )
{
	//Make a list of observables not to integrate
	vector<string> dontIntegrate;

	//	If the Observables haven't been checked, check them.
	if( !obs_check )
	{
		vector<string> PDF_params = functionToWrap->GetPrototypeDataPoint();
		vector<string> DataPoint_params;
		if( NewDataPoint != NULL ) DataPoint_params = NewDataPoint->GetAllNames();
		vector<string>::iterator data_i = DataPoint_params.begin();

		vector<string> not_found_list;
		vector<string> not_integrable = functionToWrap->GetDoNotIntegrateList();
		for( ; data_i != DataPoint_params.end(); ++data_i )
		{
			if( StringProcessing::VectorContains( &PDF_params, &(*data_i) ) == -1 )
			{
				not_found_list.push_back( *data_i );
			}
		}
		checked_list = StringProcessing::CombineUniques( not_found_list, not_integrable );
		obs_check = true;
	}

	//	Observables have been checked, set the dontIntegrate list to be the same as the checked list
	dontIntegrate = StringProcessing::CombineUniques( checked_list, input );
	return dontIntegrate;
}

double RapidFitIntegrator::TestIntegral( DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary )
{
	NewDataPoint->SetPhaseSpaceBoundary( NewBoundary );
	//Make a list of observables not to integrate
	vector<string> dontIntegrate = functionToWrap->GetDoNotIntegrateList();

	double testIntegral = 0.;
	try
	{
		testIntegral = functionToWrap->Integral( NewDataPoint, NewBoundary );
	}
	catch(...)
	{
		cerr << "RapidFitIntegrator: PDF " << functionToWrap->GetLabel() << " Failed to Integrate Correctly" << endl;
		throw(-762);
	}
	double numericalIntegral = 0.;

	int NumberCombinations = NewBoundary->GetNumberCombinations();

	if( NumberCombinations == 1 || NumberCombinations == 0 )
	{
		numericalIntegral = this->NumericallyIntegrateDataPoint( NewDataPoint, NewBoundary, dontIntegrate );
	}
	else
	{
		numericalIntegral = this->NumericallyIntegratePhaseSpace( NewBoundary, dontIntegrate );
	}

	//Check if the function has an integrate method
	if( testIntegral > 0.0 )
	{
		//Check if the function's integrate method agrees with the numerical integral
		//Trust the function's integration
		ratioOfIntegrals = numericalIntegral/testIntegral;
		functionCanIntegrate = true;
		functionToWrap->SetCache( testIntegral, NewDataPoint, NewBoundary );

		haveTestedIntegral = true;
		return testIntegral;
	}
	haveTestedIntegral = true;
	return this->NumericallyIntegrateDataPoint( NewDataPoint, NewBoundary, dontIntegrate );
}

double RapidFitIntegrator::GetRatioOfIntegrals()
{
	return ratioOfIntegrals;
}

IPDF * RapidFitIntegrator::GetPDF() const
{
	return functionToWrap;
}

void RapidFitIntegrator::SetPDF( IPDF* input )
{
	functionToWrap = input;
}

double RapidFitIntegrator::OneDimentionIntegral( IPDF* functionToWrap, IntegratorOneDim * oneDimensionIntegrator, const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary,
		ComponentRef* componentIndex, vector<string> doIntegrate, vector<string> dontIntegrate, bool haveTestedIntegral, DebugClass* debug )
{
	(void) haveTestedIntegral;
	IntegratorFunction* quickFunction = new IntegratorFunction( functionToWrap, NewDataPoint, doIntegrate, dontIntegrate, NewBoundary, componentIndex );
	if( debug != NULL )
	{
		quickFunction->SetDebug( debug );
		if( debug->DebugThisClass( "RapidFitIntegrator" ) )
		{
			cout << "RapidFitIntegrator: 1D Setting Up Constraints" << endl;
		}
	}
	//cout << "here" << endl;
	//Find the observable range to integrate over
	IConstraint * newConstraint = NULL;
	try
	{
		newConstraint = NewBoundary->GetConstraint( doIntegrate[0] );
	}
	catch(...)
	{
		cerr << "RapidFitIntegrator: Could NOT find required Constraint " << doIntegrate[0] << " in Data PhaseSpaceBoundary" << endl;
		cerr << "RapidFitIntegrator: Please Fix this by Adding the Constraint to your PhaseSpace for PDF: " << functionToWrap->GetLabel() << endl;
		cerr << endl;
		exit(-8737);
	}
	//IConstraint * newConstraint = NewBoundary->GetConstraint( doIntegrate[0] );
	double minimum = newConstraint->GetMinimum();
	double maximum = newConstraint->GetMaximum();

	//cout << "here2" << endl;
	//cout << "1D Integration" << endl;
	//Do a 1D integration
	oneDimensionIntegrator->SetFunction(*quickFunction);

	//cout << "Performing 1DInt" << endl;
	double output = oneDimensionIntegrator->Integral( minimum, maximum );
	//cout << output << endl;

	//cout << "here3" << endl;
	delete quickFunction;

	return output;
}

double RapidFitIntegrator::PseudoRandomNumberIntegral( IPDF* functionToWrap, const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary,
		ComponentRef* componentIndex, vector<string> doIntegrate, vector<string> dontIntegrate, unsigned int GSLFixedPoints )
{
#ifdef __RAPIDFIT_USE_GSL

	//Make arrays of the observable ranges to integrate over
	double* minima = new double[ doIntegrate.size() ];
	double* maxima = new double[ doIntegrate.size() ];

	for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
	{
		IConstraint * newConstraint = NULL;
		try
		{
			newConstraint = NewBoundary->GetConstraint( doIntegrate[observableIndex] );
		}
		catch(...)
		{
			cerr << "RapidFitIntegrator: Could NOT find required Constraint " << doIntegrate[observableIndex] << " in Data PhaseSpaceBoundary" << endl;
			cerr << "RapidFitIntegrator: Please Fix this by Adding the Constraint to your PhaseSpace for PDF: " << functionToWrap->GetLabel() << endl;
			cerr << endl;
			exit(-8737);
		}
		minima[observableIndex] = (double)newConstraint->GetMinimum();
		maxima[observableIndex] = (double)newConstraint->GetMaximum();
	}

	unsigned int npoint = GSLFixedPoints;
	vector<double> * integrationPoints = new std::vector<double>[doIntegrate.size()];

	//gsl_qrng * q = gsl_qrng_alloc( gsl_qrng_sobol, (unsigned)(doIntegrate.size()) );
	gsl_qrng * q = gsl_qrng_alloc( gsl_qrng_niederreiter_2, (unsigned)(doIntegrate.size()) );
	for (unsigned int i = 0; i < npoint; i++)
	{
		double v[doIntegrate.size()];
		gsl_qrng_get(q, v);
		for ( unsigned int j = 0; j < doIntegrate.size(); j++)
		{
			integrationPoints[j].push_back(v[j]);
		}
	}
	gsl_qrng_free(q);

	vector<double> minima_v, maxima_v;
	for( unsigned int i=0; i< doIntegrate.size(); ++i )
	{
		minima_v.push_back( minima[i] );
		maxima_v.push_back( maxima[i] );
	}
	IntegratorFunction* quickFunction = new IntegratorFunction( functionToWrap, NewDataPoint, doIntegrate, dontIntegrate, NewBoundary, componentIndex, minima_v, maxima_v );

	//cout << endl;
	double result = 0.;
	for (unsigned int i = 0; i < integrationPoints[0].size(); ++i)
	{
		double* point = new double[ doIntegrate.size() ];
		for ( unsigned int j = 0; j < doIntegrate.size(); j++)
		{
			//cout << doIntegrate[j] << " " << maxima[j] << " " << minima[j] << " " << integrationPoints[j][i] << endl;
			point[j] = integrationPoints[j][i]*(maxima[j]-minima[j])+minima[j];
		}
		//cout << i << "\r" << flush;
		result += quickFunction->DoEval( point );
		delete[] point;
	}
	double factor=1.;
	for( unsigned int i=0; i< (unsigned)doIntegrate.size(); ++i )
	{
		factor *= maxima[i]-minima[i];
	}
	result /= ( (double)integrationPoints[0].size() / factor );

	//cout << "GSL :D" << endl;

	delete[] minima; delete[] maxima;
	//delete[] integrationPoints;
	delete quickFunction;
	return result;
#else
	(void) functionToWrap; (void) NewDataPoint; (void) NewBoundary; (void) componentIndex; (void) doIntegrate; (void) dontIntegrate;
	return -1.;
#endif
}

vector<DataPoint*> RapidFitIntegrator::initGSLDataPoints( unsigned int number, vector<double> maxima, vector<double> minima, DataPoint* templateDataPoint, vector<string> doIntegrate )
{
	vector<DataPoint*> doEval_points;
	pthread_mutex_lock( &GSL_DATAPOINT_MANAGEMENT_THREADLOCK );

#ifdef __RAPIDFIT_USE_GSL

	unsigned int nDim = (unsigned) minima.size();

	unsigned int npoint = number;
	vector<double> * integrationPoints = new vector<double>[ nDim ];

	//pthread_mutex_lock( &gsl_mutex );
	cout << "Allocating GSL Integration Tool. nDim " << nDim << endl;
	gsl_qrng * q = NULL;
	try
	{
		//q = gsl_qrng_alloc( gsl_qrng_sobol, (unsigned)nDim );
		q = gsl_qrng_alloc( gsl_qrng_niederreiter_2, (unsigned)nDim );
	}
	catch(...)
	{
		cout << "Can't Allocate Integration Tool for GSL Integral." << endl;
		cout << " Dim: " << nDim << endl;
		exit(-742);
	}

	if( q == NULL )
	{
		cout << "Can't Allocate Integration Tool for GSL Integral." << endl;
		cout << " Dim: " << nDim << endl;
		exit(-741);
	}

	double* v = new double[ nDim ];
	for (unsigned int i = 0; i < npoint; i++)
	{
		gsl_qrng_get( q, v );
		for( unsigned int j = 0; j < nDim; j++)
		{
			integrationPoints[j].push_back( v[j] );
		}
	}
	delete v;
	gsl_qrng_free(q);
	//cout << "Freed GSL Integration Tool" << endl;
	//pthread_mutex_unlock( &gsl_mutex );

	/*vector<double> minima_v, maxima_v;
	  for( unsigned int i=0; i< nDim; ++i )
	  {
	  minima_v.push_back( minima[i] );
	  maxima_v.push_back( maxima[i] );
	  }*/
	//cout << "Constructing Functions" << endl;
	//IntegratorFunction* quickFunction = new IntegratorFunction( functionToWrap, NewDataPoint, doIntegrate, dontIntegrate, NewBoundary, componentIndex, minima_v, maxima_v );

	for (unsigned int i = 0; i < integrationPoints[0].size(); ++i)
	{
		double* point = new double[ nDim ];
		for ( unsigned int j = 0; j < nDim; j++)
		{
			//cout << doIntegrate[j] << " " << maxima[j] << " " << minima[j] << " " << integrationPoints[j][i] << endl;
			point[j] = integrationPoints[j][i]*(maxima[j]-minima[j])+minima[j];
		}

		DataPoint* thisPoint = new DataPoint( *templateDataPoint );

		for( unsigned int k=0; k< doIntegrate.size(); ++k )
		{
			thisPoint->SetObservable( doIntegrate[k], point[k], "noUnitsHere" );
		}

		doEval_points.push_back( thisPoint );
		//result += quickFunction->DoEval( point );
		//delete[] point;
	}

	delete[] integrationPoints;
#endif

	(void) maxima; (void) minima; (void) number;
	pthread_mutex_unlock( &GSL_DATAPOINT_MANAGEMENT_THREADLOCK );
	return doEval_points;
}

vector<DataPoint*> RapidFitIntegrator::getGSLIntegrationPoints( unsigned int number, vector<double> maxima, vector<double> minima, DataPoint* templateDataPoint, vector<string> doIntegrate )
{
	pthread_mutex_lock( &GSL_DATAPOINT_GET_THREADLOCK );
	if( ( number != _global_doEval_points.size() ) || ( ( _global_range_minima != minima ) || ( _global_range_maxima != maxima ) ) )
	{
		clearGSLIntegrationPoints();
		_global_doEval_points = initGSLDataPoints( number, maxima, minima, templateDataPoint, doIntegrate );
		_global_range_minima = minima;
		_global_range_maxima = maxima;
		_global_observable_names = doIntegrate;
	}
	pthread_mutex_unlock( &GSL_DATAPOINT_GET_THREADLOCK );
	return _global_doEval_points;
}

void RapidFitIntegrator::clearGSLIntegrationPoints()
{
	pthread_mutex_lock( &GSL_DATAPOINT_MANAGEMENT_THREADLOCK );
	while( !_global_doEval_points.empty() )
	{
		if( _global_doEval_points.back() != NULL ) delete _global_doEval_points.back();
	}
	_global_doEval_points.pop_back();
	_global_range_minima.clear();
	_global_range_maxima.clear();
	_global_observable_names.clear();
	pthread_mutex_unlock( &GSL_DATAPOINT_MANAGEMENT_THREADLOCK );
}

double RapidFitIntegrator::PseudoRandomNumberIntegralThreaded( IPDF* functionToWrap, const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary,
		ComponentRef* componentIndex, vector<string> doIntegrate, vector<string> dontIntegrate, unsigned int num_threads, unsigned int GSLFixedPoints, DebugClass* debug )
{
#ifdef __RAPIDFIT_USE_GSL

	(void) dontIntegrate;
	/*
	   cout << endl << "Using: " << num_threads << endl;
	   cout << "Do:" << endl;
	   for( unsigned int i=0; i< doIntegrate.size();++i )
	   {
	   cout << doIntegrate[i] << endl;
	   }
	   cout << "Dont:" << endl;
	   for( unsigned int i=0; i< dontIntegrate.size();++i )
	   {
	   cout << dontIntegrate[i] << endl;
	   }
	   */

	//Make arrays of the observable ranges to integrate over
	double* minima = new double[ doIntegrate.size() ];
	double* maxima = new double[ doIntegrate.size() ];

	if( debug != NULL )
	{
		if( debug->DebugThisClass( "RapidFitIntegrator" ) )
		{
			cout << "RapidFitIntegrator: Starting to use GSL PseudoRandomNumberThreaded :D" << endl;
		}
	}

	for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
	{
		IConstraint * newConstraint = NULL;
		try
		{
			newConstraint = NewBoundary->GetConstraint( doIntegrate[observableIndex] );
		}
		catch(...)
		{
			cerr << "RapidFitIntegrator: Could NOT find required Constraint " << doIntegrate[observableIndex] << " in Data PhaseSpaceBoundary" << endl;
			cerr << "RapidFitIntegrator: Please Fix this by Adding the Constraint to your PhaseSpace for PDF: " << functionToWrap->GetLabel() << endl;
			cerr << endl;
			exit(-8737);
		}
		minima[observableIndex] = (double)newConstraint->GetMinimum();
		maxima[observableIndex] = (double)newConstraint->GetMaximum();
	}

	if( debug != NULL )
	{
		if( debug->DebugThisClass( "RapidFitIntegrator" ) )
		{
			cout << "RapidFitIntegrator: Doing GSL stuff..." << endl;
		}
	}

	vector<double> minima_v, maxima_v;
	for( unsigned int i=0; i< doIntegrate.size(); ++i )
	{
		minima_v.push_back( minima[i] );
		maxima_v.push_back( maxima[i] );
	}

	//	the whole new/delete here is because I'm too lazy to lookup how to cast away the const from a pointer...

	DataPoint* templateDataPoint = new DataPoint( *NewDataPoint );

	vector<DataPoint*> doEval_points = RapidFitIntegrator::getGSLIntegrationPoints( GSLFixedPoints, maxima_v, minima_v, templateDataPoint, doIntegrate );

	delete templateDataPoint;

	PhaseSpaceBoundary* thisBound = new PhaseSpaceBoundary( *NewBoundary );

	IDataSet* thisDataSet = new MemoryDataSet( thisBound, doEval_points );

	delete thisBound;

	ThreadingConfig* thisConfig = new ThreadingConfig();
	thisConfig->numThreads = num_threads;
	thisConfig->MultiThreadingInstance = "pthreads";
	thisConfig->wantedComponent = new ComponentRef( *componentIndex );

	vector<double>* thisSet = MultiThreadedFunctions::ParallelEvaulate( functionToWrap, thisDataSet, thisConfig );

	delete[] minima;
	delete[] maxima;
	delete thisConfig->wantedComponent;
	delete thisConfig;

	double result=0.;
	for( unsigned int i=0; i< thisSet->size(); ++i )	result+=thisSet->at( i );

	delete thisSet;

	double factor=1.;
	for( unsigned int i=0; i< (unsigned)doIntegrate.size(); ++i )
	{
		factor *= maxima[i]-minima[i];
	}

	result /= ( (double)doEval_points.size() / factor );

	return result;
#else
	(void) functionToWrap; (void) NewDataPoint; (void) NewBoundary; (void) componentIndex; (void) doIntegrate; (void) dontIntegrate; (void) num_threads;
	return -99999.;
#endif
}

double RapidFitIntegrator::MultiDimentionIntegral( IPDF* functionToWrap, AdaptiveIntegratorMultiDim* multiDimensionIntegrator, const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary,
		ComponentRef* componentIndex, vector<string> doIntegrate, vector<string> dontIntegrate, bool haveTestedIntegral, DebugClass* debug )
{
	(void) haveTestedIntegral;
	//Make arrays of the observable ranges to integrate over
	double* minima = new double[ doIntegrate.size() ];
	double* maxima = new double[ doIntegrate.size() ];

	for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
	{
		IConstraint * newConstraint = NULL;
		try
		{
			newConstraint = NewBoundary->GetConstraint( doIntegrate[observableIndex] );
		}
		catch(...)
		{
			cerr << "RapidFitIntegrator: Could NOT find required Constraint " << doIntegrate[observableIndex] << " in Data PhaseSpaceBoundary" << endl;
			cerr << "RapidFitIntegrator: Please Fix this by Adding the Constraint to your PhaseSpace for PDF: " << functionToWrap->GetLabel() << endl;
			cerr << endl;
			exit(-8737);
		}
		minima[observableIndex] = (double)newConstraint->GetMinimum();
		maxima[observableIndex] = (double)newConstraint->GetMaximum();
	}

	//Do a 2-15D integration

	vector<double> minima_v, maxima_v;
	for( unsigned int i=0; i< doIntegrate.size(); ++i )
	{
		minima_v.push_back( minima[i] );
		maxima_v.push_back( maxima[i] );
	}

	IntegratorFunction* quickFunction = new IntegratorFunction( functionToWrap, NewDataPoint, doIntegrate, dontIntegrate, NewBoundary, componentIndex, minima_v, maxima_v );
	if( debug != NULL ) quickFunction->SetDebug( debug );

	multiDimensionIntegrator->SetFunction( *quickFunction );

	double output =  multiDimensionIntegrator->Integral( minima, maxima );

	delete[] minima; delete[] maxima;
	delete quickFunction;

	return output;
}

//	Integrate over all possible discrete combinations within the phase-space
double RapidFitIntegrator::NumericallyIntegratePhaseSpace( PhaseSpaceBoundary* NewBoundary, vector<string> DontIntegrateThese, ComponentRef* componentIndex )
{
	return this->DoNumericalIntegral( NULL, NewBoundary, DontIntegrateThese, componentIndex, false );
}

//	Integrate over the given phase space using this given DataPoint. i.e. ONLY 1 UNIQUE discrete combination
double RapidFitIntegrator::NumericallyIntegrateDataPoint( DataPoint* NewDataPoint, PhaseSpaceBoundary* NewBoundary, vector<string> DontIntegrateThese, ComponentRef* componentIndex )
{
	return this->DoNumericalIntegral( NewDataPoint, NewBoundary, DontIntegrateThese, componentIndex, true );
}

//Actually perform the numerical integration
double RapidFitIntegrator::DoNumericalIntegral( const DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary, const vector<string> DontIntegrateThese, ComponentRef* componentIndex, const bool IntegrateDataPoint )
{
	//Make lists of observables to integrate and not to integrate
	vector<string> doIntegrate, dontIntegrate;
	StatisticsFunctions::DoDontIntegrateLists( functionToWrap, NewBoundary, &DontIntegrateThese, doIntegrate, dontIntegrate );

	dontIntegrate = StringProcessing::CombineUniques( dontIntegrate, DontIntegrateThese );
	dontIntegrate = this->DontNumericallyIntegrateList( NewDataPoint, dontIntegrate );

	vector<string> CANNOT_INTEGRATE_LIST = NewBoundary->GetDiscreteNames();
	dontIntegrate = StringProcessing::CombineUniques( dontIntegrate, CANNOT_INTEGRATE_LIST );

	vector<string> safeDoIntegrate;
	for( unsigned int i=0; i< doIntegrate.size(); ++i )
	{
		if( StringProcessing::VectorContains( &dontIntegrate, &(doIntegrate[i]) ) == -1 )
		{
			safeDoIntegrate.push_back( doIntegrate[i] );
		}
	}
	doIntegrate = safeDoIntegrate;

	vector<DataPoint*> DiscreteIntegrals;

	vector<string> required = functionToWrap->GetPrototypeDataPoint();
	bool isFixed=true;
	for( unsigned int i=0; i< required.size(); ++i )
	{
		isFixed = isFixed && NewBoundary->GetConstraint( required[i] )->IsDiscrete();
	}
	if( isFixed )
	{
		if( NewDataPoint != NULL )
		{
			DataPoint* thisDataPoint = new DataPoint( *NewDataPoint );
			thisDataPoint->SetPhaseSpaceBoundary( NewBoundary );
			double returnVal = functionToWrap->Evaluate( thisDataPoint );
			delete thisDataPoint;
			return returnVal;
		}
		else
		{
			DiscreteIntegrals = NewBoundary->GetDiscreteCombinations();
			double returnVal=0.;
			for( unsigned int i=0; i< DiscreteIntegrals.size(); ++i )
			{
				returnVal+=functionToWrap->Evaluate( DiscreteIntegrals[i] );
			}
			while( !DiscreteIntegrals.empty() )
			{
				if( DiscreteIntegrals.back() != NULL ) delete DiscreteIntegrals.back();
				DiscreteIntegrals.pop_back();
			}
			return returnVal;
		}
	}

	if( IntegrateDataPoint )
	{
		DiscreteIntegrals.push_back( new DataPoint(*NewDataPoint) );
		DiscreteIntegrals.back()->SetPhaseSpaceBoundary( NewBoundary );
	}
	else
	{
		DiscreteIntegrals = NewBoundary->GetDiscreteCombinations();
	}

	double output_val = 0.;

	//If there are no observables left to integrate over, just evaluate the function
	if( doIntegrate.empty() || doIntegrate.size() == 0 )
	{
		for( vector<DataPoint*>::iterator dataPoint_i = DiscreteIntegrals.begin(); dataPoint_i != DiscreteIntegrals.end(); ++dataPoint_i )
		{
			try{
				output_val += functionToWrap->Integral( *dataPoint_i, NewBoundary );
			}
			catch(...)
			{
				cerr << "Analytical Integral Fell over!" << endl;
				return -9999.;
			}
		}
	}
	else
	{
		for( vector<DataPoint*>::iterator dataPoint_i = DiscreteIntegrals.begin(); dataPoint_i != DiscreteIntegrals.end(); ++dataPoint_i )
		{
			double numericalIntegral = 0.;
			//Chose the one dimensional or multi-dimensional method
			if( doIntegrate.size() == 1 )
			{
				if( debug != NULL )
				{
					if( debug->DebugThisClass( "RapidFitIntegrator" ) )
					{
						cout << "RapidFitIntegrator: One Dimensional Integral" << endl;
					}
				}
				pthread_mutex_lock( &one_dim_lock );
				numericalIntegral += this->OneDimentionIntegral( functionToWrap, oneDimensionIntegrator, *dataPoint_i, NewBoundary, componentIndex, doIntegrate, dontIntegrate, debug );
				//cout << "ret: " << numericalIntegral << endl;
				pthread_mutex_unlock( &one_dim_lock );
			}
			else
			{
				if( debug != NULL )
				{
					if( debug->DebugThisClass( "RapidFitIntegrator" ) )
					{
						cout << "RapidFitIntegrator: Multi Dimensional Integral" << endl;
					}
				}
				if( !pseudoRandomIntegration )
				{
					pthread_mutex_lock( &multi_dim_lock );
					numericalIntegral += this->MultiDimentionIntegral( functionToWrap, multiDimensionIntegrator, *dataPoint_i, NewBoundary, componentIndex, doIntegrate, dontIntegrate, debug );
					pthread_mutex_unlock( &multi_dim_lock );
				}
				else
				{
					if( debug != NULL )
					{
						if( debug->DebugThisClass( "RapidFitIntegrator" ) )
						{
							cout << "RapidFitIntegrator: Using GSL PseudoRandomNumber :D" << endl;
						}
					}
					//numericalIntegral += this->PseudoRandomNumberIntegral( functionToWrap, *dataPoint_i, NewBoundary, componentIndex, doIntegrate, dontIntegrate, GSLFixedPoints );
					numericalIntegral += this->PseudoRandomNumberIntegralThreaded( functionToWrap, *dataPoint_i, NewBoundary, componentIndex, doIntegrate, dontIntegrate, num_threads, GSLFixedPoints, debug );
					if( debug != NULL )
					{
						if( debug->DebugThisClass( "RapidFitIntegrator" ) )
						{
							cout << "RapidFitIntegrator: Finished: " << numericalIntegral << endl;
						}
					}

					if( numericalIntegral <= -99999. )
					{
						pthread_mutex_lock( &check_settings_lock );
						cout << "Calculated a -ve Integral: " << numericalIntegral << ". Did you Compile with the GSL options Enabled with gsl available?" << endl;
						//cout << endl;	exit(-2356);
						cout << "Reverting to non-GSL integration" << endl;

						if( doIntegrate.size() == 1 )
						{
							pthread_mutex_lock( &one_dim_lock );
							numericalIntegral = this->OneDimentionIntegral( functionToWrap, oneDimensionIntegrator, *dataPoint_i, NewBoundary, componentIndex, doIntegrate, dontIntegrate, debug );
							pthread_mutex_unlock( &one_dim_lock );
						}
						else
						{
							pthread_mutex_lock( &multi_dim_lock );
							numericalIntegral = this->MultiDimentionIntegral( functionToWrap, multiDimensionIntegrator, *dataPoint_i, NewBoundary, componentIndex, doIntegrate, dontIntegrate, debug );
							pthread_mutex_unlock( &multi_dim_lock );
						}
						this->SetUseGSLIntegrator( false );
						pthread_mutex_unlock( &check_settings_lock );
					}
				}
			}

			if( !haveTestedIntegral && !functionToWrap->GetNumericalNormalisation() )
			{
				double testIntegral = functionToWrap->Integral( *dataPoint_i, NewBoundary );
				cout << "Integration Test: numerical : analytical  " << setw(7) << numericalIntegral << " : " << testIntegral;
				string description = NewBoundary->DiscreteDescription( *dataPoint_i );
				description = description.substr(0, description.size()-2);
				cout << "  "  << description << "  " << functionToWrap->GetLabel() << endl;
			}

			output_val += numericalIntegral;
		}
		//cout << "output: " << output_val << endl;
	}

	while( !DiscreteIntegrals.empty() )
	{
		if( DiscreteIntegrals.back() != NULL ) delete DiscreteIntegrals.back();
		DiscreteIntegrals.pop_back();
	}

	//cout << "ret2: " << output_val << endl;
	return output_val;
}

//Return the integral over all observables except one
double RapidFitIntegrator::ProjectObservable( DataPoint* NewDataPoint, PhaseSpaceBoundary * NewBoundary, string ProjectThis, ComponentRef* Component )
{
	//Make the list of observables not to integrate
	vector<string> dontIntegrate = functionToWrap->GetDoNotIntegrateList();
	double value = -1.;

	vector<string> allIntegrable = functionToWrap->GetPrototypeDataPoint();

	vector<string> testedIntegrable;

	for( unsigned int i=0; i< allIntegrable.size(); ++i )
	{
		if( StringProcessing::VectorContains( &dontIntegrate, &(allIntegrable[i]) ) == -1 )
		{
			testedIntegrable.push_back( allIntegrable[i] );
		}
	}

	vector<string> CANNOT_INTEGRATE_LIST = NewBoundary->GetDiscreteNames();
	dontIntegrate = StringProcessing::CombineUniques( dontIntegrate, CANNOT_INTEGRATE_LIST );

	dontIntegrate.push_back(ProjectThis);

	vector<string> doIntegrate;

	for( unsigned int i=0; i< allIntegrable.size(); ++i )
	{
		if( StringProcessing::VectorContains( &dontIntegrate, &(allIntegrable[i]) ) == -1 )
		{
			doIntegrate.push_back( allIntegrable[i] );
		}
	}

	if( debug != NULL )
	{
		if( debug->DebugThisClass( "RapidFitIntegrator" ) )
		{
			cout << endl << "Dont Integrate:" << endl;
			for( unsigned int i=0; i< dontIntegrate.size(); ++i )
			{
				cout << dontIntegrate[i] << "  ";
			}
			cout << endl << "Do Integrate:" << endl;
			for( unsigned int i=0; i< allIntegrable.size(); ++i )
			{
				cout << allIntegrable[i] << "  ";
			}
			cout << endl << "Left to Integrate:" << endl;
			for( unsigned int i=0; i< doIntegrate.size(); ++i )
			{
				cout << doIntegrate[i] << " ";
			}
			cout << endl;
		}
	}

	if( testedIntegrable.size() <= 1 || doIntegrate.empty() )
	{
		if( doIntegrate.empty() )
		{
			try{
				value = functionToWrap->EvaluateComponent( NewDataPoint, Component );
			}
			catch(...)
			{
				cerr << "Cannot Evalate PDF AT:" << endl;
				NewDataPoint->Print();
				return 0.;
			}
		}
		else if( testedIntegrable.size() == 1 )
		{
			if( testedIntegrable[0] == ProjectThis )
			{
				try{
					value = functionToWrap->EvaluateComponent( NewDataPoint, Component );
				}
				catch(...)
				{
					cerr << "Cannot Evalate PDF AT:" << endl;
					NewDataPoint->Print();
					return 0.;
				}
			}
			else
			{
				cerr << "This PDF only knows how to Integrate: " << testedIntegrable[0] << " CAN NOT INTEGRATE: " << ProjectThis << endl << endl;
				//throw(-342);
				return 0.;
			}
		}
		else
		{
			cerr << "Unknown Error" << endl;
			return 0.;
		}
	}
	else
	{
		value = this->NumericallyIntegrateDataPoint( NewDataPoint, NewBoundary, dontIntegrate, Component );
	}
	return value;
}

void RapidFitIntegrator::ForceTestStatus( bool input )
{
	haveTestedIntegral = input;
}

void RapidFitIntegrator::SetDebug( DebugClass* input_debug )
{
	if( debug != NULL ) delete debug;
	debug = new DebugClass( *input_debug );
	//if( functionToWrap != NULL ) functionToWrap->SetDebug( input_debug );
}

