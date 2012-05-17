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
//	System Headers
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <float.h>
#include <pthread.h>
#include <exception>

pthread_mutex_t multi_mutex;
pthread_mutex_t multi_mutex2;
pthread_mutex_t multi_mutex3;

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

using namespace::std;

//1% tolerance on integral value
//const double INTEGRAL_PRECISION_THRESHOLD = 0.01;

//Constructor with correct argument
RapidFitIntegrator::RapidFitIntegrator( IPDF * InputFunction, bool ForceNumerical ) :
	ratioOfIntegrals(-1.), fastIntegrator(NULL), functionToWrap(InputFunction), multiDimensionIntegrator(NULL), oneDimensionIntegrator(NULL),
	functionCanIntegrate(false), haveTestedIntegral(false),
	RapidFitIntegratorNumerical( ForceNumerical ), obs_check(false), checked_list()
{
	multiDimensionIntegrator = new AdaptiveIntegratorMultiDim();
#if ROOT_VERSION_CODE > ROOT_VERSION(5,28,0)
	multiDimensionIntegrator->SetAbsTolerance( 1E-9 );      //      Absolute error for things such as plots
	multiDimensionIntegrator->SetRelTolerance( 1E-6 );
	multiDimensionIntegrator->SetMaxPts( 1000000 );		//	These are the defaults, and it's unlikely you will be able to realistically push the integral without using "double double"'s
#endif

	ROOT::Math::IntegrationOneDim::Type type = ROOT::Math::IntegrationOneDim::kGAUSS;
	oneDimensionIntegrator = new IntegratorOneDim(type);
}

RapidFitIntegrator::RapidFitIntegrator( const RapidFitIntegrator& input ) : ratioOfIntegrals( input.ratioOfIntegrals ),
	fastIntegrator( NULL ), functionToWrap( input.functionToWrap ), multiDimensionIntegrator( NULL ), oneDimensionIntegrator( NULL ),
	functionCanIntegrate( input.functionCanIntegrate ), haveTestedIntegral( input.haveTestedIntegral ),
	RapidFitIntegratorNumerical( input.RapidFitIntegratorNumerical ), obs_check( input.obs_check ), checked_list( input.checked_list )
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
		multiDimensionIntegrator->SetAbsTolerance( 1E-9 );      //      Absolute error for things such as plots
		multiDimensionIntegrator->SetRelTolerance( 1E-6 );
		multiDimensionIntegrator->SetMaxPts( 1000000 );
#endif
	}
	if( input.oneDimensionIntegrator != NULL )
	{
		ROOT::Math::IntegrationOneDim::Type type = ROOT::Math::IntegrationOneDim::kGAUSS;
		(void) type;	//	For the global state of ROOT
		oneDimensionIntegrator = new IntegratorOneDim( );//*(input.oneDimensionIntegrator) );
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
}

//Return the integral over all observables
double RapidFitIntegrator::Integral( DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary )
{
	if( functionToWrap == NULL )
	{
		cerr << "WHAT ARE YOU DOING TO ME, YOUR TRYING TO INTEGRATE OVER A NULL OBJECT!!!" << endl;
		exit(-69);
	}

	double PDF_test_result = functionToWrap->Integral( NewDataPoint, NewBoundary );

	//cout << "PDF: " << PDF_test_result << endl;

	functionCanIntegrate = ! functionToWrap->GetNumericalNormalisation();

	if( functionCanIntegrate == false )
	{
		//	PDF doesn't know how to Normalise no need to check PDF!
		haveTestedIntegral = true;
	}

	bool cacheEnabled = functionToWrap->GetCachingEnabled();

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
	//Make a list of observables not to integrate
	vector<string> dontIntegrate = functionToWrap->GetDoNotIntegrateList();

	double testIntegral = functionToWrap->Integral( NewDataPoint, NewBoundary );
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

double RapidFitIntegrator::OneDimentionIntegral( const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary, ComponentRef* componentIndex, vector<string> doIntegrate, vector<string> dontIntegrate )
{
	IntegratorFunction* quickFunction = new IntegratorFunction( functionToWrap, NewDataPoint, doIntegrate, dontIntegrate, NewBoundary, componentIndex );
	//Find the observable range to integrate over
	IConstraint * newConstraint = NewBoundary->GetConstraint( doIntegrate[0] );
	double minimum = newConstraint->GetMinimum();
	double maximum = newConstraint->GetMaximum();

	//cout << "1D Integration" << endl;
	//Do a 1D integration
	oneDimensionIntegrator->SetFunction(*quickFunction);

	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;
	cout_bak = cout.rdbuf();
	cerr_bak = cerr.rdbuf();
	clog_bak = clog.rdbuf();
	cout.rdbuf(0);
	cerr.rdbuf(0);
	clog.rdbuf(0);

	double output = oneDimensionIntegrator->Integral( minimum, maximum );

	cout.rdbuf( cout_bak );
	cerr.rdbuf( cerr_bak );
	clog.rdbuf( clog_bak );

	delete quickFunction;

	return output;
}

double RapidFitIntegrator::MultiDimentionIntegral( IPDF* functionToWrap, AdaptiveIntegratorMultiDim* multiDimensionIntegrator, const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary,
		ComponentRef* componentIndex, vector<string> doIntegrate, vector<string> dontIntegrate )
{
	//Make arrays of the observable ranges to integrate over
	double* minima = new double[ doIntegrate.size() ];
	double* maxima = new double[ doIntegrate.size() ];

	for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
	{
		IConstraint * newConstraint = NewBoundary->GetConstraint( doIntegrate[observableIndex] );
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

	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;
	cout_bak = cout.rdbuf();
	cerr_bak = cerr.rdbuf();
	clog_bak = clog.rdbuf();
	//cout.rdbuf(0);
	//cerr.rdbuf(0);
	//clog.rdbuf(0);
	multiDimensionIntegrator->SetFunction( *quickFunction );
	double output =  multiDimensionIntegrator->Integral( minima, maxima );

	delete minima; delete maxima;
	delete quickFunction;

	cout.rdbuf( cout_bak );
	cerr.rdbuf( cerr_bak );
	clog.rdbuf( clog_bak );

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
double RapidFitIntegrator::DoNumericalIntegral( const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary, const vector<string> DontIntegrateThese, ComponentRef* componentIndex, const bool IntegrateDataPoint )
{
	//Make lists of observables to integrate and not to integrate
	vector<string> doIntegrate, dontIntegrate;
	StatisticsFunctions::DoDontIntegrateLists( functionToWrap, NewBoundary, &DontIntegrateThese, doIntegrate, dontIntegrate );

	dontIntegrate = StringProcessing::CombineUniques( dontIntegrate, DontIntegrateThese );
	dontIntegrate = this->DontNumericallyIntegrateList( NewDataPoint, dontIntegrate );

	vector<DataPoint*> DiscreteIntegrals;

	if( IntegrateDataPoint )
	{
		DiscreteIntegrals.push_back( new DataPoint(*NewDataPoint) );
	}
	else
	{
		DiscreteIntegrals = NewBoundary->GetDiscreteCombinations();
	}

	double output_val = 0.;

	//If there are no observables left to integrate over, just evaluate the function
	if( doIntegrate.size() == 0 )
	{
		for( vector<DataPoint*>::iterator dataPoint_i = DiscreteIntegrals.begin(); dataPoint_i != DiscreteIntegrals.end(); ++dataPoint_i )
		{
			output_val += functionToWrap->Integral( *dataPoint_i, const_cast<PhaseSpaceBoundary*>(NewBoundary) );
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
				numericalIntegral += this->OneDimentionIntegral( *dataPoint_i, NewBoundary, componentIndex, doIntegrate, dontIntegrate );
			}
			else
			{
				numericalIntegral += this->MultiDimentionIntegral( functionToWrap, multiDimensionIntegrator, *dataPoint_i, NewBoundary, componentIndex, doIntegrate, dontIntegrate );
			}

			if( !haveTestedIntegral )
			{
				double testIntegral = functionToWrap->Integral( *dataPoint_i, const_cast<PhaseSpaceBoundary*>(NewBoundary) );
				cout << "Integration Test: numerical : analytical  " << setw(7) << numericalIntegral << " : " << testIntegral;
				cout << "  " << NewBoundary->DiscreteDescription( *dataPoint_i );
			}

			output_val += numericalIntegral;
		}
	}

	while( !DiscreteIntegrals.empty() )
	{
		if( DiscreteIntegrals.back() != NULL ) delete DiscreteIntegrals.back();
		DiscreteIntegrals.pop_back();
	}

	return output_val;
}

//Return the integral over all observables except one
double RapidFitIntegrator::ProjectObservable( DataPoint* NewDataPoint, PhaseSpaceBoundary * NewBoundary, string ProjectThis, ComponentRef* Component )
{
	//Make the list of observables not to integrate
	vector<string> dontIntegrate = functionToWrap->GetDoNotIntegrateList();
	dontIntegrate.push_back(ProjectThis);
	double value = -1.;

	vector<string> allIntegrable = functionToWrap->GetPrototypeDataPoint();

	if( allIntegrable.size() == 1 )
	{
		if( allIntegrable[0] == ProjectThis )
		{
			value = functionToWrap->EvaluateComponent( NewDataPoint, Component );
		}
		else
		{
			cerr << "This PDF only knows about: " << allIntegrable[0] << " CAN NOT INTEGRATE: " << ProjectThis << endl << endl;
			exit(-342);
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


