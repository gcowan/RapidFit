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
//	System Headers
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <float.h>
#include <pthread.h>
#include <exception>

#ifdef __RAPIDFIT_USE_GSL
// GSL for MC integration
#include <gsl/gsl_qrng.h>
#endif

pthread_mutex_t multi_mutex;
pthread_mutex_t multi_mutex2;
pthread_mutex_t multi_mutex3;

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

using namespace::std;

//1% tolerance on integral value
//const double INTEGRAL_PRECISION_THRESHOLD = 0.01;

//Constructor with correct argument
RapidFitIntegrator::RapidFitIntegrator( IPDF * InputFunction, bool ForceNumerical, bool UsePseudoRandomIntegration ) :
	ratioOfIntegrals(-1.), fastIntegrator(NULL), functionToWrap(InputFunction), multiDimensionIntegrator(NULL), oneDimensionIntegrator(NULL),
	functionCanIntegrate(false), haveTestedIntegral(false),
	RapidFitIntegratorNumerical( ForceNumerical ), obs_check(false), checked_list(), debug(new DebugClass(false) ), pseudoRandomIntegration( UsePseudoRandomIntegration )
{
	multiDimensionIntegrator = new AdaptiveIntegratorMultiDim();
#if ROOT_VERSION_CODE > ROOT_VERSION(5,28,0)
	multiDimensionIntegrator->SetAbsTolerance( 1E-9 );      //      Absolute error for things such as plots
	multiDimensionIntegrator->SetRelTolerance( 1E-9 );
	multiDimensionIntegrator->SetMaxPts( 1000000 );		//	These are the defaults, and it's unlikely you will be able to realistically push the integral without using "double double"'s
#endif

	ROOT::Math::IntegrationOneDim::Type type = ROOT::Math::IntegrationOneDim::kGAUSS;
	oneDimensionIntegrator = new IntegratorOneDim(type);
}

RapidFitIntegrator::RapidFitIntegrator( const RapidFitIntegrator& input ) : ratioOfIntegrals( input.ratioOfIntegrals ),
	fastIntegrator( NULL ), functionToWrap( input.functionToWrap ), multiDimensionIntegrator( NULL ), oneDimensionIntegrator( NULL ),
	pseudoRandomIntegration(input.pseudoRandomIntegration),
	functionCanIntegrate( input.functionCanIntegrate ), haveTestedIntegral( true ),//input.haveTestedIntegral ),
	RapidFitIntegratorNumerical( input.RapidFitIntegratorNumerical ), obs_check( input.obs_check ), checked_list( input.checked_list ), debug((input.debug==NULL)?NULL:new DebugClass(*input.debug))
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
		multiDimensionIntegrator->SetRelTolerance( 1E-9 );
		multiDimensionIntegrator->SetMaxPts( 1000000 );
#endif
	}
	if( input.oneDimensionIntegrator != NULL )
	{
		ROOT::Math::IntegrationOneDim::Type type = ROOT::Math::IntegrationOneDim::kGAUSS;
		(void) type;	//	For the global state of ROOT
		oneDimensionIntegrator = new IntegratorOneDim( );//*(input.oneDimensionIntegrator) );
	}

	if(input.debug!=NULL) if( !(input.debug->GetStatus()) ) debug->SetStatus(false);
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
	if( functionToWrap == NULL )
	{
		cerr << "WHAT ARE YOU DOING TO ME, YOUR TRYING TO INTEGRATE OVER A NULL OBJECT!!!" << endl;
		exit(-69);
	}

	double PDF_test_result = 0.;

	try
	{
		PDF_test_result = functionToWrap->Integral( NewDataPoint, NewBoundary );
	}
	catch(...)
	{
		cerr << "RapidFitIntegrator: PDF " << functionToWrap->GetLabel() << " has failed the basic test of not crashing, check your PhaseSpace and/or return -1 for Numerical Integration" << endl;
		throw(1280);
	}

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
	IntegratorFunction* quickFunction = new IntegratorFunction( functionToWrap, NewDataPoint, doIntegrate, dontIntegrate, NewBoundary, componentIndex );
	if( debug != NULL ) quickFunction->SetDebug( debug );
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

	//cout << "1D Integration" << endl;
	//Do a 1D integration
	oneDimensionIntegrator->SetFunction(*quickFunction);

	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;

	if( debug != NULL )
	{
		if( !debug->GetStatus() || haveTestedIntegral )
		{
			cout_bak = cout.rdbuf();
			cerr_bak = cerr.rdbuf();
			clog_bak = clog.rdbuf();
			cout.rdbuf(0);
			cerr.rdbuf(0);
			clog.rdbuf(0);
		}
	}

	double output = oneDimensionIntegrator->Integral( minimum, maximum );

	if( debug != NULL )
	{
		if( !debug->GetStatus() || haveTestedIntegral )
		{
			cout.rdbuf( cout_bak );
			cerr.rdbuf( cerr_bak );
			clog.rdbuf( clog_bak );
		}
	}

	delete quickFunction;

	return output;
}

#ifdef __RAPIDFIT_USE_GSL
double RapidFitIntegrator::PseudoRandomNumberIntegral( IPDF* functionToWrap, const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary, ComponentRef* componentIndex, 
				vector<string> doIntegrate )
{
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

	//Do a 4D integration
	int npoint = 10000;
	std::vector<double> integrationPoints1;
	std::vector<double> integrationPoints2;
	std::vector<double> integrationPoints3;
	std::vector<double> integrationPoints4;

	gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, 4);
	for (int i = 0; i < npoint; i++)
	{
		double v[4];
		gsl_qrng_get (q, v);
		integrationPoints1.push_back(v[0]);
		integrationPoints2.push_back(v[1]);
		integrationPoints3.push_back(v[2]);
		integrationPoints4.push_back(v[3]);
	}
	gsl_qrng_free (q);

	double result=0;
	DataPoint * point = new DataPoint(*NewDataPoint);
	for (unsigned int i=0; i<integrationPoints1.size();++i)
	{
		point->SetObservable("m23", 	  integrationPoints1[i]*(maxima[0]-minima[0])+minima[0], 0.0, "GeV/c^{2}"); 
		point->SetObservable("cosTheta1", integrationPoints2[i]*(maxima[1]-minima[1])+minima[1], 0.0, " "); 
		point->SetObservable("cosTheta2", integrationPoints3[i]*(maxima[2]-minima[2])+minima[2], 0.0, " "); 
		point->SetObservable("phi",       integrationPoints4[i]*(maxima[3]-minima[3])+minima[3], 0.0, "rad"); 
		result += functionToWrap->Evaluate( point );
	}
	result /= double(integrationPoints1.size());

	delete minima; delete maxima;
	delete point;
	return result;
}
#endif

double RapidFitIntegrator::MultiDimentionIntegral( IPDF* functionToWrap, AdaptiveIntegratorMultiDim* multiDimensionIntegrator, const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary,
		ComponentRef* componentIndex, vector<string> doIntegrate, vector<string> dontIntegrate, bool haveTestedIntegral, DebugClass* debug )
{
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
	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;

	if( debug != NULL )
	{
		if( !debug->GetStatus() || haveTestedIntegral )
		{
			cout_bak = cout.rdbuf();
			cerr_bak = cerr.rdbuf();
			clog_bak = clog.rdbuf();
			cout.rdbuf(0);
			cerr.rdbuf(0);
			clog.rdbuf(0);
		}
	}

	multiDimensionIntegrator->SetFunction( *quickFunction );

	double output =  multiDimensionIntegrator->Integral( minima, maxima );

	delete minima; delete maxima;
	delete quickFunction;

	if( debug != NULL )
	{
		if( !debug->GetStatus() || haveTestedIntegral )
		{
			cout.rdbuf( cout_bak );
			cerr.rdbuf( cerr_bak );
			clog.rdbuf( clog_bak );
		}
	}

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
	if( doIntegrate.empty() || doIntegrate.size() == 0 )
	{
		for( vector<DataPoint*>::iterator dataPoint_i = DiscreteIntegrals.begin(); dataPoint_i != DiscreteIntegrals.end(); ++dataPoint_i )
		{
			output_val += functionToWrap->Integral( *dataPoint_i, NewBoundary );
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
				numericalIntegral += this->OneDimentionIntegral( functionToWrap, oneDimensionIntegrator, *dataPoint_i, NewBoundary, componentIndex, doIntegrate, dontIntegrate, debug );
			}
			else
			{
				if ( !pseudoRandomIntegration )
				{
					numericalIntegral += this->MultiDimentionIntegral( functionToWrap, multiDimensionIntegrator, *dataPoint_i, NewBoundary, componentIndex, doIntegrate, dontIntegrate, debug );
				} else {
					numericalIntegral += this->PseudoRandomNumberIntegral( functionToWrap, *dataPoint_i, NewBoundary, componentIndex, doIntegrate );
				}
			}
			if( !haveTestedIntegral )
			{
				double testIntegral = functionToWrap->Integral( *dataPoint_i, NewBoundary );
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

	dontIntegrate.push_back(ProjectThis);

	if( testedIntegrable.size() == 1 )
	{
		if( testedIntegrable[0] == ProjectThis )
		{
			value = functionToWrap->EvaluateComponent( NewDataPoint, Component );
		}
		else
		{
			cerr << "This PDF only knows how to Integrate: " << testedIntegrable[0] << " CAN NOT INTEGRATE: " << ProjectThis << endl << endl;
			throw(-342);
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
	if( input_debug != NULL )
	{
		if( debug != NULL ) delete debug;
		debug = new DebugClass(*input_debug);

		if( debug->DebugThisClass( "RapidFitIntegrator" ) )
		{
			debug->SetStatus(true);
			cout << "RapidFitIntegrator: Debugging Enabled!" << endl;
		}
		else
		{
			debug->SetStatus(false);
		}
	}
	else
	{
		if( debug != NULL ) delete debug;
		debug = new DebugClass(false);
	}
	if( functionToWrap != NULL ) functionToWrap->SetDebug( input_debug );
}

