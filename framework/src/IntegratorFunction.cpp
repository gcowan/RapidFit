/**
  @class IntegratorFunction

  A wrapper to make IPDF usable by the numerical integrator

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-8
 */

//	RapidFit Headers
#include "StringProcessing.h"
#include "IntegratorFunction.h"
#include "ObservableRef.h"
#include "ClassLookUp.h"
#include "PhaseSpaceBoundary.h"
#include "ComponentRef.h"
//	System Headers
#include <iostream>
#include <cstdlib>
#include <float.h>

using namespace::std;

//	Constructor for Integrator Objects
IntegratorFunction::IntegratorFunction( IPDF * InputFunction, const DataPoint * InputPoint, vector<string> IntegrateThese, vector<string> DontIntegrateThese,
		const PhaseSpaceBoundary* inputPhaseSpaceBoundary, ComponentRef* Index, vector<double> new_lower_limit, vector<double> new_upper_limit ) :
	wrappedFunction(ClassLookUp::CopyPDF(InputFunction)), currentPoint(new DataPoint(*InputPoint) ), doIntegrate(IntegrateThese), dontIntegrate(DontIntegrateThese),
	minima(), ranges(), cache_positions(), componentIndex( Index==NULL?NULL:(new ComponentRef(*Index)) ), newDataPoint(NULL), cache_lookup(),
	lower_limit(new_lower_limit), upper_limit(new_upper_limit), generateFunc(false), integrateFunc(true), myPhaseSpaceBoundary( new PhaseSpaceBoundary( *inputPhaseSpaceBoundary ) ),
	debug(new DebugClass(false))
{
	//	Chose to use perform the lookups in the constructor to keep the Eval statement as const as possible
	vector<ObservableRef> lookups;
	for( vector<string>::iterator val_i=doIntegrate.begin(); val_i != doIntegrate.end(); ++val_i )
	{
		lookups.push_back( ObservableRef( *val_i ) );
	}
	for( vector<string>::iterator val_i=dontIntegrate.begin(); val_i != dontIntegrate.end(); ++val_i )
	{
		lookups.push_back( ObservableRef( *val_i ) );
	}
	for( vector<ObservableRef>::iterator obs_i=lookups.begin(); obs_i != lookups.end(); ++obs_i )
	{
		currentPoint->GetObservable( *obs_i );
		cache_lookup.push_back( obs_i->GetIndex() );
	}
	newDataPoint = new DataPoint( *InputPoint );
	newDataPoint->SetPhaseSpaceBoundary( myPhaseSpaceBoundary );
	currentPoint->SetPhaseSpaceBoundary( myPhaseSpaceBoundary );
}

//	Constructor with additional information needed for the Foam coordinate transform
IntegratorFunction::IntegratorFunction( IPDF * InputFunction, const DataPoint * InputPoint, vector<string> IntegrateThese, vector<string> DontIntegrateThese,
		vector<double> InputMinima, vector<double> InputRanges, const PhaseSpaceBoundary* inputPhaseSpaceBoundary ) :
	wrappedFunction(ClassLookUp::CopyPDF(InputFunction)), currentPoint(new DataPoint(*InputPoint) ), doIntegrate(IntegrateThese), dontIntegrate(DontIntegrateThese),
	minima(InputMinima), ranges(InputRanges), cache_positions(), componentIndex(NULL), newDataPoint(NULL), cache_lookup(), lower_limit(InputMinima),
	upper_limit(), generateFunc(true), integrateFunc(false), myPhaseSpaceBoundary( new PhaseSpaceBoundary( *inputPhaseSpaceBoundary) ),
	debug(new DebugClass(false))
{
	//      Chose to use perform the lookups in the constructor to keep the Eval statement as const as possible
	vector<ObservableRef> lookups;
	for( vector<string>::iterator val_i=doIntegrate.begin(); val_i != doIntegrate.end(); ++val_i )
	{
		lookups.push_back( ObservableRef( *val_i ) );
	}
	for( vector<string>::iterator val_i=dontIntegrate.begin(); val_i != dontIntegrate.end(); ++val_i )
	{
		lookups.push_back( ObservableRef( *val_i ) );
	}
	for( vector<ObservableRef>::iterator obs_i=lookups.begin(); obs_i != lookups.end(); ++obs_i )
	{
		currentPoint->GetObservable( *obs_i );
		cache_lookup.push_back( obs_i->GetIndex() );
	}
	newDataPoint = new DataPoint( *InputPoint );
	newDataPoint->SetPhaseSpaceBoundary( myPhaseSpaceBoundary );
	currentPoint->SetPhaseSpaceBoundary( myPhaseSpaceBoundary );
	vector<double>::iterator num_i=lower_limit.begin();
	vector<double>::iterator num_j=ranges.begin();
	for( ; num_i != lower_limit.end(); ++num_i, ++num_j )
	{
		upper_limit.push_back( *num_i + *num_j );
	}
}

//Destructor
IntegratorFunction::~IntegratorFunction()
{
	if( wrappedFunction != NULL ) delete wrappedFunction;
	if( currentPoint != NULL ) delete currentPoint;
	if( newDataPoint != NULL ) delete newDataPoint;
	if( myPhaseSpaceBoundary != NULL ) delete myPhaseSpaceBoundary;
	if( componentIndex != NULL ) delete componentIndex;
	if( debug != NULL ) delete debug;
}

IntegratorFunction::IntegratorFunction ( const IntegratorFunction& input ) :
	IBaseFunctionMultiDim( input ), IBaseFunctionOneDim( input ), TFoamIntegrand( input ),
	wrappedFunction( ClassLookUp::CopyPDF( input.wrappedFunction ) ), currentPoint( new DataPoint(*input.currentPoint) ), doIntegrate( input.doIntegrate ), dontIntegrate( input.dontIntegrate ),
	minima( input.minima ), ranges( input.ranges ), cache_positions( input.cache_positions ), componentIndex( NULL ),
	newDataPoint( new DataPoint(*input.newDataPoint) ), cache_lookup( input.cache_lookup ), lower_limit( input.lower_limit ), upper_limit( input.upper_limit ), generateFunc( input.generateFunc ),
	integrateFunc( input.integrateFunc ), myPhaseSpaceBoundary( new PhaseSpaceBoundary( *input.myPhaseSpaceBoundary ) ), debug( (input.debug==NULL)?NULL:new DebugClass(*input.debug) )
{
	if( input.componentIndex != NULL )
	{
		componentIndex = new ComponentRef( *(input.componentIndex) );
	}
}

//Return the IPDF inside the wrapper
IPDF * IntegratorFunction::GetWrappedFunction() const
{
	return wrappedFunction;
}

//Copy the wrapper
IntegratorFunction * IntegratorFunction::Clone() const
{
	return new IntegratorFunction( *this );//wrappedFunction, currentPoint, doIntegrate, dontIntegrate, componentIndex );
}

//Copy the wrapper
IntegratorFunction * IntegratorFunction::Clone( const char *newname ) const
{
	(void) newname;
	return new IntegratorFunction( *this );//wrappedFunction, currentPoint, doIntegrate, dontIntegrate, componentIndex );
}

//Return the number of dimensions
unsigned int IntegratorFunction::NDim() const
{
	return unsigned(int(doIntegrate.size()));
}

//Return the function value at x
double IntegratorFunction::operator()( const Double_t * x ) const
{
	return DoEval(x);
}
double IntegratorFunction::operator()( Double_t x ) const
{
	return DoEval(x);
}

//Assignment operator
IntegratorFunction& IntegratorFunction::operator=( const IntegratorFunction & NewFunction )
{
	if( &NewFunction != this )
	{
		if( this->wrappedFunction != NULL ) delete this->wrappedFunction;
		this->wrappedFunction = ClassLookUp::CopyPDF( NewFunction.wrappedFunction );
		if( this->currentPoint != NULL ) delete this->currentPoint;
		this->currentPoint = new DataPoint( *NewFunction.currentPoint );
		this->doIntegrate = NewFunction.doIntegrate;
		this->dontIntegrate = NewFunction.dontIntegrate;
		this->minima = NewFunction.minima;
		this->ranges = NewFunction.ranges;
		this->cache_positions = NewFunction.cache_positions;
		if( this->componentIndex != NULL ) delete this->componentIndex;
		this->componentIndex = NewFunction.componentIndex==NULL ? NULL : (new ComponentRef(*NewFunction.componentIndex));
		if( this->newDataPoint != NULL ) delete this->newDataPoint;
		this->newDataPoint = new DataPoint( *NewFunction.newDataPoint );
		this->cache_lookup = NewFunction.cache_lookup;
		this->lower_limit = NewFunction.lower_limit;
		this->upper_limit = NewFunction.upper_limit;
		this->generateFunc = NewFunction.generateFunc;
		this->integrateFunc = NewFunction.integrateFunc;
		if( this->myPhaseSpaceBoundary != NULL ) delete this->myPhaseSpaceBoundary;
		this->myPhaseSpaceBoundary = new PhaseSpaceBoundary( *NewFunction.myPhaseSpaceBoundary );
	}
	return *(this);
}

//Return the function value at x
double IntegratorFunction::DoEval( const Double_t * x ) const
{
	if( currentPoint->GetPhaseSpaceBoundary() == NULL ) currentPoint->SetPhaseSpaceBoundary( myPhaseSpaceBoundary );
	currentPoint->ClearPseudoObservable();
	unsigned int true_index=100;
	//Load the array into the data point
	for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
	{
		true_index = (unsigned)cache_lookup[observableIndex];
		Observable* currentObservable = currentPoint->GetObservable( true_index );
		currentObservable->SetBinNumber(-1);
		if( lower_limit.size() != 0 )           //Integrating over a fixed variable then this won't be propogated correctly
		{
			if( (double)x[observableIndex] < lower_limit[observableIndex] )
			{
				//cout << x[observableIndex] << " :\t" << lower_limit[observableIndex] << endl;
				return 0.;
			}
		}
		if( upper_limit.size() != 0 )           //Integrating over a fixed variable then this won't be propogated correctly
		{
			if( (double)x[observableIndex] > upper_limit[observableIndex] )
			{
				//cout << x[observableIndex] << " :\t" << upper_limit[observableIndex] << endl;
				return 0.;
			}
		}
		if( lower_limit.empty() )
		{

		}
		newDataPoint->SetObservable( doIntegrate[observableIndex], (double)x[observableIndex], " ", true, (int)true_index );
		newDataPoint->GetObservable( true_index )->SetBinNumber(-1);
		newDataPoint->GetObservable( true_index )->SetBkgBinNumber(-1);
	}

	//	Fixed values, no need to set them equal to themselves
	//Load values of other observables
	for (unsigned int observableIndex = 0; observableIndex < dontIntegrate.size(); ++observableIndex )
	{
		true_index = (unsigned)cache_lookup[observableIndex + (unsigned)doIntegrate.size()];
		Observable* currentObservable = currentPoint->GetObservable( true_index );
		currentObservable->SetBinNumber(-1);
		currentObservable->SetBkgBinNumber(-1);
		newDataPoint->SetObservable( dontIntegrate[observableIndex], currentObservable->GetValue(), " ", true, (int)true_index );
	}

	if( newDataPoint->GetPhaseSpaceBoundary() == NULL ) newDataPoint->SetPhaseSpaceBoundary( myPhaseSpaceBoundary );

	double result=-1;

	//cout << "IntegratorFunction::Eval: " << componentIndex << endl;

	try
	{
		//	if( ! (componentIndex.empty() || componentIndex == "") || !(componentIndex == "0") )
		if( componentIndex != NULL )
		{
			result = wrappedFunction->EvaluateComponent( newDataPoint, componentIndex );
			if( std::isnan(result) || fabs(result)>=DBL_MAX )
			{
				if( debug->DebugThisClass( "IntegratorFunction" ) )
				{
					cout << "Component Value:" << result << endl;
					newDataPoint->Print();
					myPhaseSpaceBoundary->Print();
				}
			}
		}
		else if( integrateFunc == true )
		{
			result = wrappedFunction->EvaluateForNumericIntegral( newDataPoint );
			if( std::isnan(result) || fabs(result)>=DBL_MAX )
			{
				if( debug->DebugThisClass( "IntegratorFunction" ) )
				{
					cout << "Evaluate for Numerical Integral Value:" << result << endl;
					newDataPoint->Print();
					myPhaseSpaceBoundary->Print();
				}
			}
		}
		else if( generateFunc == true )
		{
			result = wrappedFunction->EvaluateForNumericGeneration( newDataPoint );
		}
		else
		{
			result = wrappedFunction->Evaluate( newDataPoint );
		}
	}
	catch ( int e )
	{
		result = 0.;
	}
	catch (...)
	{
		result = 0.;
	}

	//if( componentIndex == NULL )
	//{
	//	if( result < 0. ) result = 0.;
	//}

	if( fabs(result) >= DBL_MAX ) result = 0.;
	if( std::isnan(result) ) result = 0.;

	//if( result == 0 ) newDataPoint->Print();
	newDataPoint->ClearPseudoObservable();
	return result;
}

double IntegratorFunction::DoEval( Double_t x ) const
{
	if ( doIntegrate.size() == 1 )
	{
		Double_t* xArray = new Double_t[1];
		xArray[0] = x;
		double return_val = DoEval(xArray);
		delete [] xArray;
		return return_val;
	}
	else
	{
		cerr << "One dimensional evaluation has been called for multidimensional function" << endl;
		return 0.0;
	}
}

Double_t IntegratorFunction::Density( Int_t ndim, Double_t * xArray )
{
	if ( ndim == int(doIntegrate.size()) )
	{
		//Coordinate transform
		double* transformedArray = new double[unsigned(ndim)];
		for ( int observableIndex = 0; observableIndex < ndim; ++observableIndex )
		{
			transformedArray[unsigned(observableIndex)] = minima[unsigned(observableIndex)] + ( ranges[unsigned(observableIndex)] * xArray[unsigned(observableIndex)] );
		}
		double return_val = DoEval(transformedArray);
		delete [] transformedArray;
		return return_val;
	}
	else
	{
		cerr << "TFoamIntegrand problem - dimension number mismatch" << endl;
		return 0.0;
	}
}

void IntegratorFunction::SetDebug( DebugClass* input_debug )
{
	if( debug != NULL ) delete debug;
	debug = new DebugClass( *input_debug );
}

DataPoint* IntegratorFunction::GetCurrentDataPoint() const
{
	return newDataPoint;
}


