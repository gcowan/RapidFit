/*
 * @class BasePDF
 *
 * Class that provides a general implementation of IPDF.
 * Can inherit from this to make a PDF without worrying about the details.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */
///	ROOT Headers
#include "TRandom3.h"
///	RapidFit Headers
#include "BasePDF.h"
#include "BasePDF_Framework.h"
#include "BasePDF_MCCaching.h"
#include "StringProcessing.h"
#include "ObservableRef.h"
#include "PhaseSpaceBoundary.h"
#include "RapidFitIntegrator.h"
///	System Headers
#include <iostream>
#include <cmath>
#include <sstream>
#include <pthread.h>

using namespace::std;

//Constructor
BasePDF::BasePDF() : BasePDF_Framework( this ), BasePDF_MCCaching(),
	numericalNormalisation(false), allParameters( vector<string>() ), allObservables(), doNotIntegrateList(), observableDistNames(), observableDistributions(),
	component_list(), requiresBoundary(false), cachingEnabled( true ), haveTestedIntegral( false ), discrete_Normalisation( false ), DiscreteCaches(new vector<double>()),
	debug_mutex(NULL), can_remove_mutex(true), debug(NULL), fixed_checked(false), isFixed(false), fixedID(0), _basePDFComponentStatus(false)
{
	debug = new DebugClass();
	component_list.push_back( "0" );
	debug_mutex = new pthread_mutex_t();
}

BasePDF::BasePDF( const BasePDF& input ) : BasePDF_Framework( input ), BasePDF_MCCaching( input ),
	numericalNormalisation( input.numericalNormalisation ), allParameters( input.allParameters.GetAllNames() ),
	allObservables( input.allObservables ), doNotIntegrateList( input.doNotIntegrateList ), observableDistNames( input.observableDistNames ), observableDistributions( input.observableDistributions ),
	component_list( input.component_list ), requiresBoundary( input.requiresBoundary ),
	cachingEnabled( input.cachingEnabled ), haveTestedIntegral( input.haveTestedIntegral ),
	discrete_Normalisation( input.discrete_Normalisation ), DiscreteCaches(NULL),
	debug_mutex(input.debug_mutex), can_remove_mutex(false), debug(NULL), fixed_checked(input.fixed_checked), isFixed(input.isFixed), fixedID(input.fixedID),
	_basePDFComponentStatus(input._basePDFComponentStatus)
{
	allParameters.SetPhysicsParameters( &(input.allParameters) );
	DiscreteCaches = new vector<double>( input.DiscreteCaches->size() );
	for( vector<double>::iterator cache_i = DiscreteCaches->begin(); cache_i != DiscreteCaches->end(); ++cache_i )
	{
		(*cache_i) = -1;
	}

	debug = (input.debug==NULL)?NULL:new DebugClass(*input.debug);

	this->GetPDFIntegrator()->SetPDF( this );
	this->GetPDFIntegrator()->SetDebug( debug );
}

//Destructor
BasePDF::~BasePDF()
{
	//	If we weren't using ROOT this would a) be safe and b) be advisible, but ROOT segfaults...
	//if( seed_function != NULL ) delete seed_function;
	if( DiscreteCaches != NULL ) delete DiscreteCaches;
	if( debug_mutex != NULL && can_remove_mutex == true ) delete debug_mutex;

	if( debug != NULL ) delete debug;
}

void BasePDF::SetComponentStatus( const bool input )
{
	this->ReallySetComponentStatus( input );
}

void BasePDF::ReallySetComponentStatus( const bool input )
{
	_basePDFComponentStatus = input;
}

bool BasePDF::GetComponentStatus() const
{
	return this->GetComponentStatus();
}

bool BasePDF::ReallyGetComponentStatus() const
{
	return _basePDFComponentStatus;
}

void BasePDF::TurnCachingOff()
{
	this->ReallyTurnCachingOff();
}

void BasePDF::ReallyTurnCachingOff()
{
	cachingEnabled=false;
}

bool BasePDF::CacheValid( DataPoint* InputPoint, PhaseSpaceBoundary* InputBoundary )
{
	if( this->IsDebuggingON() )
	{
		PDF_THREAD_LOCK
			cout << "BasePDF: Checking Cache" << endl;
		PDF_THREAD_UNLOCK
	}

	//	Force the cache to be not valid when caching is NOT enabled
	if( !cachingEnabled ) return false;

	//	Care about the input
	this->CheckCaches( InputBoundary );

	unsigned int cacheIndex = InputBoundary->GetDiscreteIndex( InputPoint );

	double norm = DiscreteCaches->at(cacheIndex);

	if( norm > 0 ) return true;
	else return false;
}

bool BasePDF::GetCachingEnabled() const
{
	return cachingEnabled;
}

void BasePDF::SetCachingEnabled( bool input )
{
	cachingEnabled = input;
}

bool BasePDF::GetNumericalNormalisation() const
{
	return numericalNormalisation;
}

//Externally Set whether the PDF wants to use Numerical Normalisation
void BasePDF::SetNumericalNormalisation( bool input )
{
	numericalNormalisation = input;
}

void BasePDF::SetCache( double input, DataPoint* InputPoint, PhaseSpaceBoundary* InputBoundary )
{
	if( this->IsDebuggingON() )
	{
		PDF_THREAD_LOCK
			cout << "BasePDF: Setting Cache. Discrete DataSet " << InputBoundary->GetDiscreteIndex( InputPoint ) << " Value: " << input << endl;
		PDF_THREAD_UNLOCK
	}

	this->CheckCaches( InputBoundary );
	/*!
	 * Get the Number Corresponding to the Unique Combination that this DataPoint is in
	 */
	unsigned int thisIndex = InputBoundary->GetDiscreteIndex( InputPoint );

	DiscreteCaches->at(thisIndex) = input;

}

void BasePDF::CheckCaches( PhaseSpaceBoundary* InputBoundary )
{
	/*!
	 * If the internal Caches don't match re-allocate the whole array to be -1
	 */
	unsigned int totalCombinations = (unsigned)InputBoundary->GetNumberCombinations();
	//cout << "I have: " << DiscreteCaches->size() << " I want : " << totalCombinations << endl;
	if( DiscreteCaches->size() != totalCombinations )
	{
		//cout << "Resizing!" << endl;
		DiscreteCaches->resize( totalCombinations );
		//cout << "Unsetting!" << endl;
		this->UnsetCache();
		//cout << "Zeroth element: " << DiscreteCaches->at(0) << endl;
	}
}

void BasePDF::UnsetCache()
{
	for(vector<double>::iterator this_i = DiscreteCaches->begin(); this_i != DiscreteCaches->end(); ++this_i )
	{
		*this_i = -1.;
	}
}


void BasePDF::UpdatePhysicsParameters( ParameterSet* Input )
{
	if( allParameters.GetAllNames().size() != 0 )
	{
		//  Invalidate the cache
		if( !(allParameters == (*Input)) ) this->UnsetCache();
		allParameters.SetPhysicsParameters( Input );
	}
	else
	{
		allParameters.AddPhysicsParameters( Input );
		this->UnsetCache();
	}

	this->SetPhysicsParameters( Input );
}

//Set the function parameters
bool BasePDF::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	(void) NewParameterSet;
	return true;
}

double BasePDF::GetCache( DataPoint* InputPoint, PhaseSpaceBoundary* NewBoundary )
{
	if( this->IsDebuggingON() )
	{
		PDF_THREAD_LOCK
			cout << "BasePDF: Requesting Cache For DataSet " << NewBoundary->GetDiscreteIndex( InputPoint ) << " of " << DiscreteCaches->size() << endl;
		PDF_THREAD_UNLOCK
	}

	this->CheckCaches( NewBoundary );

	unsigned int thisIndex = NewBoundary->GetDiscreteIndex( InputPoint );

	double thisVal = DiscreteCaches->at(thisIndex);

	if( this->IsDebuggingON() )
	{
		PDF_THREAD_LOCK
			cout << "BasePDF: Cache Value is " << thisVal << endl;
		PDF_THREAD_UNLOCK
	}


	if( thisVal > 0 )
	{
		return thisVal;
	}
	else
	{
		return -1.;
	}
}

bool BasePDF::CheckFixed( PhaseSpaceBoundary* NewBoundary )
{
	if( !( fixed_checked && (fixedID == NewBoundary->GetID()) ) )
	{
		fixedID = NewBoundary->GetID();
		vector<string> required = this->GetPrototypeDataPoint();
		isFixed = true;
		for( unsigned int i=0; i< required.size(); ++i )
		{
			isFixed = isFixed && NewBoundary->GetConstraint( required[i] )->IsDiscrete();
		}
		fixed_checked = true;
	}
	return isFixed;
}

//Return the integral of the function over the given boundary
double BasePDF::Integral(DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary)
{
	double thisCache = this->GetCache( NewDataPoint, NewBoundary );

	if( this->IsDebuggingON() )
	{
		if( cachingEnabled )
		{
			PDF_THREAD_LOCK
				cout << "BasePDF: Caching Enabled!" << endl;
			PDF_THREAD_UNLOCK
		}
		else
		{
			PDF_THREAD_LOCK
				cout << "BasePDF: Caching Disabled!" << endl;
			PDF_THREAD_UNLOCK
		}
	}

	if( this->CheckFixed( NewBoundary ) )
	{
		if( this->IsDebuggingON() )
		{
			PDF_THREAD_LOCK
				cout << "BasePDF: FixedPhaseSpace. Integral == Evaluate " << endl;
			PDF_THREAD_UNLOCK
		}
		if( this->CacheValid( NewDataPoint, NewBoundary ) )
		{
			return thisCache;
		}
		else
		{
			double thisEval = this->Evaluate( NewDataPoint );
			if( cachingEnabled )
			{
				this->SetCache( thisEval, NewDataPoint, NewBoundary );
			}
			return thisEval;
		}
	}

	if( this->IsDebuggingON() )
	{
		PDF_THREAD_LOCK
			cout << "BasePDF: thisCache " << thisCache << endl;
		PDF_THREAD_UNLOCK
	}

	// If the PDF has been configured to rely on Numerical Normalisation we can't calculate the Normalisation here at all
	// Just return what has been cached so far
	if( numericalNormalisation )
	{
		if( this->IsDebuggingON() )
		{
			PDF_THREAD_LOCK
				cout << "BasePDF: using Numerical Normalisation" << endl;
			PDF_THREAD_UNLOCK
		}

		if( this->CacheValid( NewDataPoint, NewBoundary ) )
		{
			return thisCache;
		}
		else
		{
			if( this->IsDebuggingON() )
			{
				PDF_THREAD_LOCK
					cout << "BasePDF: Calculating Numerical Integral ";
				PDF_THREAD_UNLOCK
			}

			double thisNumericalIntegral =-1.;

			try
			{
				NewDataPoint->SetPhaseSpaceBoundary( NewBoundary );
				thisNumericalIntegral = this->GetPDFIntegrator()->NumericallyIntegrateDataPoint( NewDataPoint, NewBoundary, this->GetDoNotIntegrateList() );
			}
			catch(...)
			{
				//	It's legitimate to fall over due to this throwing as we don't know what to do afterwards
				cerr << "CANNOT ANALYTICALLY OR NUMERICALLY INTEGRATE THIS PDF IT FELL OVER!!\t" << this->GetLabel() << endl;
				exit(-834);
			}
			if( this->IsDebuggingON() )
			{
				PDF_THREAD_LOCK
					cout << thisNumericalIntegral << endl;
				PDF_THREAD_UNLOCK
			}

			if( thisNumericalIntegral < 0 )
			{
				//	Component Integrals may be 0 by definition and the PDF can potentially be constructed in a non-trivial way
				//	We MUST allow for this possibility
				if( !( this->ReallyGetComponentStatus() ) )
				{
					cerr << "CANNOT ANALYTICALLY OR NUMERICALLY INTEGRATE THIS PDF!!\t" << this->GetLabel() << endl;
					exit(-833);
				}
			}

			if( cachingEnabled )
			{
				this->SetCache( thisNumericalIntegral, NewDataPoint, NewBoundary );
			}
			return thisNumericalIntegral;
		}
	}


	//	Once we have checked which analytical form to use:
	if( haveTestedIntegral )
	{
		if( this->IsDebuggingON() )
		{
			PDF_THREAD_LOCK
				cout << "BasePDF: Have Previously Tested Analytical Integral and it is GOOD :D" << endl;
			PDF_THREAD_UNLOCK
		}

		if( this->CacheValid( NewDataPoint, NewBoundary ) )
		{
			if( this->IsDebuggingON() )
			{
				PDF_THREAD_LOCK
					cout << "BasePDF: Using Cache " << thisCache << endl;
				PDF_THREAD_UNLOCK
			}
			return thisCache;
		}
		else
		{
			if( this->IsDebuggingON() )
			{
				PDF_THREAD_LOCK
					cout << "BasePDF: Calculating Integral Value" << endl;
				PDF_THREAD_UNLOCK
			}

			if( discrete_Normalisation )
			{
				if( this->IsDebuggingON() )
				{
					PDF_THREAD_LOCK
						cout << "BasePDF: Using Normalisation( DataPoint*, PhaseSpaceBoundary* )" << endl;
					PDF_THREAD_UNLOCK
				}

				double norm = this->Normalisation( NewDataPoint, NewBoundary );
				if( this->GetCachingEnabled() )
				{
					this->UnsetCache();
				}
				this->TurnCachingOff();
				//this->SetCache( norm, NewDataPoint, NewBoundary );
				return norm;
			}
			else
			{
				if( this->IsDebuggingON() )
				{
					PDF_THREAD_LOCK
						cout << "BasePDF: Using Normalisation( PhaseSpaceBoundary* )" << endl;
					PDF_THREAD_UNLOCK
				}

				double norm = this->Normalisation(NewBoundary);
				if( this->GetCachingEnabled() ) this->SetCachingEnabled( true );
				this->SetCache( norm, NewDataPoint, NewBoundary );
				return norm;
			}
		}
	}

	if( this->IsDebuggingON() )
	{
		PDF_THREAD_LOCK
			cout << "BasePDF: Testing Analytical Normalisation" << endl;
		PDF_THREAD_UNLOCK
	}
	//	Here we don't know if we're using an e_by_e or using a cache, let's test them


	//	Checking for event independent norm
	double Norm = this->Normalisation(NewBoundary);
	if( Norm >= 0. )		//	Valid
	{
		haveTestedIntegral = true;
		numericalNormalisation = false;
		discrete_Normalisation = false;
		cachingEnabled = true;
		this->SetCache( Norm, NewDataPoint, NewBoundary );
		if( this->IsDebuggingON() )
		{
			PDF_THREAD_LOCK
				cout << "BasePDF: PDF can use Normalisation( PhaseSpaceBoundary* )" << endl;
			PDF_THREAD_UNLOCK
		}
		return Norm;
	}

	Norm = this->Normalisation(NewDataPoint, NewBoundary);
	if( Norm >= 0. )		//	Valid
	{
		haveTestedIntegral = true;
		numericalNormalisation = false;
		discrete_Normalisation = true;
		if( cachingEnabled ) this->SetCache( Norm, NewDataPoint, NewBoundary );
		if( this->IsDebuggingON() )
		{
			PDF_THREAD_LOCK
				cout << "BasePDF: PDF can use Normalisation( DataPoint*, PhaseSpaceBoundary* )" << endl;
			PDF_THREAD_UNLOCK
		}
		return Norm;
	}

	//	PDF does NOT know how to Normalise!!
	numericalNormalisation = true;
	discrete_Normalisation = true;

	if( this->IsDebuggingON() )
	{
		PDF_THREAD_LOCK
			cout << "BasePDF: PDF requires Numerical Normalisation" << endl;
		PDF_THREAD_UNLOCK
	}

	double thisNumericalIntegral = this->GetPDFIntegrator()->NumericallyIntegrateDataPoint( NewDataPoint, NewBoundary, this->GetDoNotIntegrateList() );

	if( thisNumericalIntegral < 0 )
	{
		//      Component Integrals may be 0 by definition and the PDF can potentially be constructed in a non-trivial way
		//      We MUST allow for this possibility
		if( !( this->ReallyGetComponentStatus() ) )
		{
			cerr << "CANNOT ANALYTICALLY OR NUMERICALLY INTEGRATE THIS PDF!!\t" << this->GetLabel() << endl;
			cerr << "Integral: " << thisNumericalIntegral << endl;
			exit(-832);
		}
	}

	if( cachingEnabled )
	{
		this->SetCache( thisNumericalIntegral, NewDataPoint, NewBoundary );
	}

	return thisNumericalIntegral;
}

//Do the integration
double BasePDF::Normalisation(PhaseSpaceBoundary * NewBoundary)
{
	(void)NewBoundary;
	//Just a default value
	return -1.0;
}

//Do the integration
double BasePDF::Normalisation(DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary)
{
	(void)NewDataPoint;
	(void)NewBoundary;
	//Just a default value
	return -1.0;
}

//Calculate the function value
double BasePDF::Evaluate(DataPoint * NewDataPoint)
{
	(void)NewDataPoint;
	//Just a default value
	return  -1.0;
}

//Return the function value at the given point for generation
double BasePDF::EvaluateForNumericGeneration( DataPoint* NewDataPoint )
{
	double returnable=0.;
	try
	{
		returnable = this->Evaluate( NewDataPoint );
	}
	catch(...)
	{
		PDF_THREAD_LOCK
			cerr << "BasePDF: Failed to Correctly Generate a Sensible Value here:" << endl;
		NewDataPoint->Print();
		cerr << "BasePDF: Returning 0.!!!" << endl;
		PDF_THREAD_UNLOCK
			return 0.;
	}
	return returnable;
}

//Calculate the function value for numerical integration
double BasePDF::EvaluateForNumericIntegral( DataPoint * NewDataPoint )
{
	double returnable=0.;
	try
	{
		returnable = this->Evaluate( NewDataPoint );
	}
	catch(...)
	{
		PDF_THREAD_LOCK
			cerr << "BasePDF: Failed to Correctly Integrate a Sensible Value here:" << endl;
		NewDataPoint->Print();
		cerr << "BasePDF: Returning 0.!!!" << endl;
		PDF_THREAD_UNLOCK
			return 0.;
	}
	return returnable;
}

double BasePDF::EvaluateTimeOnly( DataPoint* )
{
	return -1.;
}

ParameterSet* BasePDF::GetPhysicsParameters()
{
	return &allParameters;
}

//Return a prototype data point
vector<string> BasePDF::GetPrototypeDataPoint()
{
	return allObservables;
}

//Return a prototype set of physics parameters
vector<string> BasePDF::GetPrototypeParameterSet()
{
	return allParameters.GetAllNames();
}

//Return a list of parameters not to be integrated
vector<string> BasePDF::GetDoNotIntegrateList()
{
	return doNotIntegrateList;
}

vector<string> BasePDF::PDFComponents()
{
	return component_list;
}

//Calculate the function value
double BasePDF::EvaluateComponent( DataPoint* NewDataPoint, ComponentRef* componentIndex )
{
	(void) componentIndex;
	//Just assume a single component equal to the ordinary evaluate method
	return this->EvaluateForNumericIntegral( NewDataPoint );
}

string BasePDF::GetComponentName( ComponentRef* input )
{
	if( input == NULL ) return this->GetName();
	else
	{
		if( input->getComponentName() == "0" ) return this->GetName();
		else return string( this->GetName() + "::" + input->getComponentName() );
	}
}

