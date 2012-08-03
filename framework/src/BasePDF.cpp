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
#include "StringProcessing.h"
#include "ObservableRef.h"
#include "PhaseSpaceBoundary.h"
///	System Headers
#include <iostream>
#include <cmath>
#include <sstream>
#include <pthread.h>

using namespace::std;

//Constructor
BasePDF::BasePDF() : numericalNormalisation(false), allParameters( vector<string>() ), allObservables(), doNotIntegrateList(), observableDistNames(), observableDistributions(),
	component_list(), cached_files(), hasCachedMCGenerator(false), seed_function(NULL), seed_num(0), PDFName("Base"), PDFLabel("Base"), copy_object( NULL ), requiresBoundary(false),
	do_i_control_the_cache(false), cachingEnabled( true ), haveTestedIntegral( false ), thisConfig(NULL), discrete_Normalisation( false ), DiscreteCaches(new vector<double>()),
	debug_mutex(NULL), can_remove_mutex(true), debug(new DebugClass(false) )
{
	component_list.push_back( "0" );
	debug_mutex = new pthread_mutex_t();
}

BasePDF::BasePDF( const BasePDF& input ) :
	numericalNormalisation( input.numericalNormalisation ), allParameters( input.allParameters.GetAllNames() ),
	allObservables( input.allObservables ), doNotIntegrateList( input.doNotIntegrateList ), observableDistNames( input.observableDistNames ), observableDistributions( input.observableDistributions ),
	component_list( input.component_list ), cached_files( input.cached_files ), hasCachedMCGenerator( input.hasCachedMCGenerator ), requiresBoundary( input.requiresBoundary ),
	seed_function( input.seed_function ), seed_num( input.seed_num ), PDFName( input.PDFName ), PDFLabel( input.PDFLabel ), copy_object( input.copy_object ),
	do_i_control_the_cache( input.do_i_control_the_cache ), cachingEnabled( input.cachingEnabled ), haveTestedIntegral( input.haveTestedIntegral ),
	thisConfig(NULL), discrete_Normalisation( input.discrete_Normalisation ), DiscreteCaches(NULL),
	debug_mutex(input.debug_mutex), can_remove_mutex(false), debug((input.debug==NULL)?NULL:new DebugClass(*input.debug))
{
	allParameters.SetPhysicsParameters( &(input.allParameters) );
	DiscreteCaches = new vector<double>( input.DiscreteCaches->size() );
	for( vector<double>::iterator cache_i = DiscreteCaches->begin(); cache_i != DiscreteCaches->end(); ++cache_i )
	{
		(*cache_i) = -1;
	}
	if( !input.debug->GetStatus() ) debug->SetStatus(false);
}

//Destructor
BasePDF::~BasePDF()
{
	Remove_Cache();
	if( thisConfig != NULL ) delete thisConfig;
	//	If we weren't using ROOT this would a) be safe and b) be advisible, but ROOT segfaults...
	//if( seed_function != NULL ) delete seed_function;
	if( DiscreteCaches != NULL ) delete DiscreteCaches;
	if( debug_mutex != NULL && can_remove_mutex == true ) delete debug_mutex;

	if( debug != NULL ) delete debug;
}

void BasePDF::SetCopyConstructor( const IPDF* input ) const
{
	copy_object = (CopyPDF_t*) input;
}

CopyPDF_t* BasePDF::GetCopyConstructor() const
{
	return copy_object;
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
		allParameters.SetPhysicsParameters( Input );
	}
	else
	{
		allParameters.AddPhysicsParameters( Input );
	}

	//Invalidate the cache
	this->UnsetCache();

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
	//cout << "I have : " << DiscreteCaches->size() << " caches." << endl;
	this->CheckCaches( NewBoundary );

	//cout << "I now have: " << DiscreteCaches->size() << " caches." << endl;

	unsigned int thisIndex = NewBoundary->GetDiscreteIndex( InputPoint );

	//cout << "I want this: " << thisIndex << " Index" << endl;

	double thisVal = DiscreteCaches->at(thisIndex);

	//cout << "I get this: " << thisVal << endl;

	if( thisVal > 0 )
	{
		return thisVal;
	}
	else
	{
		return -1.;
	}
}

//Return the integral of the function over the given boundary
double BasePDF::Integral(DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary)
{
	double thisCache = this->GetCache( NewDataPoint, NewBoundary );

	//cout << "thisCache: " << thisCache << endl;

	// If the PDF has been configured to rely on Numerical Normalisation we can't calculate the Normalisation here at all
	// Just return what has been cached so far
	if( numericalNormalisation && cachingEnabled )
	{
		return thisCache;
	}


	//	Once we have checked which analytical form to use:
	if( haveTestedIntegral )
	{
		//cout << "Have Tested" << endl;
		if( this->CacheValid( NewDataPoint, NewBoundary ) )
		{
			//cout << this->GetLabel() << "\t:\t" << thisCache << endl;
			return thisCache;
		}
		else
		{
			if( discrete_Normalisation )
			{
				double norm = this->Normalisation(NewDataPoint, NewBoundary);
				this->SetCache( norm, NewDataPoint, NewBoundary );
				return norm;
			}
			else
			{
				double norm = this->Normalisation(NewBoundary);
				this->SetCache( norm, NewDataPoint, NewBoundary );
				return norm;
			}
		}
	}

	//cout << "Performing Test: " << endl;
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
		//cout << "This Val1: " << Norm << endl;
		//cout << "Using Analytical Integration1 For: " << this->GetName() << endl;
		return Norm;
	}

	Norm = this->Normalisation(NewDataPoint, NewBoundary);
	if( Norm >= 0. )		//	Valid
	{
		haveTestedIntegral = true;
		numericalNormalisation = false;
		discrete_Normalisation = true;
		this->SetCache( Norm, NewDataPoint, NewBoundary );
		//cout << "This Val2: " << Norm << endl;
		//cout << "Using Analytical Integration2 For: " << this->GetName() << endl;
		return Norm;
	}

	//	PDF does NOT know how to Normalise!!
	numericalNormalisation = true;
	discrete_Normalisation = true;

	//cout << "Using Numerical Integration For:" << this->GetName() << endl;
	//	We can't produce a valid PDF Normalisation
	return -1.;
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
	return this->Evaluate( NewDataPoint );
}

//Calculate the function value for numerical integration
double BasePDF::EvaluateForNumericIntegral( DataPoint * NewDataPoint )
{
	return this->Evaluate( NewDataPoint );
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

void BasePDF::SetObservableDistribution( string inputObservable, IPDF* Observable_1D_dist )
{
	if( StringProcessing::VectorContains( &doNotIntegrateList, &inputObservable ) == -1 )
	{
		doNotIntegrateList.push_back( inputObservable );
	}
	observableDistNames.push_back( inputObservable );
	Observable_1D_dist->SetRandomFunction( this->GetRandomFunction() );
	observableDistributions.push_back( Observable_1D_dist );
}

IPDF* BasePDF::GetObservableDistribution( string inputObservable )
{
	int location = StringProcessing::VectorContains( &doNotIntegrateList, &inputObservable );
	if( location == -1 ) return NULL;
	else return observableDistributions[ (unsigned)location ];
}

//  Get a pointer to the seed function
//  Using a pointer so we have one seed per normal study
TRandom3 * BasePDF::GetRandomFunction() const
{
	if( seed_function == NULL ) seed_function = new TRandom3(0);
	return seed_function;
}

//  Set the Random Generator to be some externally defined instance
void BasePDF::SetRandomFunction( TRandom3 * new_function )
{
	if( seed_function != NULL ) delete seed_function;
	seed_function = new_function;
}

//  Seed the Random Number Generator correctly
void BasePDF::SetRandomFunction( int new_seed )
{
	seed_num = new_seed;
	if( seed_function != NULL ) delete seed_function;
	seed_function = new TRandom3( (unsigned)new_seed );
}

//  Return the numerical seed
int BasePDF::GetSeedNum() const
{
	return seed_num;
}

//	Set the Status of a cache for the MC generator associated with this PDF
void BasePDF::SetMCCacheStatus( bool newStatus)
{
	if( (newStatus == false) && hasCachedMCGenerator )
	{
		//cout << GET_Name() << ":\tRemoving Cache" << endl;
		this->Remove_Cache();
	}
	hasCachedMCGenerator = newStatus;
	do_i_control_the_cache = true;
}

void BasePDF::Remove_Cache()
{
	if( do_i_control_the_cache == true )
	{
		while( !cached_files.empty() )
		{
			cached_files.back().append( ".root" );
			remove ( cached_files.back().c_str() );
			cached_files.pop_back();
		}
	}
	hasCachedMCGenerator = false;
}

void BasePDF::Can_Remove_Cache( bool input )
{
	do_i_control_the_cache = input;
}

//	Get the Status of the MC generator for this PDF
bool BasePDF::GetMCCacheStatus() const
{
	return hasCachedMCGenerator;
}

vector<string> BasePDF::GetMCCacheNames() const
{
	return cached_files;
}

void BasePDF::AddCacheObject( string obj_name )
{
	cached_files.push_back( obj_name );
}

string BasePDF::GetName() const
{
	return PDFName;
}

void BasePDF::SetName( string input )
{
	PDFName = input;
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

string BasePDF::GetLabel() const
{
	return PDFLabel;
}

void BasePDF::SetLabel( string input )
{
	PDFLabel = input;
}

string BasePDF::XML() const
{
	stringstream xml;
	xml << "<PDF>";
	xml << this->GetName();
	xml << "</PDF>" << endl;
	return xml.str();
}


void BasePDF::SetConfigurator( PDFConfigurator* config )
{
	if( thisConfig != NULL ) delete thisConfig;
	thisConfig = new PDFConfigurator( *config );
}

PDFConfigurator* BasePDF::GetConfigurator() const
{
	return thisConfig;
}

pthread_mutex_t* BasePDF::DebugMutex() const
{
	return debug_mutex;
}

void BasePDF::SetDebugMutex( pthread_mutex_t* Input, bool can_remove )
{
	can_remove_mutex = can_remove;
	if( debug_mutex != NULL && can_remove_mutex ) delete debug_mutex;
	debug_mutex = Input;
}

void BasePDF::SetDebug( DebugClass* input_debug )
{
	if( input_debug != NULL )
	{
		if( debug != NULL  ) delete debug;
		debug = new DebugClass(*input_debug);

		if( debug->DebugThisClass("BasePDF") )
		{
			debug->SetStatus(true);
			cout << "BasePDF: Debugging Enabled!" << endl;
		}
		else if( debug->DebugThisClass("PDF") )
		{
			debug->SetStatus(true);
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


