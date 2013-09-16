/*!
 * @interface IPDF
 *
 * @brief Common interface for all PDFs
 *
 * For more details and the implementation of most of these functions see BasePDF
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef IPDF_H
#define IPDF_H

///	ROOT Headers
#include "TString.h"
#include "TRandom3.h"
///	RapidFit Headers
#include "DataPoint.h"
#include "PDFConfigurator.h"
#include "PhaseSpaceBoundary.h"
#include "ParameterSet.h"
#include "ComponentRef.h"
#include "DebugClass.h"
#include "RapidFitIntegrator.h"
#include "RapidFitIntegratorConfig.h"
///	System Headers
#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#ifdef __CINT__ 
#undef __GNUC__ 
#define _SYS__SELECT_H_
struct pthread_mutex_t;
#undef __SYS__SELECT_H_
#define __GNUC__
#endif

using namespace::std;

class IPDF;
class PDFConfigurator;
class RapidFitIntegrator;
class RapidFitIntegratorConfig;

/*!
 * @brief typedef for the class-factory objects which actually create the new class instances in memory
 */
typedef IPDF* CreatePDF_t( PDFConfigurator* );
/*!
 * @brief typedef for the class-factory objects which actually copy a new class instances into memory
 */
typedef IPDF* CopyPDF_t( const IPDF& );

class IPDF
{
	public:
		/*!
		 * Virtual Destructor
		 */
		virtual ~IPDF()
		{
		};

		/*!
		 * Interface Function:
		 * Turns all PDF caching OFF
		 */
		virtual void TurnCachingOff() = 0;

		/*!
		 * Interface Function:
		 * Is the Normalisation Cache Valid for this Discrete Combination?
		 */
		virtual bool CacheValid( DataPoint*, PhaseSpaceBoundary* ) = 0;

		/*!
		 * Interface Function:
		 * Is it valid to Cache the Normalisation?
		 */
		virtual bool GetCachingEnabled() const = 0;

		virtual void SetComponentStatus( const bool input ) = 0;

		virtual bool GetComponentStatus() const = 0;

		/*!
		 * Interface Function:
		 * Change the behaviour of Caching in this PDF.
		 */
		virtual void SetCachingEnabled( bool ) = 0;

		/*!
		 *
		 */
		virtual RapidFitIntegrator* GetPDFIntegrator() const = 0;

		virtual void SetUpIntegrator( const RapidFitIntegratorConfig* thisConfig ) = 0;

		/*!
		 * Interface Function:
		 * Does the PDF want to use Numerical Normalisation
		 */
		virtual bool GetNumericalNormalisation() const = 0;

		/*!
		 * Interface Function:
		 * Externally Set whether the PDF wants to use Numerical Normalisation
		 */
		virtual void SetNumericalNormalisation( bool ) = 0;

		/*!
		 * Interface Function:
		 * Force the Normalisation Cache relating to this DataPoint to be set
		 */
		virtual void SetCache( double, DataPoint*, PhaseSpaceBoundary* ) = 0;

		/*!
		 * Interface Function:
		 * Externally invalidate the Normalisation Cache
		 */
		virtual void UnsetCache() = 0;

		/*!
		 * Interface Function:
		 * Externally update the PDFs in the PDF
		 */
		virtual void UpdatePhysicsParameters( ParameterSet* ) = 0;

		/*!
		 * Interface Function:
		 * Return the integral of the function over the given boundary
		 */
		virtual double Integral( DataPoint*, PhaseSpaceBoundary* ) = 0;

		/*!
		 * Interface Function:
		 * Return the function value at the given point
		 */
		virtual double Evaluate( DataPoint* ) = 0;

		/*!
		 * Interface Function:
		 * Return the function value at the given point for generation
		 */
		virtual double EvaluateForNumericGeneration( DataPoint* ) = 0;

		/*!
		 * Interface Function
		 * Return the function value at the given point for use by numeric integral
		 */
		virtual double EvaluateForNumericIntegral( DataPoint* ) = 0;

		virtual double EvaluateTimeOnly( DataPoint* ) = 0;

		/*!
		 * Interface Function:
		 * Get the function parameters
		 */
		virtual ParameterSet* GetPhysicsParameters() = 0;

		/*!
		 * Interface Function:
		 * Return a prototype data point
		 */
		virtual vector<string> GetPrototypeDataPoint() = 0;

		/*!
		 * Interface Function:
		 * Return a prototype set of physics parameters
		 */
		virtual vector<string> GetPrototypeParameterSet() = 0;

		/*!
		 * Interface Function:
		 * Return a list of parameters not to be integrated
		 */
		virtual vector<string> GetDoNotIntegrateList() = 0;

		/*!
		 * Interface Function:
		 * Provide a custom distribution for a continuous observable
		 */
		virtual void SetObservableDistribution( string, IPDF* ) = 0;

		/*!
		 * Interface Function:
		 * Get a custom distribution for a requested observable (returns NULL when non initialized)
		 */
		virtual IPDF* GetObservableDistribution( string ) = 0;

		/*!
		 * Interface Function:
		 * Set the virtual cache status
		 */
		virtual void SetMCCacheStatus( bool ) = 0;

		/*!
		 * Interface Function:
		 * Get the virtual cache status
		 */
		virtual bool GetMCCacheStatus() const = 0;

		/*!
		 * Interface Function:
		 * Get the names of the MC cache files (without the .root extensions)
		 */
		virtual vector<string> GetMCCacheNames() const = 0;

		/*!
		 * Interface Function:
		 * Start a new TRandom3 instance with a given seed value
		 */
		virtual void SetRandomFunction( int ) = 0;

		/*!
		 * Interface Function:
		 * Replace TRandom3 instance with a new function
		 */
		virtual void SetRandomFunction( TRandom3* ) = 0;

		/*!
		 * Interface Function:
		 * Get the seed number used to initialise the random number gen
		 */
		virtual int GetSeedNum() const = 0;

		/*!
		 * Interface Function:
		 * Remove the cached files
		 */
		virtual void Remove_Cache() = 0;

		/*!
		 * Interface Function:
		 * Set whether this PDF can remove the MC cache object, i.e. is this instance the owner
		 */
		virtual void Can_Remove_Cache( bool ) = 0;

		/*!
		 * Interface Function:
		 * Add a virtual cache object
		 */
		virtual void AddCacheObject( string ) = 0;

		/*!
		 * Interface Function:
		 * Get the Random function stored in this PDF
		 */
		virtual TRandom3* GetRandomFunction() const = 0;

		/*!
		 * Interface Function:
		 * Get the Name of the PDF
		 */
		virtual string GetName() const = 0;

		/*!
		 * Interface Function:
		 * Return a list of PDF components addresses in string format
		 */
		virtual vector<string> PDFComponents() = 0;

		virtual string GetComponentName( ComponentRef* ) = 0;

		/*!
		 * Interface Function:
		 * Return the function value at the given point
		 */
		virtual double EvaluateComponent( DataPoint*, ComponentRef* ) = 0;

		/*!
		 * Interface Function:
		 * Get the user defined label for this PDF
		 */
		virtual string GetLabel() const = 0;

		/*!
		 * Interface Function:
		 * Set the user defined label for this PDF
		 */
		virtual void SetLabel( string ) = 0;

		/*!
		 * Interface Function:
		 * Can the PDF be safely copied through it's copy constructor?
		 */
		virtual bool IsCopyConstructorSafe() const = 0;

		/*!
		 * Interface Function:
		 * Give the PDF a pointer to the template of it's copy constructor object
		 */
		virtual void SetCopyConstructor( CopyPDF_t* ) const = 0;

		/*!
		 * Interface Function:
		 * Get the pointer to the PDF copy constructor
		 */
		virtual CopyPDF_t* GetCopyConstructor() const = 0;

		/*!
		 * Interface Function:
		 * Return the required XML for this PDF
		 */
		virtual string XML() const = 0;

		/*!
		 * Interface Function:
		 * Sets the internal copy of the the PDF configuration
		 */
		virtual void SetConfigurator( PDFConfigurator* ) = 0;

		/*!
		 * Interface Function:
		 * Returns a pointer to the internal Object which contains the PDF configuration
		 */
		virtual PDFConfigurator* GetConfigurator() const = 0;

		virtual pthread_mutex_t* DebugMutex() const = 0;

		virtual void SetDebugMutex( pthread_mutex_t* Input, bool =true ) = 0;

		virtual void SetDebug( DebugClass* input_debug ) = 0;

		virtual void Print() const = 0;

		/*!
		 * Interface Function:
		 * Get the Name of the PDF
		 */
		virtual void SetName( string ) = 0;

		virtual void ChangePhaseSpace( PhaseSpaceBoundary * InputBoundary ) = 0;

	protected:
		/*!
		 * Default Constructor
		 */
		IPDF() {};

		/*!
		 * Interface Function:
		 * Set if the PDF be safely copied through it's copy constructor?
		 */
		virtual void SetCopyConstructorSafe( bool = true ) = 0;

		virtual bool SetPhysicsParameters( ParameterSet* Input ) = 0;

		/*!
		 * Interface Function:
		 * Return the Integral over the whole PhaseSpace
		 */
		virtual double Normalisation( PhaseSpaceBoundary* ) = 0;

		/*!
		 * Interface Function:
		 * Return the Integral of this DataPoint in this PhaseSpace
		 */
		virtual double Normalisation( DataPoint*, PhaseSpaceBoundary* ) = 0;
};

#ifndef __CINT__

/*!
 * @brief Macro for adding a class lookup and copy function instance to the main function index
 *
 * This is VERY easy to maintain from a PDF developers point of view as it avoids
 * having to re-write these 2 macros for each PDF by hand
 */

#define SIMPLE_PDF_CREATOR()\
	extern "C" IPDF* CreatePDF_##X() {\
		return (IPDF*) new X();\
	}\
extern "C" IPDF* CopyPDF_##X( const IPDF& input ) { \
	return (IPDF*) new X( (X&) input ); \
}

/*!
 * @brief Macro for adding a class lookup and copy function instance to the main function index
 *
 * This allows for any standard Configuration functions inherited from each PDF that must be called after the object has been initialized to be called here
 *
 * This allows for objects to be correctly configured without the PDF developer having to care about initializing the PDFs
 *
 * @return This Returns the PDF that has been constructed
 */
#define PDF_CREATOR( X ) \
	extern "C" IPDF* CreatePDF_##X( PDFConfigurator* config ) { \
		IPDF* thisObject = (IPDF*) new X( config );\
		thisObject->SetConfigurator( config );\
		thisObject->SetName( #X );\
		thisObject->SetLabel( #X );\
		return thisObject;\
	} \
extern "C" IPDF* CopyPDF_##X( const IPDF& input ) { \
	IPDF* returnable = (IPDF*) new X( (X&) input );\
	returnable->SetName( #X );\
	return returnable;\
}


/*!
 * @brief some common thread locking commands
 */
#define PDF_THREAD_LOCK\
	bool hasOwnership=false;\
if( this->DebugMutex() != NULL )\
{\
	hasOwnership=true;\
	pthread_mutex_lock( this->DebugMutex() );\
}\


#define PDF_THREAD_UNLOCK\
	if( this->DebugMutex() != NULL && hasOwnership ) pthread_mutex_unlock( this->DebugMutex() );

#endif

#endif

