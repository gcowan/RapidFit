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
#ifndef IPDF_H_FRAMEWORK
#define IPDF_H_FRAMEWORK

///	ROOT Headers
#include "TString.h"
#include "TRandom3.h"
///	RapidFit Headers
#include "PhaseSpaceBoundary.h"
#include "DebugClass.h"
#include "PDFConfigurator.h"
#include "RapidFitIntegrator.h"
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

/*!
 *  * @brief typedef for the class-factory objects which actually create the new class instances in memory
 *   */
typedef IPDF* CreatePDF_t( PDFConfigurator* );
/*!
 *  * @brief typedef for the class-factory objects which actually copy a new class instances into memory
 *   */
typedef IPDF* CopyPDF_t( const IPDF& );


class IPDF_Framework
{
	public:
		/*!
		 * Virtual Destructor
		 */
		virtual ~IPDF_Framework()
		{
		};

		/*!
		 *
		 */
		virtual RapidFitIntegrator* GetPDFIntegrator() const = 0;

		virtual void SetUpIntegrator( const RapidFitIntegratorConfig* thisConfig ) = 0;

		/*!
		 * Interface Function:
		 * Get the Name of the PDF
		 */
		virtual string GetName() const = 0;

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
		IPDF_Framework() {};

};

#ifndef __CINT__

/*!
 *  * @brief Macro for adding a class lookup and copy function instance to the main function index
 *   *
 *    * This is VERY easy to maintain from a PDF developers point of view as it avoids
 *     * having to re-write these 2 macros for each PDF by hand
 *      */

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
 *  * @brief some common thread locking commands
 *   */
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

