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
#ifndef IPDF_H_NORMALISATION_CACHING
#define IPDF_H_NORMALISATION_CACHING

///	ROOT Headers
#include "TString.h"
#include "TRandom3.h"
///	RapidFit Headers
#include "DataPoint.h"
#include "PhaseSpaceBoundary.h"
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

class IPDF_NormalisationCaching
{
	public:

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

	protected:
		/*!
		 * Default Constructor
		 */
		IPDF_NormalisationCaching() {};

};

#endif

