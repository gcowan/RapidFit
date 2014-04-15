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
#ifndef IPDF_H_MCCaching
#define IPDF_H_MCCaching

///	ROOT Headers
#include "TString.h"
///	RapidFit Headers
///	System Headers
#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

using namespace::std;

class IPDF_MCCaching
{
	public:
		/*!
		 * Virtual Destructor
		 */
		virtual ~IPDF_MCCaching()
		{
		};

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

	protected:
		/*!
		 * Default Constructor
		 */
		IPDF_MCCaching() {};

};

#endif

