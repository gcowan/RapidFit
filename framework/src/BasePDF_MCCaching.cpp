/*
 * @class BasePDF
 *
 * Class that provides a general implementation of IPDF.
 * Can inherit from this to make a PDF without worrying about the details.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */
///	RapidFit Headers
#include "IPDF_MCCaching.h"
#include "BasePDF_MCCaching.h"
/*#include "StringProcessing.h"
#include "ObservableRef.h"
#include "PhaseSpaceBoundary.h"
#include "RapidFitIntegrator.h"*/
///	System Headers
#include <iostream>
#include <cmath>
#include <sstream>
#include <pthread.h>

using namespace::std;

//Constructor
BasePDF_MCCaching::BasePDF_MCCaching() : IPDF_MCCaching(),
	cached_files(), hasCachedMCGenerator(false), do_i_control_the_cache(false)
{
}

BasePDF_MCCaching::BasePDF_MCCaching( const BasePDF_MCCaching& input ) : IPDF_MCCaching(),
	cached_files(input.cached_files), hasCachedMCGenerator(input.hasCachedMCGenerator), do_i_control_the_cache(false)
{
}

//Destructor
BasePDF_MCCaching::~BasePDF_MCCaching()
{
	this->Remove_Cache();
}

//	Set the Status of a cache for the MC generator associated with this PDF
void BasePDF_MCCaching::SetMCCacheStatus( bool newStatus)
{
	if( (newStatus == false) && hasCachedMCGenerator )
	{
		//cout << GET_Name() << ":\tRemoving Cache" << endl;
		this->Remove_Cache();
	}
	hasCachedMCGenerator = newStatus;
	do_i_control_the_cache = true;
}

void BasePDF_MCCaching::Remove_Cache()
{
	if( do_i_control_the_cache == true )
	{
		while( !cached_files.empty() )
		{
			cached_files.back().append( ".root" );
			remove( cached_files.back().c_str() );
			cached_files.pop_back();
		}
	}
	hasCachedMCGenerator = false;
}

void BasePDF_MCCaching::Can_Remove_Cache( bool input )
{
	do_i_control_the_cache = input;
}

//	Get the Status of the MC generator for this PDF
bool BasePDF_MCCaching::GetMCCacheStatus() const
{
	return hasCachedMCGenerator;
}

vector<string> BasePDF_MCCaching::GetMCCacheNames() const
{
	return cached_files;
}

void BasePDF_MCCaching::AddCacheObject( string obj_name )
{
	cached_files.push_back( obj_name );
}

