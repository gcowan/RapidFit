/**
        @interface IPDF

        Common interface for all PDFs

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef IPDF_H
#define IPDF_H

//	ROOT Headers
#include "TString.h"
#include "TRandom3.h"
//	RapidFit Headers
#include "DataPoint.h"
#include "PDFConfigurator.h"
#include "PhaseSpaceBoundary.h"
#include "ParameterSet.h"
//	System Headers
#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>

using namespace std;

class IPDF
{
	public:
		IPDF() {};
		virtual ~IPDF() {};

		//Indicate whether the function has been set up correctly
		virtual bool IsValid() = 0;

		//Set the function parameters
		virtual bool SetPhysicsParameters( ParameterSet* ) = 0;

		//Return the integral of the function over the given boundary
		virtual double Integral( DataPoint*, PhaseSpaceBoundary* ) = 0;

		//Return the function value at the given point
		virtual double Evaluate( DataPoint* ) = 0;

		//Return the function value at the given point for use by numeric integral
		virtual double EvaluateForNumericIntegral( DataPoint* ) = 0;
	
		//Return the components of the function value at the given point
		virtual vector<double> EvaluateComponents( DataPoint* ) = 0;

		//Return a prototype data point
		virtual vector<string> GetPrototypeDataPoint() = 0;

		//Return a prototype set of physics parameters
		virtual vector<string> GetPrototypeParameterSet() = 0;

		virtual ParameterSet* GetActualParameterSet() = 0;

		//Return a list of parameters not to be integrated
		virtual vector<string> GetDoNotIntegrateList() = 0;

		//Update the integral cache
		virtual void UpdateIntegralCache() = 0;

		virtual void SET_ID( string ) = 0;
		virtual void SET_ID( TString ) = 0;

		//	Get the virtual ID
		virtual string GET_ID() = 0;

		//	Set the virtual cache status
		virtual void SetMCCacheStatus( bool ) = 0;

		//	Get the virtual cache status
		virtual bool GetMCCacheStatus() = 0;

		//	Start a new TRandom3 instance with a given seed value
		virtual void SetRandomFunction( int ) = 0;

		//	Replace TRandom3 instance with a new function
		virtual void SetRandomFunction( TRandom3* ) = 0;

		virtual int GetSeedNum() = 0;

		//	Remove the cached files
		virtual void Remove_Cache( bool=false ) = 0;

		//	Add a virtual cache object
		virtual void AddCacheObject( string ) = 0;
		virtual void AddCacheObject( TString ) = 0;

		//	Get the Random function stored in this PDF
		virtual TRandom3* GetRandomFunction() = 0;

		virtual string GetName() = 0;
		virtual void SetName( string ) = 0;
};


//	Macro for adding a class lookup and copy function instance to the main function index
//
//	This is VERY easy to maintain from a PDF developers point of view as it avoids
//	having to re-write these 2 macros for each PDF by hand
#define PDF_CREATOR( X ) \
	extern "C" IPDF* CreatePDF_##X( PDFConfigurator* config ) { \
		return (IPDF*) new X( config ); \
	} \
	extern "C" IPDF* CopyPDF_##X( IPDF& input ) { \
		return (IPDF*) new X( (X&) input ); \
	}

//	typedef for the class-factory objects which actually create the new class instances in memory
typedef IPDF* CreatePDF_t( PDFConfigurator* );
typedef IPDF* CopyPDF_t( IPDF& );

#endif
