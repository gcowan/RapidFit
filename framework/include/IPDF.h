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
		IPDF();
		virtual ~IPDF();

		//Indicate whether the function has been set up correctly
		virtual bool IsValid() = 0;

		//Set the function parameters
		virtual bool SetPhysicsParameters( ParameterSet* ) = 0;

		//Return the integral of the function over the given boundary
		virtual double Integral( DataPoint*, PhaseSpaceBoundary* ) = 0;

		//Return the function value at the given point
		virtual double Evaluate( DataPoint* ) = 0;

		//Return the components of the function value at the given point
		virtual vector<double> EvaluateComponents( DataPoint* ) = 0;
	
		//Return a prototype data point
		virtual vector<string> GetPrototypeDataPoint() = 0;

		//Return a prototype set of physics parameters
		virtual vector<string> GetPrototypeParameterSet() = 0;

		//Return a list of parameters not to be integrated
		virtual vector<string> GetDoNotIntegrateList() = 0;

		//Update the integral cache
		virtual void UpdateIntegralCache() = 0;


		//	These should __NOT__ be virtual due to the inheritance structure
		//	These can be virtual but then this requires re-implementing
		//	the same code within Normalised, Sum, Product... PDFs

		//	Set the virtual ID
		void SET_ID( string );
		void SET_ID( TString );

		//	Get the virtual ID
		string GET_ID();

		//	Set the virtual cache status
		void SetMCCacheStatus( bool );

		//	Get the virtual cache status
		bool GetMCCacheStatus();

		//	Start a new TRandom3 instance with a given seed value
		void SetRandomFunction( int );

		//	Replace TRandom3 instance with a new function
		void SetRandomFunction( TRandom3* );

		int GetSeedNum();

		//	Remove the cached files
		void Remove_Cache();

		//	Add a virtual cache object
		void AddCacheObject( string );
		void AddCacheObject( TString );

		//	Get the Random function stored in this PDF
		TRandom3* GetRandomFunction();

	private:
		//	Uncopyable!
		IPDF ( const IPDF& );
		IPDF& operator = ( const IPDF& );
  
		//	Varibles required for caching MC
		vector<string> cached_files;
		string stored_ID;
		bool hasCachedMCGenerator;

		//	Variables required for storing the Seed
		vector<TRandom3 *> seed_function;
		vector<int> seed_num;

};

#endif
