/**
  @class BasePDF

  Class that provides a general implementation of IPDF.
  Can inherit from this to make a PDF without worrying about the details.

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */


#ifndef BASE_PDF_H
#define BASE_PDF_H

//	RapidFit Headers
#include "IPDF.h"
#include "ObservableRef.h"
#include "PDFConfigurator.h"
//	System Headers
#include <vector>
#include <cmath>

class BasePDF : public IPDF
{
	public:
		BasePDF();
		virtual ~BasePDF();

		//Indicate whether the function has been set up correctly
		virtual bool IsValid();

		//Set the function parameters
		virtual bool SetPhysicsParameters( ParameterSet* );

		//Return the integral of the function over the given boundary
		virtual double Integral( DataPoint*, PhaseSpaceBoundary* );

		//Return the function value at the given point
		virtual double Evaluate( DataPoint* );

		//Return the function value at the given point for use in numeric integral
		virtual double EvaluateForNumericIntegral( DataPoint* );

		//Return the components of the function value at the given point
		virtual vector<double> EvaluateComponents( DataPoint* );

		//Return a prototype data point
		virtual vector<string> GetPrototypeDataPoint();

		//Return a prototype set of physics parameters
		virtual vector<string> GetPrototypeParameterSet();

		virtual ParameterSet* GetActualParameterSet();

		//Return a list of parameters not to be integrated
		virtual vector<string> GetDoNotIntegrateList();

		//Update integral cache
		virtual void UpdateIntegralCache();

		//	Set the virtual ID
		virtual void SET_ID( string );
		virtual void SET_ID( TString );

		//	Get the virtual ID
		virtual string GET_ID();

		//	Set the virtual cache status
		virtual void SetMCCacheStatus( bool );

		//	Get the virtual cache status
		virtual bool GetMCCacheStatus();

		//	Start a new TRandom3 instance with a given seed value
		virtual void SetRandomFunction( int );

		//	Replace TRandom3 instance with a new function
		virtual void SetRandomFunction( TRandom3* );

		virtual int GetSeedNum();

		//	Remove the cached files
		virtual void Remove_Cache( bool=false );

		//	Add a virtual cache object
		virtual void AddCacheObject( string );
		virtual void AddCacheObject( TString );

		//	Get the Random function stored in this PDF
		virtual TRandom3* GetRandomFunction();

		virtual string GetName();
		virtual void SetName( string );
	protected:
		//Do the evaluation
		//virtual double Value(DataPoint*);

		//Do the integration
		virtual double Normalisation( PhaseSpaceBoundary* );
		virtual double Normalisation( DataPoint*, PhaseSpaceBoundary* );

		double cachedIntegral;
		bool cacheValid;
		ParameterSet allParameters;
		vector<string> allObservables;
		bool valid;
		ObservableRef observables;

	private:

		//      Varibles required for caching MC
		vector<string> cached_files;
		string stored_ID;
		bool hasCachedMCGenerator;

		//      Variables required for storing the Seed
		vector<TRandom3 *> seed_function;
		vector<int> seed_num;

		string PDFName;
};

#endif

