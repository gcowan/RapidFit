/**
  @class SumPDF

  An implementation of IPDF for adding the values of two other IPDFs

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

#pragma once
#ifndef SUM_PDF_H
#define SUM_PDF_H

//	RapidFit Headers
#include "IPDF.h"
#include "BasePDF.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class SumPDF : public BasePDF
{
	public:
		SumPDF( IPDF*, IPDF*, PhaseSpaceBoundary*, string );
		SumPDF( const SumPDF& );
		~SumPDF();

		virtual void SetUseGSLIntegrator( bool input );

		void TurnCachingOff();

		//Set the function parameters
		bool SetPhysicsParameters( ParameterSet* );

		//Return the integral of the function over the given boundary
		double Normalisation( DataPoint*, PhaseSpaceBoundary* );

		//Return the function value at the given point
		double Evaluate( DataPoint* );

		//Return a prototype data point
		virtual vector<string> GetPrototypeDataPoint();

		//Return a prototype set of physics parameters
		virtual vector<string> GetPrototypeParameterSet();

		//Return a list of parameters not to be integrated
		virtual vector<string> GetDoNotIntegrateList();

		//      Return components, component 0 by default
		double EvaluateComponent( DataPoint*, ComponentRef* );

		bool GetNumericalNormalisation() const;

		IPDF* GetFirstPDF() const;

		IPDF* GetSecondPDF() const;

		void SetCachingEnabled( bool );

		bool GetCachingEnabled() const;

		string XML() const;

		void SetDebugMutex( pthread_mutex_t* Input, bool =true );

		void SetDebug( DebugClass* input_debug );

		virtual string GetComponentName( ComponentRef* = NULL );
	private:
		//	Uncopyable!
		//SumPDF ( const SumPDF& );
		SumPDF& operator = ( const SumPDF& );
		void MakePrototypes( PhaseSpaceBoundary* );

		void TurnThisCachingOff();

		vector<string> prototypeDataPoint, prototypeParameterSet, doNotIntegrateList;
		IPDF * firstPDF;
		IPDF * secondPDF;
		double firstFraction, firstIntegralCorrection, secondIntegralCorrection;
		string fractionName;
};

#endif

