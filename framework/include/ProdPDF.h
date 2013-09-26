/*!
 * @class ProdPDF
 *
 * @brief An implementation of IPDF for multiplying the values of two other IPDFs
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie
 */

#pragma once
#ifndef PROD_PDF_H
#define PROD_PDF_H

//	RapidFit Headers
#include "IPDF.h"
#include "BasePDF.h"
#include "ComponentRef.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class ProdPDF : public BasePDF
{
	public:
		/*!
		 * @brief Correct Constructor which takes in 2 arguemnt PDF objects
		 *
		 * @param First   First fully constructed PDF object
		 *
		 * @param Second  Second fully constructed PDF object
		 */
		ProdPDF( PDFConfigurator* config );
		//ProdPDF( IPDF* First, IPDF* Second );

		/*!
		 * @brief Copy Constructor
		 */
		ProdPDF( const ProdPDF& );

		/*!
		 * @brief Destructor
		 */
		~ProdPDF();

		void SetComponentStatus( const bool input );

		bool GetComponentStatus() const;

		virtual void SetUpIntegrator( const RapidFitIntegratorConfig* thisConfig );

		vector<string> PDFComponents();

		void TurnCachingOff();

		/*!
		 * @brief Set the function parameters
		 *
		 * @param Input
		 */
		bool SetPhysicsParameters( ParameterSet* Input );

		/*!
		 * @brief Return the integral of the function over the given boundary
		 */
		double Normalisation( DataPoint*, PhaseSpaceBoundary* );

		//Return the function value at the given point
		double Evaluate( DataPoint* );
		double EvaluateForNumericIntegral( DataPoint* );

		//Return a prototype data point
		vector<string> GetPrototypeDataPoint();

		//Return a prototype set of physics parameters
		vector<string> GetPrototypeParameterSet();

		//Return a list of parameters not to be integrated
		vector<string> GetDoNotIntegrateList();

		bool GetNumericalNormalisation() const;

		//      Return components, component 0 by default
		double EvaluateComponent( DataPoint*, ComponentRef* );

		IPDF* GetFirstPDF() const;

		IPDF* GetSecondPDF() const;

		bool GetCachingEnabled() const;

		void SetCachingEnabled( bool );

		/*!
		 * @brief Provide the XML which is capable of constructing this PDF in it's current state
		 */
		string XML() const;

		void SetDebugMutex( pthread_mutex_t* Input, bool =true );

		void SetDebug( DebugClass* input_debug );

		virtual string GetComponentName( ComponentRef* = NULL );
	private:
		//	Uncopyable!
		//ProdPDF ( const ProdPDF& );
		ProdPDF& operator = ( const ProdPDF& );

		/*!
		 * @brief Called by the Constructor to deal with initializing objects now that the PDF has been initalized
		 */
		void MakePrototypes();

		void TurnThisCachingOff();

		void MakePrototypes( PhaseSpaceBoundary* );
		vector<string> prototypeDataPoint;
		vector<string> prototypeParameterSet;
		vector<string> doNotIntegrateList;
		IPDF * firstPDF;
		IPDF * secondPDF;
};

#endif

