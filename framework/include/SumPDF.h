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

		/*!
		 * @brief Constructor
		 *
		 */
		//SumPDF( IPDF*, IPDF*, PhaseSpaceBoundary*, string );
		SumPDF( PDFConfigurator* config );

		/*!
		 * @brief Copy Constructor
		 */
		SumPDF( const SumPDF& );

		/*!
		 * @brief Destructor Function
		 */
		~SumPDF();

		void SetComponentStatus( const bool input );

		bool GetComponentStatus() const;

		vector<string> PDFComponents();

		virtual void SetUpIntegrator( const RapidFitIntegratorConfig* thisConfig );

		/*!
		 * @brief Explicitly request caching to be turned OFF
		 *
		 * @return Void
		 */
		void TurnCachingOff();

		/*!
		 * @brief   Interface Function:  This function is called ONCE per call from Minuit
		 *
		 *
		 * Using this, complex variables from the ParameterSet can be calculated/cached.
		 *
		 * The derived PDF MUST, either:   1) Use the Internal ParameterSet allParameters
		 *                                 2) Overload this function but set the internal ParameterSet as appropriate
		 *
		 * @param Input    ParameterSet Normally as defined by the IMinimiser
		 *
		 * @return        true = success , false = fail... (this can probably be replaced with a Void return)
		 */
		bool SetPhysicsParameters( ParameterSet* );

		/*!
		 * @brief Protected Function for each PDF which provides a method for the PDF to analytically integrate over the whole phase space
		 *
		 * @param Input This is the PhaseSpace to be Integrated
		 *
		 * @return Returns the value of the Normalisation for the whole PhaseSpace
		 */
		double Normalisation( DataPoint*, PhaseSpaceBoundary* );

		/*!
		 * @brief   Interface Function:  Return the function value at the given point
		 *
		 * This MUST BE OVERLOADED BY derived PDF
		 *
		 * @param Input   DataPoint that should be Evaluated
		 *
		 * @return        (+ve) Returns the Likelihood of this DataPoint for this PDF (not required to be normalised)
		 */
		double Evaluate( DataPoint* );

		/*!
		 * @brief Interface Function: Return a prototype data point
		 *
		 * This returns the list of Observables which this PDF requires which can be passed to a DataPoint Constructor
		 *
		 * The derived PDF MUST, either:   1)   Put the names of all required Observables in member object allObservables
		 *                                 2)   Overload This Function to provide a vector of all required Observable Names
		 *
		 * @return        Vector of individual Names of the Observables required by the PDF
		 */
		virtual vector<string> GetPrototypeDataPoint();

		/*!
		 * @brief Interface Function: Return a prototype set of physics parameters
		 *
		 * Returns a list of the Physics Parameters this PDF requires
		 *
		 * The derived PDF MUST, either:   1)   Use allParameters for the internal ParameterSet (wraps around to ParameterSet::AllNames() )
		 *                                 2)   Overload this function to provide the list of all PhysicsParameters it requires
		 *
		 * @return        Vector of individual PhysicsParameter Names which can be passed to a ParameterSet constructor
		 */
		virtual vector<string> GetPrototypeParameterSet();

		/*!
		 * @brief Interface Function: Return a list of parameters not to be integrated by this PDF
		 *
		 * Returns a list of the Observables the PDF cannot correctly model
		 *
		 * The derived PDF MUST, either:   1)   Put the name of Observables it requires, but can't correctly integrate over in the member object doNotIntegrateList
		 *                                 2)   Overload This function to provide the list of Observables it can't integrate over
		 *
		 * @return        Vector of Observable Names that the PDF requires but doesn't fully model, e.g. mistag/time resolution (most detector influenced parameters)
		 */
		virtual vector<string> GetDoNotIntegrateList();

		/*!
		 * @brief Interface Function: Return the value of the component given by ComponentRef, component 0 by default
		 *
		 * When no ComponentRef is provided this wraps around to the Evaluate method
		 *
		 * When the Name of the ComponentRef the PDF should return the same value as the Evaluate Method.
		 * (This is not tested, but a 0'th component is ALWAYS added to PDFComponents when missing for design reasons)
		 *
		 * The derived PDF MUST, either:   1)   Do nothing and not expect to provide components
		 *                                 2)   Provide some method to evaluate the PDF for the given DataPoint(Unique Discrete Combination) for the requested CompoentRef
		 *
		 * @param InputDataPoint   This is the DataPoint (a single Unique Discrete Combination from the PhaseSpaceBoundary) being evaluated
		 * @param InputRef         This is the ComponentRef object being interrogated  (Optional)
		 *                         I make use of this and not a simple string as at a higher level I want to avoid a compound PDF performing many lookups per single call (of potentially hundreds)
		 *
		 * @return        Value of this single component evaluated by the PDF
		 */
		double EvaluateComponent( DataPoint*, ComponentRef* );

		/*!
		 * @brief   Interface Function: Does the PDF want to use Numerical Normalisation
		 *
		 * @return        true = use Numerical normalisation , false = use Analytical normalisation
		 */
		bool GetNumericalNormalisation() const;

		/*!
		 * @brief return a pointer to the first internal PDF object
		 */
		IPDF* GetFirstPDF() const;

		/*!
		 * @brief return a pointer to the second internal PDF object
		 */
		IPDF* GetSecondPDF() const;

		/*!
		 * @brief   Interface Function: Change the behaviour of the Caching in this PDF
		 *
		 * Only specialized PDFs should overload this function!
		 *
		 * @param Input   true = Enable Caching,  false Disable Caching
		 *
		 * @return Void
		 */
		void SetCachingEnabled( bool );

		/*!
		 * @brief   Interface Function: Is it valid to Cache the Numerical Integral
		 *
		 * Only specialized PDFs should overload this function!
		 *
		 * @return        true = Caching has been Enabled, false = Caching has NOT been Enabled
		 */
		bool GetCachingEnabled() const;

		/*!
		 * @brief Return the required XML for this PDF
		 *
		 * @return Returns the name of the PDF in PDF tags
		 */
		string XML() const;

		void SetDebugMutex( pthread_mutex_t* Input, bool =true );

		void SetDebug( DebugClass* input_debug );

		/*!
		 * @brief Returns the name of this Component
		 *
		 * This allows the PDF to have a Component title which doesn't match the lyteral PDF Name in the code
		 *
		 * @returns a string Name
		 */
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
		PhaseSpaceBoundary * integrationBoundary;
};

#endif

