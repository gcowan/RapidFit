/*!
 * @ingroup Configurators This Generator is capable of Constructing the FitFunction objects which are to be Minimised
 * @class FitFunctionConfiguration
 *
 * @brief Container that stores all information related to FitFunction configuration, and returns an appropriate instance of a FitFunction
 *
 * @author Benjamin M Wynne bwynne@cern.ch
*/

#pragma once
#ifndef FIT_FUNCTION_CONFIGURATION_H
#define FIT_FUNCTION_CONFIGURATION_H

///	ROOT Headers
#include "TString.h"
///	RapidFit Headers
#include "FitFunction.h"
#include "PhysicsBottle.h"
///	System Headers
#include <string>

using namespace::std;

class FitFunctionConfiguration
{
	public:
		/*!
		 * Default Constructor
		 */
		FitFunctionConfiguration();

		/*!
		 * Constructor With Function Name
		 */
		FitFunctionConfiguration( string );

		/*!
		 * Constructor With Function Name and Weight Name
		 */
		FitFunctionConfiguration( string, string );

		/*!
		 *	Destructor
		 */
		~FitFunctionConfiguration();

		/*!
		 * Get a pointer to a newly constructed FitFunction
		 */
		FitFunction * GetFitFunction();
	
		/*!
		 * Were Weights Used?
		 */
		bool GetWeightsWereUsed();

		/*!
		 * What was the Weight Observable Name
		 */
		string GetWeightName();

		/*!
		 * Setup a Trace and write the output to the given filename
		 */
		void SetupTrace( TString );

		/*!
		 * Get the Name of the Strategy Used
		 */
		string GetStrategy();

		/*!
		 * Set the Name of the Strategy to Use
		 */
		void SetStrategy( string );

		/*!
		 * Set the Number of Threads to Construct the FitFunction with
		 */
		void SetThreads( int );

		/*!
		 * Set wether The FitFunction should test the Integrator
		 */
		void SetIntegratorTest( bool );

		/*!
		 * Prints some useful Print Configuration
		 */
		void Print() const;

		/*!
		 * Return the XML configuration required to reproduce this class in string format
		 */
		string XML() const;

		void SetNormaliseWeights( bool Input );

		bool GetNormaliseWeights() const;
	private:

		string functionName, weightName;/*!	Name of the Function and Weight to use		*/
		bool hasWeight;			/*!	Should a Weight be used by the FitFunction	*/
		bool wantTrace;			/*!	Should a Trace Be Performed			*/
		TString TraceFileName;		/*!	Trace Output FileName				*/
		int traceCount;			/*!	Trace Number, to avoid overwriting the output file	*/
		int Threads;			/*!	Number of Threads to Construct FitFunction With	*/
		string Strategy;		/*!	Name of Strategy to use				*/
		bool testIntegrator;		/*!	Decision about wether to test the integrator	*/
		bool NormaliseWeights;
};

#endif

