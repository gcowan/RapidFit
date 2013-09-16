/*!
 * @class IFitFunction
 *
 * @brief Interface class for the function to minimise
 * 
 * Overload the evaluate methods and UP value for Chi2, NLL, etc.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
*/

#pragma once
#ifndef _I_FIT_FUNCTION_H
#define _I_FIT_FUNCTION_H

///	ROOT Headers
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
///	RapidFit Headers
#include "PhysicsBottle.h"
#include "RapidFitIntegrator.h"
#include "RapidFitIntegratorConfig.h"
#include "ObservableRef.h"
#include "DebugClass.h"

#include <vector>
#include <string>

using namespace::std;

class IFitFunction
{
	public:
		/*!
		 * @brief Default Detructor
		 */
		virtual ~IFitFunction() {};

		/*!
		 * @brief Setup the Trace to record all of the output from the Minimiser
		 */
		virtual void SetupTrace( const TString FileName, const int traceNum ) = 0;

		/*!
		 * @brief Set the name of the numerical integration method to use 
		 */
		virtual void SetIntegratorConfig( const RapidFitIntegratorConfig* gsl ) = 0;

		/*!
		 * @brief Set the Physics Bottle to be used
		 *
		 * @param Input This is the PhysicsBottle that the IFitFunction will use to calculate the result in Evalute
		 *
		 * @return Void
		 */
		virtual void SetPhysicsBottle( const PhysicsBottle* Input ) = 0;

		/*!
		 * @brief Get the Pointer to the Internal Physics Bottle
		 *
		 * @return Returns a pointer to the PhysicsBottle the IFitFunction is using to Evaluate
		 */
		virtual PhysicsBottle * GetPhysicsBottle() const = 0;

		/*!
		 * @brief Change/Update the ParameterSet
		 *
		 * @param Input  This is the ParameterSet that the internal ParameterSet(s) should be chaned to
		 *
		 * @return true Changed the ParameterSet with no problems, false there was an error (I don't think this can be trusted and this should be made void)
		 */
		virtual void SetParameterSet( const ParameterSet* ) = 0;

		/*!
		 * @brief Get a Pointer to the Internal Parameter Set
		 *
		 * @return Returns a pointer to the Internal ParameterSet
		 */
		virtual ParameterSet * GetParameterSet() const = 0;

		/*!
		 * @brief Evaluate the IFitFunction
		 *
		 * @return Returns a final value of the whole ParameterSet as a single double
		 */
		virtual double Evaluate() = 0;

		/*!
		 * @brief Set the Name of the Weights to use and the fact that Weights were used in the fit
		 *
		 * @param Name    This sets the name of the Weights to be used when Evaluating the DataSet
		 *
		 * @return Void
		 */
		virtual void UseEventWeights( const string Name ) = 0;

		/*!
		 * @brief Return if Weights were used as part of the analysis
		 *
		 * @return bool  true = Weights were used, false = Weights were NOT used
		 */
		virtual bool GetWeightsWereUsed() const = 0;

		virtual string GetWeightName() const = 0;
		
		/*!
		 * @brief Set the IFitFunction to use Weights squared
		 *
		 * @param Input   true this causes Weights squared to be used in the fit, false don't use Weights squared
		 *
		 * @return Void
		 */
		virtual void SetUseWeightsSquared( const bool Input ) = 0;

		/*!
		 * @brief Set the Number of threads to be used by the IFitFunction
		 *
		 * @param Input   This sets the number of threads that should be used by this IFitFunction during Evaluate
		 *
		 * @return Void
		 */
		virtual void SetThreads( const int Input ) = 0;

		/*!
		 * @brief Get the Number of threads this IFitFunction is using
		 *
		 * @return Returns the number of Threads this IFitFunction is attempting to use
		 */
		virtual int GetThreads() const = 0;

		/*!
		 * @brief Set whether any RapidFitIntegrator Objects created internally should check the PDF/Numerical Integral
		 *
		 * @param Input  Should The Integrals be tested?  true = yes  false = no
		 *
		 * @return Void
		 */
		virtual void SetIntegratorTest( const bool Input ) = 0;

		/*!
		 * @brief What size of step in the function defines the error, 0.5 for NLL, 1 for chi2... etc.
		 *
		 * @input n  This is the nsigma Error you wish to calculate
		 *
		 * @return Returns the equivalent rise in the Function value that should be used to calculate the Error Value
		 */
		virtual double UpErrorValue( const int n ) = 0;

		virtual void SetDebug( DebugClass* debug ) = 0;

		virtual unsigned int GetCallNum() = 0;

	protected:
		IFitFunction() {};

		/*!
		 * Don't Copy the class this way!
		 */
		IFitFunction ( const IFitFunction& );

		/*!
		 * Don't Copy the class this way!
		 */
		IFitFunction& operator = ( const IFitFunction& );

};

#endif

