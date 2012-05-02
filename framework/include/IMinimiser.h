/*!
 * @interface IMinimiser
 *
 * @brief Interface for all function minimisers
 *
 * @author Benjamin M Wynne bwynne@cern.ch
*/

#pragma once
#ifndef I_MINIMISER_H
#define I_MINIMISER_H

///	RapidFit Headers
#include "FitFunction.h"
#include "FitResult.h"
#include "PhysicsParameter.h"
///	System Headers
#include <string>
#include <vector>

using namespace::std;

class IMinimiser
{
	public:
		/*!
		 * @brief Virtual Destructor
		 */
		virtual ~IMinimiser(){};

		/*!
		 * @brief Interface Function:
		 *        Initialize everything in the Minimiser before the actual call to Minimise the FitFunction is called
		 *
		 * This makes up for the fact you can't template Constructors in c++ so can't guarantee that all of the Minimisers must accept a FitFunction at construction regardless of how they use it then
		 *
		 * @return Void
		 */
		virtual void SetupFit( FitFunction* ) = 0;

		/*!
		 * @brief Interface Function:
		 *        Fix the given Parameters to the given Values
		 *
		 * @param Values  These are the Values that you wish to fix the values to
		 *
		 * @pram Names    These are the Names of the parameters you wish to fix
		 *
		 * @return Void
		 */
		virtual void FixParameters( vector<double> Values, vector<string> Names ) = 0;

		/*!
		 * @brief Interface Function:
		 *        Actually Find the Minima of the FitFunction
		 *
		 * @return Void
		 */
		virtual void Minimise() = 0;

		/*!
		 * @brief Interface Function:
		 *        Get the FitResult from the Minimiser after it has minimsed a FitFunction
		 *
		 * @return Returns a single FitResult instance
		 */
		virtual FitResult * GetFitResult() = 0;

		/*!
		 * @brief Interface Function:
		 *        Get the Contours which have been generated using any internal function within the Minimiser
		 *
		 * The behaviour of Most Minimisers is to move out until it has discovered the contour and to move around the contour itself to draw the objects
		 *
		 * This DOES NOT always produce the same result as analysing the whole PhaseSpace in a grid. You were warned HERE!
		 *
		 * @param Input a vector of pairs containing the pairs of Names defining the 2D planes which the Minimiser should project a Contour over
		 */
		virtual void ContourPlots( vector< pair< string, string > > Input ) = 0;

		/*!
		 * @brief Interface Function:
		 *        Set the OutputLevel, 0 is normally standard larger is more verbose and -ve is normally quiet
		 *
		 * @param Input  This is the verbosity Level you wish to set the Minimiser to take
		 *
		 * @return Void
		 */
		virtual void SetOutputLevel( int Input ) = 0;

		/*!
		 * @param Interface Function:
		 *        Set the maximum number of Steps that the minimiser can take before it has minimised
		 *
		 * This is likely common to all Minimisers so is part of the interface, but the Minimiser isn't demanded to have to honor this request
		 *
		 * @param Input   This is the Maximum number of steps you wish to constrict the Minimiser to
		 *
		 * @return Void
		 */
		virtual void SetSteps( int Input ) = 0;

		/*!
		 * @brief Interface Function:
		 *        Set the Tolerance at which the Minimiser should decide it has converged
		 *
		 * @param Input   This is the Minimim Tolerance which is the distance from the Minima the Minimser should converge to
		 *
		 * @return Void
		 */
		virtual void SetTolerance( double Input ) = 0;

		/*!
		 * @brief Interface Function:
		 *        Provide a vector of options in string format to customize the Minimiser behaviour
		 *
		 * @param Input   These are Strings contained in the XML which are Minimiser Specific Configuration Options and so this is a common way of communicating them to the Minmiser
		 *
		 * @return Void
		 */
		virtual void SetOptions( vector<string> Input ) = 0;

		/*!
		 * @brief Interface Function:
		 *        Set the Quality of the Fit Convergence
		 *
		 * @param Input   This is the quality of the Fit you wish to use, generally the Higher the number the better the quality, This may be Minimiser Dependent
		 *
		 * @return Void
		 */
		virtual void SetQuality( int ) = 0;

		/*!
		 * @brief Interface Function:
		 *         Return a Pointer to the FitFunction used internally
		 *
		 * @return Returns a pointer to the FitFunction that was Minimised in it's Minimised State
		 */
		virtual FitFunction* GetFitFunction() = 0;
};

#endif

