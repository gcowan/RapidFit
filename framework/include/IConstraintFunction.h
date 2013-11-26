/*!
 * @class ExternalConstraint
 *
 * @brief A class that holds experimentally dervied constraints on fit parameters
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef _CONSTRAINT_FUNCTION_H
#define _CONSTRAINT_FUNCTION_H

#include "ParameterSet.h"
//	System Headers
#include <string>

using namespace::std;

class IConstraintFunction
{
	protected:
		/*!
		 * @brief Default Constructor
		 */
		IConstraintFunction() {}

	public:

		/*!
		 * Destructor
		 */
		virtual ~IConstraintFunction() {};

		/*!
		 * @brief Get the Name of the ExternalConstraint
		 *
		 * @return Returns the name of the PhysicsParameter Being Constrained
		 */
		virtual string GetName() const = 0;

		virtual void SetPhysicsParameters( const ParameterSet* ) = 0;

		virtual bool CanApply( const ParameterSet* input ) const = 0;

		virtual double GetChi2() const = 0;

		virtual void Print() const = 0;

		virtual string XML() const = 0;

		virtual bool isExternalConst() const = 0;

		virtual bool isExternalMatrix() const = 0;
};

#endif

