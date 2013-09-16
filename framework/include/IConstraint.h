/*!
 * @interface IConstraint
 *
 * @brief Interface for all constraints defining a phase space boundary
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef I_CONSTRAINT_H
#define I_CONSTRAINT_H

///	ROOT Headers
#include "TRandom3.h"
///	RapidFit Jeaders
#include "Observable.h"
///	System Headers
#include <string>
#include <vector>

using namespace::std;

class IConstraint
{
	public:
		/*!
		 * @brief Default Destructor
		 */
		virtual ~IConstraint() {};

		/*!
		 * @brief Interface Function:
		 *        Return true if an observable value is compatible with the constraint
		 *
		 * @param Input    This is the Observable which is checked against the Constraint to see if it is compatible or not
		 *
		 * @return true, The Observable is within/compatible the Constraint  false, The Observable is NOT within/compatible with the Constraint
		 */
		virtual bool CheckObservable( Observable* ) const = 0;

		/*!
		 * @brief Interface Function:
		 *        Create a Random Observable, uses gRandom so may not be reproducible
		 *
		 * @warning The Observable returned here is not in anyway looked after by the Constraint object
		 *
		 * @return This returns a new Observable which has a Random value compatible with this Constraint
		 */
		virtual Observable * CreateObservable() const = 0;

		/*!
		 * @brief Interface Function:
		 *        Create a Random Observable matching this constraint using this random number generator
		 *
		 * @warning The Observable returned here is not in anyway looked after by the Constraint object
		 *
		 * @return This returns a new Observable which has a Random value compatible with this Constraint
		 */
		virtual Observable * CreateObservable( TRandom3* ) const = 0;

		/*!
		 * @brief Interface Function:
		 *        Return the name of the Constraint 'Parameter'
		 *
		 * @return This Returns the Name of this Constraint
		 */
		virtual string GetName() const = 0;

		/*!
		 * @brief Interface Function:
		 *        Return the unit of this Constraint
		 *
		 * @returns The unit corresponding to the Constraint as a string
		 */
		virtual string GetUnit() const = 0;

		/*!
		 * @brief Interface Function:
		 *        Return the Maximum Value of this Constraint (Continuous)
		 *
		 * @return Returns the maximum Value consistent with this Constraint when it's Continuous
		 */
		virtual double GetMaximum() const = 0;

		/*!
		 * @brief Interface Function:
		 *        Return the Minimum value of this Constraint (Continuous)
		 *
		 * @return Returns the minimum Value consistent with this Constraint when it's Continuous
		 */
		virtual double GetMinimum() const = 0;

		/*!
		 * @brief Interface Function:
		 *        Return the possible values of this Constraint (Discrete)
		 *
		 * @return a vector of the possible values that this Constraint allows when it's Discrete
		 */
		virtual vector<double> GetValues() const = 0;

		/*!
		 * @brief Interface Function:
		 *        true if Discrete
		 *        false if Continuous
		 *
		 * @return true = Discrete, false = NOT Discrete
		 */
		virtual bool IsDiscrete() const = 0;

		/*!
		 * @brief Print some helpful debugging info
		 *
		 * @return Void
		 */
		virtual void Print() const = 0;

		/*!
		 * @brief Interface Function
		 *        Create a new Observable which has a value equal to the minimum of this discrete set
		 *
		 * @return Returns an Observable which has a value equal to the Mid Point of this Constraint
		 */
		virtual Observable* GetMidRangeValue() const = 0;

		/*!
		 * @brief Interface Function
		 *        Return the TF1 associated with reading this constraint in from a ROOT file
		 *
		 * ATM this is only used to construct an Observable from a ROOT File, but in theory this TF1 or a similar could be used for Observable transforms
		 *
		 * @return Returns the TF1 in a string object
		 */
		virtual string GetTF1() const = 0;

		/*!
		 * @brief Interface Function
		 *        Return a string which gives the XML required to construct this Constraint
		 *
		 * @return Returns the XML to construct this object as a flat string
		 */
		virtual string XML() const = 0;

	protected:
		/*!
		 * @brief Default Constructor
		 */
		IConstraint() {};
};

#endif

