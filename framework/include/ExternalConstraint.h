/*!
 * @class ExternalConstraint
 *
 * @brief A class that holds experimentally dervied constraints on fit parameters
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef EXTERNAL_CONSTRAINT_H
#define EXTERNAL_CONSTRAINT_H

//	System Headers
#include <string>

using namespace::std;

class ExternalConstraint
{
	public:
		/*!
		 * @brief Default Constructor
		 */
		ExternalConstraint();

		/*!
		 * @brief Constructor containing the external constraint name, value and error
		 *
		 * @param Name  Name of the Physics Parameter we're applying the External Constrint to
		 *
		 * @param Value This is CV of the external constraint being applied
		 *
		 * @param Error This is the Error on the Constraint 'width of the Gaussian'
		 */
		ExternalConstraint( string Name, double Value, double Error );

		/*!
		 * Destructor
		 */
		~ExternalConstraint();

		/*!
		 * @brief Get the Name of the ExternalConstraint
		 *
		 * @return Returns the name of the PhysicsParameter Being Constrained
		 */
		string GetName() const;

		/*!
		 * @brief Get the Value of the ExternalConstraint
		 *
		 * @return the Central Value of the Constraint
		 */
		double GetValue() const;

		/*!
		 * @brief Get the Error of the ExternalConstraint
		 *
		 * @return returns the Error on the Constraint
		 */
		double GetError() const;

		void Print() const;

		string XML() const;

	private:
		string name;		/*!	External Constraint Name		*/
		double value, error;	/*!	External Constraint Value and Error	 */
};

#endif

