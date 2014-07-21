
#ifndef RAPIDFITRANDOM_H
#define RAPIDFITRANDOM_H

#include "TRandom3.h"

class RapidFitRandom
{
	public:

		/*!
		 * @brief Start a new TRandom3 instance with a given seed value
		 *
		 * @param num   This is the new number to use as a seed for creating a new TRandom3 instance and delete the existing one
		 *
		 * @return Void
		 */
		static void SetRandomFunction( int num );

		/*!
		 * @brief Replace TRandom3 instance with a new function
		 *
		 * @param Input This is the new TRandom3 function we want to associate with this PDF
		 >                 *
		 * @return Void
		 */
		static void SetRandomFunction( TRandom3* Input );

		/*!
		 * @brief Get the seed number used to initialise the random number gen
		 *
		 * @return Returns the number used as the seed of the TRandom3 instance
		 */
		static int GetSeedNum();

		/*!
		 * @brief Get the Random function stored for Physics Sim
		 *
		 * @return pointer to the TRandom3 instance inside the centinel
		 */
		static TRandom3* GetRandomFunction();

		/*!
		 * @brief Get the Random function intended to be used for internal framework use
		 *
		 * @return pointer to the TRandom3 instance for framework inside the centinel
		 */
		static TRandom3* GetFrameworkRandomFunction();

		static int GetSeedNumFramework();

	private:

		/*!
		 * vector of a single TRandom3 object for this PDF, this allows us to have a reproducable result for a defined seed
		 * the seed_function.empty() is used to see if this is defined, should probably check for NULL pointer, but oh well
		 */
		static TRandom3 * seed_function;

		static TRandom3 * seed_function_Framework;

		/*!
		 * vector of single seed value used in the construction of the TRandom3, again should probably check for a null value, but oh well
		 */
		static int seed_num;

};

#endif

