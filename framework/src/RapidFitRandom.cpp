
#include "TRandom3.h"
#include "RapidFitRandom.h"

#include <cmath>


/*!
 * @brief Start a new TRandom3 instance with a given seed value
 *
 * @param num   This is the new number to use as a seed for creating a new TRandom3 instance and delete the existing one
 *
 * @return Void
 */
void RapidFitRandom::SetRandomFunction( int num )
{
	RapidFitRandom::seed_num = num;
	if( RapidFitRandom::seed_function != NULL ) delete RapidFitRandom::seed_function;
	RapidFitRandom::seed_function = new TRandom3( RapidFitRandom::seed_num );
}

/*!
 * @brief Replace TRandom3 instance with a new function
 *
 * @param Input This is the new TRandom3 function we want to associate with this PDF
 >                 *
 * @return Void
 */
void RapidFitRandom::SetRandomFunction( TRandom3* Input )
{
	RapidFitRandom::seed_num = (int)ceil(Input->Rndm() * 1E5);
	if( RapidFitRandom::seed_function != NULL ) delete RapidFitRandom::seed_function;
	RapidFitRandom::seed_function = new TRandom3( RapidFitRandom::seed_num );
}

/*!
 * @brief Get the seed number used to initialise the random number gen
 *
 * @return Returns the number used as the seed of the TRandom3 instance
 */
int RapidFitRandom::GetSeedNum()
{
	return RapidFitRandom::seed_num;
}

/*!
 * @brief Get the Random function stored in this PDF
 *
 * @return pointer to the TRandom3 instance inside the PDF
 */
TRandom3* RapidFitRandom::GetRandomFunction()
{
	if( RapidFitRandom::seed_function == NULL )
	{
		if( RapidFitRandom::seed_num < 0 )
		{
			RapidFitRandom::seed_function = new TRandom3( (unsigned)seed_num );
			RapidFitRandom::seed_num = (unsigned) RapidFitRandom::seed_num;
		}
		else
		{
			RapidFitRandom::seed_function = new TRandom3( seed_num );
		}
	}
	return RapidFitRandom::seed_function;
}

TRandom3* RapidFitRandom::GetFrameworkRandomFunction()
{
	if( RapidFitRandom::seed_function_Framework == NULL )
	{
		//	We Explicityly use this to avoid collisions when small objects for caching are duplicated
		//	Hence, as this is NEVER INTENDED TO BE USED FOR PHYSICS SIM it is safe to be initialized
		//	using as much entropy as we can get
		RapidFitRandom::seed_function_Framework = new TRandom3(0);
	}
	return RapidFitRandom::seed_function_Framework;
}

int RapidFitRandom::GetSeedNumFramework()
{
	return RapidFitRandom::GetFrameworkRandomFunction()->Rndm();
}

/*!
 * vector of a single TRandom3 object for this PDF, this allows us to have a reproducable result for a defined seed
 * the seed_function.empty() is used to see if this is defined, should probably check for NULL pointer, but oh well
 */
TRandom3* RapidFitRandom::seed_function = NULL;

TRandom3 * RapidFitRandom::seed_function_Framework = NULL;

/*!
 * vector of single seed value used in the construction of the TRandom3, again should probably check for a null value, but oh well
 */
int RapidFitRandom::seed_num = -1;

