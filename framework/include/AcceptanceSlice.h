
#ifndef _ACCEPTANCE_SLICE_H
#define _ACCEPTANCE_SLICE_H

#include <iostream>

class AcceptanceSlice
{
	public:
		/*!
		 * @brief Constructor
		 *
		 * @param tl Lower edge of the time Acceptance bin
		 *
		 * @param th Upper edge of the time Acceptance bin
		 *
		 * @param h Height of the time Acceptance bin
		 *
		 */
		AcceptanceSlice( double tl, double th, double h ) : _tlow(tl), _thigh(th), _height(h) {}

		/*!
		 * @brief Copy Constructor
		 */
		AcceptanceSlice( const AcceptanceSlice& input );

		/*!
		 * @brief Return the lower limit of the bin
		 */
		double tlow()   const { return _tlow ; }

		/*!
		 * @brief Return the upper limit of the bin
		 */
		double thigh()  const { return _thigh ; }

		/*!
		 * @brief Return the Height of the bin
		 */
		double height() const { return _height ; }

		/*!
		  * @brief Print the contents of the bin for debugging
		  */
		void Print() const { cout << "Min: " << _tlow << endl << "Max: " << _thigh << endl << "Height: " << _height << endl; }

	private:
		/*!
		 * @brief Don't Copy This Way
		 */
		AcceptanceSlice& operator = ( const AcceptanceSlice& );

		/*!
		 * Internal Objects to an Acceptance Slice
		 */
		double _tlow;
		double _thigh;
		double _height;
};

#endif

