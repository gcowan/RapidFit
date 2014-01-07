
#include <iostream>

class AcceptanceSlice
{
	public:
		/*!
		 * @brief Constructor
		 *
		 * @param tl
		 *
		 * @param th
		 *
		 * @param h
		 *
		 */
		AcceptanceSlice( double tl, double th, double h ) : _tlow(tl), _thigh(th), _height(h) {}

		/*!
		 * @brief Copy Constructor
		 */
		AcceptanceSlice( const AcceptanceSlice& input );

		/*!
		 * @brief
		 */
		double tlow()   const { return _tlow ; }

		/*!
		 * @brief
		 */
		double thigh()  const { return _thigh ; }

		/*!
		 * @brief
		 */
		double height() const { return _height ; }

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

