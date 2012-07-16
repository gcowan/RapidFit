#ifndef _MATHEMATICS_CLASS_H
#define _MATHEMATICS_CLASS_H

#include "TPolyMarker3D.h"
#include "TF1.h"

#include <vector>
#include <algorithm>

using namespace::std;

class Mathematics
{
	public:

		/*!
		 * @brief a standard double comparison routine
		 *
		 * @return true they're the same, false they're different
		 */
		static bool double_equals_test(Double_t first, Double_t second);

		/*!
		 * @brief This will find all of the unique coordinates stored in a 1D vector of doubles
		 * 
		 * @return Returns the unique doubles each wrapped in a vector in a larger vector
		*/
		static vector<vector<Double_t> > Unique_Coords( vector<Double_t> input );

		/*!
		 * @brief This will find all of the unique coordinates stored in a 2D vector of pairs of doubles
		 *
		 * @return Returns the unique 2D coordinates in a vector
		 */
		static vector<vector<Double_t> > Unique_Coords( vector<pair<Double_t,Double_t> > input );

		/*!
		 * @brief This will find all of the unique coordinates stored in a n-D vector of doubles
		 * 
		 * @warning UNTESTED UNTESTED UNTESTED UNTESTED
		 *          This may just work, or it may burn down your house and delete your inbox,
		 *          I have no idea what will happen, but it should work
		 *
		 * @return Returns the unique corrdinates in a vector for n=1,2,3 and in principle x
		 */
		static vector<vector<Double_t> > Unique_Coords( vector<vector<Double_t> > input );

		/*!
		 * @brief Function which returns the Unique corrdinates contained within a TPolyMarker3D* set of 3d points
		 *
		 */
		static vector<vector<Double_t> > Unique_Coords( TPolyMarker3D *pm );

		/*!
		 * @brief More correct to fit to the gamma distribution than the landau function for those that require it
		 *
		 */
		static TF1* gamma_func( int OutputLevel=-1 );

		static TF1* raw_gamma_func( int OutputLevel=-1 );

		static TF1* gammaDist_func( int OutputLevel=-1 );

		/*!
		 * @brief More correct to fit to the landau dist for some other dists from toys
		 *
		 */
		static TF1* landau_func();

		/*!
		 * @brief Algorithm for comparing a 2D point and returning a boolean comparison
		 *        used for stl::unique
		 */
		static bool Unique_2D_Double( pair<double,double> one_pair, pair<double,double> two_pair );

		/*!
		 * @brief Algorithm for comparing a 2D point with another and returning a decison on the biggest one
		 *        used for stl::sort
		 */
		static bool Sort_first_Double( pair<double,double> one_pair, pair<double,double> two_pair );

		/*!
		 * @brief
		 *
		 */
		static bool Sort_second_Double( pair<double,double> one_pair, pair<double,double> two_pair );
};

#endif

