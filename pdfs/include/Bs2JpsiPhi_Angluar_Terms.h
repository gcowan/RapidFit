#ifndef _Jpp_Angular_H
#define _Jpp_Angular_H

#include <vector>

using namespace::std;

class Bs2JpsiPhi_Angular_Terms
{
	public:

		/*!
		 * These functions Assume that the they are provided with an input of vector<double>
		 *
		 * The format of the input is:
		 *
		 * Transversity:   [0] = cosTheta , [1] = cosPsi   , [2] = phi
		 *
		 * Helicity:       [0] = cosThetaK, [1] = cosThetaL, [2] = phi
		 */

		//------------------------------------------------------
		// Angle factors for three angle PDFs  in transversity basis

		//........ for one angle tests ...................
		static double TangleFactorEven( vector<double> input );
		static double TangleFactorOdd( vector<double> input );

		//........ P Wave ..........
		//...........................
		static double TangleFactorA0A0( vector<double> input );
		//...........................
		static double TangleFactorAPAP( vector<double> input );
		//...........................
		static double TangleFactorATAT( vector<double> input );
		//...........................
		static double TangleFactorImAPAT( vector<double> input );
		//...........................
		static double TangleFactorReA0AP( vector<double> input );
		//...........................
		static double TangleFactorImA0AT( vector<double> input );

		//......  S wave  ....
		//.............................
		static double TangleFactorASAS( vector<double> input );
		//...........................
		static double TangleFactorReASAP( vector<double> input );
		//...........................
		// This appears to be different to the LHCB note, but on inspection it is not. It is the difference in sign of ImASAT <=> ImATAS
		static double TangleFactorImASAT( vector<double> input );
		//...........................
		static double TangleFactorReASA0( vector<double> input );



		//------------------------------------------------------
		// Angle factors for three angle PDFs  in helicity basis

		//........ for one angle tests ...................
		static double HangleFactorEven( vector<double> input );
		static double HangleFactorOdd( vector<double> input );


		//........ P Wave ..........
		// .................
		static double HangleFactorA0A0( vector<double> input );
		//..................
		static double HangleFactorAPAP( vector<double> input );
		//..................
		static double HangleFactorATAT( vector<double> input );
		//..................
		static double HangleFactorImAPAT( vector<double> input );
		//..................
		static double HangleFactorReA0AP( vector<double> input );
		//..................
		static double HangleFactorImA0AT( vector<double> input );

		//......  S wave  ....
		//.............................
		static double HangleFactorASAS( vector<double> input );
		//...........................
		static double HangleFactorReASAP( vector<double> input );
		//...........................
		static double HangleFactorImASAT( vector<double> input );
		//...........................
		static double HangleFactorReASA0( vector<double> input );

};

#endif

