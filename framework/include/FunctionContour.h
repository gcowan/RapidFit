/*!
 * @class FunctionContour
 *
 * @brief A class holding a contour plot retrieved from Minuit
 *
 * @author Benjamin Wynne bwynne@cern.ch
  */

#pragma once
#ifndef FUNCTION_CONTOUR_H
#define FUNCTION_CONTOUR_H

//	System Headers
#include <vector>
#include <string>

using namespace::std;

class FunctionContour
{
	public:
		FunctionContour();
		FunctionContour( string, string, int );
		~FunctionContour();

		string GetXName();
		string GetYName();
		int GetContourNumber();
		void SetPlot( int, int, double*, double* );
		void SetPlot( int, vector< pair< double, double > > );
		vector< pair< double, double > > GetPlot(int);

	private:
		string xName, yName;
		vector< vector< pair< double, double > > > allContours;
};

#endif

