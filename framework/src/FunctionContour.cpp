/**
  @class FunctionContour

  A class holding a contour plot retrieved from minuit

  @author Benjamin Wynne bwynne@cern.ch
  @date 2009-11-27
 */

//	RapidFit Headers
#include "FunctionContour.h"
//	System Headers
#include <iostream>

//Default constructor
FunctionContour::FunctionContour() : xName(), yName(), allContours()
{
}

//Constructor with correct arguments
FunctionContour::FunctionContour( string XName, string YName, int ContourNumber ) : xName(XName), yName(YName), allContours()
{
	//Initialise the contour storage
	for ( int contourIndex = 0; contourIndex < ContourNumber; ++contourIndex )
	{
		vector< pair< double, double > > empty;
		allContours.push_back(empty);
	}
}

//Destructor
FunctionContour::~FunctionContour()
{
}

//Get the name of the x-axis parameter
string FunctionContour::GetXName()
{
	return xName;
}

//Get the name of  the y-axis parameter
string FunctionContour::GetYName()
{
	return yName;
}

//Get the number of contours
int FunctionContour::GetContourNumber()
{
	return int(allContours.size());
}

//Store a contour
void FunctionContour::SetPlot( int Sigma, int NumberPoints, double * XValues, double * YValues )
{
	if ( Sigma > int(allContours.size()) || Sigma < 1 )
	{
		cerr << "Contour sigma value (" << Sigma << ") is invalid" << endl;
	}
	else
	{
		vector< pair< double, double > > newContourPlot;

		//Store each pair of points
		for ( int pointIndex = 0; pointIndex < NumberPoints; ++pointIndex )
		{
			newContourPlot.push_back( make_pair( XValues[pointIndex], YValues[pointIndex] ) );
		}

		allContours[ unsigned(Sigma - 1) ] = newContourPlot;
	}
}
void FunctionContour::SetPlot( int Sigma, vector< pair< double, double > > Contour )
{
	if ( Sigma > int(allContours.size()) || Sigma < 1 )
	{
		cerr << "Contour sigma value (" << Sigma << ") is invalid" << endl;
	}       
	else
	{
		allContours[ unsigned(Sigma - 1) ] = Contour;
	}
}

//Retrieve a contour
vector< pair< double, double > > FunctionContour::GetPlot( int Sigma )
{
	if ( Sigma > int(allContours.size()) || Sigma < 1 )
	{
		cerr << "Contour sigma value (" << Sigma << ") is invalid" << endl;
	}
	else
	{
		return allContours[ unsigned(Sigma - 1) ];
	}
	return *(new vector<pair<double, double> >);
}
