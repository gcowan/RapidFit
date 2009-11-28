/**
	@class MinimiserConfiguration

	Container that stores all information related to minimiser configuration, and returns an appropriate instance of a minimiser

	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-11-27
*/

#include "MinimiserConfiguration.h"
#include "ClassLookUp.h"

//Default constructor
MinimiserConfiguration::MinimiserConfiguration()
{
}

//Constructor for a minimiser only specified by name
MinimiserConfiguration::MinimiserConfiguration( string InputName ) : minimiserName(InputName)
{
}

//Constructor for a minimiser with requested contour plots
MinimiserConfiguration::MinimiserConfiguration( string InputName, vector< pair< string, string > > InputContours ) : minimiserName(InputName), contours(InputContours)
{
}

//Destructor
MinimiserConfiguration::~MinimiserConfiguration()
{
}

//Return an appropriate minimiser instance
IMinimiser * MinimiserConfiguration::GetMinimiser( int ParameterNumber )
{
	IMinimiser * theMinimiser = ClassLookUp::LookUpMinimiserName( minimiserName, ParameterNumber );

	//Request contour plots if they are specified
	if ( contours.size() > 0 )
	{
		theMinimiser->ContourPlots(contours);
	}

	return theMinimiser;
}
