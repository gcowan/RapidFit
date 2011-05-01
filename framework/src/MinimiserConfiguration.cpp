/**
	@class MinimiserConfiguration

	Container that stores all information related to minimiser configuration, and returns an appropriate instance of a minimiser

	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-11-27
*/

//	RapidFit Headers
#include "MinimiserConfiguration.h"
#include "ClassLookUp.h"

//Default constructor
MinimiserConfiguration::MinimiserConfiguration() : theMinimiser(), OutputLevel(), minimiserName(), contours()
{
	theMinimiser = NULL;
	OutputLevel=0;
}

//Constructor for a minimiser only specified by name
MinimiserConfiguration::MinimiserConfiguration( string InputName ) : theMinimiser(), OutputLevel(), minimiserName(InputName), contours()
{
	theMinimiser = NULL;
	OutputLevel=0;
}

//Constructor for a minimiser with requested contour plots
MinimiserConfiguration::MinimiserConfiguration( string InputName, OutputConfiguration * Formatting ) : theMinimiser(), OutputLevel(), minimiserName(InputName), contours( Formatting->GetContourPlots() )
{
	theMinimiser = NULL;
	OutputLevel=0;
}

//Destructor
MinimiserConfiguration::~MinimiserConfiguration()
{
	delete theMinimiser;
}

void MinimiserConfiguration::SetOutputLevel( int output_Level )
{
	OutputLevel=output_Level;
}

//Return an appropriate minimiser instance
IMinimiser * MinimiserConfiguration::GetMinimiser( int ParameterNumber )
{
	theMinimiser = ClassLookUp::LookUpMinimiserName( minimiserName, ParameterNumber );
	theMinimiser->SetOutputLevel(short(OutputLevel));

	//Supply the output configuration
	if ( contours.size() > 0 )
	{
		theMinimiser->ContourPlots(contours);
	}

	return theMinimiser;
}

void MinimiserConfiguration::RemoveMinimiser()
{
	delete theMinimiser;
	theMinimiser=NULL;
}
