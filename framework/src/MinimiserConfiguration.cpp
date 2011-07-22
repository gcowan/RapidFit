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
MinimiserConfiguration::MinimiserConfiguration() : theMinimiser(), OutputLevel(), minimiserName(), contours(), maxSteps(), bestTolerance(), MinimiseOptions(), MultiMini(false)
{
	theMinimiser = NULL;
	OutputLevel=0;
}

//Constructor for a minimiser only specified by name
MinimiserConfiguration::MinimiserConfiguration( string InputName ) : theMinimiser(), OutputLevel(), minimiserName(InputName), contours(), maxSteps(), bestTolerance(), MinimiseOptions(), MultiMini(false)
{
	theMinimiser = NULL;
	OutputLevel=0;
}

//Constructor for a minimiser with requested contour plots
MinimiserConfiguration::MinimiserConfiguration( string InputName, OutputConfiguration * Formatting ) : theMinimiser(), OutputLevel(), minimiserName(InputName), contours( Formatting->GetContourPlots() ), maxSteps(), bestTolerance(), MinimiseOptions(), MultiMini(false)
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

IMinimiser* MinimiserConfiguration::GetMinimiser()
{
	return theMinimiser;
}

//Return an appropriate minimiser instance
IMinimiser * MinimiserConfiguration::GetMinimiser( int ParameterNumber )
{
	if( theMinimiser != NULL && MultiMini ) {
		return theMinimiser;
	} else {
		if( theMinimiser != NULL )
		{
			delete theMinimiser;
			theMinimiser = NULL;
		}
	}

	theMinimiser = ClassLookUp::LookUpMinimiserName( minimiserName, ParameterNumber );
	theMinimiser->SetOutputLevel(short(OutputLevel));

	//Supply the output configuration
	if ( contours.size() > 0 )
	{
		theMinimiser->ContourPlots(contours);
	}
	theMinimiser->SetSteps( maxSteps );
	theMinimiser->SetTolerance( bestTolerance );
	theMinimiser->SetOptions( MinimiseOptions );
	theMinimiser->SetQuality( Quality );
	return theMinimiser;
}

void MinimiserConfiguration::RemoveMinimiser()
{
	delete theMinimiser;
	theMinimiser=NULL;
}

void MinimiserConfiguration::SetSteps( int steps )
{
	maxSteps = steps;
}

void MinimiserConfiguration::SetTolerance( double newTolerance )
{
	bestTolerance = newTolerance;
}

void MinimiserConfiguration::SetOptions( vector<string> Options )
{
	MinimiseOptions = Options;
}

void MinimiserConfiguration::SetQuality( int newQuality )
{
	Quality = newQuality;
}

void MinimiserConfiguration::SetMultiMini( bool decision )
{
	MultiMini = decision;
}
