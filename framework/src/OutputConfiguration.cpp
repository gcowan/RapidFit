/**
  @class OutputConfiguration

  A container for all information configuring RapidFit output

  @author Benjamin Wynne bwynne@cern.ch
  */

#include "OutputConfiguration.h"

//Default constructor
OutputConfiguration::OutputConfiguration()
{
}

//Constructor with correct arguments
OutputConfiguration::OutputConfiguration( vector< pair< string, string > > InputContours ) : contours(InputContours)
{
}

//Destructor
OutputConfiguration::~OutputConfiguration()
{
}

//Return the requested contour plots
vector< pair< string, string > > OutputConfiguration::GetContourPlots()
{
	return contours;
}
