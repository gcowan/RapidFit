/**
  @class OutputConfiguration

  A container for all information configuring RapidFit output

  @author Benjamin Wynne bwynne@cern.ch
  */

#include "OutputConfiguration.h"

//Default constructor
OutputConfiguration::OutputConfiguration() : doPulls(false)
{
}

//Constructor with correct arguments
OutputConfiguration::OutputConfiguration( vector< pair< string, string > > InputContours, vector<string> InputProjections, bool DoPullPlots ) : contours(InputContours), projections(InputProjections), doPulls(DoPullPlots)
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

//Return the requested projections
vector<string> OutputConfiguration::GetProjections()
{
	return projections;
}

//Return whether to do pull plots
bool OutputConfiguration::DoPullPlots()
{
	return doPulls;
}
