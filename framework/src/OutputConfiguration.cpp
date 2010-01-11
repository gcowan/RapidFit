/**
  @class OutputConfiguration

  A container for all information configuring RapidFit output

  @author Benjamin Wynne bwynne@cern.ch
 */

#include "OutputConfiguration.h"
#include "ResultFormatter.h"
#include "Plotter.h"
#include <time.h>

//Default constructor
OutputConfiguration::OutputConfiguration() : pullType("None"), makeAllPlots(false), pullFileName("pullPlots.root"), projectionFileName("projectionPlots.root"), contourFileName("contourPlots.root")
{
}

//Constructor with correct arguments
OutputConfiguration::OutputConfiguration( vector< pair< string, string > > InputContours, vector<string> InputProjections, string PullPlotType ) : contours(InputContours), projections(InputProjections), pullType(PullPlotType),
	makeAllPlots(false), pullFileName("pullPlots.root"), projectionFileName("projectionPlots.root"), contourFileName("contourPlots.root")
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
	return !( pullType == "None" );
}

//Make the requested output from a single result
void OutputConfiguration::OutputFitResult( FitResult * TheResult )
{
	//Output information aboout the fit
	ResultFormatter::LatexOutputFitResult(TheResult);
	ResultFormatter::LatexOutputCovarianceMatrix(TheResult);

	//Output any calculated contours
	ResultFormatter::PlotFitContours( TheResult, contourFileName );

	//Make any requested projections
	//for ( int projectionIndex = 0; projectionIndex < projections.size(); projectionIndex++ )
	for ( int projectionIndex = 0; projectionIndex < 1; projectionIndex++ )
	{
		PhysicsBottle * resultBottle = TheResult->GetPhysicsBottle();
		
		//Loop over all PDFs, and plot
		for ( int resultIndex = 0; resultIndex < resultBottle->NumberResults(); resultIndex++ )
		{
			Plotter * testPlotter = new Plotter( resultBottle->GetResultPDF(resultIndex), resultBottle->GetResultDataSet(resultIndex) );
			char fileNumber[100];
			sprintf( fileNumber, "fit%d.", resultIndex );

			if (makeAllPlots)
			{
				testPlotter->PlotAllObservables( fileNumber + projectionFileName );
			}
			else
			{
				testPlotter->PlotObservables( fileNumber + projectionFileName, projections );
			}
		}
	}
}

//Make the requested output from a toy study
void OutputConfiguration::OutputToyResult( ToyStudyResult * TheResult )
{
	ResultFormatter::MakePullPlots( pullType, pullFileName, TheResult );
}

//Make all possible function projections, and save them in given location
void OutputConfiguration::MakeAllPlots( string FileName )
{
	makeAllPlots = true;
	projectionFileName = FileName;
}

//Set the location to store contour plots
void OutputConfiguration::SetContourFileName( string FileName )
{
	contourFileName = FileName;
}

//Set the location to store pull plots
void OutputConfiguration::SetPullFileName( string FileName )
{
	pullFileName = FileName;
}
