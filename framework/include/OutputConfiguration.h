/**
  @class OutputConfiguration

  A container for all information configuring RapidFit output

  @author Benjamin Wynne bwynne@cern.ch
  */

#ifndef OUTPUT_CONFIGURATION_H
#define OUTPUT_CONFIGURATION_H

#include "ToyStudyResult.h"
#include <vector>
#include <string>

using namespace std;

class OutputConfiguration
{
	public:
		OutputConfiguration();
		OutputConfiguration( vector< pair< string, string > >, vector<string>, string );
		~OutputConfiguration();

		//Return configuration data
		vector< pair< string, string > > GetContourPlots();
		vector<string> GetProjections();
		bool DoPullPlots();

		//Make output requested
		void OutputFitResult( FitResult* );
		void OutputToyResult( ToyStudyResult* );

		//Change the configuration
		void MakeAllPlots(string);
		void SetPullFileName(string);
		void SetContourFileName(string);

	private:
		vector< pair< string, string > > contours;
		vector<string> projections;
		bool makeAllPlots;
		string pullFileName, projectionFileName, contourFileName, pullType;
};

#endif
