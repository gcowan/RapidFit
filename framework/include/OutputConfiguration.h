/**
  @class OutputConfiguration

  A container for all information configuring RapidFit output

  @author Benjamin Wynne bwynne@cern.ch
  */

#ifndef OUTPUT_CONFIGURATION_H
#define OUTPUT_CONFIGURATION_H

#include "LLscanResult.h"
#include "LLscanResult2D.h"
#include "ToyStudyResult.h"
#include <vector>
#include <string>

using namespace std;

class OutputConfiguration
{
	public:
		OutputConfiguration();
		OutputConfiguration( vector< pair< string, string > >, vector<string>, vector<string>, string );
		~OutputConfiguration();

		//Return configuration data
		vector< pair< string, string > > GetContourPlots();
		vector<string> GetProjections();
		vector<string> GetLLscanList() ;
		bool DoPullPlots();

		//Make output requested
		void OutputFitResult( FitResult* );
		void OutputLLscanResult( vector<LLscanResult*>  );
		void OutputLLcontourResult( vector<LLscanResult2D*>  );
		void OutputToyResult( ToyStudyResult* );

		//Change the configuration
		void MakeAllPlots(string);
		void SetPullFileName(string);
		void SetContourFileName(string);
		void SetLLscanFileName(string);
		void SetLLcontourFileName(string);
		void SetWeightsWereUsed( string ) ;

	private:
		vector< pair< string, string > > contours;
		vector<string> projections;
		vector<string> LLscanList;
		bool makeAllPlots;
		bool weightedEventsWereUsed ;
		string weightName ;
		string LLcontourFileName, LLscanFileName, pullFileName, projectionFileName, contourFileName, pullType;
};

#endif
