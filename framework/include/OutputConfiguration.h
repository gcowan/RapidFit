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
#include "ScanParam.h"

using namespace std;

class OutputConfiguration
{
	public:
		OutputConfiguration();
		OutputConfiguration( vector< pair< string, string > >, vector<string>, string, vector<ScanParam*>, vector<pair<ScanParam*, ScanParam*> > );
		~OutputConfiguration();

		//Return configuration data
		vector< pair< string, string > > GetContourPlots();
		vector<pair<string, string> > Get2DScanList( string input_type="LLscan" );
		vector<pair<string, string> > Get2DLLscanList();
		vector<string> GetProjections();
		vector<string> GetLLscanList() ;
		vector<string> GetCVScanList() ;
		vector<string> GetScanList( string="LLscan" ) ;

		ScanParam* GetScanParam( string param_name, string="LLscan" );
		pair<ScanParam*, ScanParam*> Get2DScanParams( string , string, string="LLscan", string="LLscan" );

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
		void SetInputResults( ResultParameterSet* );

	private:
		vector<double> GetRange( string );		//  Return Max, Min and Resolution
		vector<double> GetRange( ScanParam* );		//  Return Max, Min and Resolution
		pair<vector<double>, vector<double> > Get2DRange( string, string );

		vector< pair< string, string > > contours;
		vector<string> projections;
		vector<string> LLscanList;
		bool makeAllPlots;
		bool weightedEventsWereUsed ;
		string weightName ;
		string LLcontourFileName, LLscanFileName, pullFileName, projectionFileName, contourFileName, pullType;
		vector<ScanParam*> Global_Scan_List;
		vector<pair<ScanParam*, ScanParam*> > Global_2DScan_List;

		vector<ResultParameterSet* > Stored_Fit_Results;
};

#endif
