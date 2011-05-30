/**
  @class OutputConfiguration

  A container for all information configuring RapidFit output

  @author Benjamin Wynne bwynne@cern.ch
  */

#ifndef OUTPUT_CONFIGURATION_H
#define OUTPUT_CONFIGURATION_H

//	RapidFit Headers
//#include "LLscanResult.h"
//#include "LLscanResult2D.h"
#include "ToyStudyResult.h"
#include "FitResult.h"
#include "ScanParam.h"
//	System Headers
#include <vector>
#include <string>

using namespace std;

class OutputConfiguration
{
	public:
		OutputConfiguration();
		OutputConfiguration( vector< pair< string, string > >, vector<string>, string, vector<ScanParam*>, vector<pair<ScanParam*, ScanParam*> > );
		~OutputConfiguration();

		//Return configuration data
		vector< pair< string, string > > GetContourPlots();
		vector<pair<string, string> > Get2DScanList();
		vector<string> GetProjections();
		vector<string> GetScanList();

		ScanParam* GetScanParam( string param_name );
		pair<ScanParam*, ScanParam*> Get2DScanParams( string , string );

		void Clear2DScanList( );

		bool DoPullPlots();

		//Make output requested
		void OutputFitResult( FitResult* );
//		void OutputLLscanResult( vector<LLscanResult*>  );
//		void OutputLLcontourResult( vector<LLscanResult2D*>  );
		void OutputToyResult( ToyStudyResult* );

		//Change the configuration
		void MakeAllPlots( string );
		void SetPullFileName( string );
		void SetContourFileName( string );
		void SetLLscanFileName( string );
		void SetLLcontourFileName( string );
		void SetWeightsWereUsed( string ) ;
		void SetInputResults( ResultParameterSet* );
		
		void AddContour( const string, const string );
		void AddScan( const string );

	private:
		vector<double> GetRange( const string );		//  Return Max, Min and Resolution
		vector<double> GetRange( ScanParam* );		//  Return Max, Min and Resolution
		pair<vector<double>, vector<double> > Get2DRange( const string, const string );

		vector< pair< string, string > > contours;
		vector<string> projections;
		vector<string> LLscanList;
		string pullType;
		bool makeAllPlots;
		string weightName ;
		string pullFileName,projectionFileName,LLscanFileName,LLcontourFileName,contourFileName;
		bool weightedEventsWereUsed ;
		vector<ScanParam*> Global_Scan_List;
		vector<pair<ScanParam*, ScanParam*> > Global_2DScan_List;

		vector<ResultParameterSet* > Stored_Fit_Results;
};

#endif
