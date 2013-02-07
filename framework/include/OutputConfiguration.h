/*!
 * @ingroup Generators
 * @class OutputConfiguration
 *
 * A container for all information configuring RapidFit output
 *
 * @author Benjamin Wynne bwynne@cern.ch
  */

#pragma once
#ifndef OUTPUT_CONFIGURATION_H
#define OUTPUT_CONFIGURATION_H

//	RapidFit Headers
//#include "LLscanResult.h"
//#include "LLscanResult2D.h"
#include "FitResultVector.h"
#include "FitResult.h"
#include "ScanParam.h"
#include "ComponentPlotter.h"
//	System Headers
#include <vector>
#include <string>
#include <utility>

using namespace::std;

class OutputConfiguration
{
	public:
		OutputConfiguration( vector<pair<string, string> > contour, string, vector<ScanParam*>, vector<pair<ScanParam*, ScanParam*> >, vector<CompPlotter_config*> );
		~OutputConfiguration();

		//Return configuration data
		vector< pair< string, string > > GetContourPlots();
		vector<pair<string, string> > Get2DScanList();
		vector<string> GetProjections();
		vector<string> GetScanList();

		ScanParam* GetScanParam( string param_name );
		pair<ScanParam*, ScanParam*> Get2DScanParams( string , string );

		void ClearScanList();
		void Clear2DScanList();

		bool DoPullPlots();

		//Make output requested
		void OutputFitResult( FitResult* );
		void OutputToyResult( FitResultVector* );

		void OutputCompProjections( FitResult* );

		string GetPullFileName() const;

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

		void SetDebug( DebugClass* debug );
	private:
		OutputConfiguration ( const OutputConfiguration& );
		OutputConfiguration& operator= ( const OutputConfiguration& );
		vector<double> GetRange( const string );		//  Return Max, Min and Resolution
		vector<double> GetRange( ScanParam* );		//  Return Max, Min and Resolution
		pair<vector<double>, vector<double> > Get2DRange( const string, const string );

		vector< pair< string, string > > contours;
		vector<string> LLscanList;
		string pullType;
		bool makeAllPlots;
		string weightName ;
		string pullFileName,projectionFileName,LLscanFileName,LLcontourFileName,contourFileName;
		bool weightedEventsWereUsed ;
		vector<ScanParam*> Global_Scan_List;
		vector<pair<ScanParam*, ScanParam*> > Global_2DScan_List;

		vector<ResultParameterSet* > Stored_Fit_Results;

		vector<CompPlotter_config*> proj_components;

		DebugClass* debug;
};

#endif

