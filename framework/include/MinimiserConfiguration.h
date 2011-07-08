/**
	@class MinimiserConfiguration

	Container that stores all information related to minimiser configuration, and returns an appropriate instance of a minimiser

	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-11-27
*/

#ifndef MINIMISER_CONFIGURATION_H
#define MINIMISER_CONFIGURATION_H

//	RapidFit Headers
#include "IMinimiser.h"
#include "OutputConfiguration.h"

class MinimiserConfiguration
{
	public:
		MinimiserConfiguration();
		MinimiserConfiguration(string);
		MinimiserConfiguration( string, OutputConfiguration* );
		~MinimiserConfiguration();

		IMinimiser * GetMinimiser(int);
		void SetOutputLevel(int);
		void RemoveMinimiser();

		void SetSteps(int);
		void SetTolerance(double);
		void SetOptions(vector<string>);
		void SetQuality( int );

	private:
		//	Uncopyable!
		MinimiserConfiguration ( const MinimiserConfiguration& );
		MinimiserConfiguration& operator = ( const MinimiserConfiguration& );

		IMinimiser* theMinimiser;
		int OutputLevel;
		string minimiserName;
		vector< pair< string, string > > contours;
		int maxSteps;
		double bestTolerance;
		vector<string> MinimiseOptions;
		int Quality;
};

#endif
