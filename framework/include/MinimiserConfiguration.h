/**
	@class MinimiserConfiguration

	Container that stores all information related to minimiser configuration, and returns an appropriate instance of a minimiser

	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-11-27
*/

#ifndef MINIMISER_CONFIGURATION_H
#define MINIMISER_CONFIGURATION_H

#include "IMinimiser.h"

class MinimiserConfiguration
{
	public:
		MinimiserConfiguration();
		MinimiserConfiguration(string);
		MinimiserConfiguration( string, vector< pair< string, string > > );
		~MinimiserConfiguration();

		IMinimiser * GetMinimiser(int);

	private:
		string minimiserName;
		vector< pair< string, string > > contours;
};

#endif
