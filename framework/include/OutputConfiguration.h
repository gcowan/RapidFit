/**
  @class OutputConfiguration

  A container for all information configuring RapidFit output

  @author Benjamin Wynne bwynne@cern.ch
  */

#ifndef OUTPUT_CONFIGURATION_H
#define OUTPUT_CONFIGURATION_H

#include <vector>
#include <string>

using namespace std;

class OutputConfiguration
{
	public:
		OutputConfiguration();
		OutputConfiguration( vector< pair< string, string > > );
		~OutputConfiguration();

		vector< pair< string, string > > GetContourPlots();

	private:
		vector< pair< string, string > > contours;
};

#endif
