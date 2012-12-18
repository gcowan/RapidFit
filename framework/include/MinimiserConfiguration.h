/**
 * @ingroup Configurators  This Generator is capable of Constructing the Minimiser which minimises the FitFunction
 * @class MinimiserConfiguration
 *
 *Container that stores all information related to minimiser configuration, and returns an appropriate instance of a minimiser
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef MINIMISER_CONFIGURATION_H
#define MINIMISER_CONFIGURATION_H

//	RapidFit Headers
#include "IMinimiser.h"
#include "OutputConfiguration.h"
#include "DebugClass.h"

#include <string>
#include <vector>

using namespace::std;

class MinimiserConfiguration
{
	public:
		MinimiserConfiguration( string );
		MinimiserConfiguration( string, OutputConfiguration* );
		~MinimiserConfiguration();

		IMinimiser * GetMinimiser();
		IMinimiser * GetMinimiser(int);
		void SetOutputLevel(int);
		void RemoveMinimiser();

		void SetSteps(int);
		void SetTolerance(double);
		void SetOptions(vector<string>);
		void SetQuality( int );
		void SetMultiMini( bool );
		void SetNSigma( int );

		//Output some debugging info
		void Print() const;

		string XML() const;

		void SetDebug( DebugClass* debug );
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
		bool MultiMini;
		int Quality;
		int nSigma;
		DebugClass* debug;
};

#endif

