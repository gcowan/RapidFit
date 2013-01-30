/*!
 * @brief This class is intended to parse the command line arguments (argc,argv) and change the objects in RapidFitConfiguration accordingly
 */

#pragma once
#ifndef _RAPIDFIT_PARSE_COMMANDLINE_H
#define _RAPIDFIT_PARSE_COMMANDLINE_H

#include "RapidFitConfiguration.h"

#include <string>

using namespace::std;

class ParseCommandLine
{
	public:
	/*!
	 * @brief This is a static method intended to contain all of the code to parse the command line options
	 * 
	 */
	static int ParseThisCommandLine( RapidFitConfiguration& config, vector<string> input );

	/*!
	 * @brief This contains all of the Help command from RapidFitHelp
	 */
	 static void RapidFitHelp();

	/*!
	 * @brief This prints some information about RapidFitAbout
	 */
	static void RapidFitAbout( string name );
};

const int BAD_COMMAND_LINE_ARG = -39;
const int NO_XML = -41;
const int exit_RapidFit = 2;

#endif


