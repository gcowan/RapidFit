/*!
 * @class MinimiserConfiguration
 *
 * Container that stores all information related to minimiser configuration, and returns an appropriate instance of a minimiser
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

///	RapidFit Headers
#include "MinimiserConfiguration.h"
#include "ClassLookUp.h"
///	System Headers
#include <string>
#include <sstream>

using namespace::std;

//Constructor for a minimiser only specified by name
MinimiserConfiguration::MinimiserConfiguration( string InputName ) :
	theMinimiser(), OutputLevel(), minimiserName(InputName), contours(), maxSteps(), bestTolerance(), MinimiseOptions(), MultiMini(false), Quality(), debug(new DebugClass(false) ), nSigma(1)
{
	theMinimiser = NULL;
	OutputLevel=0;
}

//Constructor for a minimiser with requested contour plots
MinimiserConfiguration::MinimiserConfiguration( string InputName, OutputConfiguration * Formatting ) :
	theMinimiser(), OutputLevel(), minimiserName(InputName), nSigma(1),
	contours( Formatting != NULL ? Formatting->GetContourPlots() : vector< pair< string, string > >() ), maxSteps(), bestTolerance(), MinimiseOptions(), MultiMini(false), Quality(), debug(new DebugClass(false) )
{
	theMinimiser = NULL;
	OutputLevel=0;
}

//Destructor
MinimiserConfiguration::~MinimiserConfiguration()
{
	delete theMinimiser;
	if( debug != NULL ) delete debug;
}

void MinimiserConfiguration::SetOutputLevel( int output_Level )
{
	OutputLevel=output_Level;
}

IMinimiser* MinimiserConfiguration::GetMinimiser()
{
	return theMinimiser;
}

//Return an appropriate minimiser instance
IMinimiser * MinimiserConfiguration::GetMinimiser( int ParameterNumber )
{
	if( theMinimiser != NULL && MultiMini )
	{
		return theMinimiser;
	}
	else
	{
		if( theMinimiser != NULL )
		{
			delete theMinimiser;
			theMinimiser = NULL;
		}
	}

	theMinimiser = ClassLookUp::LookUpMinimiserName( minimiserName, ParameterNumber );
	theMinimiser->SetOutputLevel( OutputLevel );

	//Supply the output configuration
	if( contours.size() > 0 )
	{
		theMinimiser->ContourPlots(contours);
	}
	theMinimiser->SetSteps( maxSteps );
	theMinimiser->SetTolerance( bestTolerance );
	theMinimiser->SetOptions( MinimiseOptions );
	theMinimiser->SetQuality( Quality );
	theMinimiser->SetNSigma( nSigma );
	return theMinimiser;
}

void MinimiserConfiguration::RemoveMinimiser()
{
	delete theMinimiser;
	theMinimiser=NULL;
}

void MinimiserConfiguration::SetSteps( int steps )
{
	maxSteps = steps;
}

void MinimiserConfiguration::SetTolerance( double newTolerance )
{
	bestTolerance = newTolerance;
}

void MinimiserConfiguration::SetOptions( vector<string> Options )
{
	MinimiseOptions = Options;
}

void MinimiserConfiguration::SetQuality( int newQuality )
{
	Quality = newQuality;
}

void MinimiserConfiguration::SetMultiMini( bool decision )
{
	MultiMini = decision;
}

void MinimiserConfiguration::SetNSigma( int input )
{
	nSigma = input;
}

void MinimiserConfiguration::Print() const
{
	cout << "MinimiserConfiguration:" << endl;
	cout << "This is to be coded up when needed" << endl;
}

string MinimiserConfiguration::XML() const
{
	stringstream xml;

	xml << "<Minimiser>" << endl;
	xml << "\t" << "<MinimiserName>" << minimiserName << "</MinimiserName>" << endl;
	xml << "\t" << "<MaxSteps>" << maxSteps << "</MaxSteps>" << endl;
	xml << "\t" << "<GradTolerance>" << bestTolerance << "</GradTolerance>" << endl;
	xml << "\t" << "<Quality>" << Quality << "</Quality>" << endl;
	xml << "</Minimiser>" << endl;

	return xml.str();
}

void MinimiserConfiguration::SetDebug( DebugClass* input_debug )
{
	if( debug != NULL ) delete debug;
	debug = new DebugClass( *input_debug );
}

