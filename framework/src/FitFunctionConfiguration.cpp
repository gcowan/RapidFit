/*!
 * @class FitFunctionConfiguration
 *
 * Container that stores all information related to FitFunction configuration, and returns an appropriate instance of a FitFunction
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

//	RapidFit Headers
#include "FitFunctionConfiguration.h"
#include "ClassLookUp.h"
///	System Headers
#include <string>
#include <sstream>
#include <iostream>

using namespace::std;

//Constructor with only name of FitFunction
FitFunctionConfiguration::FitFunctionConfiguration( string InputName ) :
	functionName(InputName), weightName(), hasWeight(false), wantTrace(false), TraceFileName(), traceCount(0), Threads(0), Strategy(), testIntegrator(true)
{
}

//Constructor for FitFunction with event weights
FitFunctionConfiguration::FitFunctionConfiguration( string InputName, string InputWeight ) :
	functionName(InputName), weightName(InputWeight), hasWeight(true), wantTrace(false), TraceFileName(), traceCount(0), Threads(0), Strategy(), testIntegrator(true)
{
}

//Destructor
FitFunctionConfiguration::~FitFunctionConfiguration()
{
}

//Return appropriate instance of FitFunction
FitFunction * FitFunctionConfiguration::GetFitFunction()
{
	FitFunction * theFunction = ClassLookUp::LookUpFitFunctionName(functionName);

	//Use event weights if specified
	if (hasWeight)
	{
		theFunction->UseEventWeights(weightName);
	}

	if( wantTrace )
	{
		theFunction->SetupTrace( TraceFileName, traceCount );
		++traceCount;
	}

	theFunction->SetThreads( Threads );

	theFunction->SetIntegratorTest( testIntegrator );

	return theFunction;
}

//Return whether weights are being used
bool FitFunctionConfiguration::GetWeightsWereUsed()
{
	return hasWeight;
}

//Return weight name
string FitFunctionConfiguration::GetWeightName()
{
	return weightName;
}

void FitFunctionConfiguration::SetupTrace( TString FileName )
{
	wantTrace = true;
	TraceFileName = FileName;
}

void FitFunctionConfiguration::SetStrategy( string newStrategy )
{
	Strategy = newStrategy;
}

string FitFunctionConfiguration::GetStrategy()
{
	return Strategy;
}

void FitFunctionConfiguration::SetThreads( int input )
{
	Threads = input;
}

void FitFunctionConfiguration::SetIntegratorTest( bool input )
{
	testIntegrator = input;
}

void FitFunctionConfiguration::Print() const
{
	cout << "FitFunctionConfiguration::Print() Yet to be written" << endl;
}

string FitFunctionConfiguration::XML() const
{
	stringstream xml;

	xml << "<FitFunction>" << endl;
	xml << "\t" << "<FunctionName>" << functionName << "</FunctionName>" << endl;
	if(Threads>0) xml << "\t" << "<Threads>" << Threads << "</Threads>" << endl;
	xml << "\t" << "<SetIntegratorTest>";
	if( testIntegrator == true ) xml << "True</SetIntegratorTest>" << endl;
	if( testIntegrator == false ) xml << "False</SetIntegratorTest>" << endl;
	if( !Strategy.empty() ) xml << "<Strategy>" << Strategy << "</Strategy>" << endl;
	xml << "</FitFunction>" << endl;

	return xml.str();
}

