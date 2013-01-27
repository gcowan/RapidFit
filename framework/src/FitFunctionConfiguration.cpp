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
	functionName(InputName), weightName(), hasWeight(false), wantTrace(false), TraceFileName(), traceCount(0),
	Threads(0), Strategy(), testIntegrator(true), NormaliseWeights(false), gslIntegrator(false), SingleNormaliseWeights(false), alphaName("undefined"), hasAlpha(false)
{
}

//Constructor for FitFunction with event weights
FitFunctionConfiguration::FitFunctionConfiguration( string InputName, string InputWeight ) :
	functionName(InputName), weightName(InputWeight), hasWeight(true), wantTrace(false), TraceFileName(), traceCount(0),
	Threads(0), Strategy(), testIntegrator(true), NormaliseWeights(false), gslIntegrator(false), SingleNormaliseWeights(false), alphaName("undefined"), hasAlpha(false)
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
	if( hasWeight )
	{
		theFunction->UseEventWeights(weightName);
	}

	theFunction->SetGSLIntegrator(gslIntegrator);

	if( wantTrace )
	{
		theFunction->SetupTrace( TraceFileName, traceCount );
		++traceCount;
	}

	theFunction->SetThreads( Threads );

	theFunction->SetIntegratorTest( testIntegrator );

	return theFunction;
}

void FitFunctionConfiguration::SetGSLIntegrator( bool gsl )
{
	gslIntegrator = gsl;
}

//Return whether weights are being used
bool FitFunctionConfiguration::GetWeightsWereUsed() const
{
	return hasWeight;
}

//Return weight name
string FitFunctionConfiguration::GetWeightName() const
{
	return weightName;
}

bool FitFunctionConfiguration::GetUseCustomAlpha() const
{
	return hasAlpha;
}

string FitFunctionConfiguration::GetAlphaName() const
{
	return alphaName;
}

void FitFunctionConfiguration::SetAlphaName( const string input )
{
	hasAlpha=true;
	alphaName = input;
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

string FitFunctionConfiguration::GetStrategy() const
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
	if( hasWeight == true ) xml << "<WeightName>" << weightName << "</WeightName>" << endl;
	if( hasAlpha == true ) xml << "<AlphaName>" << alphaName << "</AlphaName>" << endl;
	xml << "\t" << "<SetIntegratorTest>";
	if( testIntegrator == true ) xml << "True</SetIntegratorTest>" << endl;
	if( testIntegrator == false ) xml << "False</SetIntegratorTest>" << endl;
	if( !Strategy.empty() ) xml << "<Strategy>" << Strategy << "</Strategy>" << endl;
	xml << "</FitFunction>" << endl;

	return xml.str();
}

void FitFunctionConfiguration::SetNormaliseWeights( const bool Input )
{
	NormaliseWeights = Input;
}

void FitFunctionConfiguration::SetSingleNormaliseWeights( const bool Input )
{
	SingleNormaliseWeights = Input;
}

bool FitFunctionConfiguration::GetNormaliseWeights() const
{
	return NormaliseWeights;
}

bool FitFunctionConfiguration::GetSingleNormaliseWeights() const
{
	return SingleNormaliseWeights;
}

