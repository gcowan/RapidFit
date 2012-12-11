
#include "TMath.h"
#include "Landau.h"
#include "Mathematics.h"
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

PDF_CREATOR( Landau );

Landau::Landau( PDFConfigurator* config ) :
	xName( config->getName("x") ),
	sigmaName( config->getName("sigma") ),
	mpvName( config->getName("mpv") ),
	CName( config->getName("C") )
{
	this->MakePrototypes();
}

void Landau::MakePrototypes()
{
	allObservables.push_back( string(xName) );
	vector<string> allParams; allParams.push_back( sigmaName );
	allParams.push_back( mpvName ); allParams.push_back( CName );
	allParameters = ParameterSet( allParams );
}

Landau::~Landau()
{}

double Landau::Evaluate( DataPoint* input )
{
	double sigmaVal = allParameters.GetPhysicsParameter( sigmaName )->GetValue();
	double mpvValue = allParameters.GetPhysicsParameter( mpvName )->GetValue();
	double CVal = allParameters.GetPhysicsParameter( CName )->GetValue();
	double xVal = input->GetObservable( xName )->GetValue();
	double numerator = CVal * TMath::Landau( xVal, mpvValue, sigmaVal, true ) / sqrt(sigmaVal);
	//double numerator = 351 *  TMath::Gaus(xVal-mpvValue,sigmaVal,sqrt(sigmaVal),true);
	if( numerator < 1E-9 ) numerator=1E-9;
	return numerator;
}

double Landau::Normalisation( PhaseSpaceBoundary* range )
{
	(void) range;
	return -1.;
}

