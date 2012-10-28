
#include "SingleGauss.h"
#include "Mathematics.h"
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

PDF_CREATOR( SingleGauss );

SingleGauss::SingleGauss( PDFConfigurator* config ) :
	xName( config->getName("x") ),
	sigmaName( config->getName("sigma") ),
	centerName( config->getName("center") ),
	centerValue(0.), sigma_den(0.), Norm(0.)
{
	this->MakePrototypes();
}

void SingleGauss::MakePrototypes()
{
	allObservables.push_back( xName );
	vector<string> allParam; allParam.push_back( sigmaName ); allParam.push_back( centerName );
	allParameters = ParameterSet( allParam );
}

SingleGauss::~SingleGauss()
{}

bool SingleGauss::SetPhysicsParameters( ParameterSet* input )
{
	allParameters = *input;
	centerValue = allParameters.GetPhysicsParameter( centerName )->GetValue();
	double sigmaValue = allParameters.GetPhysicsParameter( sigmaName )->GetValue();
	sigma_den = 1./(2.*sigmaValue*sigmaValue);
	Norm = sigmaValue*sqrt(2.*Mathematics::Pi());
	return true;
}

double SingleGauss::Evaluate( DataPoint* input )
{
	double xVal = input->GetObservable( xName )->GetValue();
	xVal -= centerValue;
	double numerator = exp( -(xVal*xVal) * sigma_den );
	return numerator;
}

double SingleGauss::Normalisation( PhaseSpaceBoundary* range )
{
	(void) range;
	return Norm;
}

