
#include "StudentT.h"
#include "Mathematics.h"
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

PDF_CREATOR( StudentT );

StudentT::StudentT( PDFConfigurator* config ) :
	massName( config->getName("mass") ),
	sName( config->getName("s") ),
	muName( config->getName("mu") ),
	nName( config->getName("n") ),
	sValue(0.), muValue(0.), nValue(0.)
{
	this->MakePrototypes();
}

void StudentT::MakePrototypes()
{
	allObservables.push_back( massName );
	vector<string> allParam;
    allParam.push_back( sName );
    allParam.push_back( muName );
    allParam.push_back( nName );
	allParameters = ParameterSet( allParam );
}

StudentT::~StudentT()
{}

bool StudentT::SetPhysicsParameters( ParameterSet* input )
{
	allParameters = *input;
	muValue = allParameters.GetPhysicsParameter( muName )->GetValue();
	sValue  = allParameters.GetPhysicsParameter( sName )->GetValue();
	nValue  = allParameters.GetPhysicsParameter( nName )->GetValue();
	return true;
}

double StudentT::Evaluate( DataPoint* input )
{
    double mass = input->GetObservable( massName )->GetValue();
    const double massmu2 = (mass - muValue)*(mass - muValue);
    const double ns2 = nValue*sValue*sValue;
    const double factor = pow((1 + (massmu2/ns2)), 0.5*(nValue+1));
    const double Z = sqrt(TMath::Pi()*ns2)*TMath::Gamma(0.5*nValue)/TMath::Gamma(0.5*(nValue+1));
    return 1./(factor*Z);
}

double StudentT::Normalisation( PhaseSpaceBoundary* range )
{
  (void) range;
  return 1.0;
}

