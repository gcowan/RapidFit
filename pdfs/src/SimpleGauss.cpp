
#include "SimpleGauss.h"
#include "Mathematics.h"
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

PDF_CREATOR( SimpleGauss );

SimpleGauss::SimpleGauss( PDFConfigurator* config ) :
  xName( "x" ), sigmaName( "sigma" )
{
  this->MakePrototypes();
}

void SimpleGauss::MakePrototypes()
{
  allObservables.push_back( string(xName) );
  allParameters = ParameterSet( vector<string>( 1, string(sigmaName)) );
}

SimpleGauss::~SimpleGauss()
{}

double SimpleGauss::Evaluate( DataPoint* input )
{
  double sigmaVal = allParameters.GetPhysicsParameter( sigmaName )->GetValue();
  double xVal = input->GetObservable( xName )->GetValue();
  double numerator = exp(-(xVal*xVal)/(2.*sigmaVal*sigmaVal));
  return numerator;
}

double SimpleGauss::Normalisation( PhaseSpaceBoundary* range )
{
  double sigmaVal = allParameters.GetPhysicsParameter( sigmaName )->GetValue();
  double denominator = (sigmaVal*sqrt(2.*Mathematics::Pi()));
  return denominator;
}

