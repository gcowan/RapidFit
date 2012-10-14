
#include "SimpleGauss2D.h"
#include "Mathematics.h"
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

PDF_CREATOR( SimpleGauss2D );

SimpleGauss2D::SimpleGauss2D( PDFConfigurator* config ) :
  xName( "x" ), yName( "y" ), sigmaName( "sigma" ), sigma2Name( "sigma2" )
{
  this->MakePrototypes();
}

void SimpleGauss2D::MakePrototypes()
{
  allObservables.push_back( string(xName) );
  allObservables.push_back( string(yName) );
  vector<string> Params; Params.push_back( sigmaName ); Params.push_back( sigma2Name );
  allParameters = ParameterSet( Params );
}

SimpleGauss2D::~SimpleGauss2D()
{}

double SimpleGauss2D::Evaluate( DataPoint* input )
{
  double sigmaVal = allParameters.GetPhysicsParameter( sigmaName )->GetValue();
  double sigma2Val = allParameters.GetPhysicsParameter( sigma2Name )->GetValue();
  double xVal = input->GetObservable( xName )->GetValue();
  double yVal = input->GetObservable( yName )->GetValue();
  double numerator = exp(-(xVal*xVal)/(2.*sigmaVal*sigmaVal));
  numerator *= exp(-(yVal*yVal)/(2.*sigma2Val*sigma2Val));
  return numerator;
}

double SimpleGauss2D::Normalisation( PhaseSpaceBoundary* range )
{
  double sigmaVal = allParameters.GetPhysicsParameter( sigmaName )->GetValue();
  double sigma2Val = allParameters.GetPhysicsParameter( sigma2Name )->GetValue();
  double denominator =(sigmaVal*sqrt(2.*Mathematics::Pi()));
  denominator *= (sigma2Val*sqrt(2.*Mathematics::Pi()));
  return denominator;
}

