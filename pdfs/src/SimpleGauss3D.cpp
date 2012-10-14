
#include "SimpleGauss3D.h"
#include "Mathematics.h"
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

PDF_CREATOR( SimpleGauss3D );

SimpleGauss3D::SimpleGauss3D( PDFConfigurator* config ) :
  xName( "x" ), yName( "y" ), zName( "z" ),
  sigmaName( "sigma" ), sigma2Name( "sigma2" ), sigma3Name( "sigma3" )
{
  this->MakePrototypes();
}

void SimpleGauss3D::MakePrototypes()
{
  allObservables.push_back( string(xName) );
  allObservables.push_back( string(yName) );
  allObservables.push_back( string(zName) );
  vector<string> Params;
  Params.push_back( sigmaName ); Params.push_back( sigma2Name );
  Params.push_back( sigma3Name );
  allParameters = ParameterSet( Params );
}

SimpleGauss3D::~SimpleGauss3D()
{}

double SimpleGauss3D::Evaluate( DataPoint* input )
{
  double sigmaVal = allParameters.GetPhysicsParameter( sigmaName )->GetValue();
  double sigma2Val = allParameters.GetPhysicsParameter( sigma2Name )->GetValue();
  double sigma3Val = allParameters.GetPhysicsParameter( sigma3Name )->GetValue();
  double xVal = input->GetObservable( xName )->GetValue();
  double yVal = input->GetObservable( yName )->GetValue();
  double zVal = input->GetObservable( zName )->GetValue();
  double numerator = exp(-(xVal*xVal)/(2.*sigmaVal*sigmaVal));
  numerator *= exp(-(yVal*yVal)/(2.*sigma2Val*sigma2Val));
  numerator *= exp(-(zVal*zVal)/(2.*sigma3Val*sigma3Val));
  if( numerator <= 0. ) return 1E-10;
  else return numerator;
}

double SimpleGauss3D::Normalisation( PhaseSpaceBoundary* range )
{
  double sigmaVal = allParameters.GetPhysicsParameter( sigmaName )->GetValue();
  double sigma2Val = allParameters.GetPhysicsParameter( sigma2Name )->GetValue();
  double sigma3Val = allParameters.GetPhysicsParameter( sigma3Name )->GetValue();
  double denominator =(sigmaVal*sqrt(2.*Mathematics::Pi()));
  denominator *= (sigma2Val*sqrt(2.*Mathematics::Pi()));
  denominator *= (sigma3Val*sqrt(2.*Mathematics::Pi()));
  return denominator;
}

