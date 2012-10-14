
#include "SimpleDoubleGauss.h"
#include "Mathematics.h"
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

PDF_CREATOR( SimpleDoubleGauss );

SimpleDoubleGauss::SimpleDoubleGauss( PDFConfigurator* config ) :
  xName( "x" ), fracName( "f_sig1" ), sigma1Name( "sigma" ), sigma2Name( "sigma2" ),
  plotComponents(false), componentIndex(0)
{
  plotComponents = config->isTrue( "PlotComponents" );
  this->MakePrototypes();
}

void SimpleDoubleGauss::MakePrototypes()
{
  allObservables.push_back( xName );
  vector<string> ParamNames;
  ParamNames.push_back( sigma1Name ); ParamNames.push_back( sigma2Name );
  ParamNames.push_back( fracName );
  allParameters = ParameterSet( ParamNames );
}

SimpleDoubleGauss::~SimpleDoubleGauss()
{
}

double SimpleDoubleGauss::Evaluate( DataPoint* input )
{
  //  Read Physics Parameters
  double sigma1Val = allParameters.GetPhysicsParameter( sigma1Name )->GetValue();
  double sigma2Val = allParameters.GetPhysicsParameter( sigma2Name )->GetValue();
  double f = allParameters.GetPhysicsParameter( fracName )->GetValue();
  //  Read Observables
  double xVal = input->GetObservable( xName )->GetValue();

  // Perform Calculation
  double numerator = f * exp(-(xVal*xVal)/(2.*sigma1Val*sigma1Val));
  numerator +=  (1.-f) * exp(-(xVal*xVal)/(2.*sigma2Val*sigma2Val));

  return numerator;
}

double SimpleDoubleGauss::Normalisation( PhaseSpaceBoundary* range )
{
  //  Read Physics Parameters
  double sigma1Val = allParameters.GetPhysicsParameter( sigma1Name )->GetValue();
  double sigma2Val = allParameters.GetPhysicsParameter( sigma2Name )->GetValue();
  double f = allParameters.GetPhysicsParameter( fracName )->GetValue();

  // Perform Calculation
  double denominator =  f * (sigma1Val*sqrt(2.*Mathematics::Pi()));
  denominator += (1. - f) * (sigma2Val*sqrt(2.*Mathematics::Pi()));

  return denominator;
}

vector<string> SimpleDoubleGauss::PDFComponents()
{
  vector<string> ComponentNames;
  //  The RapidFit convention is that component "0" corresponds to the total of the whole PDF
  ComponentNames.push_back( "0" );

  //  By Default we don't always want to plot the components to save wasting CPU.
  //  This can be configured by passing the relevant Configuration Option to the PDF in the XML
  if( plotComponents )
  {
    //  We do Want to Plot the Components this time so Add them to the list of Components this PDF knows about
    ComponentNames.push_back( "sigma1" );
    ComponentNames.push_back( "sigma2" );
  }
  return ComponentNames;
}

double SimpleDoubleGauss::EvaluateComponent( DataPoint* input, ComponentRef* compRef )
{
  //  Check to see if we already know what component we are trying to Evaluate
  //  if we DO    know what the componentIndex this is the Number is NOT -1
  //  if we DON'T know what the componentIndex is the Number is -1
  if( compRef->getComponentNumber() != -1 )
  {
    //  We know what the corresponding componentIndex is so lets use it
    componentIndex = compRef->getComponentNumber();
  }
  else
  {
      //  This is the Name of the Component that we want to plot
      string Name = compRef->getComponentName();

      if( Name == "sigma1" )
      {
        //  Found this Component, save the result for the future
        componentIndex = 1;
        compRef->setComponentNumber( 1 );
      }
      else if( Name == "sigma2" )
      {
        //  Found this Component, save the result for the future
        componentIndex = 2;
        compRef->setComponentNumber( 2 );
      }
      else
      {
        //  Either wanted component "0" or didn't know what to do
        //  As good practice provide something in the case we don't know what to do
        componentIndex = 0;
        compRef->setComponentNumber( 0 );
      }
  }

  //  componentIndex has now been set corresponding to which component we want to Evaluate
  //  This code could re-use Evaluate to reduce the code

  //  Read Physics Parameters
  double sigma1Val = allParameters.GetPhysicsParameter( sigma1Name )->GetValue();
  double sigma2Val = allParameters.GetPhysicsParameter( sigma2Name )->GetValue();
  double f = allParameters.GetPhysicsParameter( fracName )->GetValue();
  //  Read Observables
  double xVal = input->GetObservable( xName )->GetValue();

  // Perform Calculation
  double numerator = 0.;

  //  To save time in calculations only calculate what we want depending on what component we are after
  //  It is perfectly acceptable to use an if/else chain here but switch is faster

  switch( componentIndex )
  {
    case 1:
      numerator = f * exp(-(xVal*xVal)/(2.*sigma1Val*sigma1Val));
      break;
    case 2:
      numerator =  (1.-f) * exp(-(xVal*xVal)/(2.*sigma2Val*sigma2Val));
      break;
    default:
      numerator =        f * exp(-(xVal*xVal)/(2.*sigma1Val*sigma1Val));
      numerator +=  (1.-f) * exp(-(xVal*xVal)/(2.*sigma2Val*sigma2Val));
      break;
  }

  return numerator;
}

