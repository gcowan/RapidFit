
//	PDF Header File where methods are defined
#include "OptimisedDoubleGauss.h"
//	RapidFit Headers
#include "Mathematics.h"
//	STL Headers
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

//  Advertise this PDF as being available to RapidFit
PDF_CREATOR( OptimisedDoubleGauss );


//  PDF constructor
OptimisedDoubleGauss::OptimisedDoubleGauss( PDFConfigurator* config ) :
  //  Notice that the BasePDF is initialized by this PDF but PDFs are only constructed/copied from within ClassLookUp
  //  ClassLookUp assigns this PDFConriguator to this PDF 'behind the scenes' and is available on request for derrived objects
  BasePDF(),
  //  Constructing objects here in the 'initialization' phase is recommended as it's quicker
	  xName( config->getName("x") )
	, fracName( config->getName("f_sig1") )
	, sigma1Name( config->getName("sigma") )
	, sigma2Name( config->getName("sigma2") )
	, centerName( config->getName("center") ),
  plotComponents(false), componentIndex(0), sigma1_denom(0.), sigma2_denom(0.), f(0.), f2(0.), denominator(0.), center(0.)
{
  //  Work out of we want to plot Components when we do Projections
  plotComponents = config->isTrue( "PlotComponents" );
  //  Setup the Observables and ParameterSet
  this->MakePrototypes();
}

// Initialize the Observables and ParameterSet
void OptimisedDoubleGauss::MakePrototypes()
{
  //  Setup the Observables we expect to see
  allObservables.push_back( xName );

  //  Setup the PhysicsParameters we expect to find
  vector<string> ParamNames;
  ParamNames.push_back( sigma1Name ); ParamNames.push_back( sigma2Name );
  ParamNames.push_back( fracName ); ParamNames.push_back( centerName );
  allParameters = ParameterSet( ParamNames );
}

//  Copy Constructor
//  This is NOT REQUIRED for this PDF but is here for reference with some notation
OptimisedDoubleGauss::OptimisedDoubleGauss( const OptimisedDoubleGauss& input ) :
  BasePDF( input ), xName( input.xName ), fracName( input.fracName ), sigma1Name( input.sigma1Name ), sigma2Name( input.sigma2Name ),
  plotComponents( input.plotComponents ), componentIndex( input.componentIndex ), sigma1_denom( input.sigma1_denom ),
  sigma2_denom( input.sigma2_denom ), f( input.f ), f2( input.f2 ), denominator( input.denominator ),
  centerName( input.centerName ), center( input.center )
{
  //  Copy BY HAND all objects which are initialized with the new pointer here unless you know a direct copy by reference is safe.
  //  You can either define an object as global and give copied/derrived PDFs no ability to remove the object
  //  or, you can create an explicit new instance of each object per-pdf and make sure that each PDF removes the object in the destructor
}

//  Default Destructor, we don't NEED this but it's nice to have it
OptimisedDoubleGauss::~OptimisedDoubleGauss()
{
  //  Remove ANY objects created with 'new' and belonging to this PDF here
}

//  This function is calculated ONCE per Minuit Call
//  Whatever calculations we do here can be saved and reduce the time on the CPU
//  100 Minuit calls for 10,000 DataPoints calls this ONLY 100 times, keep expensive calculations here where possible
bool OptimisedDoubleGauss::SetPhysicsParameters( ParameterSet* input )
{
  //  Change our internal ParameterSet to match the given ParameterSet from Minuit
  allParameters = *input;

  //  Read the Physics Values
  double sigma1Value = allParameters.GetPhysicsParameter( sigma1Name )->GetValue();
  double sigma2Value = allParameters.GetPhysicsParameter( sigma2Name )->GetValue();
  f = allParameters.GetPhysicsParameter( fracName )->GetValue();
  f2 = 1. - f;

  center = allParameters.GetPhysicsParameter( centerName )->GetValue();

  // Divisions are VERY expensive computationally, cache results where possible
  sigma1_denom = 1. / (2.*sigma1Value*sigma1Value);
  sigma2_denom = 1. / (2.*sigma2Value*sigma2Value);

  // Whole Denominator for this PDF only depends on the ParameterSet so can be pre-Calculated once per Minuit call
  denominator = ( f*sigma1Value + (1. - f)*sigma2Value )*sqrt(2.*Mathematics::Pi());

  //  This is only here because we don't want to break existing PDFs, nothing more
  return true;
}

//  This function is called once per DataPoint for every single call from Minuit
//  100 Minuit calls for 10,000 DataPoints = 1,000,000 calls so reduce the amount of maths in this part of your PDF!
double OptimisedDoubleGauss::Evaluate( DataPoint* input )
{
  //  Read Observables
  double xVal = input->GetObservable( xName )->GetValue() - center;

  //  PreCalculate this to save CPU
  double xVal_sq = xVal*xVal;

  // Simple PDF contains all of the calculations in these few lines
  switch( componentIndex )
  {
    //  First Gaussian
    case 1:
      return f * exp(- xVal_sq * sigma1_denom );

    //  Second Gaussian
    case 2:
      return f2 * exp(- xVal_sq * sigma2_denom );

    //  Total of all of the PDF
    default:
      return f * exp(- xVal_sq * sigma1_denom ) + f2 * exp(- xVal_sq * sigma2_denom );
  }
}

//  This function is called once per DataPoint for every single call from Minuit
//  100 Minuit calls for 10,000 DataPoints = 1,000,000 calls so reduce the amount of maths in this part of your PDF!
double OptimisedDoubleGauss::Normalisation( PhaseSpaceBoundary* range )
{
  //  Because we don't actually use range in this function gcc complains
  //  We assume that the PDF is defined in such a PhaseSpace that the range doesn't effect the normalisation.
  //  When this is NOT true the PDF is NOT valid. But you will see this in your projections
  (void) range;

  return denominator;
}

vector<string> OptimisedDoubleGauss::PDFComponents()
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

//  This function is potentially called a LOT of times so try to reduce the total amount of code here per-call
double OptimisedDoubleGauss::EvaluateComponent( DataPoint* input, ComponentRef* compRef )
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

  return this->Evaluate( input );
}

