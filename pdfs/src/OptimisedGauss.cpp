
//	PDF Header File where methods are defined
#include "OptimisedGauss.h"
//	RapidFit Headers
#include "Mathematics.h"
//	STL Headers
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

//  Advertise this PDF as being available to RapidFit
PDF_CREATOR( OptimisedGauss );


//  PDF constructor
OptimisedGauss::OptimisedGauss( PDFConfigurator* config ) :
  //  Notice that the BasePDF is initialized by this PDF but PDFs are only constructed/copied from within ClassLookUp
  //  ClassLookUp assigns this PDFConriguator to this PDF 'behind the scenes' and is available on request for derrived objects
  BasePDF(),
  //  Constructing objects here in the 'initialization' phase is recommended as it's quicker
  xName( config->getName("x") ), sigmaName( config->getName("sigma") ), centreName( config->getName("centre") ),
  sigma_denom(0.), denominator(0.), centre(0.)
{
  //  Setup the Observables and ParameterSet
  this->MakePrototypes();
}

// Initialize the Observables and ParameterSet
void OptimisedGauss::MakePrototypes()
{
  //  Setup the Observables we expect to see
  allObservables.push_back( xName );

  //  Setup the PhysicsParameters we expect to find
  vector<string> ParamNames;
  ParamNames.push_back( sigmaName ); ParamNames.push_back( centreName );
  allParameters = ParameterSet( ParamNames );
}

//  Copy Constructor
//  This is NOT REQUIRED for this PDF but is here for reference with some notation
OptimisedGauss::OptimisedGauss( const OptimisedGauss& input ) :
  BasePDF( input ), xName( input.xName ), sigmaName( input.sigmaName ), sigma_denom( input.sigma_denom ),
  denominator( input.denominator ), centreName( input.centreName ), centre( input.centre )
{
  //  Copy BY HAND all objects which are initialized with the new pointer here unless you know a direct copy by reference is safe.
  //  You can either define an object as global and give copied/derrived PDFs no ability to remove the object
  //  or, you can create an explicit new instance of each object per-pdf and make sure that each PDF removes the object in the destructor
}

//  Default Destructor, we don't NEED this but it's nice to have it
OptimisedGauss::~OptimisedGauss()
{
  //  Remove ANY objects created with 'new' and belonging to this PDF here
}

//  This function is calculated ONCE per Minuit Call
//  Whatever calculations we do here can be saved and reduce the time on the CPU
//  100 Minuit calls for 10,000 DataPoints calls this ONLY 100 times, keep expensive calculations here where possible
bool OptimisedGauss::SetPhysicsParameters( ParameterSet* input )
{
  //  Change our internal ParameterSet to match the given ParameterSet from Minuit
  allParameters = *input;

  //  Read the Physics Values
  double sigmaValue = allParameters.GetPhysicsParameter( sigmaName )->GetValue();

  centre = allParameters.GetPhysicsParameter( centreName )->GetValue();

  // Divisions are VERY expensive computationally, cache results where possible
  sigma_denom = 1. / (2.*sigmaValue*sigmaValue);

  // Whole Denominator for this PDF only depends on the ParameterSet so can be pre-Calculated once per Minuit call
  denominator = sigmaValue * sqrt(2.*Mathematics::Pi());

  //  This is only here because we don't want to break existing PDFs, nothing more
  return true;
}

//  This function is called once per DataPoint for every single call from Minuit
//  100 Minuit calls for 10,000 DataPoints = 1,000,000 calls so reduce the amount of maths in this part of your PDF!
double OptimisedGauss::Evaluate( DataPoint* input )
{
  //  Read Observables
  double xVal = input->GetObservable( xName )->GetValue() - centre;

  //  PreCalculate this to save CPU
  double xVal_sq = xVal*xVal;

  return exp(- xVal_sq * sigma_denom );
}

//  This function is called once per DataPoint for every single call from Minuit
//  100 Minuit calls for 10,000 DataPoints = 1,000,000 calls so reduce the amount of maths in this part of your PDF!
double OptimisedGauss::Normalisation( PhaseSpaceBoundary* range )
{
  //  Because we don't actually use range in this function gcc complains
  //  We assume that the PDF is defined in such a PhaseSpace that the range doesn't effect the normalisation.
  //  When this is NOT true the PDF is NOT valid. But you will see this in your projections
  (void) range;

  return denominator;
}

