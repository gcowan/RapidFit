/**
  @class PerEventResModelWithOffset

  A class for holding a sliced propertime acceptance

  @author Pete Clarke
  @data 2011-06-07
  */


#include "PerEventResModelWithOffset.h"
#include "StringProcessing.h"
#include "Mathematics.h"

#include <stdio.h>
#include <vector>
#include <string>

using namespace::std;

RESMODEL_CREATOR( PerEventResModelWithOffset );

//............................................
// Constructor 
PerEventResModelWithOffset::PerEventResModelWithOffset( PDFConfigurator* configurator, bool quiet ) :
	resScaleName		( configurator->getName("timeResolutionScale") ),
	resOffsetName		( configurator->getName("timeResolutionOffset") ),
	eventResolutionName	( configurator->getName("eventResolution") ),
	numberComponents( 1 ), isCacheValid(false)
{
	if( !quiet) cout << "PerEventResModelWithOffset:: Instance created " << endl ;
}


//..........................
//This method allows the instance to add the parameters it needs to the list
void PerEventResModelWithOffset::addParameters( vector<string> & parameterNames )
{

	parameterNames.push_back( resScaleName );
	parameterNames.push_back( resOffsetName );
	return;
}

//..........................
//To take the current value of a parameter into the instance
void PerEventResModelWithOffset::setParameters( ParameterSet & parameters )
{
	isCacheValid =(resScale - parameters.GetPhysicsParameter( resScaleName )->GetValue()) <= numeric_limits<double>::epsilon();
        isCacheValid = isCacheValid && ((resOffset - parameters.GetPhysicsParameter( resOffsetName )->GetValue()) <= numeric_limits<double>::epsilon());
	resScale = parameters.GetPhysicsParameter( resScaleName )->GetValue();
	resOffset = parameters.GetPhysicsParameter( resOffsetName )->GetValue();
	return;
}

//..........................
//This method allows the instance to add the specific observables it needs to the list
void PerEventResModelWithOffset::addObservables( vector<string> & observableNames )
{
	observableNames.push_back( eventResolutionName );
	return;
}

//..........................
//To take the current value of an obserable into the instance
void PerEventResModelWithOffset::setObservables( DataPoint * measurement )
{
	eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();
	return;
}

//..........................
//To take the current value of an obserable into the instance
bool PerEventResModelWithOffset::isPerEvent( ) {  return true ; }


//..............................
// Primitive Functions
double PerEventResModelWithOffset::Exp( double time, double gamma ) {
	return Mathematics::Exp( time, gamma, (eventResolution*resScale)+resOffset  ) ;
}

double PerEventResModelWithOffset::ExpInt( double tlow, double thigh, double gamma ) {
	return Mathematics::ExpInt( tlow, thigh, gamma, (eventResolution*resScale)+resOffset ) ;
}

double PerEventResModelWithOffset::ExpSin( double time, double gamma, double dms ) {
	return Mathematics::ExpSin( time, gamma, dms, (eventResolution*resScale)+resOffset) ;
}
double PerEventResModelWithOffset::ExpSinInt( double tlow, double thigh, double gamma, double dms ) {
	return Mathematics::ExpSinInt( tlow, thigh, gamma, dms, (eventResolution*resScale)+resOffset ) ;
}

double PerEventResModelWithOffset::ExpCos( double time, double gamma, double dms ) {
	return Mathematics::ExpCos( time, gamma, dms, (eventResolution*resScale)+resOffset) ;
}
double PerEventResModelWithOffset::ExpCosInt( double tlow, double thigh, double gamma, double dms ) {
	//cout << " tlow" << tlow << "   thigh  "  << thigh << "   gamma  "  << gamma << "  dms  "  << dms  << "    res  " << eventResolution*resScale << endl;
	return Mathematics::ExpCosInt( tlow, thigh, gamma, dms, (eventResolution*resScale)+resOffset) ;
}

unsigned int PerEventResModelWithOffset::numComponents()
{
	return numberComponents;
}

void PerEventResModelWithOffset::requestComponent( unsigned int wanted )
{
	(void) wanted;
	return;
}

double PerEventResModelWithOffset::GetFraction( unsigned int input )
{
	(void) input;
	return 1.;
}

bool PerEventResModelWithOffset::CacheValid() const
{
	return isCacheValid;
}

pair<double,double> PerEventResModelWithOffset::ExpCosSin( double time, double gamma, double dms )
{
	return Mathematics::ExpCosSin( time, gamma, dms, (eventResolution*resScale)+resOffset);
}

pair<double,double> PerEventResModelWithOffset::ExpCosSinInt( double tlow, double thigh, double gamma, double dms )
{
	return Mathematics::ExpCosSinInt( tlow, thigh, gamma, dms, (eventResolution*resScale)+resOffset);
}

