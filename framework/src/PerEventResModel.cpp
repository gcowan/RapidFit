/**
  @class PerEventResModel

  A class for holding a sliced propertime acceptance

  @author Pete Clarke
  @data 2011-06-07
  */


#include "PerEventResModel.h"
#include "StringProcessing.h"
#include "Mathematics.h"

#include <stdio.h>
#include <vector>
#include <string>

using namespace::std;


//............................................
// Constructor 
PerEventResModel::PerEventResModel( PDFConfigurator* configurator, bool quiet ) :
	resScaleName		( configurator->getName("timeResolutionScale") ),
	eventResolutionName	( configurator->getName("eventResolution") ),
	numberComponents( 1 )
{
	if( !quiet) cout << "PerEventResModel:: Instance created " << endl ;
}


//..........................
//This method allows the instance to add the parameters it needs to the list
void PerEventResModel::addParameters( vector<string> & parameterNames )
{
	parameterNames.push_back( resScaleName );
	return;
}

//..........................
//To take the current value of a parameter into the instance
void PerEventResModel::setParameters( ParameterSet & parameters )
{
	resScale = parameters.GetPhysicsParameter( resScaleName )->GetValue();
	return;
}

//..........................
//This method allows the instance to add the specific observables it needs to the list
void PerEventResModel::addObservables( vector<string> & observableNames )
{
	observableNames.push_back( eventResolutionName );
	return;
}

//..........................
//To take the current value of an obserable into the instance
void PerEventResModel::setObservables( DataPoint * measurement )
{
	eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();
	return;
}

//..........................
//To take the current value of an obserable into the instance
bool PerEventResModel::isPerEvent( ) {  return true ; }


//..............................
// Primitive Functions
double PerEventResModel::Exp( double time, double gamma ) {
	return Mathematics::Exp( time, gamma, eventResolution*resScale  ) ;
}

double PerEventResModel::ExpInt( double tlow, double thigh, double gamma ) {
	return Mathematics::ExpInt( tlow, thigh, gamma, eventResolution*resScale ) ;
}

double PerEventResModel::ExpSin( double time, double gamma, double dms ) {
	return Mathematics::ExpSin( time, gamma, dms, eventResolution*resScale) ;
}
double PerEventResModel::ExpSinInt( double tlow, double thigh, double gamma, double dms ) {
	return Mathematics::ExpSinInt( tlow, thigh, gamma, dms, eventResolution*resScale) ;
}

double PerEventResModel::ExpCos( double time, double gamma, double dms ) {
	return Mathematics::ExpCos( time, gamma, dms, eventResolution*resScale) ;
}
double PerEventResModel::ExpCosInt( double tlow, double thigh, double gamma, double dms ) {
	//cout << " tlow" << tlow << "   thigh  "  << thigh << "   gamma  "  << gamma << "  dms  "  << dms  << "    res  " << eventResolution*resScale << endl;
	return Mathematics::ExpCosInt( tlow, thigh, gamma, dms, eventResolution*resScale) ;
}

unsigned int PerEventResModel::numComponents()
{
	return numberComponents;
}

void PerEventResModel::requestComponent( unsigned int wanted )
{
	(void) wanted;
	return;
}

double PerEventResModel::GetFraction( unsigned int input )
{
	(void) input;
	return 1.;
}

