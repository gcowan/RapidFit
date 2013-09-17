/**
  @class FixedResolutionModel

  A class for holding a sliced propertime acceptance

  @author Pete Clarke
  @data 2011-06-07
  */


#include "FixedResolutionModel.h"
#include "StringProcessing.h"
#include "Mathematics.h"


#include <stdio.h>
#include <vector>
#include <string>

using namespace::std;


//............................................
// Constructor 
FixedResolutionModel::FixedResolutionModel( PDFConfigurator* configurator, bool quiet ) :
	resScaleName		( configurator->getName( "timeResolutionScale" ) ),
	eventResolutionName	( configurator->getName( "eventResolution" ) ),
	numberComponents( 1 ), wantedComponent( 1 )
{
	if( !quiet) cout << "FixedResolutionModel:: Instance created " << endl ;
}


//..........................
//This method allows the instance to add the parameters it needs to the list
void FixedResolutionModel::addParameters( vector<string> & parameterNames )
{
	parameterNames.push_back( resScaleName );
	parameterNames.push_back( eventResolutionName );
	return;
}

//..........................
//To take the current value of a parameter into the instance
void FixedResolutionModel::setParameters( ParameterSet & parameters )
{
	eventResolution = parameters.GetPhysicsParameter( eventResolutionName )->GetValue();
	resScale = parameters.GetPhysicsParameter( resScaleName )->GetValue();
	return;
}

//..........................
//This method allows the instance to add the specific observables it needs to the list
void FixedResolutionModel::addObservables( vector<string> & observableNames )
{
	(void) observableNames;
	//observableNames.push_back( eventResolutionName );
	return;
}

//..........................
//To take the current value of an obserable into the instance
void FixedResolutionModel::setObservables( DataPoint * measurement )
{
	(void) measurement;
	//eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();
	return;
}

//..........................
//To take the current value of an obserable into the instance
bool FixedResolutionModel::isPerEvent( ) {  return false; }

//..............................
// Primitive Functions
double FixedResolutionModel::Exp( double time, double gamma ) {
	return Mathematics::Exp( time, gamma, this->GetThisScale() );
}

double FixedResolutionModel::ExpInt( double tlow, double thigh, double gamma ) {
	return Mathematics::ExpInt( tlow, thigh, gamma, this->GetThisScale() );
}

double FixedResolutionModel::ExpSin( double time, double gamma, double dms ) {
	return Mathematics::ExpSin( time, gamma, dms, this->GetThisScale() );
}
double FixedResolutionModel::ExpSinInt( double tlow, double thigh, double gamma, double dms ) {
	return Mathematics::ExpSinInt( tlow, thigh, gamma, dms, this->GetThisScale() );
}

double FixedResolutionModel::ExpCos( double time, double gamma, double dms ) {
	return Mathematics::ExpCos( time, gamma, dms, this->GetThisScale() );
}
double FixedResolutionModel::ExpCosInt( double tlow, double thigh, double gamma, double dms ) {
	//cout << " tlow" << tlow << "   thigh  "  << thigh << "   gamma  "  << gamma << "  dms  "  << dms  << "    res  " << eventResolution*resScale << endl;
	return Mathematics::ExpCosInt( tlow, thigh, gamma, dms, this->GetThisScale() );
}

double FixedResolutionModel::GetThisScale()
{
	double thisRes = eventResolution;
	thisRes *= resScale;
	return thisRes;
}

unsigned int FixedResolutionModel::numComponents()
{
	return numberComponents;
}

void FixedResolutionModel::requestComponent( unsigned int wanted )
{
	wantedComponent = wanted;
	return;
}

double FixedResolutionModel::GetFraction( unsigned int input )
{
	(void) input;
	return 1.;
}

