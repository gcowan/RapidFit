/**
  @class DoubleFixedResModel

  A class for holding a sliced propertime acceptance

  @author Pete Clarke
  @data 2011-06-07
  */


#include "DoubleFixedResModel.h"
#include "StringProcessing.h"
#include "Mathematics.h"


#include <stdio.h>
#include <vector>
#include <string>

using namespace::std;

RESMODEL_CREATOR( DoubleFixedResModel ); 

//............................................
// Constructor 
DoubleFixedResModel::DoubleFixedResModel( PDFConfigurator* configurator, bool quiet ) :
	Resolution1Name			( configurator->getName( "timeResolution1" ) ),
	Resolution2Name			( configurator->getName( "timeResolution2" ) ),
	Resolution2FractionName		( configurator->getName( "timeResolution2Fraction" ) ),
	ResolutionScaleName		( configurator->getName( "timeResolutionScale" ) ),
	numberComponents( 2 ), wantedComponent( 1 )
{
	if( !quiet) cout << "DoubleFixedResModel:: Instance created " << endl ;
}


//..........................
//This method allows the instance to add the parameters it needs to the list
void DoubleFixedResModel::addParameters( vector<string> & parameterNames )
{
	parameterNames.push_back( Resolution1Name );
	parameterNames.push_back( Resolution2Name );
	parameterNames.push_back( Resolution2FractionName );
	parameterNames.push_back( ResolutionScaleName );
	return;
}

//..........................
//To take the current value of a parameter into the instance
void DoubleFixedResModel::setParameters( ParameterSet & parameters )
{
	ResolutionScale = parameters.GetPhysicsParameter( ResolutionScaleName )->GetValue();
	Resolution1 = parameters.GetPhysicsParameter( Resolution1Name )->GetValue();
	Resolution2 = parameters.GetPhysicsParameter( Resolution2Name )->GetValue();
	Resolution2Fraction = parameters.GetPhysicsParameter( Resolution2FractionName )->GetValue();
	return;
}

//..........................
//This method allows the instance to add the specific observables it needs to the list
void DoubleFixedResModel::addObservables( vector<string> & observableNames )
{
	(void) observableNames;
	//observableNames.push_back( eventResolutionName );
	return;
}

//..........................
//To take the current value of an obserable into the instance
void DoubleFixedResModel::setObservables( DataPoint * measurement )
{
	(void) measurement;
	//eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();
	return;
}

//..........................
//To take the current value of an obserable into the instance
bool DoubleFixedResModel::isPerEvent( ) {  return false; }

//..............................
// Primitive Functions
double DoubleFixedResModel::Exp( double time, double gamma ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::Exp( time, gamma, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::Exp( time, gamma, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double DoubleFixedResModel::ExpInt( double tlow, double thigh, double gamma ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpInt( tlow, thigh, gamma, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpInt( tlow, thigh, gamma, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double DoubleFixedResModel::ExpSin( double time, double gamma, double dms ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpSin( time, gamma, dms, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpSin( time, gamma, dms, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}
double DoubleFixedResModel::ExpSinInt( double tlow, double thigh, double gamma, double dms ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpSinInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpSinInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double DoubleFixedResModel::ExpCos( double time, double gamma, double dms ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpCos( time, gamma, dms, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpCos( time, gamma, dms, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double DoubleFixedResModel::ExpCosInt( double tlow, double thigh, double gamma, double dms ) {
	//cout << " tlow" << tlow << "   thigh  "  << thigh << "   gamma  "  << gamma << "  dms  "  << dms  << "    res  " << eventResolution*resScale << endl;
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpCosInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpCosInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double DoubleFixedResModel::GetThisScale()
{
	double thisRes = ResolutionScale;
	if( wantedComponent == 1 ) { thisRes *= Resolution1; }
	else if( wantedComponent == 2 ) { thisRes *= Resolution2; }
	else { thisRes *= Resolution1; }
	return thisRes;
}

unsigned int DoubleFixedResModel::numComponents()
{
	return numberComponents;
}

void DoubleFixedResModel::requestComponent( unsigned int wanted )
{
	wantedComponent = wanted;
	return;
}

double DoubleFixedResModel::GetFraction( unsigned int input )
{
	if( input == 2 )
	{
		return Resolution2Fraction;
	}
	else if( input == 1 )
	{
		return (1.-Resolution2Fraction);
	}
	else return 0.;
}

