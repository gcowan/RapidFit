/**
  @class TripleFixedResModel

  A class for holding a sliced propertime acceptance

  @author Pete Clarke
  @data 2011-06-07
  */


#include "TripleFixedResModel.h"
#include "StringProcessing.h"
#include "Mathematics.h"


#include <stdio.h>
#include <vector>
#include <string>

using namespace::std;


//............................................
// Constructor 
TripleFixedResModel::TripleFixedResModel( PDFConfigurator* configurator, bool quiet ) :
	resScaleName		( configurator->getName( "timeResolutionScale" ) ),
	resScale2Name            ( configurator->getName( "timeResolutionScale2" ) ),
	resScale3Name            ( configurator->getName( "timeResolutionScale3" ) ),
	timeResFrac2Name            ( configurator->getName( "timeResFraction2" ) ),
	timeResFrac3Name            ( configurator->getName( "timeResFraction3" ) ),
	eventResolutionName	( configurator->getName( "eventResolution" ) ),
	numberComponents( 3 ), wantedComponent( 1 )
{
	if( !quiet) cout << "TripleFixedResModel:: Instance created " << endl ;
}


//..........................
//This method allows the instance to add the parameters it needs to the list
void TripleFixedResModel::addParameters( vector<string> & parameterNames )
{
	parameterNames.push_back( resScaleName );
	parameterNames.push_back( resScale2Name );
	parameterNames.push_back( resScale3Name );
	parameterNames.push_back( timeResFrac2Name );
	parameterNames.push_back( timeResFrac3Name );
	parameterNames.push_back( eventResolutionName );
	return;
}

//..........................
//To take the current value of a parameter into the instance
void TripleFixedResModel::setParameters( ParameterSet & parameters )
{
	eventResolution = parameters.GetPhysicsParameter( eventResolutionName )->GetValue();
	resScale = parameters.GetPhysicsParameter( resScaleName )->GetValue();
	resScale2 = parameters.GetPhysicsParameter( resScale2Name )->GetValue();
	resScale3 = parameters.GetPhysicsParameter( resScale3Name )->GetValue();
	resFrac2 = parameters.GetPhysicsParameter( timeResFrac2Name )->GetValue();
	resFrac3 = parameters.GetPhysicsParameter( timeResFrac3Name )->GetValue();
	return;
}

//..........................
//This method allows the instance to add the specific observables it needs to the list
void TripleFixedResModel::addObservables( vector<string> & observableNames )
{
	(void) observableNames;
	//observableNames.push_back( eventResolutionName );
	return;
}

//..........................
//To take the current value of an obserable into the instance
void TripleFixedResModel::setObservables( DataPoint * measurement )
{
	(void) measurement;
	//eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();
	return;
}

//..........................
//To take the current value of an obserable into the instance
bool TripleFixedResModel::isPerEvent( ) {  return false; }

//..............................
// Primitive Functions
double TripleFixedResModel::Exp( double time, double gamma ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::Exp( time, gamma, this->GetThisScale() ) * this->GetFraction(1);
	this->requestComponent( 2 );
	returnable += Mathematics::Exp( time, gamma, this->GetThisScale() ) * this->GetFraction(2);
	this->requestComponent( 3 );
	returnable += Mathematics::Exp( time, gamma, this->GetThisScale() ) * this->GetFraction(3);
	return returnable;
}

double TripleFixedResModel::ExpInt( double tlow, double thigh, double gamma ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpInt( tlow, thigh, gamma, this->GetThisScale() ) * this->GetFraction(1);
	this->requestComponent( 2 );
	returnable += Mathematics::ExpInt( tlow, thigh, gamma, this->GetThisScale() ) * this->GetFraction(2);
	this->requestComponent( 3 );
	returnable += Mathematics::ExpInt( tlow, thigh, gamma, this->GetThisScale() ) * this->GetFraction(3);
	return returnable;
}

double TripleFixedResModel::ExpSin( double time, double gamma, double dms ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpSin( time, gamma, dms, this->GetThisScale() ) * this->GetFraction(1);
	this->requestComponent( 2 );
	returnable += Mathematics::ExpSin( time, gamma, dms, this->GetThisScale() ) * this->GetFraction(2);
	this->requestComponent( 3 );
	returnable += Mathematics::ExpSin( time, gamma, dms, this->GetThisScale() ) * this->GetFraction(3);
	return returnable;
}

double TripleFixedResModel::ExpSinInt( double tlow, double thigh, double gamma, double dms ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpSinInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction(1);
	this->requestComponent( 2 );
	returnable += Mathematics::ExpSinInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction(2);
	this->requestComponent( 3 );
	returnable += Mathematics::ExpSinInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction(3);
	return returnable;
}

double TripleFixedResModel::ExpCos( double time, double gamma, double dms ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpCos( time, gamma, dms, this->GetThisScale() ) * this->GetFraction(1);
	this->requestComponent( 2 );
	returnable += Mathematics::ExpCos( time, gamma, dms, this->GetThisScale() ) * this->GetFraction(2);
	this->requestComponent( 3 );
	returnable += Mathematics::ExpCos( time, gamma, dms, this->GetThisScale() ) * this->GetFraction(3);
	return returnable;
}

double TripleFixedResModel::ExpCosInt( double tlow, double thigh, double gamma, double dms ) {
	//cout << " tlow" << tlow << "   thigh  "  << thigh << "   gamma  "  << gamma << "  dms  "  << dms  << "    res  " << eventResolution*resScale << endl;
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpCosInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction(1);
	this->requestComponent( 2 );
	returnable += Mathematics::ExpCosInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction(2);
	this->requestComponent( 3 );
	returnable += Mathematics::ExpCosInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction(3);
	return returnable;
}

double TripleFixedResModel::GetThisScale()
{
	double thisRes = eventResolution;
	if( wantedComponent == 1 ) { thisRes *= resScale; }
	else if( wantedComponent == 2 ) { thisRes *= resScale2; }
	else if( wantedComponent == 3 ) { thisRes *= resScale3; }
	else { thisRes *= resScale; }
	return thisRes;
}

unsigned int TripleFixedResModel::numComponents()
{
	return numberComponents;
}

void TripleFixedResModel::requestComponent( unsigned int wanted )
{
	wantedComponent = wanted;
	return;
}

double TripleFixedResModel::GetFraction( unsigned int input )
{
	if( input == 1 )
	{
		return ( 1. - (resFrac2+resFrac3) );
	}
	else if( input == 2 )
	{
		return resFrac2;
	}
	else if( input == 3 )
	{
		return resFrac3;
	}
	else return 0.;
}

