/**
  @class DoubleResolutionModel

  A class for holding a sliced propertime acceptance

  @author Pete Clarke
  @data 2011-06-07
  */


#include "DoubleResolutionModel.h"
#include "StringProcessing.h"
#include "Mathematics.h"


#include <stdio.h>
#include <vector>
#include <string>

using namespace::std;


//............................................
// Constructor 
DoubleResolutionModel::DoubleResolutionModel( PDFConfigurator* configurator, bool quiet ) :
	resScaleName		( configurator->getName( "timeResolutionScale" ) ),
	resScale2Name            ( configurator->getName( "timeResolutionScale2" ) ),
	timeResFracName            ( configurator->getName( "timeResFraction" ) ),
	eventResolutionName	( configurator->getName( "eventResolution" ) ),
	numberComponents( 2 ), wantedComponent( 1 ), resFrac( 1. )
{
	if( !quiet) cout << "DoubleResolutionModel:: Instance created " << endl ;
}


//..........................
//This method allows the instance to add the parameters it needs to the list
void DoubleResolutionModel::addParameters( vector<string> & parameterNames )
{
	parameterNames.push_back( timeResFracName );
	parameterNames.push_back( resScaleName );
	parameterNames.push_back( resScale2Name );
	return;
}

//..........................
//To take the current value of a parameter into the instance
void DoubleResolutionModel::setParameters( ParameterSet & parameters )
{
	resFrac = parameters.GetPhysicsParameter( timeResFracName )->GetValue();
	resScale = parameters.GetPhysicsParameter( resScaleName )->GetValue();
	resScale2 = parameters.GetPhysicsParameter( resScale2Name )->GetValue();
	return;
}

//..........................
//This method allows the instance to add the specific observables it needs to the list
void DoubleResolutionModel::addObservables( vector<string> & observableNames )
{
	observableNames.push_back( eventResolutionName );
	return;
}

//..........................
//To take the current value of an obserable into the instance
void DoubleResolutionModel::setObservables( DataPoint * measurement )
{
	eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();
	return;
}

//..........................
//To take the current value of an obserable into the instance
bool DoubleResolutionModel::isPerEvent( ) {  return true; }

//..............................
// Primitive Functions
double DoubleResolutionModel::Exp( double time, double gamma ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::Exp( time, gamma, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::Exp( time, gamma, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double DoubleResolutionModel::ExpInt( double tlow, double thigh, double gamma ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpInt( tlow, thigh, gamma, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpInt( tlow, thigh, gamma, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double DoubleResolutionModel::ExpSin( double time, double gamma, double dms ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpSin( time, gamma, dms, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpSin( time, gamma, dms, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double DoubleResolutionModel::ExpSinInt( double tlow, double thigh, double gamma, double dms ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpSinInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpSinInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double DoubleResolutionModel::ExpCos( double time, double gamma, double dms ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpCos( time, gamma, dms, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpCos( time, gamma, dms, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double DoubleResolutionModel::ExpCosInt( double tlow, double thigh, double gamma, double dms ) {
	//cout << " tlow" << tlow << "   thigh  "  << thigh << "   gamma  "  << gamma << "  dms  "  << dms  << "    res  " << eventResolution*resScale << endl;
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpCosInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpCosInt( tlow, thigh, gamma, dms, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double DoubleResolutionModel::GetThisScale()
{
        double thisRes = eventResolution;
        if( wantedComponent == 1 ) { eventResolution *= resScale; }
        else if( wantedComponent == 2 ) { eventResolution *= resScale2; }
        else { eventResolution *= resScale; }
	return thisRes;
}

unsigned int DoubleResolutionModel::numComponents()
{
	return numberComponents;
}

void DoubleResolutionModel::requestComponent( unsigned int wanted )
{
	wantedComponent = wanted;
	return;
}

double DoubleResolutionModel::GetFraction( unsigned int input )
{
	if( input == 1 )
	{
		return resFrac;
	}
	else if( input == 2 )
	{
		return (1.-resFrac);
	}
	else return 0.;
}

