/**
  @class Phis2012ResolutionModel

  A class for holding a sliced propertime acceptance

  @author Dianne Ferguson
  @data 2013-10-17
  */


#include "Phis2012ResolutionModel.h"
#include "StringProcessing.h"
#include "Mathematics.h"


#include <stdio.h>
#include <vector>
#include <string>

using namespace::std;

RESMODEL_CREATOR( Phis2012ResolutionModel ); 

//............................................
// Constructor 
Phis2012ResolutionModel::Phis2012ResolutionModel( PDFConfigurator* configurator, bool quiet ) :
	timeResFracName            ( configurator->getName( "timeResFraction" ) ),
	eventResolutionName	( configurator->getName( "eventResolution" ) ),
	sfBarOffsetName		( configurator->getName( "sfBarOffset")),
	sfBarSlopeName		( configurator->getName( "sfBarSlope")),
	sfSigmaOffsetName	( configurator->getName( "sfSigmaOffset")),
	sfSigmaSlopeName	( configurator->getName( "sfSigmaSlope")),
	sigmaBarName		( configurator->getName( "sigmaBar")),
	muName			( configurator->getName( "timeResMu" )),
	numberComponents( 2 ), wantedComponent( 1 ), resFrac( 1. )
{
	if( !quiet) cout << "DoubleResolutionModel:: Instance created " << endl ;
}


//..........................
//This method allows the instance to add the parameters it needs to the list
void Phis2012ResolutionModel::addParameters( vector<string> & parameterNames )
{
	parameterNames.push_back( timeResFracName );
	parameterNames.push_back( sfBarOffsetName );
	parameterNames.push_back( sfBarSlopeName );
	parameterNames.push_back( sfSigmaOffsetName );
	parameterNames.push_back( sfSigmaSlopeName );
	parameterNames.push_back( sigmaBarName );
	parameterNames.push_back( muName );
	return;
}

//..........................
//To take the current value of a parameter into the instance
void Phis2012ResolutionModel::setParameters( ParameterSet & parameters )
{
	resFrac = parameters.GetPhysicsParameter( timeResFracName )->GetValue();
	sfBarOffset = parameters.GetPhysicsParameter( sfBarOffsetName)->GetValue();
	sfBarSlope = parameters.GetPhysicsParameter( sfBarSlopeName)->GetValue();
	sfSigmaOffset = parameters.GetPhysicsParameter( sfSigmaOffsetName)->GetValue();
	sfSigmaSlope = parameters.GetPhysicsParameter( sfSigmaSlopeName)->GetValue();
	sigmaBar = parameters.GetPhysicsParameter( sigmaBarName)->GetValue();
        mu = parameters.GetPhysicsParameter( muName )->GetValue();
	return;
}

//..........................
//This method allows the instance to add the specific observables it needs to the list
void Phis2012ResolutionModel::addObservables( vector<string> & observableNames )
{
	observableNames.push_back( eventResolutionName );
	return;
}

//..........................
//To take the current value of an obserable into the instance
void Phis2012ResolutionModel::setObservables( DataPoint * measurement )
{
	eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();
	return;
}

//..........................
//To take the current value of an obserable into the instance
bool Phis2012ResolutionModel::isPerEvent( ) {  return true; }

//..............................
// Primitive Functions
double Phis2012ResolutionModel::Exp( double time, double gamma ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::Exp( time-mu, gamma, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::Exp( time-mu, gamma, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double Phis2012ResolutionModel::ExpInt( double tlow, double thigh, double gamma ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpInt( tlow-mu, thigh-mu, gamma, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpInt( tlow-mu, thigh-mu, gamma, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double Phis2012ResolutionModel::ExpSin( double time, double gamma, double dms ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpSin( time-mu, gamma, dms, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpSin( time-mu, gamma, dms, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double Phis2012ResolutionModel::ExpSinInt( double tlow, double thigh, double gamma, double dms ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpSinInt( tlow-mu, thigh-mu, gamma, dms, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpSinInt( tlow-mu, thigh-mu, gamma, dms, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double Phis2012ResolutionModel::ExpCos( double time, double gamma, double dms ) {
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpCos( time-mu, gamma, dms, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpCos( time-mu, gamma, dms, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double Phis2012ResolutionModel::ExpCosInt( double tlow, double thigh, double gamma, double dms ) {
	//cout << " tlow" << tlow << "   thigh  "  << thigh << "   gamma  "  << gamma << "  dms  "  << dms  << "    res  " << eventResolution*resScale << endl;
	double returnable = 0.;
	this->requestComponent( 1 );
	returnable += Mathematics::ExpCosInt( tlow-mu, thigh-mu, gamma, dms, this->GetThisScale() ) * this->GetFraction( 1 );
	this->requestComponent( 2 );
	returnable += Mathematics::ExpCosInt( tlow-mu, thigh-mu, gamma, dms, this->GetThisScale() ) * this->GetFraction( 2 );
	return returnable;
}

double Phis2012ResolutionModel::GetThisScale()
{
        double thisRes = eventResolution;
        if( wantedComponent == 1 ) {
		sfbar = sfBarOffset + (sfBarSlope)*(eventResolution - sigmaBar);
		sfsigma = sfSigmaOffset + (sfSigmaSlope)*(eventResolution - sigmaBar);
		resScale = -1.0*sqrt(resFrac/(1.0-resFrac))*sfsigma + sfbar;
		thisRes *= resScale; 
	}
        else if( wantedComponent == 2 ) {
		sfbar = sfBarOffset + (sfBarSlope)*(eventResolution - sigmaBar);
		sfsigma = sfSigmaOffset + (sfSigmaSlope)*(eventResolution - sigmaBar);
		resScale2 =sqrt((1.0-resFrac)/resFrac)*sfsigma + sfbar;
		thisRes *= resScale2; 
	}
        else { thisRes *= resScale; }

	return thisRes;
}

unsigned int Phis2012ResolutionModel::numComponents()
{
	return numberComponents;
}

void Phis2012ResolutionModel::requestComponent( unsigned int wanted )
{
	wantedComponent = wanted;
	return;
}

double Phis2012ResolutionModel::GetFraction( unsigned int input )
{
	if( input == 2 )
	{
		return resFrac;
	}
	else if( input == 1 )
	{
		return (1.-resFrac);
	}
	else return 0.;
}

