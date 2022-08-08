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

RESMODEL_CREATOR( TripleFixedResModel );
//............................................
// Constructor 
TripleFixedResModel::TripleFixedResModel( PDFConfigurator* configurator, bool quiet ) :
	Resolution1Name			( configurator->getName( "timeResolution1" ) ),
	Resolution2Name			( configurator->getName( "timeResolution2" ) ),
	Resolution3Name			( configurator->getName( "timeResolution3" ) ),
	Resolution2FractionName		( configurator->getName( "timeResolution2Fraction" ) ),
	Resolution3FractionName		( configurator->getName( "timeResolution3Fraction" ) ),
	ResolutionScaleName		( configurator->getName( "timeResolutionScale" ) ),
	numberComponents( 3 ), wantedComponent( 1 ), isCacheValid(false)
{
	if( !quiet) cout << "TripleFixedResModel:: Instance created " << endl ;
}


//..........................
//This method allows the instance to add the parameters it needs to the list
void TripleFixedResModel::addParameters( vector<string> & parameterNames )
{
	parameterNames.push_back( Resolution1Name );
	parameterNames.push_back( Resolution2Name );
	parameterNames.push_back( Resolution3Name );
	parameterNames.push_back( Resolution2FractionName );
	parameterNames.push_back( Resolution3FractionName );
	parameterNames.push_back( ResolutionScaleName );
	return;
}

//..........................
//To take the current value of a parameter into the instance
void TripleFixedResModel::setParameters( ParameterSet & parameters )
{
	if( parameters.GetPhysicsParameter( ResolutionScaleName )->isFixed() && 
		parameters.GetPhysicsParameter( Resolution1Name )->isFixed() &&
		parameters.GetPhysicsParameter( Resolution2Name )->isFixed() &&
		parameters.GetPhysicsParameter( Resolution3Name )->isFixed() &&
		parameters.GetPhysicsParameter( Resolution2FractionName )->isFixed() &&
		parameters.GetPhysicsParameter( Resolution3FractionName )->isFixed() )
	{
		isCacheValid = true;
	}
	else
	{
		isCacheValid = ( abs(ResolutionScale - parameters.GetPhysicsParameter( ResolutionScaleName )->GetValue()) <= numeric_limits<float>::epsilon() );
		isCacheValid = isCacheValid && ( abs(Resolution1 - parameters.GetPhysicsParameter( Resolution1Name )->GetValue()) <= numeric_limits<float>::epsilon() );
		isCacheValid = isCacheValid && ( abs(Resolution2 - parameters.GetPhysicsParameter( Resolution2Name )->GetValue()) <= numeric_limits<float>::epsilon() );
		isCacheValid = isCacheValid && ( abs(Resolution3 - parameters.GetPhysicsParameter( Resolution3Name )->GetValue()) <= numeric_limits<float>::epsilon() );
		isCacheValid = isCacheValid && ( abs(Resolution2Fraction - parameters.GetPhysicsParameter( Resolution2FractionName )->GetValue()) <= numeric_limits<float>::epsilon() );
		isCacheValid = isCacheValid && ( abs(Resolution3Fraction - parameters.GetPhysicsParameter( Resolution3FractionName )->GetValue()) <= numeric_limits<float>::epsilon() );
	}

	ResolutionScale = parameters.GetPhysicsParameter( ResolutionScaleName )->GetValue();
	Resolution1 = parameters.GetPhysicsParameter( Resolution1Name )->GetValue();
	Resolution2 = parameters.GetPhysicsParameter( Resolution2Name )->GetValue();
	Resolution3 = parameters.GetPhysicsParameter( Resolution3Name )->GetValue();
	Resolution2Fraction = parameters.GetPhysicsParameter( Resolution2FractionName )->GetValue();
	Resolution3Fraction = parameters.GetPhysicsParameter( Resolution3FractionName )->GetValue();
	return;
}

bool TripleFixedResModel::CacheValid() const
{
	return isCacheValid;
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
	double thisRes = ResolutionScale;
	if( wantedComponent == 1 ) { thisRes *= Resolution1; }
	else if( wantedComponent == 2 ) { thisRes *= Resolution2; }
	else if( wantedComponent == 3 ) { thisRes *= Resolution3; }
	else { thisRes *= ResolutionScale; }
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
		return ( 1. - (Resolution2Fraction+Resolution3Fraction) );
	}
	else if( input == 2 )
	{
		return Resolution2Fraction;
	}
	else if( input == 3 )
	{
		return Resolution3Fraction;
	}
	else return 0.;
}

