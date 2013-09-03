/**
  @class ResolutionModel

  A class for holding a sliced propertime acceptance

  @author Pete Clarke
  @data 2011-06-07
  */


#include "ResolutionModel.h"
#include "StringProcessing.h"
#include "Mathematics.h"


#include <stdio.h>
#include <vector>
#include <string>

using namespace::std;


//............................................
// Constructor 
ResolutionModel::ResolutionModel( PDFConfigurator* configurator, bool quiet ) :
    resScaleName		( configurator->getName("timeResolutionScale") ),
    eventResolutionName	( configurator->getName("eventResolution") )
{
    if( !quiet) cout << "ResolutionModel:: Instance created " << endl ;
}


//..........................
//This method allows the instance to add the parameters it needs to the list
void ResolutionModel::addParameters( vector<string> & parameterNames )
{
    parameterNames.push_back( resScaleName );
    return ;
}

//..........................
//To take the current value of a parameter into the instance
void ResolutionModel::setParameters( ParameterSet & parameters )
{
	resScale = parameters.GetPhysicsParameter( resScaleName )->GetValue();
    return ;
}

//..........................
//This method allows the instance to add the specific observables it needs to the list
void ResolutionModel::addObservables( vector<string> & observableNames )
{
    observableNames.push_back( eventResolutionName );
    return ;
}

//..........................
//To take the current value of an obserable into the instance
void ResolutionModel::setObservables( DataPoint * measurement )
{
    eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();
    return ;
}

//..........................
//To take the current value of an obserable into the instance
bool ResolutionModel::isPerEvent( ) {  return true ; }


//..............................
// Primitive Functions
double ResolutionModel::Exp( double time, double gamma ) {
    return Mathematics::Exp( time, gamma, eventResolution*resScale  ) ;
}

double ResolutionModel::ExpInt( double tlow, double thigh, double gamma ) {
    return Mathematics::ExpInt( tlow, thigh, gamma, eventResolution*resScale ) ;
}

double ResolutionModel::ExpSin( double time, double gamma, double dms ) {
    return Mathematics::ExpSin( time, gamma, dms, eventResolution*resScale) ;
}
double ResolutionModel::ExpSinInt( double tlow, double thigh, double gamma, double dms ) {
    return Mathematics::ExpSinInt( tlow, thigh, gamma, dms, eventResolution*resScale) ;
}

double ResolutionModel::ExpCos( double time, double gamma, double dms ) {
    return Mathematics::ExpCos( time, gamma, dms, eventResolution*resScale) ;
}
double ResolutionModel::ExpCosInt( double tlow, double thigh, double gamma, double dms ) {
    //cout << " tlow" << tlow << "   thigh  "  << thigh << "   gamma  "  << gamma << "  dms  "  << dms  << "    res  " << eventResolution*resScale << endl;
    return Mathematics::ExpCosInt( tlow, thigh, gamma, dms, eventResolution*resScale) ;
}

//..............................
// Wrappers of Primitive Functions for Robs clever stuff
double ResolutionModel::Exp_Wrapper( vector<double> input ) {
    return this->Exp( input[0], input[1]  ) ;
}

double ResolutionModel::ExpInt_Wrapper( vector<double> input ) {
    return this->ExpInt( input[0], input[1], input[2] ) ;
}

double ResolutionModel::ExpSin_Wrapper( vector<double> input ) {
    return this->ExpSin( input[0], input[1], input[2] ) ;
}
double ResolutionModel::ExpSinInt_Wrapper( vector<double> input ) {
    return this->ExpSinInt( input[0], input[1], input[2], input[3] ) ;
}

double ResolutionModel::ExpCos_Wrapper( vector<double> input ) {
    return this->ExpCos( input[0], input[1], input[2] ) ;
}
double ResolutionModel::ExpCosInt_Wrapper( vector<double> input ) {
    return this->ExpCosInt( input[0], input[1], input[2], input[3] )  ;
}



