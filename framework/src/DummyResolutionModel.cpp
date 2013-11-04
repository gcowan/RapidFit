/**
  @class DummyResolutionModel

  A class for holding a sliced propertime acceptance

  @author Pete Clarke
  @data 2011-06-07
  */


#include "DummyResolutionModel.h"
#include "StringProcessing.h"
#include "Mathematics.h"
#include "ClassLookUp.h"

#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <string>

using namespace::std;

class DummyResolutionModel;

RESMODEL_CREATOR( DummyResolutionModel )

	//............................................
	// Constructor 
DummyResolutionModel::DummyResolutionModel( PDFConfigurator* configurator, bool quiet )
{
	(void) configurator;
	if( !quiet) cout << "DummyResolutionModel:: Instance created " << endl ;
}


//..........................
//This method allows the instance to add the parameters it needs to the list
void DummyResolutionModel::addParameters( vector<string> & parameterNames )
{
	(void) parameterNames;
	return ;
}

//..........................
//To take the current value of a parameter into the instance
void DummyResolutionModel::setParameters( ParameterSet & parameters )
{
	(void) parameters;
	return ;
}

//..........................
//This method allows the instance to add the specific observables it needs to the list
void DummyResolutionModel::addObservables( vector<string> & observableNames )
{
	(void) observableNames;
	return ;
}

//..........................
//To take the current value of an obserable into the instance
void DummyResolutionModel::setObservables( DataPoint * measurement )
{
	(void) measurement;
	return;
}

//..........................
//To take the current value of an obserable into the instance
bool DummyResolutionModel::isPerEvent( ) {  return true ; }


//..............................
// Primitive Functions
double DummyResolutionModel::Exp( double time, double gamma ) {
	return Mathematics::Exp( time, gamma, 0. );
}

double DummyResolutionModel::ExpInt( double tlow, double thigh, double gamma ) {
	return Mathematics::ExpInt( tlow, thigh, gamma, 0. );
}

double DummyResolutionModel::ExpSin( double time, double gamma, double dms ) {
	return Mathematics::ExpSin( time, gamma, dms, 0. );
}
double DummyResolutionModel::ExpSinInt( double tlow, double thigh, double gamma, double dms ) {
	return Mathematics::ExpSinInt( tlow, thigh, gamma, dms, 0. );
}

double DummyResolutionModel::ExpCos( double time, double gamma, double dms ) {
	return Mathematics::ExpCos( time, gamma, dms, 0. );
}
double DummyResolutionModel::ExpCosInt( double tlow, double thigh, double gamma, double dms ) {
	return Mathematics::ExpCosInt( tlow, thigh, gamma, dms, 0. );
}

