// $Id: Exponential.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Exponential Exponential.cpp
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-11-13
 */

#include "Exponential.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"

PDF_CREATOR( Exponential );

//Constructor
Exponential::Exponential( PDFConfigurator* configurator) :

	// Physics parameters
	tauLL1Name	( configurator->getName("tau1") )
	, sigmaLL1Name	( configurator->getName("timeResolution1") )
	, sigmaLL2Name	( configurator->getName("timeResolution2") )
	, timeResLL1FracName( configurator->getName("timeResolution1Fraction") )
	// Observables
	, timeName      ( configurator->getName("time") )
	//objects used in XML
	,tauLL1(), sigmaLL(), sigmaLL1(), sigmaLL2(), timeResLL1Frac(), tlow(), thigh(), time()
{
	MakePrototypes();
}

//Make the data point and parameter set
void Exponential::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( tauLL1Name );
	parameterNames.push_back( timeResLL1FracName );
	parameterNames.push_back( sigmaLL1Name );
	parameterNames.push_back( sigmaLL2Name );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Exponential::~Exponential()
{
}

bool Exponential::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	tauLL1      = allParameters.GetPhysicsParameter( tauLL1Name )->GetValue();
	sigmaLL1    = allParameters.GetPhysicsParameter( sigmaLL1Name )->GetValue();
	sigmaLL2    = allParameters.GetPhysicsParameter( sigmaLL2Name )->GetValue();
	timeResLL1Frac = allParameters.GetPhysicsParameter( timeResLL1FracName )->GetValue();
	return isOK;
}

//Calculate the function value
double Exponential::Evaluate(DataPoint * measurement)
{
	// Observable
	time = measurement->GetObservable( timeName )->GetValue();

	if( timeResLL1Frac >= 0.9999 )
	{
		// Set the member variable for time resolution to the first value and calculate
		sigmaLL = sigmaLL1;
		return buildPDFnumerator();
	}
	else
	{
		// Set the member variable for time resolution to the first value and calculate
		sigmaLL = sigmaLL1;
		double val1 = buildPDFnumerator();
		// Set the member variable for time resolution to the second value and calculate
		sigmaLL = sigmaLL2;
		double val2 = buildPDFnumerator();
		return timeResLL1Frac*val1 + (1. - timeResLL1Frac)*val2;
	}
}

double Exponential::buildPDFnumerator()
{
	// Sum of two exponentials, using the time resolution functions

	if( tauLL1 <= 0 ) {
		cout << " In Exponential() you gave a negative or zero lifetime for tauLL1 " << endl ;
		throw(10) ;
	}
	double val = Mathematics::Exp(time, 1./tauLL1, sigmaLL);
	return val;
}

double Exponential::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	(void) measurement;
	IConstraint * timeBound = boundary->GetConstraint( timeName );
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		return -1.;
	}
	else
	{
		tlow = timeBound->GetMinimum();
		thigh = timeBound->GetMaximum();
	}

	if( timeResLL1Frac >= 0.9999 )
	{
		// Set the member variable for time resolution to the first value and calculate
		sigmaLL = sigmaLL1;
		return buildPDFdenominator();
	}
	else
	{
		// Set the member variable for time resolution to the first value and calculate
		sigmaLL = sigmaLL1;
		double val1 = buildPDFdenominator();
		// Set the member variable for time resolution to the second value and calculate
		sigmaLL = sigmaLL2;
		double val2 = buildPDFdenominator();
		return timeResLL1Frac*val1 + (1. - timeResLL1Frac)*val2;
	}
}

double Exponential::buildPDFdenominator()
{
	// Sum of two exponentials, using the time resolution functions

	if( tauLL1 <= 0 ) {
		cout << " In Exponential() you gave a negative or zero lifetime for tauLL1 " << endl ;
		throw(10) ;
	}
	double val = Mathematics::ExpInt(tlow, thigh, 1./tauLL1, sigmaLL);
	return val;
}

