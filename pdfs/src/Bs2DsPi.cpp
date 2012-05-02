/** @class Bs2DsPi Bs2DsPi.cpp
 *
 *  RapidFit PDF for Bs2DsPi
 *
 *  @author Gemma Fardell
 */

#include "Bs2DsPi.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( Bs2DsPi );

//Constructor
Bs2DsPi::Bs2DsPi( PDFConfigurator* configurator ) : 
	// Physics parameters
	  gammaName     ( configurator->getName("gamma") )
	, deltaGammaName( configurator->getName("deltaGamma") )
	, deltaMName    ( configurator->getName("deltaM") )

	// Observables
	, timeName	( configurator->getName("time") )
	, tagName	( configurator->getName("tag") )
	, mistagName	( configurator->getName("mistag") )
	//, timeres	( configurator->getName("resolution") )
	//, normalisationCacheValid(false)
	//, evaluationCacheValid(false)
	, timeconstraintName( configurator->getName("time") )
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2DsPi::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( tagName );
	allObservables.push_back( mistagName );

	// Need to think about additional parameters like
	// event-by-event propertime resolution and acceptance.
	// This will require event-by-event PDF normalisation,
	// but we are already doing this for tagging.

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName );
	parameterNames.push_back( deltaGammaName );
	parameterNames.push_back( deltaMName );
	allParameters = ParameterSet(parameterNames);
}

//Destructor
Bs2DsPi::~Bs2DsPi()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bs2DsPi::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	//normalisationCacheValid = false;
	//evaluationCacheValid = false;
	return allParameters.SetPhysicsParameters(NewParameterSet);
}

//Calculate the function value
double Bs2DsPi::Evaluate(DataPoint * measurement)
{
	double expLT, expHT, t3;
	getTimeDependentFuncs(  expLT, expHT, t3, measurement );	

	// Now need to know the tag and the mistag
	int q = (int)measurement->GetObservable( tagName )->GetValue();
	double omega = measurement->GetObservable( mistagName )->GetValue();	
	double D  = 1.0 - 2.0 * omega;


	return (0.25 * ( expLT + expHT + q * 2.0 * t3 * D ) ); //Normalisation from dunietz
}


double Bs2DsPi::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	double expLTInt, expHTInt, t3Int;
	getTimeDependentFuncsInt(  expLTInt, expHTInt, t3Int, boundary );	

	// Now need to know the tag and the mistag
	int q = (int)measurement->GetObservable( tagName )->GetValue();
	double omega = measurement->GetObservable( mistagName )->GetValue();	
	double D  = 1.0 - 2.0 * omega;


	return (0.25 * ( expHTInt + expLTInt + q * 2.0 * t3Int * D ) ); //Normalisation from dunietz
}

vector<string> Bs2DsPi::GetDoNotIntegrateList()
{
	vector<string> returnList;
	returnList.push_back(mistagName);
	return returnList;
}

void Bs2DsPi::getTimeDependentFuncs(  double & expLT, double & expHT, double & t3, DataPoint * measurement)
{

	// Observable
	double time = measurement->GetObservable( timeName )->GetValue();

	double gamma, deltaGamma, deltaMs;
	getPhysicsParameters( gamma, deltaGamma, deltaMs);

	double gamma_l = gamma + deltaGamma / 2.;
	double gamma_h = gamma - deltaGamma / 2.;

	if( time < 0.0 ){
		expLT = 0.0;
		expHT = 0.0;
		t3 =0.0;
	}
	else {
		expLT = exp( -gamma_l*time );	
		expHT = exp( -gamma_h*time );

		double expGT = exp( -gamma*time );
		double cosDeltaMsT = cos( deltaMs*time );

		t3 = expGT * cosDeltaMsT;
	}

	return;
}


void Bs2DsPi::getTimeDependentFuncsInt(  double & expLTInt, double & expHTInt, double & t3Int, PhaseSpaceBoundary * boundary)
{
	double tmin = 0.;
	double tmax = 0.;
	IConstraint * timeBound = boundary->GetConstraint( timeconstraintName );
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		return;
	}
	else
	{
		tmin = timeBound->GetMinimum();
		tmax = timeBound->GetMaximum();
	}

	double gamma, deltaGamma, deltaMs;
	getPhysicsParameters( gamma, deltaGamma, deltaMs);

	double gamma_l = gamma + deltaGamma / 2.;
	double gamma_h = gamma - deltaGamma / 2.;

	expLTInt = -1./gamma_l * ( exp(-gamma_l*tmax) - exp(-gamma_l*tmin));	
	expHTInt = -1./gamma_h * ( exp(-gamma_h*tmax) - exp(-gamma_h*tmin));	
	
	double denominator = deltaMs*deltaMs + gamma*gamma;
	double tmaxTerm = exp(-gamma*tmax) * (deltaMs * sin(deltaMs*tmax) - gamma * cos(deltaMs*tmax));
	double tminTerm = exp(-gamma*tmin) * (deltaMs * sin(deltaMs*tmin) - gamma * cos(deltaMs*tmin));
	t3Int = (tmaxTerm - tminTerm)/denominator;

	return;
}


void Bs2DsPi::getPhysicsParameters( double & gamma
		, double & deltaGamma
		, double & deltaM
		)
{
	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
	deltaGamma = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	deltaM     = allParameters.GetPhysicsParameter( deltaMName )->GetValue();

	return;
}
