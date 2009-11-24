/** @class RaPDF_Bs2DsPi RaPDF_Bs2DsPi.cpp
 *
 *  RapidFit PDF for Bs2DsPi
 *
 *  @author Gemma Fardell
 */

#include "RaPDF_Bs2DsPi.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

//Constructor
RaPDF_Bs2DsPi::RaPDF_Bs2DsPi() : 
	// Physics parameters
	  gammaName     ( "gamma" )
	, deltaGammaName( "deltaGamma" )
	, deltaMName    ( "deltaM")

	// Observables
	, timeName	( "time" )
	, tagName	( "tag" )
	, mistagName	( "mistag" )
	//, timeres	( "resolution" )
	//, normalisationCacheValid(false)
	//, evaluationCacheValid(false)
{
	MakePrototypes();
}

//Make the data point and parameter set
void RaPDF_Bs2DsPi::MakePrototypes()
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
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
RaPDF_Bs2DsPi::~RaPDF_Bs2DsPi()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool RaPDF_Bs2DsPi::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	//normalisationCacheValid = false;
	//evaluationCacheValid = false;
	return allParameters.SetPhysicsParameters(NewParameterSet);
}

//Calculate the function value
double RaPDF_Bs2DsPi::Evaluate(DataPoint * measurement)
{
	double expLT, expHT, t3;
	getTimeDependentFuncs(  expLT, expHT, t3, measurement );	
		
	// Now need to know the tag and the mistag
	int q = (int)measurement->GetObservable( tagName )->GetValue();
	double omega = measurement->GetObservable( mistagName )->GetValue();	
	double D  = 1.0 - 2.0 * omega;
  	
  	
	return (0.25 * ( expLT + expHT + q * 2.0 * t3 * D ) ); //Normalisation from dunietz
}


double RaPDF_Bs2DsPi::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	double expLTInt, expHTInt, t3Int;
	getTimeDependentFuncsInt(  expLTInt, expHTInt, t3Int, boundary );	
		
	// Now need to know the tag and the mistag
	int q = (int)measurement->GetObservable( tagName )->GetValue();
	double omega = measurement->GetObservable( mistagName )->GetValue();	
	double D  = 1.0 - 2.0 * omega;
  	
  	
	return (0.25 * ( expHTInt + expLTInt + q * 2.0 * t3Int * D ) ); //Normalisation from dunietz
}

vector<string> RaPDF_Bs2DsPi::GetDoNotIntegrateList()
{
	vector<string> returnList;
	returnList.push_back(mistagName);
	return returnList;
}

void RaPDF_Bs2DsPi::getTimeDependentFuncs(  double & expLT, double & expHT, double & t3, DataPoint * measurement)
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


void RaPDF_Bs2DsPi::getTimeDependentFuncsInt(  double & expLTInt, double & expHTInt, double & t3Int, PhaseSpaceBoundary * boundary)
{

	double tlow = 0.;
	double thigh = 0.;
	IConstraint * timeBound = boundary->GetConstraint("time");
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		return;
	}
	else
	{
		tlow = timeBound->GetMinimum();
		thigh = timeBound->GetMaximum();
	}

	// Observable
        
	double gamma, deltaGamma, deltaMs;
	getPhysicsParameters( gamma, deltaGamma, deltaMs);
	
	double gamma_l = gamma + deltaGamma / 2.;
	double gamma_h = gamma - deltaGamma / 2.;
	
	expLTInt = 1.0 / gamma_l;	
	expHTInt = 1.0 / gamma_h;	
	t3Int = gamma / (gamma*gamma + deltaMs*deltaMs);
									
	return;
}


void RaPDF_Bs2DsPi::getPhysicsParameters( double & gamma
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
