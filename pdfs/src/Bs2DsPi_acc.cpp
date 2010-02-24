/** @class Bs2DsPi_mistagParameter Bs2DsPi_mistagParameter.cpp
 *
 *  RapidFit PDF for Bs2DsPi
 *
 * Modified by Pete to add mistag as a fit parameter and time timeRes
 *
 *  @author Gemma Fardell
 */

#include "Bs2DsPi_acc.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"

//Constructor
Bs2DsPi_acc::Bs2DsPi_acc() : 

// Physics parameters
	  gammaName     ( "gamma" )
	, deltaGammaName( "deltaGamma" )
	, deltaMName    ( "deltaM")
    , mistagName	( "mistag" )
    , timeresName	( "timeresDsPi" )

	// Observables
	, timeName	( "time" )
	, tagName	( "tag" )

	//, normalisationCacheValid(false)
	//, evaluationCacheValid(false)
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2DsPi_acc::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( tagName );
   
	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName );
	parameterNames.push_back( deltaGammaName );
	parameterNames.push_back( deltaMName );
	parameterNames.push_back( mistagName );
	parameterNames.push_back( timeresName );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Bs2DsPi_acc::~Bs2DsPi_acc()
{
}

//.....................................................................
//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bs2DsPi_acc::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	//normalisationCacheValid = false;
	return allParameters.SetPhysicsParameters(NewParameterSet);
}

//...............................................................
//Calculate the PDF value

double Bs2DsPi_acc::Evaluate(DataPoint * measurement)
{
	// Get physics parameters and observables
	getPhysicsParameters( );
	getObservables( measurement ) ;
				
  	double D  = 1.0 - 2.0 * mistag;
  		
	return (0.25 * acc(time) * ( expL() + expH() + tag * 2.0 * expCos() * D ) ); //Normalisation from dunietz
}




//.....................................................................
//Calculate normalisation of PDF

double Bs2DsPi_acc::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	// Get physics parameters and observables
	getPhysicsParameters( );
	getObservables( measurement ) ;

	// Get time integration boundaries
	IConstraint * timeBound = boundary->GetConstraint("time");
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		return -999.;
	}
	else
	{
		tlow = timeBound->GetMinimum();
		thigh = timeBound->GetMaximum();
		if( thigh <  tlow ) return -999.0 ;
	}
	
	double D  = 1.0 - 2.0 * mistag;
  	
  	return -1; //can't deal with acceptance yet
	//return (0.25 * ( expHint() + expLint() + tag * 2.0 * expCosInt() * D ) ); //Normalisation from dunietz
}

//.......................................................
vector<string> Bs2DsPi_acc::GetDoNotIntegrateList()
{
	vector<string> returnList;
	return returnList;
}



//........................................................
// Get parameters

void Bs2DsPi_acc::getPhysicsParameters( )
{
	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
    deltaGamma = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	deltaM     = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
	mistag     = allParameters.GetPhysicsParameter( mistagName )->GetValue();
	timeRes    = allParameters.GetPhysicsParameter( timeresName )->GetValue();
	
	return;
}

//................................................................
//Get observables

void Bs2DsPi_acc::getObservables( DataPoint* measurement)
{
	// Observables
	time = measurement->GetObservable( timeName )->GetValue();
	tag =  measurement->GetObservable( tagName )->GetValue();
	
	return;
}


//...................................................................
// Adds acceptance function

double Bs2DsPi_acc::acc(double t) const 
{
	double acc_s = 1.7;
	double acc_o = 0;
	double acc_p = 2.2;

	return pow(((t-acc_o)*acc_s),acc_p)/(1.0+ pow(((t-acc_o)*acc_s),acc_p));
}


//...................................................................
// Interface to time primitives including single gaussian timeRes
//

double Bs2DsPi_acc::expL() const 
{
	return Mathematics::Exp( time, gamma_l(), timeRes ) ;
}

double Bs2DsPi_acc::expLint( ) const 
{
	return Mathematics::ExpInt( tlow, thigh, gamma_l(), timeRes ) ;
}

double Bs2DsPi_acc::expH() const 
{
	return Mathematics::Exp( time, gamma_h(), timeRes ) ;
}

double Bs2DsPi_acc::expHint( ) const 
{
	return Mathematics::ExpInt( tlow, thigh, gamma_h(), timeRes ) ;
}


double Bs2DsPi_acc::expCos() const 
{
	return Mathematics::ExpCos( time, gambar(), deltaM, timeRes ) ;
}

double Bs2DsPi_acc::expCosInt() const 
{
	return Mathematics::ExpCosInt( tlow, thigh, gambar(), deltaM, timeRes ) ;
}

//.............................................................
//Helpers

double Bs2DsPi_acc::gamma_l() const { return gamma + ( deltaGamma / 2.0 ) ; }
double Bs2DsPi_acc::gamma_h() const { return gamma - ( deltaGamma / 2.0 ) ; }
double Bs2DsPi_acc::gambar() const   { return gamma ; }
