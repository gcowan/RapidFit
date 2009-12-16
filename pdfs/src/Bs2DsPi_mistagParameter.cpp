/** @class Bs2DsPi_mistagParameter Bs2DsPi_mistagParameter.cpp
 *
 *  RapidFit PDF for Bs2DsPi
 *
 * Modified by Pete to add mistag as a fit parameter and time timeRes
 *
 *  @author Gemma Fardell
 */

#include "Bs2DsPi_mistagParameter.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

//Constructor
Bs2DsPi_mistagParameter::Bs2DsPi_mistagParameter() : 

// Physics parameters
	  gammaName     ( "gamma" )
	, deltaGammaName( "deltaGamma" )
	, deltaMName    ( "deltaM")
    , mistagName	( "mistag" )
    , timeresName	( "timeres" )

	// Observables
	, timeName	( "time" )
	, tagName	( "tag" )

	//, normalisationCacheValid(false)
	//, evaluationCacheValid(false)
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2DsPi_mistagParameter::MakePrototypes()
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
Bs2DsPi_mistagParameter::~Bs2DsPi_mistagParameter()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bs2DsPi_mistagParameter::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	//normalisationCacheValid = false;
	return allParameters.SetPhysicsParameters(NewParameterSet);
}

//Calculate the function value
double Bs2DsPi_mistagParameter::Evaluate(DataPoint * measurement)
{
	// Get physics parameters and observables
	getPhysicsParameters( );
	getObservables( measurement ) ;
				
	double D  = 1.0 - 2.0 * mistag;
  	
	return (0.25 * ( expL() + expH() + tag * 2.0 * expCos() * D ) ); //Normalisation from dunietz
}


double Bs2DsPi_mistagParameter::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
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
	
	
//	double expLTInt, expHTInt, t3Int;
//	getTimeDependentFuncsInt(  expLTInt, expHTInt, t3Int, boundary );	
		

	double D  = 1.0 - 2.0 * mistag;
  	
	return (0.25 * ( expHint() + expLint() + tag * 2.0 * expCosInt() * D ) ); //Normalisation from dunietz
}

vector<string> Bs2DsPi_mistagParameter::GetDoNotIntegrateList()
{
	vector<string> returnList;
	return returnList;
}



void Bs2DsPi_mistagParameter::getTimeDependentFuncsInt(  double & expLTInt, double & expHTInt, double & t3Int, PhaseSpaceBoundary * boundary)
{


        		
	double gamma_l = gamma + deltaGamma / 2.;
	double gamma_h = gamma - deltaGamma / 2.;
	
	expLTInt = 1.0 / gamma_l;	
	expHTInt = 1.0 / gamma_h;	
	t3Int = gamma / (gamma*gamma + deltaM*deltaM);
									
	return;
}


void Bs2DsPi_mistagParameter::getPhysicsParameters( )
{
	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
    deltaGamma = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	deltaM     = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
	mistag     = allParameters.GetPhysicsParameter( mistagName )->GetValue();
	timeRes    = allParameters.GetPhysicsParameter( timeresName )->GetValue();
	
	return;
}


void Bs2DsPi_mistagParameter::getObservables( DataPoint* measurement)
{
	// Observables
	time = measurement->GetObservable( timeName )->GetValue();
	tag =  measurement->GetObservable( tagName )->GetValue();
	
	return;
}


//-----------------------------------------------------------------
// Time primitives including single gaussian timeRes
//
// Much of the single gausian resolutino code is copied from Yue Hongs pdf

double Bs2DsPi_mistagParameter::expL() const 
{
    if(timeRes>0.) {     
		
		double theExp = exp( -time*gamma_l() + timeRes*timeRes * gamma_l()*gamma_l() / 2. ) ;
		double theErfc = RooMath::erfc(  -( time - timeRes*timeRes*gamma_l() ) /sqrt(2.)/timeRes )  ;
		return theExp * theErfc  / 2.0 ;
		
		// Yue hongs code
		//double c = gamma_l() * timeRes /sqrt(2.); 
		//double u = t / timeRes / sqrt(2.);
		//return exp( c*c - gamma_l()*t ) * RooMath::erfc(c-u) / 2.;
    }	
    else {
		if( time < 0.0 ) return 0.0 ;
		return exp( -gamma_l() * time ) ;
	}
}

double Bs2DsPi_mistagParameter::expH() const 
{
    if(timeRes>0.) {     
		
		double theExp = exp( -time*gamma_h() + timeRes*timeRes * gamma_h()*gamma_h() / 2. ) ;
		double theErfc = RooMath::erfc(  -( time - timeRes*timeRes*gamma_h() ) /sqrt(2.)/timeRes )  ;
		return theExp * theErfc  / 2.0 ;
		
		//Yue Hongs code
		//double c = gamma_h() * timeRes /sqrt(2.); 
		//double u = t / timeRes / sqrt(2.);
		//return exp( c*c - gamma_h()*t ) * RooMath::erfc(c-u) / 2.;
    }	
    else {
		if( time < 0.0 ) return 0.0 ;
		return exp( -gamma_h() * time ) ;
	}
}

double Bs2DsPi_mistagParameter::expCos() const 
{
    if( timeRes > 0. ) {
		
		//double theExp = exp( -t*gambar() + timeRes*timeRes * ( gambar()*gambar() - deltaM*deltaM ) / 2. ) ;
		//double theCos = cos( deltaM * ( t - timeRes*timeRes*gambar() ) ) ;
		//double theSin = sin( deltaM * ( t - timeRes*timeRes*gambar() ) ) ;
		//RooComplex z( -( t - timeRes*timeRes*gambar() )/sqrt(2.)/timeRes,  - timeRes*deltaM/sqrt(2.) ) ;
		//double theReErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncRe(z) ) : ( 1.0 - RooMath::FastComplexErrFuncRe(z) );
		//double theImErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncIm(z) ) : ( 1.0 - RooMath::FastComplexErrFuncIm(z) );
		//return theExp * ( theCos*theReErfc - theSin*theImErfc ) / 2.0 ;
		
		//Yue Hongs code
		double c = gambar() * timeRes/sqrt(2.);
		double u = time / timeRes / sqrt(2.) ;
		double wt = deltaM / gambar() ;
		return ( evalCerfRe(wt,-u,c) + evalCerfRe(-wt,-u,c) ) / 4. ;
	}
	else {
		if( time < 0.0 ) return 0.0 ;
		return exp( -gambar() *time ) * cos( deltaM * time )  ;
	}
	
}

//-----------------------------------------------------------
// Time integrals 

double Bs2DsPi_mistagParameter::expLint( ) const {

    if( timeRes > 0. ) {
		// This is a placeholder as I havnt put the correct code in yet as i dont know it.
		// So it only works if time limits are large and start from < 0
		if( ( tlow > -2 ) || ( thigh < 10 ) ) {
			std::cerr << " Bs2DsPi_mistagParameter cannot handle tlow > -2 or thigh < 10  with resolution on" << std::endl ;
			return -999.0 ;				
		}
		return (1/gamma_l()) * ( 1.0 - exp(-gamma_l()*thigh) ) ;

	}
	else {
		if( tlow <= 0. ) return (1/gamma_l()) * ( 1.0 - exp(-gamma_l()*thigh) ) ;
		else return (1/gamma_l()) * ( exp(-gamma_l()*tlow) - exp(-gamma_l()*thigh) ) ;
	}
}

double Bs2DsPi_mistagParameter::expHint( ) const {
	
    if( timeRes > 0. ) {
		// This is a placeholder as I havnt put the correct code in yet as i dont know it.
		// So it only works if time limits are large and start from < 0
		if( ( tlow > -2 ) || ( thigh < 10 ) ) {
			std::cerr << " Bs2DsPi_mistagParameter cannot handle tlow > -2 or thigh < 10  with resolution on" << std::endl ;
			return -999.0 ;				
		}
		return (1/gamma_h()) * ( 1.0 - exp(-gamma_h()*thigh) ) ;
		
	}
	else {
		if( tlow <= 0. ) return (1/gamma_h()) * ( 1.0 - exp(-gamma_h()*thigh) ) ;
		else return (1/gamma_h()) * ( exp(-gamma_h()*tlow) - exp(-gamma_h()*thigh) ) ;
	}
}




double Bs2DsPi_mistagParameter::expCosInt(  ) const {

	double tl, th ;
	
	tl = tlow ;
	th = thigh ;
	
	if( timeRes > 0. ) {
		// This is a placeholder as I havnt put the correct code in yet as i dont know it.
		// So it only works if time limits are large and start from < 0
		if( ( tlow > -2 ) || ( thigh < 10 ) ) {
			std::cerr << " Bs2DsPi_mistagParameter cannot handle tlow > -2 or thigh < 10  with resolution on" << std::endl ;
			return -999.0 ;				
		}
	}
	
	if( tlow <= 0. ) tl = 0. ;

	return (1/(gambar()*gambar() + deltaM*deltaM)) * (
					( exp(-gambar()*tl)* (gambar()*cos(deltaM*tl) - deltaM*sin(deltaM*tl)))
				   -( exp(-gambar()*th)* (gambar()*cos(deltaM*th) - deltaM*sin(deltaM*th)))
			 );
}


//-------------------------------------------------------
// All of these functions were taken from Yue Hongs code

// Calculate exp(-u^2) cwerf(swt*c + i(u+c)), taking care of numerical instabilities
//RooComplex Bs2DsPi_mistagParameter::evalCerf(double swt, double u, double c) const {
//  RooComplex z(swt*c,u+c);
//  return (z.im()>-4.0) ? RooMath::FastComplexErrFunc(z)*exp(-u*u) : evalCerfApprox(swt,u,c) ;
///}  DIDNT APPEAR TO BE USED

// Calculate Re(exp(-u^2) cwerf(swt*c + i(u+c))), taking care of numerical instabilities
double Bs2DsPi_mistagParameter::evalCerfRe(double swt, double u, double c) const {
    RooComplex z(swt*c,u+c);
    return (z.im()>-4.0) ? RooMath::FastComplexErrFuncRe(z)*exp(-u*u) : evalCerfApprox(swt,u,c).re() ;
}

// Calculate Im(exp(-u^2) cwerf(swt*c + i(u+c))), taking care of numerical instabilities
double Bs2DsPi_mistagParameter::evalCerfIm(double swt, double u, double c) const {
    RooComplex z(swt*c,u+c);
    return (z.im()>-4.0) ? RooMath::FastComplexErrFuncIm(z)*exp(-u*u) : evalCerfApprox(swt,u,c).im() ;
}

// use the approximation: erf(z) = exp(-z*z)/(sqrt(pi)*z) to explicitly cancel the divergent exp(y*y) behaviour of
// CWERF for z = x + i y with large negative y
RooComplex Bs2DsPi_mistagParameter::evalCerfApprox(double swt, double u, double c) const
{
	static double rootpi= sqrt(atan2(0.,-1.));
	RooComplex z(swt*c,u+c);  
	RooComplex zc(u+c,-swt*c);
	RooComplex zsq= z*z;
	RooComplex v= -zsq - u*u;
	return v.exp()*(-zsq.exp()/(zc*rootpi) + 1)*2 ; //why shoule be a 2 here?
	//   return v.exp()*(-zsq.exp()/(zc*rootpi) + 1);
	
}

double Bs2DsPi_mistagParameter::gamma_l() const { return gamma + ( deltaGamma / 2.0 ) ; }
double Bs2DsPi_mistagParameter::gamma_h() const { return gamma - ( deltaGamma / 2.0 ) ; }
double Bs2DsPi_mistagParameter::gambar() const   { return gamma ; }
