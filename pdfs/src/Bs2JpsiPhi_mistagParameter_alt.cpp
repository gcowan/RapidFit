// $Id: Bs2JpsiPhi_mistagParameter_alt.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2JpsiPhi_mistagParameter_alt Bs2JpsiPhi_mistagParameter_alt.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2009-12-06
 */

#include "Bs2JpsiPhi_mistagParameter_alt.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"
#include "Mathematics.h"

//......................................
//Constructor

Bs2JpsiPhi_mistagParameter_alt::Bs2JpsiPhi_mistagParameter_alt() : 
	// Physics parameters
	  gammaName     ( "gamma" )
	, deltaGammaName( "deltaGamma" )
	, deltaMName    ( "deltaM")
	, Phi_sName     ( "Phi_s")
	, Azero_sqName  ( "Azero_sq" )
	, Aperp_sqName  ( "Aperp_sq" )
	, delta_zeroName( "delta_zero" )
	, delta_paraName( "delta_para" )
	, delta_perpName( "delta_perp" )
    , mistagName	( "mistag" )
	, res1Name	( "timeResolution1" )
	, res2Name	( "timeResolution2" )
	, res1FractionName	( "timeResolution1Fraction" )
	// Observables
	, timeName	    ( "time" )
	, cosThetaName	( "cosTheta" )
	, phiName	    ( "phi" )
	, cosPsiName	( "cosPsi" )
	, tagName	    ( "tag" )
	// Other things
	, normalisationCacheValid(false)
{
	MakePrototypes();
	
	std::cout << "Constructing alternative (PELC original) J/PsiPhi classic PDF" << std::endl ;
}

//Make the data point and parameter set
void Bs2JpsiPhi_mistagParameter_alt::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( cosThetaName );
	allObservables.push_back( phiName );
	allObservables.push_back( cosPsiName );
	allObservables.push_back( tagName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName );
	parameterNames.push_back( deltaGammaName );
	parameterNames.push_back( Aperp_sqName );
	parameterNames.push_back( Azero_sqName );
	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
	parameterNames.push_back( delta_zeroName );
	parameterNames.push_back( deltaMName );
	parameterNames.push_back( Phi_sName );
	parameterNames.push_back( mistagName );
	parameterNames.push_back( res1Name );
	parameterNames.push_back( res2Name );
	parameterNames.push_back( res1FractionName );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}


//........................................................
//Destructor
Bs2JpsiPhi_mistagParameter_alt::~Bs2JpsiPhi_mistagParameter_alt()
{
}

//........................................................
//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bs2JpsiPhi_mistagParameter_alt::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	return allParameters.SetPhysicsParameters(NewParameterSet);
}

//.........................................................
//Return a list of observables not to be integrated
vector<string> Bs2JpsiPhi_mistagParameter_alt::GetDoNotIntegrateList()
{
	vector<string> list;
	return list;
}

//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2JpsiPhi_mistagParameter_alt::Evaluate(DataPoint * measurement)
{
	
    // Get parameters into member variables
	double dummy_R0, delta_zero, delta_para, delta_perp ;	
	double res1, res2, res1Frac ;
	getPhysicsParameters( gamma_in, dgam, delta_ms, phi_s, dummy_R0, Rp, Rt, delta_zero, delta_para, delta_perp, tagFraction, res1, res2, res1Frac);
	delta1 = delta_perp -  delta_para ;    
	delta2 = delta_perp -  delta_zero ;
	
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue();
	ctheta_tr = measurement->GetObservable( cosThetaName )->GetValue();
	phi_tr      = measurement->GetObservable( phiName )->GetValue();
	ctheta_1   = measurement->GetObservable( cosPsiName )->GetValue();	
	tag = (int)measurement->GetObservable( tagName )->GetValue();

	double val1, val2 ;
	

	if(res1Frac >= 0.9999 ) {
		// Set the member variable for time resolution to the first value and calculate
		resolution = res1 ;
		val1 = this->diffXsec( );
		return val1 ;
	}
	else {
		// Set the member variable for time resolution to the first value and calculate
		resolution = res1 ;
		val1 = this->diffXsec( );
		// Set the member variable for time resolution to the second value and calculate
		resolution = res2 ;
		val2 = this->diffXsec( );
		return res1Frac*val1 + (1. - res1Frac)*val2 ;
	}
	
}


//...............................................................
//Calculate the normalisation for a given set of physics parameters and boundary

double Bs2JpsiPhi_mistagParameter_alt::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	
    // Get parameters into member variables
	double dummy_R0, delta_zero, delta_para, delta_perp ;	
	double res1, res2, res1Frac ;
	getPhysicsParameters( gamma_in, dgam, delta_ms, phi_s, dummy_R0, Rp, Rt, delta_zero, delta_para, delta_perp, tagFraction, res1, res2, res1Frac);
	delta1 = delta_perp -  delta_para ;    
	delta2 = delta_perp -  delta_zero ;
	
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue();
	ctheta_tr = measurement->GetObservable( cosThetaName )->GetValue();
	phi_tr      = measurement->GetObservable( phiName )->GetValue();
	ctheta_1   = measurement->GetObservable( cosPsiName )->GetValue();	

	// Get time boundaries into member variables
	IConstraint * timeBound = boundary->GetConstraint("time");
	if ( timeBound->GetUnit() == "NameNotFoundError" ) {
		cerr << "Bound on time not provided" << endl;
		return 0;
	}
	else {
		tlo = timeBound->GetMinimum();
		thi = timeBound->GetMaximum();
	}
	
	
	// Recalculate cached values if Physics parameters have changed
	if( ! normalisationCacheValid )  {
		for( tag = -1; tag <= 1; tag ++ ) {
            resolution =  res1 ;
			normalisationCacheValueRes1[tag+1] = this->diffXsecNorm1( );
            resolution =  res2 ;
			normalisationCacheValueRes2[tag+1] = this->diffXsecNorm1( );
		}
		normalisationCacheValid = true ;
	}	
	
	// Return normalisation value according to tag 

	tag = (int)measurement->GetObservable( tagName )->GetValue();

	return res1Frac*normalisationCacheValueRes1[tag+1] + (1. - res1Frac)*normalisationCacheValueRes2[tag+1] ;
}


//..............................................
// Extract and return physics parameters
void Bs2JpsiPhi_mistagParameter_alt::getPhysicsParameters( double & gamma
					, double & deltaGamma
					, double & deltaM
					, double & Phi_s
					, double & Azero_sq
					, double & Apara_sq
					, double & Aperp_sq
					, double & delta_zero 
					, double & delta_para
					, double & delta_perp
					, double & mistag  
					, double & res1
					, double & res2
					, double & res1Frac)
{
	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
    deltaGamma = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	deltaM     = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
	Phi_s      = allParameters.GetPhysicsParameter( Phi_sName )->GetValue();
	Azero_sq   = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
	//Apara_sq   = allParameters.GetPhysicsParameter( Apara_sqName )->GetValue();
	Aperp_sq   = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
	delta_zero = allParameters.GetPhysicsParameter( delta_zeroName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	mistag = allParameters.GetPhysicsParameter( mistagName )->GetValue();
	res1 = allParameters.GetPhysicsParameter( res1Name )->GetValue();
	res2 = allParameters.GetPhysicsParameter( res2Name )->GetValue();
	res1Frac = allParameters.GetPhysicsParameter( res1FractionName )->GetValue();

	Apara_sq = 1 - Azero_sq - Aperp_sq;

	return;
}


//....................................
//Internal helper functions

//Amplitudes Used in one angle PDF
double Bs2JpsiPhi_mistagParameter_alt::AoAo() const  { return Rt; };   
double Bs2JpsiPhi_mistagParameter_alt::AeAe() const { return 1-Rt ; };

//Amplitudes Used in three angle PDF
double Bs2JpsiPhi_mistagParameter_alt::AT() const { return sqrt(Rt) ; };
double Bs2JpsiPhi_mistagParameter_alt::AP() const { return sqrt(Rp) ; };
double Bs2JpsiPhi_mistagParameter_alt::A0() const { if( (1-Rt-Rp) < 0 ) return 0; else return sqrt(1-Rt-Rp) ; };

double Bs2JpsiPhi_mistagParameter_alt::ctrsq() const { return (ctheta_tr*ctheta_tr) ; }
double Bs2JpsiPhi_mistagParameter_alt::strsq() const { return (1.0 - ctheta_tr*ctheta_tr) ; }
double Bs2JpsiPhi_mistagParameter_alt::ct1sq() const { return (ctheta_1*ctheta_1) ; }
double Bs2JpsiPhi_mistagParameter_alt::st1sq() const { return (1.0 - ctheta_1*ctheta_1) ; }
double Bs2JpsiPhi_mistagParameter_alt::cphsq() const { return (cos(phi_tr)*cos(phi_tr)) ; }
double Bs2JpsiPhi_mistagParameter_alt::sphsq() const { return (sin(phi_tr)*sin(phi_tr)) ; }

double Bs2JpsiPhi_mistagParameter_alt::gamma_l() const { return gamma() + ( dgam / 2.0 ) ; }
double Bs2JpsiPhi_mistagParameter_alt::gamma_h() const { return gamma() - ( dgam / 2.0 ) ; }
double Bs2JpsiPhi_mistagParameter_alt::gamma() const   { return gamma_in ; }

double Bs2JpsiPhi_mistagParameter_alt::q() const { return tag ;}


//--------------------------------------------------------------------------
// Time primitives including single gaussian resolution
// These now interface to an external helper library
//
//...................................................
//Exponentials 

double Bs2JpsiPhi_mistagParameter_alt::expL() const 
{
	return Mathematics::Exp( t, gamma_l(), resolution ) ;
}

double Bs2JpsiPhi_mistagParameter_alt::expH() const 
{
	return Mathematics::Exp( t, gamma_h(), resolution ) ;
}

double Bs2JpsiPhi_mistagParameter_alt::intExpL( ) const {
	return Mathematics::ExpInt( tlo, thi, gamma_l(), resolution )  ;
}

double Bs2JpsiPhi_mistagParameter_alt::intExpH( ) const {
	return Mathematics::ExpInt( tlo, thi, gamma_h(), resolution )  ;
}


//......................................................
// Exponential x sine  and cosine

double Bs2JpsiPhi_mistagParameter_alt::expSin() const  
{
    return Mathematics::ExpSin( t, gamma(), delta_ms, resolution ) ;
}

double Bs2JpsiPhi_mistagParameter_alt::expCos() const 
{
    return Mathematics::ExpCos( t, gamma(), delta_ms, resolution ) ;
}

double Bs2JpsiPhi_mistagParameter_alt::intExpSin( ) const 
{
	return Mathematics::ExpSinInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
}

// Integral of exp( - G * t ) * cos( dm * t )  
double Bs2JpsiPhi_mistagParameter_alt::intExpCos( ) const 
{
	return Mathematics::ExpCosInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
}



//------------------------------------------------------------------------------
// These are the time factors and their analytic integrals for the one angle PDF

//..................................
double Bs2JpsiPhi_mistagParameter_alt::timeFactorEven(  )  const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = 
	( 1.0 + cos(phi_s) ) * expL( ) 
	+ ( 1.0 - cos(phi_s) ) * expH( ) 
	+ q() * ( 2.0 * sin(phi_s)   ) * expSin( ) * (1.0 - 2.0*tagFraction) ;
	return result ;
};

double Bs2JpsiPhi_mistagParameter_alt::timeFactorEvenInt(  )  const
{

	double result = 
	( 1.0 + cos(phi_s) )  * intExpL()     
	+ ( 1.0 - cos(phi_s) )  * intExpH()          
	+ q() * ( 2.0 * sin(phi_s)   ) * intExpSin( ) * (1.0 - 2.0*tagFraction) ;
	return result ;
};


//..................................
double Bs2JpsiPhi_mistagParameter_alt::timeFactorOdd(  )   const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = 
	( 1.0 - cos(phi_s) ) * expL( ) 
	+ ( 1.0 + cos(phi_s) ) * expH( ) 
	- q() * ( 2.0 * sin(phi_s)   ) * expSin( ) * (1.0 - 2.0*tagFraction) ;
	return result ;
};

double Bs2JpsiPhi_mistagParameter_alt::timeFactorOddInt(  )  const
{
	double result = 
	( 1.0 - cos(phi_s) ) * intExpL()
	+ ( 1.0 + cos(phi_s) ) * intExpH() 
	- q() * ( 2.0 * sin(phi_s)   ) * intExpSin( ) * (1.0 - 2.0*tagFraction) ;
	return result ;
};


//----------------------------------------------------------
// These are the time factors and their analytic integrals for the three angle PDF

//...........................
double Bs2JpsiPhi_mistagParameter_alt::timeFactorA0A0( )    const { return timeFactorEven( ) ; } ;      
double Bs2JpsiPhi_mistagParameter_alt::timeFactorA0A0Int( ) const { return timeFactorEvenInt( ) ; } ;

//...........................
double Bs2JpsiPhi_mistagParameter_alt::timeFactorAPAP( )    const { return timeFactorEven( ) ; } ;
double Bs2JpsiPhi_mistagParameter_alt::timeFactorAPAPInt( ) const { return timeFactorEvenInt( ) ; } ;

//...........................
double Bs2JpsiPhi_mistagParameter_alt::timeFactorATAT( )    const { return timeFactorOdd( ) ; } ;
double Bs2JpsiPhi_mistagParameter_alt::timeFactorATATInt( ) const { return timeFactorOddInt( ) ; } ;

//...........................
double Bs2JpsiPhi_mistagParameter_alt::timeFactorReA0AP( )  const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = cos(delta2-delta1) * this->timeFactorEven(  ) ;
	return result ;
} ;

double Bs2JpsiPhi_mistagParameter_alt::timeFactorReA0APInt( ) const
{
	double result = cos(delta2-delta1) * this->timeFactorEvenInt( ) ;
	return result ;
} ;

//...........................
double Bs2JpsiPhi_mistagParameter_alt::timeFactorImAPAT( ) const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = 
	q() * 2.0  * ( sin(delta1)*expCos( ) - cos(delta1)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)
	- 1.0 * ( expH( ) - expL( ) ) * cos(delta1) * sin(phi_s)  ;
	
	return result ;
} ;

double Bs2JpsiPhi_mistagParameter_alt::timeFactorImAPATInt( ) const
{
	double _tlo = tlo ;
	if(_tlo < 0.) _tlo = 0. ;
	
	double result = 
	q() * 2.0  * ( sin(delta1)*intExpCos() - cos(delta1)*cos(phi_s)*intExpSin() ) * (1.0 - 2.0*tagFraction)
	- 1.0 * ( intExpH() - intExpL() ) * cos(delta1) * sin(phi_s) ;	    
	return result ;
} ;


//...........................
double Bs2JpsiPhi_mistagParameter_alt::timeFactorImA0AT(  ) const
{
	//if( t < 0.0 ) return 0.0 ;
	double result =
	q() * 2.0  * ( sin(delta2)*expCos( ) - cos(delta2)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)	
	-1.0 * ( expH( ) - expL( ) ) * cos(delta2) * sin(phi_s) ;
	return result ;
} ;

double Bs2JpsiPhi_mistagParameter_alt::timeFactorImA0ATInt( ) const
{
	double _tlo = tlo ;
	if(_tlo < 0.) _tlo = 0. ;
	
	double result = 
	q() * 2.0  * ( sin(delta2)*intExpCos() - cos(delta2)*cos(phi_s)*intExpSin()  ) * (1.0 - 2.0*tagFraction)
	-1.0 * ( intExpH() - intExpL()  ) * cos(delta2) * sin(phi_s) ;
	return result ;
} ;


//------------------------------------------------------
// Angle factors for one angle PDF

//.................................
double Bs2JpsiPhi_mistagParameter_alt::angleFactorEven(  )  const
{
	// Note that this is normalised to 1
	double result = 3.0/8.0 * (1.0 + ctrsq() ) ;
	return result ;
};

//.................................
double Bs2JpsiPhi_mistagParameter_alt::angleFactorOdd(  )   const
{
	// Note that this is normalised to 1
	double result = 3.0/4.0 * (1.0 - ctrsq() ) ;
	return result ;
};


//------------------------------------------------------
// Angle factors for three angle PDFs


//...........................
double Bs2JpsiPhi_mistagParameter_alt::angleFactorA0A0(  ) const
{
	// Normalised to  1	
	double result = 2.0 * ct1sq() * (1.0 - strsq()*cphsq() ) * (9.0/32.0/TMath::Pi());
	return result ;	
};

//...........................
double Bs2JpsiPhi_mistagParameter_alt::angleFactorAPAP(  ) const
{
	// Normalised to  1
	double result =  st1sq() * (1.0 - strsq()*sphsq() ) * (9.0/32.0/TMath::Pi());
	return result ;	
};

//...........................
double Bs2JpsiPhi_mistagParameter_alt::angleFactorATAT(  ) const
{
	// Normalised to  1
	double result = st1sq() * strsq() * (9.0/32.0/TMath::Pi());
	return result ;
	
};

//...........................
double Bs2JpsiPhi_mistagParameter_alt::angleFactorReA0AP( ) const
{
	// Normalised to  0
	double theta_1 = acos(ctheta_1) ;	
	double result =    sin(2.0*theta_1) * strsq() * sin(2.0*phi_tr) / sqrt(2.0) * (9.0/32.0/TMath::Pi());
	return result ;	
};

//...........................
double Bs2JpsiPhi_mistagParameter_alt::angleFactorImAPAT(  ) const
{
	// Normalised to  0
	double theta_tr = acos(ctheta_tr) ;		
	double result =   -1.0 *  st1sq() * sin(2.0*theta_tr) * sin(phi_tr) * (9.0/32.0/TMath::Pi()) ;
	return result ;	
};

//...........................
double Bs2JpsiPhi_mistagParameter_alt::angleFactorImA0AT(  ) const
{
	// Normalised to  0
	double theta_tr = acos(ctheta_tr) ;		
	double theta_1 = acos(ctheta_1) ;		
	double result =  +1.0*   sin(2.0*theta_1) * sin(2.0*theta_tr) * cos(phi_tr) / sqrt(2.0) * (9.0/32.0/TMath::Pi());
	return result ;	
};


//-------------------------------------------------------------
// Putting it all together to make up the differential cross sections.


//...................................
// Diff cross sections

double Bs2JpsiPhi_mistagParameter_alt::diffXsec(  )  const
{   
	double xsec = 
	0.5 * A0()*A0() * timeFactorA0A0(  ) * angleFactorA0A0( ) +
	0.5 * AP()*AP() * timeFactorAPAP(  ) * angleFactorAPAP( ) +
	0.5 * AT()*AT() * timeFactorATAT(  ) * angleFactorATAT( ) +
	0.5 * A0()*AP() * timeFactorReA0AP(  ) * angleFactorReA0AP( ) +
	0.5 * AP()*AT() * timeFactorImAPAT(  ) * angleFactorImAPAT( ) +
	0.5 * A0()*AT() * timeFactorImA0AT(  ) * angleFactorImA0AT( ) ;
	
	return xsec ;
};

double  Bs2JpsiPhi_mistagParameter_alt::diffXsecOne(  ) const
{
	double result = 
	0.5 * AeAe() * timeFactorEven(  ) * angleFactorEven(  )  +
	0.5 * AoAo() * timeFactorOdd(  )  * angleFactorOdd(  ) ;
	return result ;
};

//...................................
// Integral over all variables: t + angles

double Bs2JpsiPhi_mistagParameter_alt::diffXsecNorm1(  ) const
{      
	double norm = 
	0.5 * A0()*A0() * timeFactorA0A0Int(  ) +    // Angle factors normalised to 1
	0.5 * AP()*AP() * timeFactorAPAPInt(  ) +
	0.5 * AT()*AT() * timeFactorATATInt(  ) ;
	
	return norm ;
};

double Bs2JpsiPhi_mistagParameter_alt::diffXsecOneNorm1(  ) const
{      
	double norm = 
	0.5 * AeAe() * timeFactorEvenInt(  )  +    // Angle factors normalised to 1
	0.5 * AoAo() * timeFactorOddInt(  )   ;
	return norm ;
};


//...................................
// Integral over angles only 3 

double Bs2JpsiPhi_mistagParameter_alt::diffXsecNorm2(  ) const
{          
	double norm = 
	0.5 * A0()*A0() * timeFactorA0A0(  ) +    // Angle factors normalised to 1
	0.5 * AP()*AP() * timeFactorAPAP(  ) +
	0.5 * AT()*AT() * timeFactorATAT(  ) ;
	
	return norm ;
};

double Bs2JpsiPhi_mistagParameter_alt::diffXsecOneNorm2(  ) const
{          
	double norm = 
	0.5 * AeAe() * timeFactorEven(  )  +     // Angle factors normalised to 1
	0.5 * AoAo() * timeFactorOdd(  )   ;
	return norm ;
};


