// $Id: Bs2JpsiPhi_mistagObservable_alt.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2JpsiPhi_mistagObservable_alt Bs2JpsiPhi_mistagObservable_alt.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2009-12-06
 */

#include "Bs2JpsiPhi_mistagObservable_alt.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"
#include "Mathematics.h"

//......................................
//Constructor

Bs2JpsiPhi_mistagObservable_alt::Bs2JpsiPhi_mistagObservable_alt() : 
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
	// Detector parameters
	, res1Name	( "timeResolution1" )
	, res2Name	( "timeResolution2" )
	, res1FractionName	( "timeResolution1Fraction" )
	, timeOffsetName	( "timeOffset" )
	// Angular acceptance factors
	, angAccI1Name ( "angAccI1" )
	, angAccI2Name ( "angAccI2" )
	, angAccI3Name ( "angAccI3" )
	, angAccI4Name ( "angAccI4" )
	, angAccI5Name ( "angAccI5" )
	, angAccI6Name ( "angAccI6" )
	// Observables
	, timeName	    ( "time" )
	, cosThetaName	( "cosTheta" )
	, phiName	    ( "phi" )
	, cosPsiName	( "cosPsi" )
	, tagName	    ( "tag" )
	, mistagName	( "mistag" )
	// Other things
	, normalisationCacheValid(false)
{
	MakePrototypes();
	
	std::cout << "Constructing signal PDF: Bs2JpsiPhi_mistagObservable_alt " << std::endl ;
}

//Make the data point and parameter set
void Bs2JpsiPhi_mistagObservable_alt::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( cosThetaName );
	allObservables.push_back( phiName );
	allObservables.push_back( cosPsiName );
	allObservables.push_back( tagName );
	allObservables.push_back( mistagName );

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
	parameterNames.push_back( res1Name );
	parameterNames.push_back( res2Name );
	parameterNames.push_back( res1FractionName );
	parameterNames.push_back( timeOffsetName );
	parameterNames.push_back( angAccI1Name );
	parameterNames.push_back( angAccI2Name );
	parameterNames.push_back( angAccI3Name );
	parameterNames.push_back( angAccI4Name );
	parameterNames.push_back( angAccI5Name );
	parameterNames.push_back( angAccI6Name );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}


//........................................................
//Destructor
Bs2JpsiPhi_mistagObservable_alt::~Bs2JpsiPhi_mistagObservable_alt()
{
}

//........................................................
//Set the physics parameters into member variables
//Indicate that the cache is no longer valid

bool Bs2JpsiPhi_mistagObservable_alt::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	
	bool result = allParameters.SetPhysicsParameters(NewParameterSet);
	
	/// Some gymnastics here to match xml parameters to my original pdf parameters 
	
	// Physics parameters. 
	gamma_in  = allParameters.GetPhysicsParameter( gammaName )->GetValue();
    dgam      = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	delta_ms  = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
	phi_s     = allParameters.GetPhysicsParameter( Phi_sName )->GetValue();

	Azero_sq = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
	if( (Azero_sq < 0.) || (Azero_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_mistagObservable_alt::SetPhysicsParameters: Azero_sq <0 or >1 but left as is" <<  endl ;	}	
	Aperp_sq = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
	if( (Aperp_sq < 0.) || (Aperp_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_mistagObservable_alt::SetPhysicsParameters: Aperp_sq <0 or >1 but left as is" <<  endl ;	}	

	Apara_sq = (1. - Azero_sq - Aperp_sq ) ;
	if( Apara_sq < 0. ) {
		cout << "Warning in Bs2JpsiPhi_mistagObservable_alt::SetPhysicsParameters: derived parameter Apara_sq <0  and so set to zero" <<  endl ;
		Apara_sq = 0. ;
	}	
	
	Rt = Aperp_sq ;  // for backwards compatibility with 1 angle fits
	
	double _delta_zero = allParameters.GetPhysicsParameter( delta_zeroName )->GetValue();
	double _delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	double _delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	delta1 = _delta_perp -  _delta_para ;    
	delta2 = _delta_perp -  _delta_zero ;
	
	// Detector parameters
	resolution1         = allParameters.GetPhysicsParameter( res1Name )->GetValue();
	resolution2         = allParameters.GetPhysicsParameter( res2Name )->GetValue();
	resolution1Fraction = allParameters.GetPhysicsParameter( res1FractionName )->GetValue();
	timeOffset          = allParameters.GetPhysicsParameter( timeOffsetName )->GetValue();
	
	// Angular acceptance factors
	angAccI1 = allParameters.GetPhysicsParameter( angAccI1Name )->GetValue();
	angAccI2 = allParameters.GetPhysicsParameter( angAccI2Name )->GetValue();
	angAccI3 = allParameters.GetPhysicsParameter( angAccI3Name )->GetValue();
	angAccI4 = allParameters.GetPhysicsParameter( angAccI4Name )->GetValue();
	angAccI5 = allParameters.GetPhysicsParameter( angAccI5Name )->GetValue();
	angAccI6 = allParameters.GetPhysicsParameter( angAccI6Name )->GetValue();

	return result;
}

//.........................................................
//Return a list of observables not to be integrated
vector<string> Bs2JpsiPhi_mistagObservable_alt::GetDoNotIntegrateList()
{
	vector<string> list;
//PELC
	list.push_back(mistagName) ;
	return list;
}

//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2JpsiPhi_mistagObservable_alt::Evaluate(DataPoint * measurement)
{
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset ;
	ctheta_tr = measurement->GetObservable( cosThetaName )->GetValue();
	phi_tr      = measurement->GetObservable( phiName )->GetValue();
	ctheta_1   = measurement->GetObservable( cosPsiName )->GetValue();	
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	tagFraction = measurement->GetObservable( mistagName )->GetValue();
	//tagFraction= 0.5;  //PELC

	
	double val1, val2 ;
	double returnValue ;
	
	if(resolution1Fraction >= 0.9999 ) {

		// Set the member variable for time resolution to the first value and calculate
		resolution = resolution1 ;

		// Pre calculate the time primitives
		expL_stored = Mathematics::Exp( t, gamma_l(), resolution ) ;
		expH_stored = Mathematics::Exp( t, gamma_h(), resolution ) ;
		expSin_stored = Mathematics::ExpSin( t, gamma(), delta_ms, resolution ) ;
		expCos_stored = Mathematics::ExpCos( t, gamma(), delta_ms, resolution ) ;
				
		val1 = this->diffXsec( );
		returnValue = val1 ;

	}
	else {

		// Set the member variable for time resolution to the first value and calculate
		resolution = resolution1 ;

		// Pre calculate the time primitives
		expL_stored = Mathematics::Exp( t, gamma_l(), resolution ) ;
		expH_stored = Mathematics::Exp( t, gamma_h(), resolution ) ;
		expSin_stored = Mathematics::ExpSin( t, gamma(), delta_ms, resolution ) ;
		expCos_stored = Mathematics::ExpCos( t, gamma(), delta_ms, resolution ) ;
				
		val1 = this->diffXsec( );

		// Set the member variable for time resolution to the second value and calculate
		resolution = resolution2 ;

		// Pre calculate the time primitives
		expL_stored = Mathematics::Exp( t, gamma_l(), resolution ) ;
		expH_stored = Mathematics::Exp( t, gamma_h(), resolution ) ;
		expSin_stored = Mathematics::ExpSin( t, gamma(), delta_ms, resolution ) ;
		expCos_stored = Mathematics::ExpCos( t, gamma(), delta_ms, resolution ) ;		
		
		val2 = this->diffXsec( );
		
		returnValue = resolution1Fraction*val1 + (1. - resolution1Fraction)*val2 ;
				
	}
	

	if(  ( (returnValue <= 0.) && (t>0.) ) || isnan(returnValue) ) {
		cout << endl ;
		cout << " Bs2JpsiPhi_mistagObservable_alt::evaluate() returns <=0 or nan :" << returnValue << endl ;
		cout << "   gamma " << gamma() ;
		cout << "   gl    " << gamma_l() ;
		cout << "   gh    " << gamma_h() ;
		cout << "   AT    " << AT() ;
		cout << "   AP    " << AP() ;
		cout << "   A0    " << A0() << endl ;
		cout << " For event with: " << endl ;
		cout << "   time " << t << endl ;
		if( isnan(returnValue) ) exit(1) ;
	}
	
	return returnValue ;
	
}


//...............................................................
//Calculate the normalisation for a given set of physics parameters and boundary

double Bs2JpsiPhi_mistagObservable_alt::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	
	
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset;
	ctheta_tr = measurement->GetObservable( cosThetaName )->GetValue();
	phi_tr      = measurement->GetObservable( phiName )->GetValue();
	ctheta_1   = measurement->GetObservable( cosPsiName )->GetValue();	
	tagFraction = measurement->GetObservable( mistagName )->GetValue();
	//tagFraction= 0.5;  //PELC

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
	// Must do this for each of the two resolutions.
//PELC
// I dont think you can cache any more as normalisation depends upon the mistag which now changes per event.
	if( true /*! normalisationCacheValid*/ )  {
		for( tag = -1; tag <= 1; tag ++ ) {
            resolution =  resolution1 ;
			normalisationCacheValueRes1[tag+1] = this->diffXsecNorm1( );
            resolution =  resolution2 ;
			normalisationCacheValueRes2[tag+1] = this->diffXsecNorm1( );
		}
		normalisationCacheValid = true ;
	}	
	
	// Return normalisation value according to tag 

	tag = (int)measurement->GetObservable( tagName )->GetValue();

	double returnValue  = resolution1Fraction*normalisationCacheValueRes1[tag+1] + (1. - resolution1Fraction)*normalisationCacheValueRes2[tag+1] ;
	
	if( (returnValue <= 0.) || isnan(returnValue) ) {
		cout << " Bs2JpsiPhi_mistagObservable_alt::Normalisation() returns <=0 or nan " << endl ;
		cout << " gamma " << gamma() ;
		cout << " gl    " << gamma_l() ;
		cout << " gh    " << gamma_h() ;
		cout << " AT    " << AT() ;
		cout << " AP    " << AP() ;
		cout << " A0    " << A0() ;
		exit(1) ;
	}
	
	return returnValue ;
}


//....................................
//Internal helper functions

//Amplitudes Used in one angle PDF
double Bs2JpsiPhi_mistagObservable_alt::AoAo() const  { return Rt; };   
double Bs2JpsiPhi_mistagObservable_alt::AeAe() const { return 1-Rt ; };

//Amplitudes Used in three angle PDF
double Bs2JpsiPhi_mistagObservable_alt::AT() const { 
	if( Aperp_sq <= 0. ) return 0. ;
	else return sqrt(Aperp_sq) ; 
};
double Bs2JpsiPhi_mistagObservable_alt::AP() const { 
	if( Apara_sq <= 0. ) return 0. ;
	else return sqrt(Apara_sq) ; 
};
double Bs2JpsiPhi_mistagObservable_alt::A0() const { 
	if( Azero_sq <= 0. ) return 0. ;
	else return sqrt(Azero_sq) ; 
};

double Bs2JpsiPhi_mistagObservable_alt::ctrsq() const { return (ctheta_tr*ctheta_tr) ; }
double Bs2JpsiPhi_mistagObservable_alt::strsq() const { return (1.0 - ctheta_tr*ctheta_tr) ; }
double Bs2JpsiPhi_mistagObservable_alt::ct1sq() const { return (ctheta_1*ctheta_1) ; }
double Bs2JpsiPhi_mistagObservable_alt::st1sq() const { return (1.0 - ctheta_1*ctheta_1) ; }
double Bs2JpsiPhi_mistagObservable_alt::cphsq() const { return (cos(phi_tr)*cos(phi_tr)) ; }
double Bs2JpsiPhi_mistagObservable_alt::sphsq() const { return (sin(phi_tr)*sin(phi_tr)) ; }

double Bs2JpsiPhi_mistagObservable_alt::gamma_l() const { 
	double gl = gamma() + ( dgam / 2.0 ) ;
	if( gl <= 0. ) {
		cout << " In Bs2JpsiPhi_mistagObservable_alt : gamma_l() < 0 so setting it to 0.0000001 " << endl ;
		return 0.0000001 ;
	}
	else
		return gl ; 
}

double Bs2JpsiPhi_mistagObservable_alt::gamma_h() const { 
	double gh = gamma() - ( dgam / 2.0 ) ;
	if( gh <= 0. ) {
		cout << " In Bs2JpsiPhi_mistagObservable_alt : gamma_h() < 0 so setting it to 0.0000001 " << endl ;
		return 0.0000001 ;
	}
	else
		return gh ; 
}

double Bs2JpsiPhi_mistagObservable_alt::gamma() const   { return gamma_in ; }

double Bs2JpsiPhi_mistagObservable_alt::q() const { return tag ;}


//--------------------------------------------------------------------------
// Time primitives including single gaussian resolution
// These now interface to an external helper library
//
//...................................................
//Exponentials 

double Bs2JpsiPhi_mistagObservable_alt::expL() const 
{
	return expL_stored ; //Mathematics::Exp( t, gamma_l(), resolution ) ;
}

double Bs2JpsiPhi_mistagObservable_alt::expH() const 
{
	return expH_stored ; //Mathematics::Exp( t, gamma_h(), resolution ) ;
}

double Bs2JpsiPhi_mistagObservable_alt::intExpL( ) const {
	return Mathematics::ExpInt( tlo, thi, gamma_l(), resolution )  ;
}

double Bs2JpsiPhi_mistagObservable_alt::intExpH( ) const {
	return Mathematics::ExpInt( tlo, thi, gamma_h(), resolution )  ;
}


//......................................................
// Exponential x sine  and cosine

double Bs2JpsiPhi_mistagObservable_alt::expSin() const  
{
    return expSin_stored ; // Mathematics::ExpSin( t, gamma(), delta_ms, resolution ) ;
}

double Bs2JpsiPhi_mistagObservable_alt::expCos() const 
{
    return expCos_stored ; // Mathematics::ExpCos( t, gamma(), delta_ms, resolution ) ;
}

double Bs2JpsiPhi_mistagObservable_alt::intExpSin( ) const 
{
	return Mathematics::ExpSinInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
}

// Integral of exp( - G * t ) * cos( dm * t )  
double Bs2JpsiPhi_mistagObservable_alt::intExpCos( ) const 
{
	return Mathematics::ExpCosInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
}



//------------------------------------------------------------------------------
// These are the time factors and their analytic integrals for the one angle PDF

//..................................
double Bs2JpsiPhi_mistagObservable_alt::timeFactorEven(  )  const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = 
	( 1.0 + cos(phi_s) ) * expL( ) 
	+ ( 1.0 - cos(phi_s) ) * expH( ) 
	+ q() * ( 2.0 * sin(phi_s)   ) * expSin( ) * (1.0 - 2.0*tagFraction) ;
	return result ;
};

double Bs2JpsiPhi_mistagObservable_alt::timeFactorEvenInt(  )  const
{

	double result = 
	( 1.0 + cos(phi_s) )  * intExpL()     
	+ ( 1.0 - cos(phi_s) )  * intExpH()          
	+ q() * ( 2.0 * sin(phi_s)   ) * intExpSin( ) * (1.0 - 2.0*tagFraction) ;
	return result ;
};


//..................................
double Bs2JpsiPhi_mistagObservable_alt::timeFactorOdd(  )   const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = 
	( 1.0 - cos(phi_s) ) * expL( ) 
	+ ( 1.0 + cos(phi_s) ) * expH( ) 
	- q() * ( 2.0 * sin(phi_s)   ) * expSin( ) * (1.0 - 2.0*tagFraction) ;
	return result ;
};

double Bs2JpsiPhi_mistagObservable_alt::timeFactorOddInt(  )  const
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
double Bs2JpsiPhi_mistagObservable_alt::timeFactorA0A0( )    const { return timeFactorEven( ) ; } ;      
double Bs2JpsiPhi_mistagObservable_alt::timeFactorA0A0Int( ) const { return timeFactorEvenInt( ) ; } ;

//...........................
double Bs2JpsiPhi_mistagObservable_alt::timeFactorAPAP( )    const { return timeFactorEven( ) ; } ;
double Bs2JpsiPhi_mistagObservable_alt::timeFactorAPAPInt( ) const { return timeFactorEvenInt( ) ; } ;

//...........................
double Bs2JpsiPhi_mistagObservable_alt::timeFactorATAT( )    const { return timeFactorOdd( ) ; } ;
double Bs2JpsiPhi_mistagObservable_alt::timeFactorATATInt( ) const { return timeFactorOddInt( ) ; } ;

//...........................
double Bs2JpsiPhi_mistagObservable_alt::timeFactorReA0AP( )  const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = cos(delta2-delta1) * this->timeFactorEven(  ) ;
	return result ;
} ;

double Bs2JpsiPhi_mistagObservable_alt::timeFactorReA0APInt( ) const
{
	double result = cos(delta2-delta1) * this->timeFactorEvenInt( ) ;
	return result ;
} ;

//...........................
double Bs2JpsiPhi_mistagObservable_alt::timeFactorImAPAT( ) const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = 
	q() * 2.0  * ( sin(delta1)*expCos( ) - cos(delta1)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)
	- 1.0 * ( expH( ) - expL( ) ) * cos(delta1) * sin(phi_s)  ;
	
	return result ;
} ;

double Bs2JpsiPhi_mistagObservable_alt::timeFactorImAPATInt( ) const
{
	double _tlo = tlo ;
	if(_tlo < 0.) _tlo = 0. ;
	
	double result = 
	q() * 2.0  * ( sin(delta1)*intExpCos() - cos(delta1)*cos(phi_s)*intExpSin() ) * (1.0 - 2.0*tagFraction)
	- 1.0 * ( intExpH() - intExpL() ) * cos(delta1) * sin(phi_s) ;	    
	return result ;
} ;


//...........................
double Bs2JpsiPhi_mistagObservable_alt::timeFactorImA0AT(  ) const
{
	//if( t < 0.0 ) return 0.0 ;
	double result =
	q() * 2.0  * ( sin(delta2)*expCos( ) - cos(delta2)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)	
	-1.0 * ( expH( ) - expL( ) ) * cos(delta2) * sin(phi_s) ;
	return result ;
} ;

double Bs2JpsiPhi_mistagObservable_alt::timeFactorImA0ATInt( ) const
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
double Bs2JpsiPhi_mistagObservable_alt::angleFactorEven(  )  const
{
	// Note that this is normalised to 1
	double result = 3.0/8.0 * (1.0 + ctrsq() ) ;
	return result ;
};

//.................................
double Bs2JpsiPhi_mistagObservable_alt::angleFactorOdd(  )   const
{
	// Note that this is normalised to 1
	double result = 3.0/4.0 * (1.0 - ctrsq() ) ;
	return result ;
};


//------------------------------------------------------
// Angle factors for three angle PDFs


//...........................
double Bs2JpsiPhi_mistagObservable_alt::angleFactorA0A0(  ) const
{
	// Normalised to  1	
	double result = 2.0 * ct1sq() * (1.0 - strsq()*cphsq() ) * (9.0/32.0/TMath::Pi());
	return result ;	
};

//...........................
double Bs2JpsiPhi_mistagObservable_alt::angleFactorAPAP(  ) const
{
	// Normalised to  1
	double result =  st1sq() * (1.0 - strsq()*sphsq() ) * (9.0/32.0/TMath::Pi());
	return result ;	
};

//...........................
double Bs2JpsiPhi_mistagObservable_alt::angleFactorATAT(  ) const
{
	// Normalised to  1
	double result = st1sq() * strsq() * (9.0/32.0/TMath::Pi());
	return result ;
	
};

//...........................
double Bs2JpsiPhi_mistagObservable_alt::angleFactorReA0AP( ) const
{
	// Normalised to  0
	double theta_1 = acos(ctheta_1) ;	
	double result =    sin(2.0*theta_1) * strsq() * sin(2.0*phi_tr) / sqrt(2.0) * (9.0/32.0/TMath::Pi());
	return result ;	
};

//...........................
double Bs2JpsiPhi_mistagObservable_alt::angleFactorImAPAT(  ) const
{
	// Normalised to  0
	double theta_tr = acos(ctheta_tr) ;		
	double result =   -1.0 *  st1sq() * sin(2.0*theta_tr) * sin(phi_tr) * (9.0/32.0/TMath::Pi()) ;
	return result ;	
};

//...........................
double Bs2JpsiPhi_mistagObservable_alt::angleFactorImA0AT(  ) const
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

double Bs2JpsiPhi_mistagObservable_alt::diffXsec(  )  const
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

double  Bs2JpsiPhi_mistagObservable_alt::diffXsecOne(  ) const
{
	double result = 
	0.5 * AeAe() * timeFactorEven(  ) * angleFactorEven(  )  +
	0.5 * AoAo() * timeFactorOdd(  )  * angleFactorOdd(  ) ;
	return result ;
};

//...................................
// Integral over all variables: t + angles

double Bs2JpsiPhi_mistagObservable_alt::diffXsecNorm1(  ) const
{      
	double reference =  32.0*TMath::Pi()/9.0 ;
	
	double norm = 
	0.5 * A0()*A0() * timeFactorA0A0Int(  ) * angAccI1   +  
	0.5 * AP()*AP() * timeFactorAPAPInt(  ) * angAccI2   +  
	0.5 * AT()*AT() * timeFactorATATInt(  ) * angAccI3   +  

	0.5 * A0()*AP() * timeFactorReA0APInt(  ) * angAccI5 +  
	0.5 * AP()*AT() * timeFactorImAPATInt(  ) * angAccI4 +  
	0.5 * A0()*AT() * timeFactorImA0ATInt(  ) * angAccI6 ;  	
	// In the canonical PDF, the ApAt term is number 4!
	return norm ;
};

//  13.432', '14.151', '14.210
double Bs2JpsiPhi_mistagObservable_alt::diffXsecOneNorm1(  ) const
{      
	double norm = 
	0.5 * AeAe() * timeFactorEvenInt(  )  +    // Angle factors normalised to 1
	0.5 * AoAo() * timeFactorOddInt(  )   ;
	return norm ;
};


//...................................
// Integral over angles only 3 

double Bs2JpsiPhi_mistagObservable_alt::diffXsecNorm2(  ) const
{          
	double norm = 
	0.5 * A0()*A0() * timeFactorA0A0(  ) +    // Angle factors normalised to 1
	0.5 * AP()*AP() * timeFactorAPAP(  ) +
	0.5 * AT()*AT() * timeFactorATAT(  ) ;
	
	return norm ;
};

double Bs2JpsiPhi_mistagObservable_alt::diffXsecOneNorm2(  ) const
{          
	double norm = 
	0.5 * AeAe() * timeFactorEven(  )  +     // Angle factors normalised to 1
	0.5 * AoAo() * timeFactorOdd(  )   ;
	return norm ;
};


