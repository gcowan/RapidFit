// $Id: Bs2JpsiPhi_SignalAlt_MP_v1.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2JpsiPhi_SignalAlt_MP_v1 Bs2JpsiPhi_SignalAlt_MP_v1.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2011-01-28
 */

#include "Bs2JpsiPhi_SignalAlt_MP_v1.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"
#include "Mathematics.h"

#include <float.h>

#define DEBUGFLAG true

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

//......................................
//Constructor

Bs2JpsiPhi_SignalAlt_MP_v1::Bs2JpsiPhi_SignalAlt_MP_v1() : 
	Bs2JpsiPhi_SignalAlt_BaseClass()
	, normalisationCacheValid(false)
{
	MakePrototypes();
	std::cout << "Constructing PDF: Bs2JpsiPhi_SignalAlt_MP_v1 " << std::endl ;
	_mistagP1=1.0 ;     // These three are for backwards compatibility
	_mistagP0=0.0 ;      // These three are for backwards compatibility
	_mistagSetPoint=0.0 ; // These three are for backwards compatibility
}


//.......................................
//Make the data point and parameter set
void Bs2JpsiPhi_SignalAlt_MP_v1::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( cosThetaName );
	allObservables.push_back( phiName );
	allObservables.push_back( cosPsiName );
	allObservables.push_back( tagName );
	allObservables.push_back( timeAcceptanceCategoryName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName );
	parameterNames.push_back( deltaGammaName );
	parameterNames.push_back( Aperp_sqName );
	parameterNames.push_back( Azero_sqName );
	parameterNames.push_back( As_sqName );
	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
	parameterNames.push_back( delta_zeroName );
	parameterNames.push_back( delta_sName );
	parameterNames.push_back( deltaMName );
	parameterNames.push_back( Phi_sName );
	parameterNames.push_back( mistagName );
	parameterNames.push_back( res1FractionName );
	parameterNames.push_back( res1Name );
	parameterNames.push_back( res2Name );
	parameterNames.push_back( timeOffsetName );
	parameterNames.push_back( angAccI1Name );
	parameterNames.push_back( angAccI2Name );
	parameterNames.push_back( angAccI3Name );
	parameterNames.push_back( angAccI4Name );
	parameterNames.push_back( angAccI5Name );
	parameterNames.push_back( angAccI6Name );
	parameterNames.push_back( angAccI7Name );
	parameterNames.push_back( angAccI8Name );
	parameterNames.push_back( angAccI9Name );
	parameterNames.push_back( angAccI10Name );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}


//........................................................
//Destructor
Bs2JpsiPhi_SignalAlt_MP_v1::~Bs2JpsiPhi_SignalAlt_MP_v1()
{
}

//........................................................
//Set the physics parameters into member variables
//Indicate that the cache is no longer valid

bool Bs2JpsiPhi_SignalAlt_MP_v1::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	
	bool result = allParameters.SetPhysicsParameters(NewParameterSet);

	// Physics parameters. 
	_gamma  = allParameters.GetPhysicsParameter( gammaName )->GetValue();
	dgam  = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	
	Azero_sq = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
	if( (Azero_sq < 0.) || (Azero_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MP_v1::SetPhysicsParameters: Azero_sq <0 or >1 but left as is" <<  endl ;	}	
	Aperp_sq = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
	if( (Aperp_sq < 0.) || (Aperp_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MP_v1::SetPhysicsParameters: Aperp_sq <0 or >1 but left as is" <<  endl ;	}	
	As_sq = allParameters.GetPhysicsParameter( As_sqName )->GetValue();
	if( (As_sq < 0.) || (As_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MP_v1::SetPhysicsParameters: As_sq <0 or >1 but left as is" <<  endl ;	}	

	Apara_sq = (1. - Azero_sq - Aperp_sq  - As_sq) ;
	if( Apara_sq < 0. ) {
		cout << "Warning in Bs2JpsiPhi_SignalAlt_MP_v1::SetPhysicsParameters: derived parameter Apara_sq <0  and so set to zero" <<  endl ;
		Apara_sq = 0. ;
	}	
		
	delta_zero = allParameters.GetPhysicsParameter( delta_zeroName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	delta_s	   = allParameters.GetPhysicsParameter( delta_sName )->GetValue();
	delta1 = delta_perp -  delta_para ;    
	delta2 = delta_perp -  delta_zero ;

	_mistag         = allParameters.GetPhysicsParameter( mistagName )->GetValue();

	delta_ms  = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
	phi_s     = allParameters.GetPhysicsParameter( Phi_sName )->GetValue();
	_cosphis = cos(phi_s) ;
	_sinphis = sin(phi_s) ;
	
	// Resolution parameters
	resolution1Fraction = allParameters.GetPhysicsParameter( res1FractionName )->GetValue();
	resolution1         = allParameters.GetPhysicsParameter( res1Name )->GetValue();
	resolution2         = allParameters.GetPhysicsParameter( res2Name )->GetValue();
	timeOffset          = allParameters.GetPhysicsParameter( timeOffsetName )->GetValue();
	
	// Angular acceptance factors
	angAccI1 = allParameters.GetPhysicsParameter( angAccI1Name )->GetValue();
	angAccI2 = allParameters.GetPhysicsParameter( angAccI2Name )->GetValue();
	angAccI3 = allParameters.GetPhysicsParameter( angAccI3Name )->GetValue();
	angAccI4 = allParameters.GetPhysicsParameter( angAccI4Name )->GetValue();
	angAccI5 = allParameters.GetPhysicsParameter( angAccI5Name )->GetValue();
	angAccI6 = allParameters.GetPhysicsParameter( angAccI6Name )->GetValue();
	angAccI7 = allParameters.GetPhysicsParameter( angAccI7Name )->GetValue();
	angAccI8 = allParameters.GetPhysicsParameter( angAccI8Name )->GetValue();
	angAccI9 = allParameters.GetPhysicsParameter( angAccI9Name )->GetValue();
	angAccI10 = allParameters.GetPhysicsParameter( angAccI10Name )->GetValue();
	
	// Do a test to ensure user is not using upper time acceptance wrongly
	if( ((fabs(resolution1-0.0)>DOUBLE_TOLERANCE) || (fabs(resolution2-0.0)>DOUBLE_TOLERANCE) || (fabs(mistag()-0.5)>DOUBLE_TOLERANCE) || (fabs(phi_s-0.0)>DOUBLE_TOLERANCE)) && useUpperTimeAcceptance() )
	{
		cout << " You appear to be trying to use the upper time acceptance but are using either resolution or are doing a tagged fit" << endl ;
		cout << " This is not possible at present" << endl ;
		cout << " Resolution1 : " << resolution1 << endl ;
		cout << " Resolution2 : " << resolution2 << endl ;
		cout << " Mistag : " << mistag() << endl ;
		cout << " Phi_s : " << phi_s <<  endl ;
		throw(10);
	}
	
	return result;
}

//.........................................................
//Return a list of observables not to be integrated
vector<string> Bs2JpsiPhi_SignalAlt_MP_v1::GetDoNotIntegrateList()
{
	vector<string> list;
	return list;
}

//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2JpsiPhi_SignalAlt_MP_v1::Evaluate(DataPoint * measurement) 
{
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset;
	ctheta_tr = measurement->GetObservable( cosThetaName )->GetValue();
	phi_tr      = measurement->GetObservable( phiName )->GetValue();
	ctheta_1   = measurement->GetObservable( cosPsiName )->GetValue();	
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	timeAcceptanceCategory = (int)measurement->GetObservable( timeAcceptanceCategoryName )->GetValue();
	
	double val1, val2 ;
	double returnValue ;
	
	if(resolution1Fraction >= 0.9999 ) {
		// Set the member variable for time resolution to the first value and calculate
		resolution = resolution1 ;
		returnValue = this->diffXsec( );
	}
	else {
		// Set the member variable for time resolution to the first value and calculate
		resolution = resolution1 ;
		val1 = this->diffXsec( );
		// Set the member variable for time resolution to the second value and calculate
		resolution = resolution2 ;
		val2 = this->diffXsec( );
		
		returnValue = resolution1Fraction*val1 + (1. - resolution1Fraction)*val2 ;				
	}
	
	//conditions to throw exception
	bool c1 = isnan(returnValue) ;
	bool c2 = ((resolution1>0.)||(resolution2>0.)) && (returnValue <= 0.) ;
	bool c3 = ((fabs(resolution1-0.)<DOUBLE_TOLERANCE)&&((fabs(resolution2-0.)<DOUBLE_TOLERANCE))) && (returnValue <= 0.) && (t>0.) ;

	if( DEBUGFLAG && (c1 || c2 || c3)  ) {
		cout << endl ;
		cout << " Bs2JpsiPhi_SignalAlt_MP_v1::evaluate() returns <=0 or nan :" << returnValue << endl ;
		cout << "   gamma " << gamma() << endl ;
		cout << "   gl    " << gamma_l() << endl ;
		cout << "   gh    " << gamma_h()  << endl;
		cout << "   AT^2    " << AT()*AT() << endl;
		cout << "   AP^2    " << AP()*AP() << endl;
		cout << "   A0^2    " << A0()*A0() << endl ;
		cout << "   AS^2    " << AS()*AS() << endl ;
		cout << "   ATOTAL  " << AS()*AS()+A0()*A0()+AP()*AP()+AT()*AT() << endl ;
		cout << "   delta_ms       " << delta_ms << endl ;
		cout << "   mistag    " << _mistag << endl ;
		cout << " For event with:  " << endl ;
		cout << "   time      " << t << endl ;
		cout << "   ctheta_tr " << ctheta_tr << endl ;
		cout << "   ctheta_1 " << ctheta_1 << endl ;
		cout << "   phi_tr " << phi_tr << endl ;
		if( isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
	
	if( useLowerTimeAcceptance() ) return returnValue * timeAcceptance.acceptance(t);
	else return returnValue ;

	
}


//...............................................................
//Calculate the normalisation for a given set of physics parameters and boundary

double Bs2JpsiPhi_SignalAlt_MP_v1::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary) 
{
		
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset;
	ctheta_tr = measurement->GetObservable( cosThetaName )->GetValue();
	phi_tr      = measurement->GetObservable( phiName )->GetValue();
	ctheta_1   = measurement->GetObservable( cosPsiName )->GetValue();	
	timeAcceptanceCategory = (int)measurement->GetObservable( timeAcceptanceCategoryName )->GetValue();
	
	// Get time boundaries into member variables
	IConstraint * timeBound = boundary->GetConstraint( timeConstraintName );
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
	if( ! normalisationCacheValid )  {
		for( tag = -1; tag <= 1; ++tag ) {
			if(resolution1Fraction >= 0.9999 ){
				resolution =  resolution1 ;
				normalisationCacheValueRes1[tag+1] = this->diffXsecCompositeNorm1( );
			}
			else {
				resolution =  resolution1 ;
				normalisationCacheValueRes1[tag+1] = this->diffXsecCompositeNorm1( );
				resolution =  resolution2 ;
				normalisationCacheValueRes2[tag+1] = this->diffXsecCompositeNorm1( );
			}
		}
		normalisationCacheValid = true ;
	}	
	
	// calculate return value according to tag 
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	double returnValue  ;
	if(resolution1Fraction >= 0.9999 )
	{
		returnValue = normalisationCacheValueRes1[tag+1] ;
	}
	else
	{
		returnValue = resolution1Fraction*normalisationCacheValueRes1[tag+1] + (1. - resolution1Fraction)*normalisationCacheValueRes2[tag+1] ;
	}
	
	if( (returnValue <= 0.) || isnan(returnValue) ) {
		cout << " Bs2JpsiPhi_SignalAlt_MP_v1::Normalisation() returns <=0 or nan " << returnValue << endl ;
		cout << " gamma " << gamma() ;
		cout << " gl    " << gamma_l() ;
		cout << " gh    " << gamma_h() ;
		cout << " AT    " << AT() ;
		cout << " AP    " << AP() ;
		cout << " A0    " << A0() ;
		cout << " AS    " << A0() ;
		throw 10 ;
	}
	
	return returnValue ;
}

