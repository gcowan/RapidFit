// $Id: Bs2Jpsifzero_SignalAlt_MP_dev.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2Jpsifzero_SignalAlt_MP_dev Bs2Jpsifzero_SignalAlt_MP_dev.cpp
 *
 *  RapidFit PDF for Bs2Jpsifzero
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2011-01-28
 */

#include "TMath.h"
#include <cmath>

#include "Bs2Jpsifzero_SignalAlt_MP_dev.h"
#include <iostream>
#include "math.h"
#include "Mathematics.h"

#include <float.h>

#define DEBUGFLAG true

PDF_CREATOR( Bs2Jpsifzero_SignalAlt_MP_dev );

//#define DOUBLE_TOLERANCE DBL_MIN
//#define DOUBLE_TOLERANCE 1E-6

//......................................
//Constructor(s)
//...........
// New with configurator
Bs2Jpsifzero_SignalAlt_MP_dev::Bs2Jpsifzero_SignalAlt_MP_dev( PDFConfigurator* config) : 
Bs2Jpsifzero_SignalAlt_BaseClass_dev(config)
{
	MakePrototypes();	
	std::cout << "Constructing PDF: Bs2Jpsifzero_SignalAlt_MP_dev " << std::endl ;
}


//.......................................
//Make the data point and parameter set
void Bs2Jpsifzero_SignalAlt_MP_dev::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( tagName );
	//X allObservables.push_back( timeAcceptanceCategoryName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName );
	parameterNames.push_back( deltaGammaName );
	parameterNames.push_back( deltaMName );

	if( _useCosAndSin ) {
		parameterNames.push_back( cosphisName );
		parameterNames.push_back( sinphisName );
	}
	else{
		parameterNames.push_back( Phi_sName );
	}
	
	parameterNames.push_back( mistagName );
	parameterNames.push_back( mistagP1Name );
	parameterNames.push_back( mistagP0Name );
	parameterNames.push_back( mistagSetPointName );

	parameterNames.push_back( resScaleName );
	parameterNames.push_back( res1Name );
	parameterNames.push_back( res2Name );
	parameterNames.push_back( res3Name );
	parameterNames.push_back( res2FractionName );
	parameterNames.push_back( res3FractionName );
	parameterNames.push_back( timeOffsetName );
	
	allParameters = ParameterSet(parameterNames);
}


//........................................................
//Destructor
Bs2Jpsifzero_SignalAlt_MP_dev::~Bs2Jpsifzero_SignalAlt_MP_dev()
{
}

//........................................................
//Set the physics parameters into member variables
//Indicate that the cache is no longer valid

bool Bs2Jpsifzero_SignalAlt_MP_dev::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	
	bool result = allParameters.SetPhysicsParameters(NewParameterSet);
	
	// Physics parameters. 
	_gamma  = allParameters.GetPhysicsParameter( gammaName )->GetValue();
    dgam      = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	
				

	delta_ms  = allParameters.GetPhysicsParameter( deltaMName )->GetValue();

	if(_useCosAndSin){
		_cosphis = allParameters.GetPhysicsParameter( cosphisName )->GetValue();
		_sinphis = allParameters.GetPhysicsParameter( sinphisName )->GetValue();
	}
	else{
		phi_s     = allParameters.GetPhysicsParameter( Phi_sName )->GetValue();
		_cosphis = cos(phi_s) ;
		_sinphis = sin(phi_s) ;
	}

	// mistag parameters
	_mistag			= allParameters.GetPhysicsParameter( mistagName )->GetValue();
	_mistagP1		= allParameters.GetPhysicsParameter( mistagP1Name )->GetValue();
	_mistagP0		= allParameters.GetPhysicsParameter( mistagP0Name )->GetValue();
	_mistagSetPoint = allParameters.GetPhysicsParameter( mistagSetPointName )->GetValue();
		
	// Detector parameters
	resolutionScale		= allParameters.GetPhysicsParameter( resScaleName )->GetValue();
	resolution1         = allParameters.GetPhysicsParameter( res1Name )->GetValue();
	resolution2         = allParameters.GetPhysicsParameter( res2Name )->GetValue();
	resolution3         = allParameters.GetPhysicsParameter( res3Name )->GetValue();
	resolution2Fraction = allParameters.GetPhysicsParameter( res2FractionName )->GetValue();
	resolution3Fraction = allParameters.GetPhysicsParameter( res3FractionName )->GetValue();
	timeOffset          = allParameters.GetPhysicsParameter( timeOffsetName )->GetValue();
	
	
	return result;
}

//.........................................................
//Return a list of observables not to be integrated
vector<string> Bs2Jpsifzero_SignalAlt_MP_dev::GetDoNotIntegrateList()
{
	vector<string> list;
	return list;
}

//.............................................................
//Calculate the PDF value for a given set of observables for use by numeric integral

double Bs2Jpsifzero_SignalAlt_MP_dev::EvaluateForNumericIntegral(DataPoint * measurement) 
{
	if( _numericIntegralTimeOnly ) return this->EvaluateTimeOnly(measurement) ;

	else return this->Evaluate(measurement) ;
}

//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2Jpsifzero_SignalAlt_MP_dev::Evaluate(DataPoint * measurement) 
{
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset ;
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	
	double val1=0. , val2=0., val3=0. ;
	double returnValue ;
	
	double resolution1Fraction = 1. - resolution2Fraction - resolution3Fraction ;
	
	if( resolutionScale <= 0. ) {
		resolution = 0. ;
		returnValue = this->diffXsec( );
	}
	else {		
		if(resolution1Fraction > 0 ) {
			resolution = resolution1 * resolutionScale ;
			val1 = this->diffXsec( );
		}
		if(resolution2Fraction > 0 ) {
			resolution = resolution2 * resolutionScale ;
			val2 = this->diffXsec( );
		}
		if(resolution3Fraction > 0 ) {
			resolution = resolution3 * resolutionScale ;
			val3 = this->diffXsec( );
		}
		returnValue = resolution1Fraction*val1 + resolution2Fraction*val2 + resolution3Fraction*val3 ;				
	}
		
	//conditions to throw exception
	bool c1 = std::isnan(returnValue) ;
	bool c2 = (resolutionScale> 0.) && (returnValue <= 0.) ;
	bool c3 = (resolutionScale<=0.) && (t>0.) && (returnValue <= 0.)  ;	
	if( DEBUGFLAG && (c1 || c2 || c3)  ) {
		cout << endl ;
		cout << " Bs2Jpsifzero_SignalAlt_MP_dev::evaluate() returns <=0 or nan :" << returnValue << endl ;
		cout << "   gamma " << gamma() << endl ;
		cout << "   gl    " << gamma_l() << endl ;
		cout << "   gh    " << gamma_h()  << endl;
		cout << "   AT^2    " << AT()*AT() << endl;
		cout << "   AP^2    " << AP()*AP() << endl;
		cout << "   A0^2    " << A0()*A0() << endl ;
		cout << "   AS^2    " << AS()*AS() << endl ;
		cout << "   ATOTAL  " << AS()*AS()+A0()*A0()+AP()*AP()+AT()*AT() << endl ;
		cout << "   delta_ms       " << delta_ms << endl ;
		cout << "   mistag         " << mistag() << endl ;
		cout << "   mistagP1       " << _mistagP1 << endl ;
		cout << "   mistagP0       " << _mistagP0 << endl ;
		cout << "   mistagSetPoint " << _mistagSetPoint << endl ;
		cout << "   resolutionScale " << resolutionScale << endl; 
		cout << "   resolution1Fraction " << resolution1Fraction << endl; 
		cout << "   resolution2Fraction " << resolution2Fraction << endl; 
		cout << "   resolution3Fraction " << resolution3Fraction << endl; 
		cout << "   val1 " << val1 << endl; 
		cout << "   val2 " << val2 << endl; 
		cout << "   val3 " << val3 << endl; 
		cout << " For event with:  " << endl ;
		cout << "   time      " << t << endl ;
		cout << "   ctheta_tr " << ctheta_tr << endl ;
		cout << "   ctheta_1 " << ctheta_1 << endl ;
		cout << "   phi_tr " << phi_tr << endl ;
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
	
	return returnValue ;	
}


//.............................................................
//Calculate the PDF value for a given time, but integrated over angles

double Bs2Jpsifzero_SignalAlt_MP_dev::EvaluateTimeOnly(DataPoint * measurement) 
{
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset ;
	//ctheta_tr = measurement->GetObservable( cosThetaName )->GetValue();
	//phi_tr      = measurement->GetObservable( phiName )->GetValue();
	//ctheta_1   = measurement->GetObservable( cosPsiName )->GetValue();	
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	
	double val1=0. , val2=0., val3=0. ;
	double returnValue ;
	
	double resolution1Fraction = 1. - resolution2Fraction - resolution3Fraction ;
	
	if( resolutionScale <= 0. ) {
		resolution = 0. ;
		returnValue = this->diffXsecTimeOnly( );
	}
	else {				
		if(resolution1Fraction > 0 ) {
			resolution = resolution1 * resolutionScale ;
			val1 = this->diffXsecTimeOnly( );
		}
		if(resolution2Fraction > 0 ) {
			resolution = resolution2 * resolutionScale ;
			val2 = this->diffXsecTimeOnly( );
		}
		if(resolution3Fraction > 0 ) {
			resolution = resolution3 * resolutionScale ;
			val3 = this->diffXsecTimeOnly( );
		}
		returnValue = resolution1Fraction*val1 + resolution2Fraction*val2 + resolution3Fraction*val3 ;	
	}
	
	
	//conditions to throw exception
	bool c1 = std::isnan(returnValue) ;
	bool c2 = (resolutionScale> 0.) && (returnValue <= 0.) ;
	bool c3 = (resolutionScale<=0.) && (t>0.) && (returnValue <= 0.)  ;	
	if( DEBUGFLAG && (c1 || c2 || c3)  ) {
		cout << endl ;
		cout << " Bs2Jpsifzero_SignalAlt_MP_dev::EvaluateTimeOnly() returns <=0 or nan :" << returnValue << endl ;
		cout << "   gamma " << gamma() << endl ;
		cout << "   gl    " << gamma_l() << endl ;
		cout << "   gh    " << gamma_h()  << endl;
		cout << "   AT^2    " << AT()*AT() << endl;
		cout << "   AP^2    " << AP()*AP() << endl;
		cout << "   A0^2    " << A0()*A0() << endl ;
		cout << "   AS^2    " << AS()*AS() << endl ;
		cout << "   ATOTAL  " << AS()*AS()+A0()*A0()+AP()*AP()+AT()*AT() << endl ;
		cout << "   delta_ms       " << delta_ms << endl ;
		cout << "   mistag         " << mistag() << endl ;
		cout << "   mistagP1       " << _mistagP1 << endl ;
		cout << "   mistagP0       " << _mistagP0 << endl ;
		cout << "   mistagSetPoint " << _mistagSetPoint << endl ;
		cout << " For event with:  " << endl ;
		cout << "   time      " << t << endl ;
		cout << "   ctheta_tr " << ctheta_tr << endl ;
		cout << "   ctheta_1 " << ctheta_1 << endl ;
		cout << "   phi_tr " << phi_tr << endl ;
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
	
	return returnValue ;	
}


//...............................................................
//Calculate the normalisation for a given set of physics parameters and boundary

double Bs2Jpsifzero_SignalAlt_MP_dev::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary) 
{
		
	if( _numericIntegralForce ) return -1. ;
	
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset;
	//X timeAcceptanceCategory = (int)measurement->GetObservable( timeAcceptanceCategoryName )->GetValue();
	
	// Get time boundaries into member variables
	IConstraint * timeBound = boundary->GetConstraint( timeConstraintName );
	if ( timeBound->GetUnit() == "NameNotFoundError" ) {
		cerr << "Bound on time not provided" << endl;
		exit(1);
	}
	else {
		tlo = timeBound->GetMinimum();
		thi = timeBound->GetMaximum();
	}
	
	// Recalculate cached values if Physics parameters have changed
	// Must do this for each of the two resolutions.
	if( ! normalisationCacheValid )  {
		for( tag = -1; tag <= 1; ++tag ) {

			double val1=0. , val2=0., val3=0. ;
			double result ;
			
			double resolution1Fraction = 1. - resolution2Fraction - resolution3Fraction ;
			
			if( resolutionScale <= 0. ) {
				resolution = 0. ;
				result = this->diffXsecCompositeNorm1( );
			}
			else {						
				if(resolution1Fraction > 0 ) {
					resolution = resolution1 * resolutionScale ;
					val1 = this->diffXsecCompositeNorm1( );
				}
				if(resolution2Fraction > 0 ) {
					resolution = resolution2 * resolutionScale ;
					val2 = this->diffXsecCompositeNorm1( );
				}
				if(resolution3Fraction > 0 ) {
					resolution = resolution3 * resolutionScale ;
					val3 = this->diffXsecCompositeNorm1( );
				}
				result = resolution1Fraction*val1 + resolution2Fraction*val2 + resolution3Fraction*val3 ;	
			}
			normalisationCacheValue[tag+1] = result ;
			
		}
		normalisationCacheValid = true ;
	}	
	
	// calculate return value according to tag 
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	double returnValue = normalisationCacheValue[tag+1] ;

	//conditions to throw exception
	bool c1 = std::isnan(returnValue)  ;
	bool c2 = (returnValue <= 0.) ;	
	if( DEBUGFLAG && (c1 || c2 ) ) {
		cout << endl ;
		cout << " Bs2Jpsifzero_SignalAlt_MP_dev::Normaisation() returns <=0 or nan :" << returnValue << endl ;
		cout << "   gamma " << gamma() << endl ;
		cout << "   gl    " << gamma_l() << endl ;
		cout << "   gh    " << gamma_h()  << endl;
		cout << "   AT^2    " << AT()*AT() << endl;
		cout << "   AP^2    " << AP()*AP() << endl;
		cout << "   A0^2    " << A0()*A0() << endl ;
		cout << "   AS^2    " << AS()*AS() << endl ;
		cout << "   ATOTAL  " << AS()*AS()+A0()*A0()+AP()*AP()+AT()*AT() << endl ;
		cout << "   delta_ms       " << delta_ms << endl ;
		cout << "   mistag         " << mistag() << endl ;
		cout << "   mistagP1       " << _mistagP1 << endl ;
		cout << "   mistagP0       " << _mistagP0 << endl ;
		cout << "   mistagSetPoint " << _mistagSetPoint << endl ;
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}

	return returnValue ;
}

