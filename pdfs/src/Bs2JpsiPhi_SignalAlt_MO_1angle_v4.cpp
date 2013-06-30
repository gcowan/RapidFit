// $Id: Bs2JpsiPhi_SignalAlt_MO_1angle_v4.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2JpsiPhi_SignalAlt_MO_1angle_v4 Bs2JpsiPhi_SignalAlt_MO_1angle_v4.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-13
 */

#include "TMath.h"
#include <cmath>

#include "Bs2JpsiPhi_SignalAlt_MO_1angle_v4.h"
#include <iostream>
#include "math.h"
#include "Mathematics.h"

#define DEBUGFLAG true

PDF_CREATOR( Bs2JpsiPhi_SignalAlt_MO_1angle_v4 );

//......................................
//Constructor(s)
//New one with configurator
Bs2JpsiPhi_SignalAlt_MO_1angle_v4::Bs2JpsiPhi_SignalAlt_MO_1angle_v4(PDFConfigurator* configurator) : 
Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4(configurator)
{
	MakePrototypes();	
	std::cout << "Constructing PDF: Bs2JpsiPhi_SignalAlt_MO_1angle_v4 " << std::endl ;
}


//......................................
//Make the data point and parameter set
void Bs2JpsiPhi_SignalAlt_MO_1angle_v4::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( cosThetaName );
	allObservables.push_back( tagName );
	allObservables.push_back( mistagName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName );
	parameterNames.push_back( deltaGammaName );
	//parameterNames.push_back( Aeven_sqName );
	parameterNames.push_back( Aodd_sqName );
	parameterNames.push_back( As_sqName );
	parameterNames.push_back( deltaMName );

	if( _useCosAndSin ) {
		parameterNames.push_back( cosphisName );
		parameterNames.push_back( sinphisName );
	}
	else{
		parameterNames.push_back( Phi_sName );
	}

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
	
	parameterNames.push_back( angAccEvenName );
	parameterNames.push_back( angAccOddName );
	allParameters = ParameterSet(parameterNames);
}


//........................................................
//Destructor
Bs2JpsiPhi_SignalAlt_MO_1angle_v4::~Bs2JpsiPhi_SignalAlt_MO_1angle_v4() { }


//........................................................
//Set the physics parameters into member variables

bool Bs2JpsiPhi_SignalAlt_MO_1angle_v4::SetPhysicsParameters( ParameterSet * NewParameterSet ) 
{
	
	bool result = Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::SetPhysicsParameters( NewParameterSet ) ;
	
	
	return result;
}

//.........................................................
//Return a list of observables not to be integrated
vector<string> Bs2JpsiPhi_SignalAlt_MO_1angle_v4::GetDoNotIntegrateList()
{
	vector<string> list;
	list.push_back(mistagName) ;
	if( _numericIntegralTimeOnly ) {
		list.push_back( cosThetaName );
	}	
	return list;
}


//.............................................................
//Calculate the PDF value for a given set of observables for use by numeric integral

double Bs2JpsiPhi_SignalAlt_MO_1angle_v4::EvaluateForNumericIntegral(DataPoint * measurement) 
{
	if( _numericIntegralTimeOnly ) return this->EvaluateTimeOnly(measurement) ;
	
	else return this->Evaluate(measurement) ;
}


//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2JpsiPhi_SignalAlt_MO_1angle_v4::Evaluate(DataPoint * measurement)
{
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset ;
	ctheta_tr = measurement->GetObservable( cosThetaName )->GetValue();
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	_mistag = measurement->GetObservable( mistagName )->GetValue();

	//Cache amplitues and angles terms used in cross section
	this->CacheAmplitudesAndAngles() ;
	
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
	if( DEBUGFLAG && (c1 || c2 || c3)  ) 
	{
		cout << endl ;
		cout << " Bs2JpsiPhi_SignalAlt_MO_1angle_v4::evaluate() returns <=0 or nan :" << returnValue << endl ;
		cout << "   gamma " << gamma() << endl ;
		cout << "   gl    " << gamma_l() << endl ;
		cout << "   gh    " << gamma_h()  << endl;
		cout << "   Aeven^2    " << Aeven_sq << endl;
		cout << "   Aodd^2    " << Aodd_sq << endl;
		cout << "   AS^2    " << As_sq << endl ;
		cout << "   ATOTAL  " <<  Aeven_sq+Aodd_sq+As_sq << endl ;
		cout << "   delta_ms       " << delta_ms << endl ;
		cout << "   mistag         " << mistag() << endl ;
		cout << "   mistagP1       " << _mistagP1 << endl ;
		cout << "   mistagP0       " << _mistagP0 << endl ;
		cout << "   mistagSetPoint " << _mistagSetPoint << endl ;
		cout << " For event with:  " << endl ;
		cout << "   time      " << t << endl ;
		cout << "   ctheta_tr " << ctheta_tr << endl ;
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
		
	return returnValue ;	
}


//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2JpsiPhi_SignalAlt_MO_1angle_v4::EvaluateTimeOnly(DataPoint * measurement)
{
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset ;
	//ctheta_tr = measurement->GetObservable( cosThetaName )->GetValue();
	//phi_tr      = measurement->GetObservable( phiName )->GetValue();
	//ctheta_1   = measurement->GetObservable( cosPsiName )->GetValue();
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	_mistag = measurement->GetObservable( mistagName )->GetValue();
	
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
		cout << " Bs2JpsiPhi_SignalAlt_MO_1angle_v4::evaluate() returns <=0 or nan :" << returnValue << endl ;
		cout << "   gamma " << gamma() << endl ;
		cout << "   gl    " << gamma_l() << endl ;
		cout << "   gh    " << gamma_h()  << endl;
		cout << "   Aeven^2    " << Aeven_sq << endl;
		cout << "   Aodd^2    " << Aodd_sq << endl;
		cout << "   AS^2    " << As_sq << endl ;
		cout << "   ATOTAL  " <<  Aeven_sq+Aodd_sq+As_sq << endl ;
		cout << "   delta_ms       " << delta_ms << endl ;
		cout << "   mistag         " << mistag() << endl ;
		cout << "   mistagP1       " << _mistagP1 << endl ;
		cout << "   mistagP0       " << _mistagP0 << endl ;
		cout << "   mistagSetPoint " << _mistagSetPoint << endl ;
		cout << " For event with:  " << endl ;
		cout << "   time      " << t << endl ;
		cout << "   ctheta_tr " << ctheta_tr << endl ;
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
	
	
	return returnValue ;
	
}


//...............................................................
//Calculate the normalisation for a given set of physics parameters and boundary

double Bs2JpsiPhi_SignalAlt_MO_1angle_v4::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{

	if( _numericIntegralForce ) return -1. ;

	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset;
	ctheta_tr = measurement->GetObservable( cosThetaName )->GetValue();
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	_mistag = measurement->GetObservable( mistagName )->GetValue() ;
	
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
	
	//First job for any new set of parameters is to Cache the time integrals
	if( ! timeIntegralCacheValid ) {
		CacheTimeIntegrals() ;
		timeIntegralCacheValid = true ;
	}
	
	
	//If this is an untagged event and the result has been cached, then it can be used
	// Otherwise must calculate the normalisation

	double returnValue ;
	
	if( (tag==0) && normalisationCacheValid ) {
		returnValue = normalisationCacheUntagged ;
	}
	
	else {
			
		double val1=0. , val2=0., val3=0. ;		
		double resolution1Fraction = 1. - resolution2Fraction - resolution3Fraction ;

		if( resolutionScale <= 0. ) {
			resolution = 0. ;
			returnValue = this->diffXsecCompositeNorm1( 0 );
		}
		else {						
			if(resolution1Fraction > 0 ) {
				resolution = resolution1 * resolutionScale ;
				val1 = this->diffXsecCompositeNorm1( 1 );
			}
			if(resolution2Fraction > 0 ) {
				resolution = resolution2 * resolutionScale ;
				val2 = this->diffXsecCompositeNorm1( 2 );
			}
			if(resolution3Fraction > 0 ) {
				resolution = resolution3 * resolutionScale ;
				val3 = this->diffXsecCompositeNorm1( 3 );
			}
			returnValue = resolution1Fraction*val1 + resolution2Fraction*val2 + resolution3Fraction*val3 ;	
		}
	}
	
	// If this is an untagged event then the normaisation is invariant and so can be cached
	if( (tag==0) && !normalisationCacheValid )  {
		normalisationCacheUntagged = returnValue ;
		normalisationCacheValid = true ;
	}
	
	// Conditions to throw exception
	bool c1 = std::isnan(returnValue)  ;
	bool c2 = (returnValue <= 0.) ;	
	if( DEBUGFLAG && (c1 || c2 ) ) {
		cout << endl ;
		cout << " Bs2JpsiPhi_SignalAlt_MO_1angle_v4::normalisation() returns <=0 or nan :" << returnValue << endl ;
		cout << "   gamma " << gamma() << endl ;
		cout << "   gl    " << gamma_l() << endl ;
		cout << "   gh    " << gamma_h()  << endl;
		cout << "   Aeven^2    " << Aeven_sq << endl;
		cout << "   Aodd^2    " << Aodd_sq << endl;
		cout << "   AS^2    " << As_sq << endl ;
		cout << "   ATOTAL  " <<  Aeven_sq+Aodd_sq+As_sq << endl ;
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



