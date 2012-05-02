// $Id: Bs2JpsiPhiPromptBkg_tripleGaussian.cpp,v 1.2 2009/11/13 15:31:52 gcowan Exp $
/** @class Bs2JpsiPhiPromptBkg_tripleGaussian Bs2JpsiPhiPromptBkg_tripleGaussian.cpp
 *
 *  PDF for Bs2JpsiPhi prompt background
 *
 *  @author Pete Clarke 
 *  @date 2010 01 22
 */
#include "Bs2JpsiPhiPromptBkg_tripleGaussian.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( Bs2JpsiPhiPromptBkg_tripleGaussian );

//Constructor
Bs2JpsiPhiPromptBkg_tripleGaussian::Bs2JpsiPhiPromptBkg_tripleGaussian( PDFConfigurator* configurator ) : 
	// Physics parameters
	  frac_sigmaPr1Name	( configurator->getName("frac_sigmaPr1") )
	, frac_sigmaPr23Name	( configurator->getName("frac_sigmaPr23" ))
	, sigmaPr1Name	( configurator->getName("sigmaPr1") )
	, sigmaPr2Name	( configurator->getName("sigmaPr2") )
	, sigmaPr3Name	( configurator->getName("sigmaPr3") )
	// Observables
	, timeName	( configurator->getName("time") )	
	, timeconstraintName( configurator->getName("time") )
{
	MakePrototypes();
	cout << "Constructing Bs2JpsiPhi prompt background with triple Gaussian" << endl;
}

//Make the data point and parameter set
void Bs2JpsiPhiPromptBkg_tripleGaussian::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( frac_sigmaPr1Name );
	parameterNames.push_back( frac_sigmaPr23Name );
	parameterNames.push_back( sigmaPr1Name );
	parameterNames.push_back( sigmaPr2Name );
	parameterNames.push_back( sigmaPr3Name );
	allParameters = ParameterSet(parameterNames);
}

//Destructor
Bs2JpsiPhiPromptBkg_tripleGaussian::~Bs2JpsiPhiPromptBkg_tripleGaussian()
{
}


//Calculate the function value
double Bs2JpsiPhiPromptBkg_tripleGaussian::Evaluate(DataPoint * measurement)
{
	//Retrieve Observables
	double time = measurement->GetObservable( timeName )->GetValue();

	//Retrieve parameters
	double sigmaPr1  = allParameters.GetPhysicsParameter( sigmaPr1Name )->GetValue();
	double sigmaPr2 =  allParameters.GetPhysicsParameter( sigmaPr2Name )->GetValue();
	double sigmaPr3 =  allParameters.GetPhysicsParameter( sigmaPr3Name )->GetValue();	
	double frac1    = allParameters.GetPhysicsParameter( frac_sigmaPr1Name )->GetValue();
	double frac23   = allParameters.GetPhysicsParameter( frac_sigmaPr23Name )->GetValue();
	
	//Reality Checks
	if( frac1 < 0. || frac1 > 1.001 )  cout << "Bs2JpsiPhiPromptBkg_tripleGaussian::Evaluate() : frac1 = " << frac1 << endl ;
	if( frac23 < 0. || frac23 > 1.001 )  cout << "Bs2JpsiPhiPromptBkg_tripleGaussian::Evaluate() : frac23 = " << frac23 << endl ;
	if( sigmaPr1 <= 0. ) cout << "Bs2JpsiPhiPromptBkg_tripleGaussian::Evaluate() : sigmaPr1 < 0 : " << sigmaPr1 << endl ;
	if( sigmaPr2 <= 0. ) cout << "Bs2JpsiPhiPromptBkg_tripleGaussian::Evaluate() : sigmaPr2 < 0 : " << sigmaPr2 << endl ;
	if( sigmaPr3 <= 0. ) cout << "Bs2JpsiPhiPromptBkg_tripleGaussian::Evaluate() : sigmaPr3 < 0 : " << sigmaPr3 << endl ;
	
	//Calculate PDF
	double gauss1=0, gauss2=0, gauss3=0 ;
	
	if( sigmaPr1 > 0. ) {
		double timeNorm = 1./( sigmaPr1 * sqrt( 2.*TMath::Pi() ) );
		gauss1    = timeNorm * exp( -time*time / ( 2. * sigmaPr1 * sigmaPr1 ) );
	}

	if( sigmaPr2 > 0. ) {
		double timeNorm = 1./( sigmaPr2 * sqrt( 2.*TMath::Pi() ) );
		gauss2    = timeNorm * exp( -time*time / ( 2. * sigmaPr2 * sigmaPr2 ) );
	}
	
	if( sigmaPr3 > 0. ) {
		double timeNorm = 1./( sigmaPr3 * sqrt( 2.*TMath::Pi() ) );
		gauss3    = timeNorm * exp( -time*time / ( 2. * sigmaPr3 * sigmaPr3 ) );
	}
		
	return frac1*gauss1 + (1.-frac1)*( frac23*gauss2 + (1.-frac23)*gauss3 ) ;
}

// int(1/(sigmaPr*sqrt(2*Pi))*exp(-(time)^2/(2*sigmaPr*sigmaPr)),time=tmin..tmax);
//                                                                   1/2                  1/2
//                                                                  2    tmin            2    tmax
//                                                         -1/2 erf(---------) + 1/2 erf(---------)
//                                                                  2 sigmaPr            2 sigmaPr
double Bs2JpsiPhiPromptBkg_tripleGaussian::Normalisation(PhaseSpaceBoundary * boundary)
{
	double tmin = 0.;
	double tmax = 0.;
	IConstraint * timeBound = boundary->GetConstraint(timeconstraintName);
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		exit(1) ;
	}
	else
	{
		tmin = timeBound->GetMinimum();
		tmax = timeBound->GetMaximum();
	}

	//Extract parameters
	double sigmaPr1 = allParameters.GetPhysicsParameter( sigmaPr1Name )->GetValue();
	double sigmaPr2 = allParameters.GetPhysicsParameter( sigmaPr2Name )->GetValue();
	double sigmaPr3 = allParameters.GetPhysicsParameter( sigmaPr3Name )->GetValue();
	double frac1 = allParameters.GetPhysicsParameter( frac_sigmaPr1Name )->GetValue();
	double frac23 = allParameters.GetPhysicsParameter( frac_sigmaPr23Name )->GetValue();
	
	//Reality checks
	if( sigmaPr1 <= 0. ) {cout << "Bs2JpsiPhiPromptBkg_tripleGaussian::Normalisation() : sigmaPr1 < 0 : " << sigmaPr1 << endl ;exit(1); }
	if( sigmaPr2 <= 0. ) {cout << "Bs2JpsiPhiPromptBkg_tripleGaussian::Normalisation() : sigmaPr2 < 0 : " << sigmaPr2 << endl ; exit(1); }
	if( sigmaPr3 <= 0. ) {cout << "Bs2JpsiPhiPromptBkg_tripleGaussian::Normalisation() : sigmaPr3 < 0 : " << sigmaPr3 << endl ; exit(1); }
	
	double val1 = 0.5 * ( erf( tmax/(sqrt(2.)*sigmaPr1) ) - erf( tmin/(sqrt(2.)*sigmaPr1 )) );
	double val2 = 0.5 * ( erf( tmax/(sqrt(2.)*sigmaPr2) ) - erf( tmin/(sqrt(2.)*sigmaPr2 )) );
	double val3 = 0.5 * ( erf( tmax/(sqrt(2.)*sigmaPr3) ) - erf( tmin/(sqrt(2.)*sigmaPr3 )) );
	
	return frac1*val1 + (1.-frac1)*(frac23*val2 + (1.-frac23)*val3) ;
}
