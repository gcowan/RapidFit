// $Id: Bs2JpsiPhiPromptBkg_withTimeResDouble.cpp,v 1.2 2009/11/13 15:31:52 gcowan Exp $
/** @class Bs2JpsiPhiPromptBkg_withTimeResDouble Bs2JpsiPhiPromptBkg_withTimeResDouble.cpp
 *
 *  PDF for Bs2JpsiPhi prompt background
 *
 *  @author Greig A Cowan greig.cowan@cern.ch
 *  @date 2009-11-12
 */
#include "Bs2JpsiPhiPromptBkg_withTimeResDouble.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( Bs2JpsiPhiPromptBkg_withTimeResDouble );

//Constructor
Bs2JpsiPhiPromptBkg_withTimeResDouble::Bs2JpsiPhiPromptBkg_withTimeResDouble(PDFConfigurator* configurator) : 
	// Physics parameters
	  frac_sigmaPrName	( configurator->getName("frac_sigmaPr") )
        , sigmaPrName		( configurator->getName("sigmaPr") )
        , sigmaPr2Name		( configurator->getName("sigmaPr2") )
        // Observables
        ,  timeName		( configurator->getName("time") )
	,  timeconstraintName	( configurator->getName("time") )
{
	MakePrototypes();
	cout << "Constructing Bs2JpsiPhi prompt background with double time resolution" << endl;
}

//Make the data point and parameter set
void Bs2JpsiPhiPromptBkg_withTimeResDouble::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( frac_sigmaPrName );
	parameterNames.push_back( sigmaPrName );
	parameterNames.push_back( sigmaPr2Name );
	allParameters = ParameterSet(parameterNames);
}

//Destructor
Bs2JpsiPhiPromptBkg_withTimeResDouble::~Bs2JpsiPhiPromptBkg_withTimeResDouble()
{
}


//Calculate the function value
double Bs2JpsiPhiPromptBkg_withTimeResDouble::Evaluate(DataPoint * measurement)
{
	double time = measurement->GetObservable( timeName )->GetValue();

	double sigmaPr = allParameters.GetPhysicsParameter( sigmaPrName )->GetValue();
	double timeNorm = 1./( sigmaPr * sqrt( 2.*TMath::Pi() ) );
	double gauss    = exp( -time*time / ( 2. * sigmaPr * sigmaPr ) );

	double sigmaPr2 = allParameters.GetPhysicsParameter( sigmaPr2Name )->GetValue();
	double timeNorm2 = 1./( sigmaPr2 * sqrt( 2.*TMath::Pi() ) );
	double gauss2    = exp( -time*time / ( 2. * sigmaPr2 * sigmaPr2 ) );
	
	double frac = allParameters.GetPhysicsParameter( frac_sigmaPrName )->GetValue();
	
	return frac*(timeNorm * gauss) + (1. - frac)*(timeNorm2 * gauss2) ;
}

// int(1/(sigmaPr*sqrt(2*Pi))*exp(-(time)^2/(2*sigmaPr*sigmaPr)),time=tmin..tmax);
//                                                                   1/2                  1/2
//                                                                  2    tmin            2    tmax
//                                                         -1/2 erf(---------) + 1/2 erf(---------)
//                                                                  2 sigmaPr            2 sigmaPr
double Bs2JpsiPhiPromptBkg_withTimeResDouble::Normalisation(PhaseSpaceBoundary * boundary)
{
        double tmin = 0.;
        double tmax = 0.;
        IConstraint * timeBound = boundary->GetConstraint( timeconstraintName );
        if ( timeBound->GetUnit() == "NameNotFoundError" )
        {
                cerr << "Bound on time not provided" << endl;
                return -1.;
        }
        else
        {
                tmin = timeBound->GetMinimum();
                tmax = timeBound->GetMaximum();
        }

		double sigmaPr = allParameters.GetPhysicsParameter( sigmaPrName )->GetValue();
		double val = 0.5 * ( erf( tmax/(sqrt(2.)*sigmaPr) ) - erf( tmin/(sqrt(2.)*sigmaPr )) );

		double sigmaPr2 = allParameters.GetPhysicsParameter( sigmaPr2Name )->GetValue();
		double val2 = 0.5 * ( erf( tmax/(sqrt(2.)*sigmaPr2) ) - erf( tmin/(sqrt(2.)*sigmaPr2 )) );

		double frac = allParameters.GetPhysicsParameter( frac_sigmaPrName )->GetValue();

		return frac*val + (1.-frac)*val2 ;
}
