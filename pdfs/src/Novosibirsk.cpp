// $Id: Novosibirsk.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Novosibirsk Novosibirsk.cpp
 *
 *  RapidFit PDF for Bs mass
 *
 *
 *  @author Pete
 *  @date 2011-07-30
 */

#include "Novosibirsk.h"
#include <iostream>
#include "TMath.h"
#include "RooMath.h"
#include "math.h"

PDF_CREATOR( Novosibirsk );

//Constructor
Novosibirsk::Novosibirsk(PDFConfigurator* configurator) :
	// Physics parameters
	BasePDF()
	, widthName	( configurator->getName("width") )
	, peakName	( configurator->getName("peak") )
	, tailName	( configurator->getName("tail") )
	// Observables
	, xName	( configurator->getName("x") )
{
	MakePrototypes();
}

//Make the data point and parameter set
void Novosibirsk::MakePrototypes()
{
	// Observables
	allObservables.push_back( xName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( widthName );
	parameterNames.push_back( peakName );
	parameterNames.push_back( tailName );
	allParameters = ParameterSet(parameterNames);
}

Novosibirsk::Novosibirsk( const Novosibirsk &copy ) :
        BasePDF( (BasePDF)copy )
        , widthName ( copy.widthName )
        , peakName ( copy.peakName )
        , tailName ( copy.tailName )
        , xName ( copy.xName )
	, width ( copy.width )
	, peak ( copy.peak )
	, tail ( copy.tail )
{
}

//Destructor
Novosibirsk::~Novosibirsk()
{
}

//Calculate the function value
double Novosibirsk::Evaluate(DataPoint * measurement)
{
	// Get the physics parameters
	width = allParameters.GetPhysicsParameter( widthName )->GetValue();
	peak  = allParameters.GetPhysicsParameter( peakName )->GetValue();
	tail  = allParameters.GetPhysicsParameter( tailName )->GetValue();

	// Get the observable
	double x = measurement->GetObservable( xName )->GetValue();


	double qa=0,qb=0,qc=0,qx=0,qy=0;

	if(TMath::Abs(tail) < 1.e-7)
		qc = 0.5*TMath::Power(((x-peak)/width),2);
	else {
		qa = tail*sqrt(log(4.));
		qb = sinh(qa)/qa;
		qx = (x-peak)/width*qb;
		qy = 1.+tail*qx;

		//---- Cutting curve from right side

		if( qy > 1.E-7)
			qc = 0.5*(TMath::Power((log(qy)/tail),2) + tail*tail);
		else
			qc = 15.0;
	}

	//cout << qa << " " <<qb << " " <<qc << " " <<qx << " " <<qy<< endl;

	//---- Normalize the result

	return exp(-qc);

}

double Novosibirsk::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void) boundary;
	
	// Get the physics parameters
	width = allParameters.GetPhysicsParameter( widthName )->GetValue();
	peak  = allParameters.GetPhysicsParameter( peakName )->GetValue();
	tail  = allParameters.GetPhysicsParameter( tailName )->GetValue();

	double tailLog4 = tail * sqrt(log(4.));
	double widthLog4 = width * sqrt(log(4.));
	double _sinh = sinh(tailLog4);
	double _csch = 1./sinh(tailLog4);
	
	double val(0.);

	double xhigh = 0.12;
	double xlow = 0.;
	
	//cout << - width*tail*sqrt(TMath::Pi()*log(2.))*_csch*erf( ( tail*tail - log( (xhigh - peak)*_sinh/widthLog4 + 1. ) ) / (sqrt(2.)*tail) ) << endl;
	//cout << + width*tail*sqrt(TMath::Pi()*log(2.))*_csch*erf( ( tail*tail - log( (xlow  - peak)*_sinh/widthLog4 + 1. ) ) / (sqrt(2.)*tail) ) << endl;

	val = - width*tail*sqrt(TMath::Pi()*log(2.))*_csch*erf( ( tail*tail - log( (xhigh - peak)*_sinh/widthLog4 + 1. ) ) / (sqrt(2.)*tail) );
	      //+ width*tail*sqrt(TMath::Pi()*log(2.))*_csch*erf( ( tail*tail - log( (xlow  - peak)*_sinh/widthLog4 + 1. ) ) / (sqrt(2.)*tail) );

	return 2.*val;
}

