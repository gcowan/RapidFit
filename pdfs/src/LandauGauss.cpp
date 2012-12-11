
#include "TMath.h"
#include "LandauGauss.h"
#include "Mathematics.h"
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

PDF_CREATOR( LandauGauss );

LandauGauss::LandauGauss( PDFConfigurator* config ) :
	xName( config->getName("x") ),
	landauSigmaName( config->getName("LandauSigma") ),
	gaussSigmaName( config->getName("GaussSigma") ),
	mpvName( config->getName("mpv") ),
	CName( config->getName("C") )
{
	this->MakePrototypes();
}

void LandauGauss::MakePrototypes()
{
	allObservables.push_back( string(xName) );
	vector<string> allParams;
	allParams.push_back( landauSigmaName ); allParams.push_back( gaussSigmaName );
	allParams.push_back( mpvName ); allParams.push_back( CName );
	allParameters = ParameterSet( allParams );
}

LandauGauss::~LandauGauss()
{}

double LandauGauss::Evaluate( DataPoint* input )
{
	double LandauSigmaVal = allParameters.GetPhysicsParameter( landauSigmaName )->GetValue();
	double GaussSigmaVal = allParameters.GetPhysicsParameter( gaussSigmaName )->GetValue();
	double mpvValue = allParameters.GetPhysicsParameter( mpvName )->GetValue();
	double CVal = allParameters.GetPhysicsParameter( CName )->GetValue();
	double xVal = input->GetObservable( xName )->GetValue();

	//double numerator = CVal * TMath::Landau( xVal, mpvValue, sigmaVal, true );
	//double numerator = 351 *  TMath::Gaus(xVal-mpvValue,sigmaVal,sqrt(sigmaVal),true);

	// Range of convolution integral
	double xlow  = xVal - 5. * GaussSigmaVal;
	double xhigh = xVal + 5. * GaussSigmaVal;

	double stepSize = (xhigh-xlow) / 1000.;

	//	Convolution integral of Landau and Gaussian by sum
	double landau_part=0.;

	double sum=0.;

	double running_xlow=xlow-0.5*stepSize;
	double running_xhigh=xhigh+0.5*stepSize;

	for( unsigned int i=1; i<= 500; ++i )
	{
		running_xlow += stepSize;
		landau_part = TMath::Landau( running_xlow, mpvValue, LandauSigmaVal, true );
		sum += landau_part * TMath::Gaus( xVal, running_xlow, GaussSigmaVal );

		running_xhigh -= stepSize;
		landau_part = TMath::Landau( running_xhigh, mpvValue, LandauSigmaVal, true );
		sum += landau_part * TMath::Gaus( xVal, running_xhigh, GaussSigmaVal );
	}

	double numerator = CVal * stepSize * sum / ( GaussSigmaVal * LandauSigmaVal );

	if( numerator < 1E-9 ) numerator=1E-9;

	return numerator;
}

double LandauGauss::Normalisation( PhaseSpaceBoundary* range )
{
	(void) range;
	return 1.;
}

