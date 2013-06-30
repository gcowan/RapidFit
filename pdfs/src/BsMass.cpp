// $Id: BsMass.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/*! @class BsMass BsMass.cpp
 *
 *  RapidFit PDF for Bs mass
 *
 *
 *  @author Pete
 *  @date 2011-07-30
 */

#include "TMath.h"
#include <cmath>

#include "BsMass.h"
#include "Mathematics.h"

#include <iostream>

using namespace::std;

PDF_CREATOR( BsMass );

//Constructor
BsMass::BsMass(PDFConfigurator* configurator) :
	// Physics parameters
	f_sig_m1Name	( configurator->getName("f_sig_m1") )
	, sigma_m1Name	( configurator->getName("sigma_m1") )
	, sigma_m2Name	( configurator->getName("sigma_m2") )
	, ratio_21Name	( configurator->getName("ratio_21") )
	, m_BsName		( configurator->getName("m_Bs") )
	// Observables
	, recoMassName	( configurator->getName("mass") )
	// Internal Parameters
	, mlow(-1.), mhigh(-1.)
	, componentIndex(0)
	, _useSig1Sig2(false)
	, denom_sigmam1_root2(0.), denom_sigmam2_root2(0.)
	, numer_factor1(0.), numer_factor2(0.)
	, exp1_denom(0.), exp2_denom(0.)
{
	cout << "Constructing BsMass" << endl;
	if( configurator->isTrue("UseSig1Sig2") )   { _useSig1Sig2 = true ; }

	MakePrototypes();

	plotComponents = configurator->isTrue( "PlotComponents" );
}

//Make the data point and parameter set
void BsMass::MakePrototypes()
{
	// Observables
	allObservables.push_back( recoMassName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( f_sig_m1Name );
	parameterNames.push_back( sigma_m1Name );

	if( _useSig1Sig2 ) parameterNames.push_back( sigma_m2Name );
	else		   parameterNames.push_back( ratio_21Name );

	parameterNames.push_back( m_BsName );
	allParameters = ParameterSet(parameterNames);
}

//Destructor
BsMass::~BsMass()
{
}

vector<string> BsMass::PDFComponents()
{
	vector<string> components;

	if( plotComponents )
	{
		components.push_back( "Gaussian1" );
		components.push_back( "Gaussian2" );
	}

	return components;
}

double BsMass::EvaluateComponent( DataPoint* measurement, ComponentRef* Component )
{
	componentIndex = Component->getComponentNumber();
	if( componentIndex == -1 )
	{
		string ComponentName = Component->getComponentName();
		if( ComponentName.compare( "Gaussian1" ) == 0 )
		{
			Component->setComponentNumber( 1 );
			componentIndex = 1;
		}
		else if( ComponentName.compare( "Gaussian2" ) == 0 )
		{
			Component->setComponentNumber( 2 );
			componentIndex = 2;
		}
		else
		{
			Component->setComponentNumber( 0 );
			componentIndex = 0;
		}
	}

	double return_value = this->Evaluate( measurement );
	componentIndex = 0;

	return return_value;
}

bool BsMass::SetPhysicsParameters( ParameterSet* Input )
{
	allParameters.SetPhysicsParameters( Input );

	double sigma_m1 = allParameters.GetPhysicsParameter( sigma_m1Name )->GetValue();

	double sigma_m2 = 0.0 ;
	if( _useSig1Sig2 ) {
		sigma_m2 = allParameters.GetPhysicsParameter( sigma_m2Name )->GetValue();
	}
	else {
		double ratio_21 = allParameters.GetPhysicsParameter( ratio_21Name )->GetValue();
		sigma_m2 = sigma_m1 * ratio_21 ;
	}

	denom_sigmam1_root2 = 1./(sigma_m1*Mathematics::SQRT_2());
	//PELC -> bug denom_sigmam2_root2 = 1./(sigma_m1*Mathematics::SQRT_2());
	denom_sigmam2_root2 = 1./(sigma_m2*Mathematics::SQRT_2());

	numer_factor1 = 1./(sigma_m1*sqrt(2.*Mathematics::Pi()));
	numer_factor2 = 1./(sigma_m2*sqrt(2.*Mathematics::Pi()));

	exp1_denom = 1. / ( 2. * sigma_m1 * sigma_m1 );
	exp2_denom = 1. / ( 2. * sigma_m2 * sigma_m2 );

	return true;
}

//Calculate the function value
double BsMass::Evaluate(DataPoint * measurement)
{
	// Get the physics parameters
	double f_sig_m1  = allParameters.GetPhysicsParameter( f_sig_m1Name )->GetValue();

	double m_Bs = allParameters.GetPhysicsParameter( m_BsName )->GetValue();

	// Get the observable
	double mass = measurement->GetObservable( recoMassName )->GetValue();

	// Temp way to initialise these - this means the first event of the first iteration is wrong
	// This also means it can never work for Toys
	if( (mlow < 0.) || (mhigh <0.))
	{
		mlow = measurement->GetPhaseSpaceBoundary()->GetConstraint( recoMassName )->GetMinimum();
		mhigh = measurement->GetPhaseSpaceBoundary()->GetConstraint( recoMassName )->GetMaximum();
	}

	double mhi_diff = mhigh-m_Bs;
	double mlo_diff = mlow-m_Bs;

	double s1_erf_factor = 0.5*( erf( mhi_diff * denom_sigmam1_root2 ) - erf( mlo_diff * denom_sigmam1_root2 ) );
	double s2_erf_factor = 0.5*( erf( mhi_diff * denom_sigmam2_root2 ) - erf( mlo_diff * denom_sigmam2_root2 ) );
	double returnValue = 0;

	if( f_sig_m1 >= 0.99999 ) // || ratio_21 < 1E-5 )
	{
		double factor1 = numer_factor1 / s1_erf_factor;
		double deltaM = mass - m_Bs;
		double deltaMsq = deltaM*deltaM;
		double exp1 = exp( -deltaMsq * exp1_denom );
		switch( componentIndex )
		{
			case 2:
				returnValue = 0;
				break;
			default:
				returnValue = factor1 * exp1;
				break;
		}
	}
	else
	{
		double factor1 = numer_factor1  / s1_erf_factor;
		double factor2 = numer_factor2  / s2_erf_factor;
		double deltaM = mass - m_Bs;
		double deltaMsq = deltaM*deltaM;
		double exp1 = exp( -deltaMsq * exp1_denom );
		double exp2 = exp( -deltaMsq * exp2_denom );
		switch( componentIndex )
		{
			case 1:
				returnValue = f_sig_m1 * factor1 * exp1;
				break;
			case 2:
				returnValue = (1. - f_sig_m1) * factor2 * exp2;
				break;
			default:
				returnValue = f_sig_m1 * factor1 * exp1 + (1. - f_sig_m1) * factor2 * exp2;
				break;
		}
	}

	if( std::isnan(returnValue) )
	{
		double sigma_m1 = allParameters.GetPhysicsParameter( sigma_m1Name )->GetValue();

		double sigma_m2 = 0.0;
		if( _useSig1Sig2 )
		{
			sigma_m2 = allParameters.GetPhysicsParameter( sigma_m2Name )->GetValue();
		}
		else
		{
			double ratio_21 = allParameters.GetPhysicsParameter( ratio_21Name )->GetValue();
			sigma_m2 = sigma_m1 * ratio_21 ;
		}

		PDF_THREAD_LOCK
		measurement->Print();
		measurement->GetPhaseSpaceBoundary()->Print();
		allParameters.Print();
		cout << "mlow: " << mlow << "  mhigh: " << mhigh << endl;
		cout << "sigma_m1: " << sigma_m1 << endl;
		cout << "PhaseSpace: " << measurement->GetPhaseSpaceBoundary() << endl;
		cout << "s1_erf_factor: " << erf((mhigh-m_Bs)/(sigma_m1*sqrt(2.))) << "-" << erf((mlow-m_Bs)/(sigma_m1*sqrt(2.))) << "   s2_erf_factor: " << s2_erf_factor << endl;
		cout << "factor1: " <<  "1./(" << sigma_m1 << " * " << sqrt(2.*Mathematics::Pi()) << ")" << endl;
		cout << "factor2: " << 1./(sigma_m2*sqrt(2.*Mathematics::Pi())) << endl;
		PDF_THREAD_UNLOCK
	}

	if( returnValue <=1E-15) returnValue=1E-15;

	return returnValue;
}


// Normalisation
double BsMass::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void) boundary;
	// Assumes that the mass integration limits are +/- Infinity
	// So take sufficiently large mass window.

	/*
	//These need to be put into the member variables in order to have then available in the evaluate method
	IConstraint * massBound = boundary->GetConstraint( recoMassName );
	if( massBound->GetUnit() == "NameNotFoundError" )
	{
	cerr << "Bound on mass not provided" << endl;
	return -1.;
	}
	else
	{
	mlow = massBound->GetMinimum();
	mhigh = massBound->GetMaximum();
	}
	 */
	double returnValue = 1.;//erf(mhigh) - erf(mlow);

	return returnValue;
}

