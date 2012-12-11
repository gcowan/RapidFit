// $Id: LogNormalDistribution.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class LogNormalDistribution LogNormalDistribution.cpp
 *
 *  RapidFit PDF for Bs mass
 *
 *
 *  @author Pete
 *  @date 2011-07-30
 */

#include "LogNormalDistribution.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( LogNormalDistribution );

//Constructor
LogNormalDistribution::LogNormalDistribution(PDFConfigurator* configurator) :
	// Physics parameters
	LNsigma1Name	( configurator->getName("LNsigma1" ) ),
	LNtheta1Name	( configurator->getName("LNtheta1" ) ),
	LNm1Name	( configurator->getName("LNm1" ) ),
	LNsigma2Name	( configurator->getName("LNsigma2" ) ),
	LNtheta2Name	( configurator->getName("LNtheta2" ) ),
	LNm2Name	( configurator->getName("LNm2" ) ),
	LNfName		( configurator->getName("LNf1" ) ),
	// Observables
	LNxName( configurator->getName("LNx") ),
	plotComponents(false)
{
	
	std::cout << "Constructing PDF: LogNormalDistribution " << std::endl ;

	plotComponents = configurator->isTrue( "PlotComponents" );

	MakePrototypes();
}

//Make the data point and parameter set
void LogNormalDistribution::MakePrototypes()
{
	// Observables
	allObservables.push_back( LNxName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( LNsigma1Name );
	parameterNames.push_back( LNtheta1Name );
	parameterNames.push_back( LNm1Name );
	parameterNames.push_back( LNsigma2Name );
	parameterNames.push_back( LNtheta2Name );
	parameterNames.push_back( LNm2Name );
	parameterNames.push_back( LNfName );
	allParameters = ParameterSet(parameterNames);
}

//Destructor
LogNormalDistribution::~LogNormalDistribution()
{
}


bool LogNormalDistribution::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  	sigma1 = allParameters.GetPhysicsParameter( LNsigma1Name )->GetValue();
  	theta1 = allParameters.GetPhysicsParameter( LNtheta1Name )->GetValue();
  	m1     = allParameters.GetPhysicsParameter( LNm1Name )->GetValue();
  	sigma2 = allParameters.GetPhysicsParameter( LNsigma2Name )->GetValue();
  	theta2 = allParameters.GetPhysicsParameter( LNtheta2Name )->GetValue();
  	m2     = allParameters.GetPhysicsParameter( LNm2Name )->GetValue();
  	f      = allParameters.GetPhysicsParameter( LNfName )->GetValue();
	
	return isOK;
}



//Calculate the function value
double LogNormalDistribution::Evaluate(DataPoint * measurement)
{
	// Get the observable
	double x = measurement->GetObservable( LNxName )->GetValue();
	
	double returnValue = 0.000000000001;

	if ( x > theta1 && x > theta2 ) returnValue = f*TMath::LogNormal( x, sigma1, theta1, m1 ) + (1.-f)*TMath::LogNormal( x, sigma2, theta2, m2 );

  	return returnValue;
}


// Normalisation
double LogNormalDistribution::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void)boundary;
	// Assumes that the mass integration limits are +/- Infinity
	// So take sufficiently large mass window.
	return 1.0;
}

vector<string> LogNormalDistribution::PDFComponents()
{
        vector<string> components;

        if( plotComponents )
        {
                components.push_back( "Log1" );
                components.push_back( "Log2" );
        }

        return components;
}

double LogNormalDistribution::EvaluateComponent( DataPoint* input, ComponentRef* Component )
{
        int componentIndex = Component->getComponentNumber();
        if( componentIndex == -1 )
        {
		double x = input->GetObservable( LNxName )->GetValue();
                string ComponentName = Component->getComponentName();
                if( ComponentName.compare( "Log1" ) == 0 )
                {
                        Component->setComponentNumber( 1 );
			return f*TMath::LogNormal( x, sigma1, theta1, m1 );
		}
		else if( ComponentName.compare( "Log2" ) == 0 )
		{
			Component->setComponentNumber( 2 );
			return (1.-f)*TMath::LogNormal( x, sigma2, theta2, m2 );
		}
		else
		{
			Component->setComponentNumber( 0 );
			return this->Evaluate( input );
		}
	}
	else if( componentIndex == 1 )
	{
		double x = input->GetObservable( LNxName )->GetValue();
		return f*TMath::LogNormal( x, sigma1, theta1, m1 );
	}
	else if( componentIndex == 2 )
	{
		double x = input->GetObservable( LNxName )->GetValue();
		return (1.-f)*TMath::LogNormal( x, sigma2, theta2, m2 );
	}
	else
	{
		return this->Evaluate( input );
	}
}

