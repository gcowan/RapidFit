/*! @class Sinusoid Sinusoid.cpp
 *
 *  RapidFit Sinusoid for Constructing new PDFs
 *
 *
 *  @author Rob
 *  @date 20xx-yy-zz
 */

#include "Sinusoid.h"
#include <iostream>

using namespace::std;

//	Uncomment this when you have a working PDF
PDF_CREATOR( Sinusoid );

//Constructor
Sinusoid::Sinusoid(PDFConfigurator* configurator) :
	// Physics parameters
	AmplitudeName		( configurator->getName("Amplitude") ),
	PhaseName		( configurator->getName("Phase") ),
	FreqName		( configurator->getName("Frequency") ),
	// Observables
	TimeName		( configurator->getName("time") ),
	plotComponents(false), useCos(false), componentIndex(0),
	A(0.), f(0.), phase(0.)
{
	this->MakePrototypes();

	plotComponents = configurator->isTrue( "PlotComponents" );
	useCos = configurator->isTrue( "UseCos" );
}

//Destructor
Sinusoid::~Sinusoid()
{ 
}

//Make the data point and parameter set
void Sinusoid::MakePrototypes()
{
	this->MakePrototypeDataPoint();
	this->MakePrototypeParameterSet();
}

void Sinusoid::MakePrototypeDataPoint()
{
	// Observables
	allObservables.push_back( TimeName );
}

void Sinusoid::MakePrototypeParameterSet()
{
	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( AmplitudeName );
	parameterNames.push_back( PhaseName );
	parameterNames.push_back( FreqName );
	allParameters = ParameterSet(parameterNames);
}

vector<string> Sinusoid::PDFComponents()
{
	vector<string> components;

	if( plotComponents )
	{
		components.push_back( "Sinusoid" );
	}

	return components;
}

bool Sinusoid::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	allParameters.SetPhysicsParameters( NewParameterSet );

	A = allParameters.GetPhysicsParameter( AmplitudeName )->GetValue();
	phase = allParameters.GetPhysicsParameter( PhaseName )->GetValue();
	f = allParameters.GetPhysicsParameter( FreqName )->GetValue();

	return true;	//	I would like to get rid of this return!
}

double Sinusoid::EvaluateComponent( DataPoint* measurement, ComponentRef* Component )
{
	componentIndex = Component->getComponentNumber();

	if( componentIndex == -1 )
	{
		string ComponentName = Component->getComponentName();
		if( ComponentName.compare( "Sinusoid" ) == 0 )
		{
			Component->setComponentNumber( 1 );
			componentIndex = 1;
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

//Calculate the function value
double Sinusoid::Evaluate(DataPoint * measurement)
{
	double t = measurement->GetObservable( TimeName )->GetValue();

	double returnValue=0.;

	switch( componentIndex )
	{
		case 2:
			break;
		default:
			if( !useCos )
			{
				returnValue+=A*(1.+sin(f*t+phase));
			}
			else
			{
				returnValue+=A*(1.+cos(f*t+phase));
			}
			break;
	}

	//cout << returnValue << endl;
	//if( measurement->GetPhaseSpaceBoundary() !=  NULL ) cout << this->Normalisation( measurement->GetPhaseSpaceBoundary() ) << endl; exit(0);
	return returnValue;
}

// Normalisation
double Sinusoid::Normalisation(PhaseSpaceBoundary * boundary)
{
	double hi = boundary->GetConstraint( TimeName )->GetMaximum();
	double lo = boundary->GetConstraint( TimeName )->GetMinimum();

	double box = fabs(hi-lo)*A;

	double sinus=0.;

	if( !useCos )
	{
		sinus = -(A/f) * ( cos( hi*f+phase ) - cos( lo*f+phase ) );
	}
	else
	{
		sinus = (A/f) * ( sin( hi*f+phase ) - sin( lo*f+phase) );
	}

	double returnValue=box+sinus;

	return returnValue;
}

