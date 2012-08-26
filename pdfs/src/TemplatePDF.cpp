/*! @class TemplatePDF TemplatePDF.cpp
 *
 *  RapidFit TemplatePDF for Constructing new PDFs
 *
 *
 *  @author Rob
 *  @date 20xx-yy-zz
 */

#include "TemplatePDF.h"
#include <iostream>

using namespace::std;

//	Uncomment this when you have a working PDF
//PDF_CREATOR( TemplatePDF );

//Constructor
TemplatePDF::TemplatePDF(PDFConfigurator* configurator) :
	// Physics parameters
	ParameterName		( configurator->getName("Parameter1") )
	// Observables
	, ObservableName	( configurator->getName("Observable1") )
{
	this->MakePrototypes();

	this->plotComponents = configurator->isTrue( "PlotComponents" );
}

//Destructor
TemplatePDF::~TemplatePDF()
{ 
}

//Make the data point and parameter set
void TemplatePDF::MakePrototypes()
{
	this->MakePrototypeDataPoint();
	this->MakePrototypeParameterSet();
}

void TemplatePDF::MakePrototypeDataPoint()
{
	// Observables
	allObservables.push_back( ObservableName );
}

void TemplatePDF::MakePrototypeParameterSet()
{
	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( ParameterName );
	allParameters = ParameterSet(parameterNames);
}

vector<string> TemplatePDF::PDFComponents()
{
	vector<string> components;

	if( plotComponents )
	{
		components.push_back( "Component1" );
		components.push_back( "Component2" );
	}

	return components;
}

bool TemplatePDF::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	allParameters.SetPhysicsParameters( NewParameterSet );
}

double TemplatePDF::EvaluateComponent( DataPoint* measurement, ComponentRef* Component )
{
	componentIndex = Component->getComponentNumber();

	if( componentIndex == -1 )
	{
		string ComponentName = Component->getComponentName();
		if( ComponentName.compare( "Component1" ) == 0 )
		{
			Component->setComponentNumber( 1 );
			componentIndex = 1;
		}
		else if( ComponentName.compare( "Component2" ) == 0 )
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

//Calculate the function value
double TemplatePDF::Evaluate(DataPoint * measurement)
{
	double somePhysicsVal = allParameters.GetPhysicsParameter( ParameterName )->GetValue();
	(void) somePhysicsVal;
	double someDataPointValue = measurement->GetObservable( ObservableName )->GetValue();
	(void) someDataPointValue;

	double returnValue=1.;

	switch( componentIndex )
	{
		case 1:
			returnValue = 0.25*returnValue;
			break;
		case 2:
			returnValue = 0.75*returnValue;
			break;
		default:
			break;
	}

	return returnValue;
}

// Normalisation
double TemplatePDF::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void) boundary;
	double returnValue = 1.;

	//double returnValue = -1.;	//	For Numerical Integration

	return returnValue;
}

