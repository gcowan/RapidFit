// $Id: Exponential.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Exponential Exponential.cpp
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-11-13
 */

#include "Exponential.h"
#include "Mathematics.h"
#include "TimeAccRes.h"
#include <iostream>
#include <cmath>

using namespace::std;

PDF_CREATOR( Exponential );

//Constructor
Exponential::Exponential( PDFConfigurator* configurator) :
	// Physics parameters
	tauName	( configurator->getName("tau") )
	// Observables
	, timeName      ( configurator->getName("time") )
	, timeConst	( configurator->getName("time") )
	//objects used in XML
	, tau(), gamma()
{
	resolutionModel = new TimeAccRes( configurator );

	this->MakePrototypes();
}

Exponential::Exponential( const Exponential &copy ) :
	BasePDF( (BasePDF)copy )
	, tauName ( copy.tauName )
	, timeName ( copy.timeName )
	, timeConst ( copy.timeConst )
	, resolutionModel( NULL )
	, tau ( copy.tau )
	, gamma ( copy.gamma )
{
	resolutionModel = new TimeAccRes( *( (TimeAccRes*)copy.resolutionModel) );
}

//Make the data point and parameter set
void Exponential::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	resolutionModel->addObservables( allObservables );
	resolutionModel->addObservables( doNotIntegrateList );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( tauName );

	resolutionModel->addParameters( parameterNames );

	allParameters = ParameterSet(parameterNames);
}

//Destructor
Exponential::~Exponential()
{
	if( resolutionModel != NULL ) delete resolutionModel;
}

bool Exponential::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);

	tau = allParameters.GetPhysicsParameter( tauName )->GetValue();
	gamma = 1./tau;

	return isOK;
}

/*
double Exponential::EvaluateComponent( DataPoint* input, ComponentRef* thisRef )
{
	(void) thisRef;		//	This PDF has no concept of a component, it is JUST an exponential
	Observable* timeObs = input->GetObservable( timeName );
	if( _useSteppedProjection )
	{
		if( timeAcc->GetIsSorted() )
		{
			unsigned int binNum = timeAcc->findSliceNum( timeObs, timeOffset );
			double t_low = timeAcc->getSlice( binNum )->tlow();
			double t_high = timeAcc->getSlice( binNum )->thigh();
			double frac=0.5;
			//frac = ( exp(-gamma*t_low) - exp(-gamma*(t_low+0.5*(t_high-t_low))) ) / ( exp(-gamma*t_low) - exp(-gamma*t_high) );
			frac = ( exp(-gamma*(t_low+0.5*(t_high-t_low))) - exp(-gamma*t_low) ) / ( exp(-gamma*t_high) - exp(-gamma*t_low) );
			double mean = t_low+(1.-frac)*( t_high-t_low );
			DataPoint* thisPoint = new DataPoint( *input );
			Observable* newTimeObs = new Observable( timeName, mean, " " );
			thisPoint->SetObservable( timeName, newTimeObs );
			double returnable = this->Evaluate( thisPoint );
			delete newTimeObs;
			delete thisPoint;
			return returnable;
		}
		else
		{
			cerr << "Acceptance Histogram has to be in an ordered format (histogram-like) in order to do this!" << endl;
			cerr << "Disabling the Stepped Projections" << endl;
			_useSteppedProjection = false;
			return this->Evaluate( input );
		}
	}
	else
	{
		return this->Evaluate( input );
	}
	return 0.;
}
*/

//Calculate the function value
double Exponential::Evaluate(DataPoint * measurement)
{
	// Observable
	Observable* timeObs = measurement->GetObservable( timeName );
	time = timeObs->GetValue();

	return resolutionModel->Exp( time, gamma );
}

double Exponential::Normalisation( DataPoint * measurement, PhaseSpaceBoundary * boundary )
{
	(void) measurement;
	IConstraint* timeC = boundary->GetConstraint( timeConst );
	double tlow = timeC->GetMinimum();
	double thigh = timeC->GetMaximum();

	return resolutionModel->ExpInt( tlow, thigh, gamma );
}

