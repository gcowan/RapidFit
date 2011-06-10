/**
        @class FitFunction

        Parent class for the function to minimise
	Overload the evaluate methods and UP value for Chi2, NLL, etc.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "FitFunction.h"

//Default constructor
FitFunction::FitFunction() : allData(), allIntegrators(), testDouble(), useWeights(false), weightObservableName()
{
}

//Destructor
FitFunction::~FitFunction()
{
	//cout << "Hello from FitFunction destructor" << endl;
}

//Set the physics bottle to fit with
void FitFunction::SetPhysicsBottle( PhysicsBottle * NewBottle )
{
	NewBottle->Finalise();
	allData = NewBottle;

	//Initialise the integrators
	for ( int resultIndex = 0; resultIndex < NewBottle->NumberResults(); ++resultIndex )
	{
		RapidFitIntegrator * resultIntegrator = new RapidFitIntegrator( NewBottle->GetResultPDF(resultIndex) );
		allIntegrators.push_back( resultIntegrator );
	}
}

//Return the physics bottle
PhysicsBottle* FitFunction::GetPhysicsBottle()
{
	return allData;
}

// Get and set the fit parameters
bool FitFunction::SetParameterSet( ParameterSet * NewParameters )
{
	return allData->SetParameterSet(NewParameters);
}
ParameterSet * FitFunction::GetParameterSet()
{
	return allData->GetParameterSet();
}

//Return the value to minimise
double FitFunction::Evaluate()
{
	double minimiseValue = 0.0;

	//Calculate the function value for each PDF-DataSet pair
	for (int resultIndex = 0; resultIndex < allData->NumberResults(); ++resultIndex)
	{
		minimiseValue += EvaluateDataSet( allData->GetResultPDF( resultIndex ), allData->GetResultDataSet( resultIndex ), allIntegrators[unsigned(resultIndex)] );
	}

	//Calculate the value of each constraint
	vector< ConstraintFunction* > constraints = allData->GetConstraints();
	for (unsigned int constraintIndex = 0; constraintIndex < constraints.size(); ++constraintIndex )
	{
		minimiseValue += constraints[constraintIndex]->Evaluate( allData->GetParameterSet() );
	}

	return minimiseValue;
}

//Return the value to minimise for a given PDF/DataSet pair
double FitFunction::EvaluateDataSet( IPDF * TestPDF, IDataSet * TestDataSet, RapidFitIntegrator * ResultIntegrator )
{
	//	Stupid gcc
	(void)TestPDF;
	(void)TestDataSet;
	(void)ResultIntegrator;

	return 1.0;
}

//Return the Up value for error calculation
double FitFunction::UpErrorValue( int Sigma )
{
	//	Stupid gcc
	(void)Sigma;
	return 1.0;
}

//Finalise the PhysicsBottle
void FitFunction::Finalise()
{
	allData->Finalise();
}

//Set the FitFunction to use per-event weights
void FitFunction::UseEventWeights( string WeightName )
{
	useWeights = true;
	weightObservableName = WeightName;
}
