/**
        @class FitFunction

        Parent class for the function to minimise
	Overload the evaluate methods and UP value for Chi2, NLL, etc.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#include "FitFunction.h"

//Default constructor
FitFunction::FitFunction() : useWeights(false)
{
}

//Destructor
FitFunction::~FitFunction()
{
}

//Set the physics bottle to fit with
void FitFunction::SetPhysicsBottle( PhysicsBottle * NewBottle )
{
	NewBottle->Finalise();
	allData = NewBottle;

	//Identify parameters that directly affect the fit value
	vector<string> allNames = NewBottle->GetParameterSet()->GetAllNames();
	for ( int nameIterator = 0; nameIterator < allNames.size(); nameIterator++ )
	{
		string type = NewBottle->GetParameterSet()->GetPhysicsParameter( allNames[nameIterator] )->GetType();
		if ( type == "GaussianConstrained" )
		{
			interestingParameters.push_back( allNames[nameIterator] );
		}
	}

	//Initialise the integrators
	for ( int resultIndex = 0; resultIndex < NewBottle->NumberResults(); resultIndex++ )
	{
		RapidFitIntegrator * resultIntegrator = new RapidFitIntegrator( NewBottle->GetResultPDF(resultIndex) );//, NewBottle->GetResultDataSet(resultIndex), NewBottle->GetParameterSet() );
		allIntegrators.push_back( resultIntegrator );
	}
}

//Return the physics bottle
PhysicsBottle * FitFunction::GetPhysicsBottle()
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
	for (int resultIndex = 0; resultIndex < allData->NumberResults(); resultIndex++)
	{
		minimiseValue += EvaluateDataSet( allData->GetResultPDF( resultIndex ), allData->GetResultDataSet( resultIndex ), allIntegrators[resultIndex] );
	}
	minimiseValue += EvaluateParameterSet( allData->GetParameterSet(), interestingParameters );

	return minimiseValue;
}

//Return the value to minimise for a given PDF/DataSet result
double FitFunction::EvaluateDataSet( IPDF * TestPDF, IDataSet * TestDataSet, RapidFitIntegrator * ResultIntegrator )
{
	return 1.0;
}

//Return the value to minimise for a ParameterSet
double FitFunction::EvaluateParameterSet( ParameterSet * TestParameterSet, vector<string> InterestingParameters )
{
	return 1.0;
}

//Return the Up value for error calculation
double FitFunction::UpErrorValue( int Sigma )
{
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
