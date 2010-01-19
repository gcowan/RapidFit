/**
  @class RapidFitIntegrator

  A numerical integrator to be used for projecting the PDF, and as a fallback if the PDF does not provide its own normalisation
  This class uses two integrator classes provided by root: AdaptiveIntegratorMultiDim and GaussLegendreIntegrator.
  Both of these assume the function to be integrated is well behaved.

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-8
  */

#include "RapidFitIntegrator.h"
#include "StringProcessing.h"
#include <iostream>
#include <cmath>

using namespace std;

//1% tolerance on integral value
const double INTEGRAL_PRECISION_THRESHOLD = 0.01;

//Default constructor
RapidFitIntegrator::RapidFitIntegrator()
{
}

//Constructor with correct argument
RapidFitIntegrator::RapidFitIntegrator( IPDF * InputFunction, bool ForceNumerical ) : functionCanIntegrate(false), functionCanProject(false), haveTestedIntegral(ForceNumerical), functionToWrap(InputFunction), testFast(false)
{
	multiDimensionIntegrator = new AdaptiveIntegratorMultiDim();
	ROOT::Math::IntegrationOneDim::Type type = ROOT::Math::IntegrationOneDim::kGAUSS;
	oneDimensionIntegrator = new IntegratorOneDim(type);
}


//Constructor to test FoamIntegrator
RapidFitIntegrator::RapidFitIntegrator( IPDF * InputFunction, IDataSet * InputData, ParameterSet * InputParameters, bool ForceNumerical ) : functionCanIntegrate(false), functionCanProject(false), haveTestedIntegral(ForceNumerical), functionToWrap(InputFunction), testFast(true)
{
        multiDimensionIntegrator = new AdaptiveIntegratorMultiDim();
        ROOT::Math::IntegrationOneDim::Type type = ROOT::Math::IntegrationOneDim::kGAUSS;
        oneDimensionIntegrator = new IntegratorOneDim(type);

	cumulativeError = 0.0;
	numberCalls = 0.0;
	fastIntegrator = new FoamIntegrator( InputFunction, InputData );
	//fastIntegrator = new BenIntegrator( InputFunction, InputData->GetBoundary(), InputParameters );
}

//Destructor
RapidFitIntegrator::~RapidFitIntegrator()
{
	//cout << "Cumulative error: " << cumulativeError << endl;
	//cout << "Average error: " << cumulativeError / numberCalls  << endl;

	delete multiDimensionIntegrator;
	delete oneDimensionIntegrator;
}

//Return the integral over all observables
double RapidFitIntegrator::Integral( DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary )
{
	//Make a list of observables not to integrate
	vector<string> dontIntegrate = functionToWrap->GetDoNotIntegrateList();

	//Test the integration method the function has provided
	if (haveTestedIntegral)
	{
		//If the function has been tested already, use the result
		if (functionCanIntegrate)
		{
			if ( testFast )
			{
				double foamIntegral = fastIntegrator->Integral(NewDataPoint, NewBoundary );
				double analyticalIntegral = functionToWrap->Integral( NewDataPoint, NewBoundary );
				cout << "Foam integral = " << foamIntegral << " vs analytical = " << analyticalIntegral << endl;
				//numberCalls += 1.0;
				//cumulativeError += foamIntegral - analyticalIntegral;
				return analyticalIntegral;
				//return fastIntegrator->Integral( NewDataPoint, NewBoundary );
			}
			else
			{
				return functionToWrap->Integral( NewDataPoint, NewBoundary );
			}
		}
		else
		{
			return DoNumericalIntegral( NewDataPoint, NewBoundary, dontIntegrate );
		}
	}
	else
	{
		haveTestedIntegral = true;
		double testIntegral = functionToWrap->Integral( NewDataPoint, NewBoundary );
		double numericalIntegral = DoNumericalIntegral( NewDataPoint, NewBoundary, dontIntegrate );

		//Check if the function has an integrate method
		if ( testIntegral > 0.0 )
		{
			//Check if the function's integrate method agrees with the numerical integral
			if ( abs( numericalIntegral - testIntegral) / testIntegral < INTEGRAL_PRECISION_THRESHOLD )
			{
				//Trust the function's integration
				cout << "Function provides acceptable integration method: numerical " << numericalIntegral << " vs analytical " << testIntegral << endl;
				functionCanIntegrate = true;
				return testIntegral;
			}
			else
			{
				//Use numerical integration
				cerr << "Function provides poor integration method: numerical " << numericalIntegral << " vs analytical " << testIntegral << endl;
				//functionCanIntegrate = false;
				//return numericalIntegral;
				functionCanIntegrate = true;
				return testIntegral;
			}
		}
		else
		{
			//Use numerical integration
			cout << "Function provides no integration method: numerical " << numericalIntegral << " vs analytical " << testIntegral << endl;
			functionCanIntegrate = false;
			return numericalIntegral;
		}
	}
}

//Actually perform the numerical integration
double RapidFitIntegrator::DoNumericalIntegral( DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary, vector<string> DontIntegrateThese )
{
	//Make lists of observables to integrate and not to integrate
	vector<string> observableNames = functionToWrap->GetPrototypeDataPoint();
	vector<string> doIntegrate, dontIntegrate;
	for ( int observableIndex = 0; observableIndex < observableNames.size(); observableIndex++ )
	{
		bool continuous = !( NewBoundary->GetConstraint( observableNames[observableIndex] )->IsDiscrete() );
		bool integrate = ( StringProcessing::VectorContains( &DontIntegrateThese, &(observableNames[observableIndex]) ) == -1 );

		if ( continuous && integrate )
		{
			doIntegrate.push_back( observableNames[observableIndex] );
		}
		else
		{
			dontIntegrate.push_back( observableNames[observableIndex] );
		}
	}

	//If there are no observables left to integrate over, just evaluate the function
	if ( doIntegrate.size() == 0 )
	{
		return functionToWrap->Evaluate(NewDataPoint);
	}
	else
	{
		//Make the function wrapper
		IntegratorFunction quickFunction( functionToWrap, NewDataPoint, doIntegrate, dontIntegrate );

		//Chose the one dimensional or multi-dimensional method
		if ( doIntegrate.size() == 1 )
		{
			//Find the observable range to integrate over
			IConstraint * newConstraint = NewBoundary->GetConstraint( doIntegrate[0] );
			double minimum = newConstraint->GetMinimum();
			double maximum = newConstraint->GetMaximum();

			//Do a 1D integration
			oneDimensionIntegrator->SetFunction(quickFunction);
			return oneDimensionIntegrator->Integral( minimum, maximum );
		}
		else
		{
			//Make arrays of the observable ranges to integrate over
			double minima[ doIntegrate.size() ];
			double maxima[ doIntegrate.size() ];
			for ( int observableIndex = 0; observableIndex < doIntegrate.size(); observableIndex++ )
			{
				IConstraint * newConstraint = NewBoundary->GetConstraint( doIntegrate[observableIndex] );
				minima[observableIndex] = newConstraint->GetMinimum();
				maxima[observableIndex] = newConstraint->GetMaximum();
			}

			//Do a 2-15D integration
			multiDimensionIntegrator->SetFunction(quickFunction);
			return multiDimensionIntegrator->Integral( minima, maxima );
		}
	}
}

//Return the integral over all observables except one
double RapidFitIntegrator::ProjectObservable( DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary, string ProjectThis)
{
	//Make the list of observables not to integrate
	vector<string> dontIntegrate = functionToWrap->GetDoNotIntegrateList();
	dontIntegrate.push_back(ProjectThis);

	return DoNumericalIntegral( NewDataPoint, NewBoundary, dontIntegrate );
}
