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
#include "StatisticsFunctions.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#define DOUBLE_TOLERANCE 1E-6

using namespace std;

//1% tolerance on integral value
const double INTEGRAL_PRECISION_THRESHOLD = 0.01;

//Default constructor
RapidFitIntegrator::RapidFitIntegrator()
{
}

//Constructor with correct argument
RapidFitIntegrator::RapidFitIntegrator( IPDF * InputFunction, bool ForceNumerical ) : functionToWrap(InputFunction), functionCanIntegrate(false), functionCanProject(false), haveTestedIntegral(false), forceNumerical(ForceNumerical), cacheSetUp(false)
{
	multiDimensionIntegrator = new AdaptiveIntegratorMultiDim();
	ROOT::Math::IntegrationOneDim::Type type = ROOT::Math::IntegrationOneDim::kGAUSS;
	oneDimensionIntegrator = new IntegratorOneDim(type);
}

//Destructor
RapidFitIntegrator::~RapidFitIntegrator()
{
	delete multiDimensionIntegrator;
	delete oneDimensionIntegrator;
}

//Return the integral over all observables
double RapidFitIntegrator::Integral( DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary, bool UseCache )
{
	//Make a list of observables not to integrate
	vector<string> dontIntegrate = functionToWrap->GetDoNotIntegrateList();

	//Test the integration method the function has provided
	if ( haveTestedIntegral || forceNumerical )
	{
		//If the function has been tested already, use the result
		if ( functionCanIntegrate && !forceNumerical )
		{
			return functionToWrap->Integral( NewDataPoint, NewBoundary );
		}
		else
		{
			//Check whether to use precalculated numerical integral, or to recalculate
			if (UseCache)
			{
				return GetCachedIntegral(NewDataPoint);
			}
			else
			{
				return DoNumericalIntegral( NewDataPoint, NewBoundary, dontIntegrate );
			}
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
			if ( true /*abs( numericalIntegral - testIntegral) / testIntegral < INTEGRAL_PRECISION_THRESHOLD*/ )
			{
				//Trust the function's integration
				cout << std::setprecision(3) ;
				cout << "Integration  Test: numerical : analytical   " << numericalIntegral << " :  " << testIntegral <<  "\t\t\tUsing ANALYTICAL in all cases" << endl;
				ratioOfIntegrals = 1.;
				functionCanIntegrate = true;
				return testIntegral;
			}
			else
			{
				//Use numerical integration
				cerr << "Function provides poor integration method: numerical " << numericalIntegral << " vs analytical " << testIntegral << endl;
				ratioOfIntegrals = testIntegral/numericalIntegral;
				//functionCanIntegrate = false;
				//return numericalIntegral;

				if (UseCache)
				{
					//Calculate the discrete combinations
					SetUpIntegralCache(NewBoundary);

					//Check validity of cache
					double cachedIntegral = GetCachedIntegral(NewDataPoint);
					if ( abs( cachedIntegral - numericalIntegral ) / numericalIntegral < INTEGRAL_PRECISION_THRESHOLD )
					{
						cout << "Integral caching assumptions seem to be valid: numerical " << numericalIntegral << " vs cache " << cachedIntegral << endl;
					}
					else
					{
						cerr << "Integral caching error: numerical " << numericalIntegral << " vs cache " << cachedIntegral << endl;
						throw( 13 );
//						exit(1);
					}
				}

				functionCanIntegrate = true;
				return testIntegral;
			}
		}
		else
		{
			//Use numerical integration
			cout << "Function provides no integration method: numerical " << numericalIntegral << " vs analytical " << testIntegral << endl;
			functionCanIntegrate = false;
			ratioOfIntegrals = 1.;

			if (UseCache)
			{
				//Calculate the discrete combinations
				SetUpIntegralCache(NewBoundary);

				//Check validity of cache
				double cachedIntegral = GetCachedIntegral(NewDataPoint);
				if ( abs( cachedIntegral - numericalIntegral ) / numericalIntegral < INTEGRAL_PRECISION_THRESHOLD )
				{
					cout << "Integral caching assumptions seem to be valid: numerical " << numericalIntegral << " vs cache " << cachedIntegral << endl;
				}
				else
				{
					cerr << "Integral caching error: numerical " << numericalIntegral << " vs cache " << cachedIntegral << endl;
					throw( 13 );
//					exit(1);
				}
			}

			return numericalIntegral;
		}
	}
}

double RapidFitIntegrator::GetRatioOfIntegrals()
{
	return ratioOfIntegrals;
}

IPDF * RapidFitIntegrator::GetPDF()
{
	return functionToWrap;
}

//Actually perform the numerical integration
double RapidFitIntegrator::DoNumericalIntegral( DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary, vector<string> DontIntegrateThese )
{
	//Make lists of observables to integrate and not to integrate
	vector<string> doIntegrate, dontIntegrate;
	StatisticsFunctions::DoDontIntegrateLists( functionToWrap, NewBoundary, &DontIntegrateThese, doIntegrate, dontIntegrate );

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
			double* minima = new double[ doIntegrate.size() ];
			double* maxima = new double[ doIntegrate.size() ];
			for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
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

//Cache a numerical integral value for each discrete observable combination
void RapidFitIntegrator::UpdateIntegralCache( PhaseSpaceBoundary * NewBoundary )
{
	//First update any internal RapdidFitIntegrator caches that the PDF is using
	functionToWrap->UpdateIntegralCache();

	//Make a list of observables not to integrate
	vector<string> dontIntegrate = functionToWrap->GetDoNotIntegrateList();

	if ( cachedIntegrals.size() > 0 ) cachedIntegrals.clear();

	//Don't do it unless discrete combinations cached
	if (cacheSetUp)
	{
		for (unsigned int combinationIndex = 0; combinationIndex < discreteCombinations.size(); ++combinationIndex )
		{
			//Make a sample data point
			DataPoint samplePoint( NewBoundary->GetAllNames() );

			//Load the calculated discrete combinations
			for (unsigned int discreteIndex = 0; discreteIndex < discreteNames.size(); ++discreteIndex )
			{
				string unit = NewBoundary->GetConstraint( discreteNames[discreteIndex] )->GetUnit();
				samplePoint.SetObservable( discreteNames[discreteIndex], discreteCombinations[combinationIndex][discreteIndex], 1, unit );
			}

			//Just generate some random continuous values
			for (unsigned int continuousIndex = 0; continuousIndex < continuousNames.size(); ++continuousIndex )
			{
				samplePoint.SetObservable( continuousNames[continuousIndex], NewBoundary->GetConstraint( continuousNames[continuousIndex] )->CreateObservable() );
			}

			//Integrate and store the result
			double integralToCache = DoNumericalIntegral( &samplePoint, NewBoundary, dontIntegrate );
			cachedIntegrals.push_back(integralToCache);

			/*
			//Debug
			cout << "For combination: ";
			for ( int discreteIndex = discreteNames.size() - 1; discreteIndex >= 0; --discreteIndex )
			{
				cout << discreteCombinations[combinationIndex][discreteIndex] << ", ";
			}
			cout << "get integral: " << integralToCache << endl;
			*/
		}
	}
}

//Return the correct cached integral value
double RapidFitIntegrator::GetCachedIntegral( DataPoint * NewDataPoint )
{
	//Don't do it unless cache done
	if (cacheSetUp)
	{
		//Generate the discrete observables, and select the correct Foam generator to use
		int combinationIndex = 0;
		int incrementValue = 1;
		//cout << "Discrete values are: ";
		for ( int discreteIndex = int(discreteNames.size()) - 1; discreteIndex >= 0; --discreteIndex )
		{
			//Retrieve the discrete value
			double currentValue = NewDataPoint->GetObservable( discreteNames[discreteIndex] )->GetValue();
			//cout << currentValue << ", ";

			//Calculate the index
			for (unsigned int valueIndex = 0; valueIndex < discreteValues[discreteIndex].size(); ++valueIndex )
			{
				if ( fabs( discreteValues[discreteIndex][valueIndex] - currentValue ) < DOUBLE_TOLERANCE )
				{
					combinationIndex += ( incrementValue * valueIndex );
					incrementValue *= int(discreteValues[discreteIndex].size());
					break;
				}
			}
		}

		//Debug
		//cout << "versus: ";
		//for ( int discreteIndex = discreteNames.size() - 1; discreteIndex >= 0; --discreteIndex )
		//{
		//	cout << discreteCombinations[combinationIndex][discreteIndex] << ", ";
		//}
		//cout << endl;

		//Return the cached integral value for this index
		return cachedIntegrals[combinationIndex];
	}
	else
	{
		return -1.0;
	}
}

//Calculate the discrete combinations
void RapidFitIntegrator::SetUpIntegralCache( PhaseSpaceBoundary * NewBoundary )
{
	//Retrieve all combinations of discrete variables
	vector<string> allNames = NewBoundary->GetAllNames();
	discreteCombinations = StatisticsFunctions::DiscreteCombinations( &allNames, NewBoundary, discreteNames, continuousNames, discreteValues );

	//Now calculate the integral values to cache
	cacheSetUp = true;
	UpdateIntegralCache(NewBoundary);
}
