#include "MakeFoam.h"
#include <queue>
#include "StatisticsFunctions.h"
#include <iostream>
#include <cmath>
#include "RapidFitIntegrator.h"

using namespace std;

const int MAXIMUM_SAMPLES = 1000;
const double MAXIMUM_GRADIENT_TOLERANCE = 0.01;
const int MAXIMUM_CELLS = 1000;

//Default constructor
MakeFoam::MakeFoam()
{
}

//Constructor with correct argument
MakeFoam::MakeFoam( IPDF * InputPDF, PhaseSpaceBoundary * InputBoundary, DataPoint * InputPoint ) : integratePDF(InputPDF)
{
	//Make the container to hold the possible cells
	queue<PhaseSpaceBoundary> possibleCells;
	PhaseSpaceBoundary firstCell = *InputBoundary;
	possibleCells.push(firstCell);

	//Make a list of observables to integrate over
	vector<string> doIntegrate, dontIntegrate;
	vector<string> pdfDontIntegrate = InputPDF->GetDoNotIntegrateList();
	StatisticsFunctions::DoDontIntegrateLists( InputPDF, InputBoundary, &(pdfDontIntegrate), doIntegrate, dontIntegrate );

	//Continue until all possible cells have been examined
	while ( !possibleCells.empty() )
	{
		//Remove the next possible cell
		PhaseSpaceBoundary currentCell = possibleCells.front();
		possibleCells.pop();

		//MC sample the cell, store the DataPoints producing the highest and lowest values
		double maximumValue, minimumValue;
		DataPoint maximumPoint, minimumPoint;
		for ( int sampleIndex = 0; sampleIndex < MAXIMUM_SAMPLES; sampleIndex++ )
		{
			//Create a data point within the current cell
			DataPoint samplePoint( InputPoint->GetAllNames() );
			for ( int observableIndex = 0; observableIndex < doIntegrate.size(); observableIndex++ )
			{
				//Generate random values to explore integrable observables
				IConstraint * temporaryConstraint = currentCell.GetConstraint( doIntegrate[observableIndex] );
				samplePoint.SetObservable( doIntegrate[observableIndex], temporaryConstraint->CreateObservable() );
			}
			for ( int observableIndex = 0; observableIndex < dontIntegrate.size(); observableIndex++ )
			{
				//Use given values for unintegrable observables
				Observable * temporaryObservable = InputPoint->GetObservable( dontIntegrate[observableIndex] );
				samplePoint.SetObservable( dontIntegrate[observableIndex], temporaryObservable );
			}

			//Evaluate the function at this point
			double sampleValue = InputPDF->Evaluate( &samplePoint );

			//Update maximum and minimum
			if ( sampleIndex == 0 )
			{
				//Populate maximum and minimum
				maximumValue = sampleValue;
				minimumValue = sampleValue;
				maximumPoint = samplePoint;
				minimumPoint = samplePoint;
			}
			else
			{
				//Find maximum
				if ( sampleValue > maximumValue )
				{
					maximumValue = sampleValue;
					maximumPoint = samplePoint;
				}

				//Find minimum
				if ( sampleValue < minimumValue )
				{
					minimumValue = sampleValue;
					minimumPoint = samplePoint;
				}
			}
		}

		//Find the maximum gradient
		vector<double> midPoints;
		string maximumGradientObservable, unit;
		double maximumGradient, lowPoint, midPoint, highPoint;
		for ( int observableIndex = 0; observableIndex < doIntegrate.size(); observableIndex++ )
		{
			//Find the size of the cell in this observable
			IConstraint * temporaryConstraint = currentCell.GetConstraint( doIntegrate[observableIndex] );
			double cellMaximum = temporaryConstraint->GetMaximum();
			double cellMinimum = temporaryConstraint->GetMinimum();

			//Find the distance between the maximum and minimum points in this observable
			double pointMaximum = maximumPoint.GetObservable( doIntegrate[observableIndex] )->GetValue();
			double pointMinimum = minimumPoint.GetObservable( doIntegrate[observableIndex] )->GetValue();
			double pointRange = pointMaximum - pointMinimum;

			//Find the gradient
			//NOTE: this is the crucial expression for determining foam generation. It might well be wrong!
			//FIRST TRY: chose the observable with smallest relative distance between max and min point
			//Failed because creates infinitely thin cells
			//SECOND TRY: chose the observable with greatest absolute distance between min and max
			double gradient = abs( ( maximumValue - minimumValue ) * pointRange );//( cellMaximum - cellMinimum ) / pointRange );

			//Store the mid point
			double observableMiddle = cellMinimum + ( ( cellMaximum - cellMinimum ) / 2.0 );
			midPoints.push_back(observableMiddle);

			//Update maximum
			if ( observableIndex == 0 || gradient > maximumGradient )
			{
				maximumGradient = gradient;
				maximumGradientObservable = doIntegrate[observableIndex];
				unit = temporaryConstraint->GetUnit();
				highPoint = cellMaximum;
				midPoint = observableMiddle;
				lowPoint = cellMinimum;
			}
		}

		//If the maximum gradient is within tolerance, the cell is finished
		if ( maximumGradient < MAXIMUM_GRADIENT_TOLERANCE )
		{
			//Store the finished cell
			finishedCells.push_back(currentCell);

			//Create a data point at the center of the current cell
			DataPoint cellCenter( InputPoint->GetAllNames() );
			for ( int observableIndex = 0; observableIndex < doIntegrate.size(); observableIndex++ )
			{
				//Use the mid points for the integrable values
				Observable * temporaryObservable = cellCenter.GetObservable( doIntegrate[observableIndex] );
				temporaryObservable->SetValue( midPoints[observableIndex] );
				cellCenter.SetObservable( doIntegrate[observableIndex], temporaryObservable );
			}
			for ( int observableIndex = 0; observableIndex < dontIntegrate.size(); observableIndex++ )
			{
				//Use given values for unintegrable observables
				Observable * temporaryObservable = InputPoint->GetObservable( dontIntegrate[observableIndex] );
				cellCenter.SetObservable( dontIntegrate[observableIndex], temporaryObservable );
			}

			//Store the center point
			centerPoints.push_back(cellCenter);
		}
		else
		{
			//Create two cells to replace the current cell
			PhaseSpaceBoundary daughterCell1( currentCell.GetAllNames() );
			PhaseSpaceBoundary daughterCell2( currentCell.GetAllNames() );

			for ( int observableIndex = 0; observableIndex < doIntegrate.size(); observableIndex++ )
			{
				IConstraint * temporaryConstraint = currentCell.GetConstraint( doIntegrate[observableIndex] );
				if ( doIntegrate[observableIndex] == maximumGradientObservable )
				{
					//Split the cells on the observable with the greatest gradient
					if ( lowPoint == midPoint )
					{
						cerr << "lowPoint == midPoint" << endl;
						return;
					}
					if ( midPoint == highPoint )
					{
						cerr << "midPoint == highPoint" << endl;
						return;
					}
					daughterCell1.SetConstraint( doIntegrate[observableIndex], lowPoint, midPoint, unit );
					daughterCell2.SetConstraint( doIntegrate[observableIndex], midPoint, highPoint, unit );
				}
				else
				{
					//Copy the continuous constraint (if it can be integrated, it must be continuous)
					daughterCell1.SetConstraint( doIntegrate[observableIndex], temporaryConstraint->GetMinimum(), temporaryConstraint->GetMaximum(), temporaryConstraint->GetUnit() );
					daughterCell2.SetConstraint( doIntegrate[observableIndex], temporaryConstraint->GetMinimum(), temporaryConstraint->GetMaximum(), temporaryConstraint->GetUnit() );
				}
			}
			for ( int observableIndex = 0; observableIndex < dontIntegrate.size(); observableIndex++ )
			{
				IConstraint * temporaryConstraint = currentCell.GetConstraint( dontIntegrate[observableIndex] );
				if ( temporaryConstraint->IsDiscrete() )
				{
					//Copy the discrete constraint
					daughterCell1.SetConstraint( dontIntegrate[observableIndex], temporaryConstraint->GetValues(), temporaryConstraint->GetUnit() );
					daughterCell2.SetConstraint( dontIntegrate[observableIndex], temporaryConstraint->GetValues(), temporaryConstraint->GetUnit() );
				}
				else
				{
					//Copy the continuous constraint
					daughterCell1.SetConstraint( dontIntegrate[observableIndex], temporaryConstraint->GetMinimum(), temporaryConstraint->GetMaximum(), temporaryConstraint->GetUnit() );
					daughterCell2.SetConstraint( dontIntegrate[observableIndex], temporaryConstraint->GetMinimum(), temporaryConstraint->GetMaximum(), temporaryConstraint->GetUnit() );
				}
			}

			//Add the two new cells to the possibles
			possibleCells.push(daughterCell1);
			possibleCells.push(daughterCell2);
		}

		//Make sure you don't exceed the maximum number of cells
		if ( finishedCells.size() + possibleCells.size() >= MAXIMUM_CELLS )
		{
			//Dump out any possible cells into finished, without any further examination
			while ( !possibleCells.empty() )
			{
				//Move the cell from "possible" to "finished"
				PhaseSpaceBoundary temporaryCell = possibleCells.front();
				possibleCells.pop();
				finishedCells.push_back(temporaryCell);

				//Create a data point at the center of the current cell
				DataPoint cellCenter( InputPoint->GetAllNames() );
				for ( int observableIndex = 0; observableIndex < doIntegrate.size(); observableIndex++ )
				{
					//Calculate the cell mid point
					IConstraint * temporaryConstraint = temporaryCell.GetConstraint( doIntegrate[observableIndex] );
					double midPoint = temporaryConstraint->GetMinimum() + ( ( temporaryConstraint->GetMaximum() - temporaryConstraint->GetMinimum() ) / 2  );

					//Use the mid points for the integrable values
					Observable * temporaryObservable = cellCenter.GetObservable( doIntegrate[observableIndex] );
					temporaryObservable->SetValue(midPoint);
					cellCenter.SetObservable( doIntegrate[observableIndex], temporaryObservable );
				}
				for ( int observableIndex = 0; observableIndex < dontIntegrate.size(); observableIndex++ )
				{
					//Use given values for unintegrable observables
					Observable * temporaryObservable = InputPoint->GetObservable( dontIntegrate[observableIndex] );
					cellCenter.SetObservable( dontIntegrate[observableIndex], temporaryObservable );
				}

				//Store the center point
				centerPoints.push_back(cellCenter);
			}

			//Exit the foam loop
			break;
		}
	}

	//Now all the cells have been made!
	//Make a numerical integrator for the function
	RapidFitIntegrator cellIntegrator( InputPDF, true );

	//Prepare the sanity check
	//double wholeIntegral = cellIntegrator.Integral( InputPoint, InputBoundary );
	//double totalCellIntegral = 0.0;

	//Find the function value at the center of each cell, and the integral of the function over the cell
	for ( int cellIndex = 0; cellIndex < finishedCells.size(); cellIndex++ )
	{
		//Integrate the cell
		double integral = cellIntegrator.Integral( InputPoint, &( finishedCells[cellIndex] ) );
		cellIntegrals.push_back(integral);
		//totalCellIntegral += integral;

		//Evaluate the function at the center of the cell
		double value = InputPDF->Evaluate( &( centerPoints[cellIndex] ) );
		centerValues.push_back(value);
	}

	//cout << "Whole integral = " << wholeIntegral << endl;
	//cout << "Total cells = " << totalCellIntegral << endl;
}

//Destructor
MakeFoam::~MakeFoam()
{
}

//Output all cells to console
void MakeFoam::Debug()
{
	cout << centerPoints.size() << " cells generated" << endl << endl;
	vector<string> allNames = centerPoints[0].GetAllNames();
	vector<string>::iterator observableIterator;
	for ( observableIterator = allNames.begin(); observableIterator != allNames.end(); observableIterator++ )
	{
		cout << *observableIterator << ", ";
	}
	cout << endl;
	for ( int pointIndex = 0; pointIndex < centerPoints.size(); pointIndex++ )
	{
		for ( observableIterator = allNames.begin(); observableIterator != allNames.end(); observableIterator++ )
		{
			cout << centerPoints[pointIndex].GetObservable( *observableIterator )->GetValue() << ", ";
		}
		cout << endl;
	}

	for ( int cellIndex = 0; cellIndex < finishedCells.size(); cellIndex++ )
	{
		cout << centerValues[cellIndex] << ", " << cellIntegrals[cellIndex] << endl;
	}
}

//Integrate the function
double MakeFoam::Integral()
{
	//Don't need any input: assume the phase space is unchanged (or else the whole initialisation is invalid)
	//Assume the correct foam has already been chosen for the data point

	//Evaluate each cell center point, weight precalculated integral by the change
	double totalIntegral = 0.0;
	for ( int cellIndex = 0; cellIndex < centerPoints.size(); cellIndex++ )
	{
		//Evaluate the function at the center of the cell
		double newCenterValue = integratePDF->Evaluate( &( centerPoints[cellIndex] ) );

		//Return the precalculated integral of the cell, scaled to the new function value
		double newIntegral = cellIntegrals[cellIndex] * newCenterValue / centerValues[cellIndex];
		totalIntegral += newIntegral;
	}

	return totalIntegral;
}
