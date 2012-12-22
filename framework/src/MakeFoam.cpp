//	RapidFit Headers
#include "MakeFoam.h"
#include "StatisticsFunctions.h"
#include "RapidFitIntegrator.h"
#include "ClassLookUp.h"
//	System Headers
#include <queue>
#include <iostream>
#include <cmath>

using namespace::std;

const int MAXIMUM_SAMPLES = 10000;
const double MAXIMUM_GRADIENT_TOLERANCE = 0.05;
const int MAXIMUM_CELLS = 100;
const int HISTOGRAM_BINS = 8;

MakeFoam::MakeFoam( const MakeFoam& input ) :
	finishedCells(), centerPoints( input.centerPoints ), centerValues( input.centerValues ), cellIntegrals( input.cellIntegrals ),
	integratePDF( ClassLookUp::CopyPDF( input.integratePDF ) )
{
	for( unsigned int i=0; i< input.finishedCells.size(); ++i )
	{
		finishedCells.push_back( new PhaseSpaceBoundary( *(input.finishedCells[i]) ) );
	}
}

//Constructor with correct argument
MakeFoam::MakeFoam( IPDF * InputPDF, PhaseSpaceBoundary * InputBoundary, DataPoint * InputPoint ) : finishedCells(), centerPoints(), centerValues(), cellIntegrals(), integratePDF(InputPDF)
{
	//Make the container to hold the possible cells
	queue<PhaseSpaceBoundary*> possibleCells;
	PhaseSpaceBoundary* firstCell = InputBoundary;
	possibleCells.push(firstCell);

	//Make a list of observables to integrate over
	vector<string> doIntegrate, dontIntegrate;
	vector<string> pdfDontIntegrate = InputPDF->GetDoNotIntegrateList();
	StatisticsFunctions::DoDontIntegrateLists( InputPDF, InputBoundary, &(pdfDontIntegrate), doIntegrate, dontIntegrate );

	//Continue until all possible cells have been examined
	while ( !possibleCells.empty() )
	{
		//Remove the next possible cell
		PhaseSpaceBoundary* currentCell = possibleCells.front();
		possibleCells.pop();

		//Set up the histogram storage
		vector< vector<double> > histogramBinHeights, histogramBinMiddles, histogramBinMaxes;
		double normalisation = 0.0;
		for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
		{
			vector<double> binHeights, binMiddles, binMaxes;
			IConstraint * temporaryConstraint = currentCell->GetConstraint( doIntegrate[observableIndex] );
			double minimum = temporaryConstraint->GetMinimum();
			double delta = ( temporaryConstraint->GetMaximum() - minimum ) / (double)HISTOGRAM_BINS;

			//Loop over bins
			for (int binIndex = 0; binIndex < HISTOGRAM_BINS; ++binIndex )
			{
				binHeights.push_back(0.0);
				binMiddles.push_back( minimum + ( delta * ( binIndex + 0.5 ) ) );
				binMaxes.push_back( minimum + ( delta * ( binIndex + 1.0 ) ) );
			}

			histogramBinHeights.push_back(binHeights);
			histogramBinMiddles.push_back(binMiddles);
			histogramBinMaxes.push_back(binMaxes);
		}

		//MC sample the cell, make projections, sort of
		for (int sampleIndex = 0; sampleIndex < MAXIMUM_SAMPLES; ++sampleIndex )
		{
			//Create a data point within the current cell
			DataPoint samplePoint( InputPoint->GetAllNames() );
			for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
			{
				//Generate random values to explore integrable observables
				IConstraint * temporaryConstraint = currentCell->GetConstraint( doIntegrate[observableIndex] );
				samplePoint.SetObservable( doIntegrate[observableIndex], temporaryConstraint->CreateObservable() );
			}
			for (unsigned int observableIndex = 0; observableIndex < dontIntegrate.size(); ++observableIndex )
			{
				//Use given values for unintegrable observables
				Observable * temporaryObservable = new Observable( *(InputPoint->GetObservable( dontIntegrate[observableIndex] )) );
				samplePoint.SetObservable( dontIntegrate[observableIndex], temporaryObservable );
				delete temporaryObservable;
			}

			//Evaluate the function at this point
			double sampleValue = InputPDF->Evaluate( &samplePoint );
			normalisation += sampleValue;

			//Populate the histogram
			for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
			{
				double observableValue = samplePoint.GetObservable( doIntegrate[observableIndex] )->GetValue();

				for ( int binIndex = 0; binIndex < HISTOGRAM_BINS; ++binIndex )
				{
					if ( observableValue < histogramBinMaxes[observableIndex][unsigned(binIndex)] )
					{
						histogramBinHeights[observableIndex][unsigned(binIndex)] += sampleValue;
						break;
					}
				}
			}
		}

		//Normalise the histograms
		for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
		{
			for ( int binIndex = 0; binIndex < HISTOGRAM_BINS; ++binIndex )
			{
				histogramBinHeights[observableIndex][unsigned(binIndex)] /= normalisation;
			}
		}

		//Find the maximum gradient
		vector<double> midPoints;
		string maximumGradientObservable, unit;
		double maximumGradient=0.;
		double lowPoint=0.;
		double splitPoint=0.;
		double highPoint=0.;
		for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
		{
			//Find the size of the cell in this observable
			IConstraint * temporaryConstraint = currentCell->GetConstraint( doIntegrate[observableIndex] );
			double cellMaximum = temporaryConstraint->GetMaximum();
			double cellMinimum = temporaryConstraint->GetMinimum();

			//Store the mid point
			double observableMiddle = cellMinimum + ( ( cellMaximum - cellMinimum ) / 2.0 );
			midPoints.push_back(observableMiddle);

			for ( int binIndex = 1; binIndex < HISTOGRAM_BINS; ++binIndex )
			{
				double gradient = abs( histogramBinHeights[observableIndex][ unsigned(binIndex - 1) ] - histogramBinHeights[observableIndex][unsigned(binIndex)] );

				//Update maximum
				if ( ( observableIndex == 0 && binIndex == 1 ) || gradient > maximumGradient )
				{
					maximumGradient = gradient;
					maximumGradientObservable = doIntegrate[observableIndex];
					unit = temporaryConstraint->GetUnit();
					highPoint = cellMaximum;
					splitPoint = ( histogramBinMiddles[observableIndex][ unsigned(binIndex - 1) ] + histogramBinMiddles[observableIndex][unsigned(binIndex)] ) / 2.0;
					lowPoint = cellMinimum;
				}
			}
		}

		//If the maximum gradient is within tolerance, the cell is finished
		if ( maximumGradient < MAXIMUM_GRADIENT_TOLERANCE )
		{
			//Store the finished cell
			finishedCells.push_back(currentCell);

			//Create a data point at the center of the current cell
			DataPoint* cellCenter = new DataPoint( InputPoint->GetAllNames() );
			for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
			{
				//Use the mid points for the integrable values
				Observable * temporaryObservable = cellCenter->GetObservable( doIntegrate[observableIndex] );
				Observable* temporaryObservable2 = new Observable( temporaryObservable->GetName(), midPoints[observableIndex], temporaryObservable->GetUnit() );
				cellCenter->SetObservable( doIntegrate[observableIndex], temporaryObservable2 );
				delete temporaryObservable2;
			}
			for (unsigned int observableIndex = 0; observableIndex < dontIntegrate.size(); ++observableIndex )
			{
				//Use given values for unintegrable observables
				Observable * temporaryObservable = new Observable( *(InputPoint->GetObservable( dontIntegrate[observableIndex] )) );
				cellCenter->SetObservable( dontIntegrate[observableIndex], temporaryObservable );
				delete temporaryObservable;
			}

			//Store the center point
			centerPoints.push_back(cellCenter);
		}
		else
		{
			//Create two cells to replace the current cell
			PhaseSpaceBoundary* daughterCell1 = new PhaseSpaceBoundary( currentCell->GetAllNames() );
			PhaseSpaceBoundary* daughterCell2 = new PhaseSpaceBoundary( currentCell->GetAllNames() );

			for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
			{
				IConstraint * temporaryConstraint = currentCell->GetConstraint( doIntegrate[observableIndex] );
				if ( doIntegrate[observableIndex] == maximumGradientObservable )
				{
					//Split the cells on the observable with the greatest gradient
					daughterCell1->SetConstraint( doIntegrate[observableIndex], lowPoint, splitPoint, unit );
					daughterCell2->SetConstraint( doIntegrate[observableIndex], splitPoint, highPoint, unit );
				}
				else
				{
					//Copy the continuous constraint (if it can be integrated, it must be continuous)
					daughterCell1->SetConstraint( doIntegrate[observableIndex], temporaryConstraint->GetMinimum(), temporaryConstraint->GetMaximum(), temporaryConstraint->GetUnit() );
					daughterCell2->SetConstraint( doIntegrate[observableIndex], temporaryConstraint->GetMinimum(), temporaryConstraint->GetMaximum(), temporaryConstraint->GetUnit() );
				}
			}
			for (unsigned int observableIndex = 0; observableIndex < dontIntegrate.size(); ++observableIndex )
			{
				IConstraint * temporaryConstraint = currentCell->GetConstraint( dontIntegrate[observableIndex] );
				if ( temporaryConstraint->IsDiscrete() )
				{
					//Copy the discrete constraint
					daughterCell1->SetConstraint( dontIntegrate[observableIndex], temporaryConstraint->GetValues(), temporaryConstraint->GetUnit() );
					daughterCell2->SetConstraint( dontIntegrate[observableIndex], temporaryConstraint->GetValues(), temporaryConstraint->GetUnit() );
				}
				else
				{
					//Copy the continuous constraint
					daughterCell1->SetConstraint( dontIntegrate[observableIndex], temporaryConstraint->GetMinimum(), temporaryConstraint->GetMaximum(), temporaryConstraint->GetUnit() );
					daughterCell2->SetConstraint( dontIntegrate[observableIndex], temporaryConstraint->GetMinimum(), temporaryConstraint->GetMaximum(), temporaryConstraint->GetUnit() );
				}
			}

			//Add the two new cells to the possibles
			possibleCells.push(daughterCell1);
			possibleCells.push(daughterCell2);
		}

		//Make sure you don't exceed the maximum number of cells
		if ( int(finishedCells.size() + possibleCells.size()) >= MAXIMUM_CELLS )
		{
			cout << "MakeFoam warning: maximum cells reached with " << possibleCells.size() << " unexplored" << endl;

			//Dump out any possible cells into finished, without any further examination
			while ( !possibleCells.empty() )
			{
				//Move the cell from "possible" to "finished"
				PhaseSpaceBoundary* temporaryCell = possibleCells.front();
				possibleCells.pop();
				finishedCells.push_back(temporaryCell);

				//Create a data point at the center of the current cell
				DataPoint* cellCenter = new DataPoint( InputPoint->GetAllNames() );
				for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); ++observableIndex )
				{
					//Calculate the cell mid point
					IConstraint * temporaryConstraint = temporaryCell->GetConstraint( doIntegrate[observableIndex] );
					double midPoint = temporaryConstraint->GetMinimum() + ( ( temporaryConstraint->GetMaximum() - temporaryConstraint->GetMinimum() ) / 2  );

					//Use the mid points for the integrable values
					Observable * temporaryObservable = cellCenter->GetObservable( doIntegrate[observableIndex] );
					Observable* temporaryObservable2 = new Observable( temporaryObservable->GetName(), midPoint, temporaryObservable->GetUnit() );
					cellCenter->SetObservable( doIntegrate[observableIndex], temporaryObservable2 );
					delete temporaryObservable2;
				}
				for (unsigned int observableIndex = 0; observableIndex < dontIntegrate.size(); ++observableIndex )
				{
					//Use given values for unintegrable observables
					Observable * temporaryObservable = new Observable( *(InputPoint->GetObservable( dontIntegrate[observableIndex] )) );
					cellCenter->SetObservable( dontIntegrate[observableIndex], temporaryObservable );
					delete temporaryObservable;
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

	//Find the function value at the center of each cell, and the integral of the function over the cell
	for (unsigned int cellIndex = 0; cellIndex < finishedCells.size(); ++cellIndex )
	{
		//Integrate the cell
		double integral = cellIntegrator.Integral( InputPoint, finishedCells[cellIndex] );
		cellIntegrals.push_back(integral);

		//Evaluate the function at the center of the cell
		double value = InputPDF->Evaluate( centerPoints[cellIndex] );
		centerValues.push_back(value);
	}
}

//Destructor
MakeFoam::~MakeFoam()
{
	while( !centerPoints.empty() )
	{
		if( centerPoints.back() != NULL ) delete centerPoints.back();
		centerPoints.pop_back();
	}
	while( !finishedCells.empty() )
	{
		if( finishedCells.back() != NULL ) delete finishedCells.back();
		finishedCells.pop_back();
	}
	if( integratePDF != NULL ) delete integratePDF;
}

//Output all cells to console
void MakeFoam::Debug()
{
	cout << centerPoints.size() << " cells generated" << endl << endl;
	vector<string> allNames = centerPoints[0]->GetAllNames();
	vector<string>::iterator observableIterator;
	for ( observableIterator = allNames.begin(); observableIterator != allNames.end(); ++observableIterator )
	{
		cout << *observableIterator << ", ";
	}
	cout << endl;
	for (unsigned int pointIndex = 0; pointIndex < centerPoints.size(); ++pointIndex )
	{
		for ( observableIterator = allNames.begin(); observableIterator != allNames.end(); ++observableIterator )
		{
			cout << centerPoints[pointIndex]->GetObservable( *observableIterator )->GetValue() << ", ";
		}
		cout << endl;
	}

	for (unsigned int cellIndex = 0; cellIndex < finishedCells.size(); ++cellIndex )
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
	for (unsigned int cellIndex = 0; cellIndex < centerPoints.size(); ++cellIndex )
	{
		//Evaluate the function at the center of the cell
		double newCenterValue = integratePDF->Evaluate( centerPoints[cellIndex] );

		//Return the precalculated integral of the cell, scaled to the new function value
		double newIntegral = cellIntegrals[cellIndex] * newCenterValue / centerValues[cellIndex];
		totalIntegral += newIntegral;
	}

	return totalIntegral;
}
