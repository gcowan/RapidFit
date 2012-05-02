/**
  @class Plotter

  A class for plotting PDF projections onto histograms

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-13
  */

//	ROOT Headers
#include "TFile.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TFolder.h"
//	RapidFit Headers
#include "Plotter.h"
#include "StatisticsFunctions.h"
#include "EdStyle.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <float.h>

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

using namespace::std;

//Default constructor
Plotter::Plotter(): weightName(), plotPDF(), plotData(), pdfIntegrator(), weightsWereUsed()
{
}

//Constructor with correct arguments
Plotter::Plotter( IPDF * NewPDF, IDataSet * NewDataSet ) : weightName(), plotPDF(NewPDF), plotData(NewDataSet), pdfIntegrator( new RapidFitIntegrator(NewPDF) ), weightsWereUsed(false)
{
	pdfIntegrator->ProjectionSettings();
}

//Destructor
Plotter::~Plotter()
{
}

//Create a root file containing a projection plot over each observable
void Plotter::PlotAllObservables( string FileName )
{
	//Work out what to plot
	vector<double> combinationWeights;
	vector<string> combinationDescriptions;
	vector<DataPoint*> allCombinations = GetDiscreteCombinations( combinationWeights, combinationDescriptions );
	vector<string> doNotIntegrate = plotPDF->GetDoNotIntegrateList();

	//Make the file
	TFile * rootFile = new TFile( FileName.c_str(), "RECREATE" );

	//Loop over all (integrable) observables
	vector<string> observableNames = plotPDF->GetPrototypeDataPoint();
	vector<string>::iterator observableIterator;
	for ( observableIterator = observableNames.begin(); observableIterator != observableNames.end(); ++observableIterator )
	{
		//Check the observable can be plotted
		bool continuous = !( plotData->GetBoundary()->GetConstraint( *observableIterator )->IsDiscrete() );
		bool doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrate, &(*observableIterator) ) == -1 );
		if ( continuous && doIntegrate )
		{
			cout << "Projecting " << *observableIterator << endl;
			MakeObservablePlots( *observableIterator, allCombinations, combinationWeights, combinationDescriptions, rootFile );
		}
	}

	rootFile->Close();
}

//Create a root file containing a projection plot over one observable
void Plotter::PlotObservables( string FileName, vector<string> ObservableNames )
{
	//Work out what to plot
	vector<double> combinationWeights;
	vector<string> combinationDescriptions;
	vector<DataPoint*> allCombinations = GetDiscreteCombinations( combinationWeights, combinationDescriptions );
	vector<string> doNotIntegrate = plotPDF->GetDoNotIntegrateList();

	//Make the file
	TFile * rootFile = new TFile( FileName.c_str(), "RECREATE" );

	//Loop over all observables
	for (unsigned int observableIndex = 0; observableIndex < ObservableNames.size(); ++observableIndex )
	{
		//Check the observable can be plotted
		bool continuous = !( plotData->GetBoundary()->GetConstraint( ObservableNames[observableIndex] )->IsDiscrete() );
		bool doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrate, &( ObservableNames[observableIndex] ) ) == -1 );
		if ( continuous && doIntegrate )
		{
			//Do the plotting
			cout << "Projecting " << ObservableNames[observableIndex] << endl;
			MakeObservablePlots( ObservableNames[observableIndex], allCombinations, combinationWeights, combinationDescriptions, rootFile );
		}
	}

	rootFile->Close();
}

//Make a directory containing all the plots for a particular observable
//NOTE: this function assumes there is an open TFile already
void Plotter::MakeObservablePlots( string ObservableName, vector<DataPoint*> AllCombinations, vector<double> CombinationWeights, vector<string> CombinationDescriptions, TFile * PlotFile )
{
	//////////////////////////////////////////////////////////////////
	//Make the histogram of this observable
	EdStyle * greigFormat = new EdStyle();
	greigFormat->SetStyle();
        gStyle->SetMarkerStyle(15);
        gStyle->SetMarkerSize(Size_t(0.8));
	//Get the data needed to make the histogram
	double maximum, minimum;
	int binNumber;
	vector<double> observableValues = GetStatistics( ObservableName, minimum, maximum, binNumber );
	double binInterval = ( maximum - minimum ) / (double)binNumber;

	cout << setprecision(10) << "Min: " << minimum << " Max: " << maximum << " binNumber: " << binNumber << "  binInterval: " << binInterval << endl;

	//New code to deal with weighted events
	vector<double> observableWeights ;
	double normalisationFactor=0. ;
	if( weightsWereUsed ) {
		double dum1, dum2 ; int idum ;
		observableWeights = GetStatistics( weightName, dum1, dum2, idum );
		double sumWeights =0;
		for( unsigned int iw=0;iw<observableWeights.size();iw++) sumWeights+=observableWeights[iw];
		normalisationFactor = sumWeights/double(observableWeights.size()) ;
	}

	//Make the histogram
	string histogramName = ObservableName + "ProjectionPlot";
	string histogramTitle = "";
	TH1F * dataHistogram = new TH1F( histogramName.c_str(), histogramTitle.c_str(), binNumber, minimum, maximum );

	dataHistogram->Sumw2();

	//Loop over all data points and add them to the histogram
	for (unsigned int dataIndex = 0; dataIndex < observableValues.size(); dataIndex++)
	{
		if( weightsWereUsed ) dataHistogram->Fill( observableValues[dataIndex], observableWeights[dataIndex] / normalisationFactor );
		else dataHistogram->Fill( observableValues[dataIndex] );
	}

	//Histogram complete
	////////////////////////////////////////////////////////////////////

	//Make and enter the directory
	string directoryName = ObservableName + "Directory";
	string directoryTitle = ObservableName + " directory";
	TDirectory * combinationDirectory = PlotFile->mkdir( directoryName.c_str(), directoryTitle.c_str() );
	combinationDirectory->cd();

	////////////////////////////////////////////////////////////////////
	//Calculate the projection for each combination, and the data averaged projection

	vector<double> dataAverageProjection;
	double averageIntegral = 0.0;
	double plotInterval;
	double ratioOfIntegrals = 1.;

	int plotNumber=0;
	if (ObservableName == "time") plotNumber = 256;
	else plotNumber = 128;

	//Initialise the data averaged projection
	for (int pointIndex = 0; pointIndex < plotNumber; ++pointIndex )
	{
		dataAverageProjection.push_back(0.0);
	}

	double test=0.;

	//Loop over all discrete combinations
	for (unsigned int combinationIndex = 0; combinationIndex < AllCombinations.size(); ++combinationIndex )
	{
		cout << "Constructing PDF Integral of " << ObservableName << " " << combinationIndex+1 << " of " << AllCombinations.size() << endl;
		//Calculate the projection for this combination
		vector<double> projectionValueVector = ProjectObservable( *AllCombinations[combinationIndex], ObservableName, minimum, maximum, plotNumber, plotInterval );
		double projectionIntegral = pdfIntegrator->Integral( AllCombinations[combinationIndex], plotData->GetBoundary() );
		ratioOfIntegrals = pdfIntegrator->GetRatioOfIntegrals();

		//Update the data average values, and make the projection graph arrays
		double* projectionValueArray = new double[unsigned(plotNumber)];
		double* observableValueArray = new double[unsigned(plotNumber)];
		cout << "Constructing DataSet Integral of " << ObservableName << " " << combinationIndex+1 << " of " << AllCombinations.size() << endl;
		for (int pointIndex = 0; pointIndex < plotNumber; ++pointIndex )
		{
			cout << setprecision(3) << setw(4) << (pointIndex+1.)/plotNumber*100 << "\%\tComplete\r\r" << flush;
			//Projection graph
			projectionValueArray[unsigned(pointIndex)] = ratioOfIntegrals * projectionValueVector[unsigned(pointIndex)] * binInterval * double(observableValues.size()) / projectionIntegral;
			observableValueArray[unsigned(pointIndex)] = minimum + ( plotInterval * pointIndex );

			//Data average
			dataAverageProjection[unsigned(pointIndex)] += projectionValueVector[unsigned(pointIndex)] * CombinationWeights[combinationIndex];
		}
		cout << endl;
		test = projectionValueArray[0];
		cout << "ratioOfIntegrals: " << ratioOfIntegrals << endl;
		cout << "binInterval: " << binInterval << endl;
		cout << "observableValues.size(): " << observableValues.size() << endl;
		cout << "Integral: " << projectionIntegral << endl;
		cout << "CombinationWeights[combinationIndex]: " << CombinationWeights[combinationIndex] << endl;

		averageIntegral += projectionIntegral * CombinationWeights[combinationIndex];

		//Plot the projection
		MakePlotCanvas( ObservableName, CombinationDescriptions[combinationIndex], dataHistogram, observableValueArray, projectionValueArray, plotNumber );
		cout << "Finished combination " << combinationIndex+1 << " of " << AllCombinations.size() << endl;
	}


	//Projecting complete
	/////////////////////////////////////////////////////////////////////

	cout << "Projection for " << ObservableName << " Finished." << endl<<endl;
	//Leave the directory
	combinationDirectory->GetMotherDir()->cd();

	//Plot the data averaged projection
	double* projectionValueArray = new double[unsigned(plotNumber)];
	double* observableValueArray = new double[unsigned(plotNumber)];
	for (int pointIndex = 0; pointIndex < plotNumber; ++pointIndex )
	{
		projectionValueArray[pointIndex] = ratioOfIntegrals * dataAverageProjection[unsigned(pointIndex)] * binInterval * double(observableValues.size()) / averageIntegral;
		observableValueArray[pointIndex] = minimum + ( plotInterval * pointIndex );
	}
	cout << "projectionValueArray[0]:  was " << test << " is now: " << projectionValueArray[0] << " avgInt: " << averageIntegral << endl;
	MakePlotCanvas( ObservableName, "", dataHistogram, observableValueArray, projectionValueArray, plotNumber );
}

//Create a projection plot for a given observable
//NOTE: this function assumes there is an open TFile already
void Plotter::MakePlotCanvas( string ObservableName, string Description, TH1F * Histogram, double * ProjectionXValues, double * ProjectionYValues, int PlotNumber )
{
	//Make the graph
	EdStyle * greigFormat = new EdStyle();
	greigFormat->SetStyle();
	gStyle->SetMarkerStyle(15);
        gStyle->SetMarkerSize(Size_t(0.8));
	TMultiGraph * graph = new TMultiGraph();
	TGraphErrors * projectionGraph = new TGraphErrors( PlotNumber, ProjectionXValues, ProjectionYValues );
	graph->Add(projectionGraph);

	//Formatting
	projectionGraph->SetLineColor(2);
	projectionGraph->SetLineWidth(4);
	string xTitle = ObservableName + " (" + plotData->GetDataPoint(0)->GetObservable(ObservableName)->GetUnit() + ")";
	Histogram->GetXaxis()->SetTitle( xTitle.c_str() );
	Histogram->GetYaxis()->SetTitle( "Events" );

	string drawOptions = "C";
	//Note "same" option is not required if graph is drawn second. "A" option - requesting axes - will overwrite whatever was there before

	// Just to get a unique name for the canvas
	TRandom3 * random = new TRandom3(0);
	double ran = random->Rndm();
	ostringstream ranString;
	ranString << ran;
	//Stick both objects on one canvas
	string canvasName = ObservableName + "Projection" + Description;// + ranString.str();
	//string canvasTitle = ObservableName + " projection " + Description;
	string canvasTitle = "";
	TCanvas * bothPlots = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), 600, 600 );
 	//gStyle->SetMarkerStyle(20);
	Histogram->SetStats(0);

//	double ymin = Histogram->GetBinContent(Histogram->GetMinimumBin());
	double ymax = Histogram->GetBinContent(Histogram->GetMaximumBin());

	if (ObservableName == "time" || ObservableName == "B_s_TAU"){
		// Attempted hack to get round the fact that ROOT is rubbish!
		//TH1F * tmpHist = new TH1F("tmp", "tmp", 1, -2., -1.);
		//tmpHist->Fill(-2.);
		//tmpHist->Draw();
		//Histogram->SetStats(1);
		bothPlots->SetLogy();
		Histogram->GetXaxis()->SetRangeUser(-2., 15.); // This does not work because of STUPID ROOT!!
		Histogram->GetYaxis()->SetRangeUser(0.1, 2*ymax);
		//projectionGraph->SetMinimum(-2.); // Don't think this works either
	}
	else
	{
		Histogram->GetYaxis()->SetRangeUser(0., ymax + ymax/3.);
	}

	TPaveText *lhcb7TeVPrelimR = new TPaveText(0.65,
					   0.75,
					   0.90,
					   0.85,
					   "BRNDC");
	lhcb7TeVPrelimR->SetFillColor(0);
	lhcb7TeVPrelimR->SetTextAlign(12);
	lhcb7TeVPrelimR->SetBorderSize(0);
	lhcb7TeVPrelimR->AddText("#splitline{#splitline{LHCb}{Preliminary}}{#scale[0.7]{#sqrt{s} = 7 TeV}}");

	Histogram->Draw("E1");
	graph->Draw( drawOptions.c_str() );
	lhcb7TeVPrelimR->Draw();

	bothPlots->Update();
	bothPlots->Write();
}

//Return the values tracing the PDF projection
vector<double> Plotter::ProjectObservable( DataPoint InputPoint, string ObservableName )
{
	double maximum, minimum, plotInterval;
	int plotNumber;
	GetStatistics( ObservableName, minimum, maximum, plotNumber );
	return ProjectObservable( InputPoint, ObservableName, minimum, maximum, plotNumber, plotInterval );
}

//Return the values tracing the PDF projection
vector<double> Plotter::ProjectObservable( DataPoint InputPoint, string ObservableName, double Minimum, double Maximum, int PlotNumber, double & PlotInterval )
{
	//cout << "Projecting Observable:\t" << ObservableName << "\tmin: " << Minimum << "\tmax: " << Maximum << "\tnum: " << PlotNumber << "\tInterval: " << PlotInterval << endl;
	PlotInterval = ( Maximum - Minimum ) / (double)( PlotNumber - 1 );

	//Find the value of the observable projection at each data point
	vector<double> pointValues;
	//cout << "looping over " << PlotNumber << " points" << endl;
	for (int pointIndex = 0; pointIndex < PlotNumber; ++pointIndex )
	{
		cout << setprecision(3) << setw(4) << (pointIndex+1.)/PlotNumber*100 << "\%\tComplete\r\r" << flush;
		//cout << "Integration " << pointIndex + 1 << " of " << PlotNumber << endl;
		double observableValue = Minimum + ( PlotInterval * pointIndex );
		Observable* temp1 = InputPoint.GetObservable(ObservableName);
		Observable* temp2 = new Observable( temp1->GetName(), observableValue, temp1->GetError(), temp1->GetUnit() );
		InputPoint.SetObservable( temp1->GetName(), temp2 );
		delete temp2;
		//cout << "Adding pointValue from pdfIntegrator->ProjectObeservable()" << endl;
		pointValues.push_back( pdfIntegrator->ProjectObservable( &InputPoint, plotData->GetBoundary(), ObservableName ) );
	}
	cout << endl;
	return pointValues;
}

//Run the statistics functions
vector<double> Plotter::GetStatistics( string ObservableName, double & Minimum, double & Maximum, int & BinNumber )
{
	//Make a vector of the values found for the observable
	vector<double> observableValues;
	for ( int dataIndex = 0; dataIndex < plotData->GetDataNumber(); ++dataIndex )
	{
		observableValues.push_back( plotData->GetDataPoint(dataIndex)->GetObservable(ObservableName)->GetValue() );
	}

	//Find out number and range of data points needed
	Maximum = StatisticsFunctions::Maximum(observableValues);
	Minimum = StatisticsFunctions::Minimum(observableValues);
	BinNumber = StatisticsFunctions::OptimumBinNumber(observableValues);
	if ( ObservableName == "time" ) BinNumber = (int)BinNumber/2;
	BinNumber = 100;
	return observableValues;
}

//Return a list of data points
//Each should take the data average value of each continuous observable
//Each should represent one combination of possible discrete values
vector<DataPoint*> Plotter::GetDiscreteCombinations( vector<double> & DataPointWeights, vector<string> & DataPointDescriptions )
{
	//Calculate all possiblle combinations of discrete observables
	vector< vector<double> > discreteValues;
	vector<string> discreteNames, continuousNames;
	vector<string> allNames = plotPDF->GetPrototypeDataPoint();
	vector< vector<double> > discreteCombinations = StatisticsFunctions::DiscreteCombinations( &allNames, plotData->GetBoundary(), discreteNames, continuousNames, discreteValues );

	//Initialise the data averaging
	vector<double> continuousSums;
	vector<long> combinationCounts;
	for (unsigned int continuousIndex = 0; continuousIndex < continuousNames.size(); ++continuousIndex )
	{
		continuousSums.push_back(0.0);
	}
	for (unsigned int combinationIndex = 0; combinationIndex < discreteCombinations.size(); ++combinationIndex )
	{
		combinationCounts.push_back(0);
	}

	//Examine the data set. Find the average value for each continuous observable, and the weight for each discrete combination
	for (int dataIndex = 0; dataIndex < plotData->GetDataNumber(); ++dataIndex )
	{
		DataPoint * readDataPoint = plotData->GetDataPoint(dataIndex);

		//Sum the continuous values, in preparation for taking the average
		for (unsigned int continuousIndex = 0; continuousIndex < continuousNames.size(); ++continuousIndex )
		{
			continuousSums[continuousIndex] += readDataPoint->GetObservable( continuousNames[continuousIndex] )->GetValue();
		}

		//Calculate the index for the discrete combination, and increment the corresponding count
		int combinationIndex = 0;
		int incrementValue = 1;
		for (int discreteIndex = int(discreteNames.size() - 1); discreteIndex >= 0; discreteIndex-- )
		{
			double currentValue = readDataPoint->GetObservable( discreteNames[unsigned(discreteIndex)] )->GetValue();

			for (unsigned int valueIndex = 0; valueIndex < discreteValues[unsigned(discreteIndex)].size(); ++valueIndex )
			{
				if ( fabs( discreteValues[unsigned(discreteIndex)][valueIndex] - currentValue ) < DOUBLE_TOLERANCE )
				{
					combinationIndex += ( incrementValue * int(valueIndex) );
					incrementValue *= int(discreteValues[unsigned(discreteIndex)].size());
					break;
				}
			}
		}
		++combinationCounts[unsigned(combinationIndex)];
	}

	//Calculate averages and weights
	vector<double> combinationWeights;
	double dataNumber = (double)plotData->GetDataNumber();
	for (unsigned int continuousIndex = 0; continuousIndex < continuousNames.size(); ++continuousIndex )
	{
		continuousSums[continuousIndex] /= dataNumber;
	}

	for (unsigned int combinationIndex = 0; combinationIndex < discreteCombinations.size(); ++combinationIndex )
	{
		combinationWeights.push_back( (double)combinationCounts[combinationIndex] / dataNumber );
	}

	//Create the data points to return
	vector<DataPoint*> newDataPoints;
	vector<string> allDescriptions;
	DataPoint* templateDataPoint = new DataPoint( *( plotData->GetDataPoint(0) ) );
	for (unsigned int continuousIndex = 0; continuousIndex < continuousNames.size(); ++continuousIndex )
	{
		Observable * newValue = templateDataPoint->GetObservable( continuousNames[continuousIndex] );
		Observable* newValue2 = new Observable( newValue->GetName(), continuousSums[continuousIndex], newValue->GetError(), newValue->GetUnit() );
		templateDataPoint->SetObservable( continuousNames[continuousIndex], newValue2 );
		delete newValue2;
	}
	for (unsigned int combinationIndex = 0; combinationIndex < discreteCombinations.size(); ++combinationIndex )
	{
		string description = "(";

		//Output the discrete values for this combination
		for (unsigned int discreteIndex = 0; discreteIndex < discreteNames.size(); ++discreteIndex )
		{
			//Set the data point
			Observable * newValue = templateDataPoint->GetObservable( discreteNames[discreteIndex] );
			Observable* newValue2 = new Observable( newValue->GetName(), discreteCombinations[combinationIndex][discreteIndex], newValue->GetError(), newValue->GetUnit() );
			templateDataPoint->SetObservable( discreteNames[discreteIndex], newValue2 );
			delete newValue2;

			//Make the description
			char value[100];
			sprintf( value, "%f", discreteCombinations[combinationIndex][discreteIndex] );
			string addToDescription;
			if ( discreteIndex == discreteNames.size() - 1 )
			{
				addToDescription = discreteNames[discreteIndex] + "=" + value;
			}
			else
			{
				addToDescription = discreteNames[discreteIndex] + "=" + value + "_";
			}
			description.append(addToDescription);
		}

		description.append(")");
		allDescriptions.push_back(description);
		newDataPoints.push_back(templateDataPoint);
	}

	//Output the results
	DataPointDescriptions = allDescriptions;
	DataPointWeights = combinationWeights;
	return newDataPoints;
}


//Setting knowledge that weights were uses
void Plotter::SetWeightsWereUsed( string _weightName )
{
	weightsWereUsed = true ;
	weightName = _weightName ;
}
