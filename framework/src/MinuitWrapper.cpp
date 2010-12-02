/**
  @class MinuitWrapper

  A wrapper to integrate Minuit with RapidFit

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

#include "MinuitWrapper.h"
#include <iostream>
#include "Rtypes.h"
#include <limits>
#include "StringProcessing.h"
#include "FunctionContour.h"
#include <ctime>

const double MAXIMUM_MINIMISATION_STEPS = 100000.0;//800.0;
const double FINAL_GRADIENT_TOLERANCE = 0.01;//0.001;
const double STEP_SIZE = 0.01;
FitFunction * MinuitWrapper::function = 0;

//Default constructor
MinuitWrapper::MinuitWrapper()
{
	minuit = new TMinuit( 1 );
}

//Constructor with correct argument
MinuitWrapper::MinuitWrapper( int NumberParameters )
{
	minuit = new TMinuit( NumberParameters );
}

//Destructor
MinuitWrapper::~MinuitWrapper()
{
	delete minuit;
}

//Use Migrad to minimise the given function
void MinuitWrapper::Minimise( FitFunction * NewFunction )
{
	function = NewFunction;
	int errorFlag = 0;
	double arguments[2] = {0.0, 0.0};

	//Set the function
	minuit->SetFCN( Function );

	//Store the parameters
	ParameterSet * newParameters = NewFunction->GetParameterSet();
	vector<string> allNames = newParameters->GetAllNames();
	for (int nameIndex = 0; nameIndex < allNames.size(); nameIndex++)
	{
		PhysicsParameter * newParameter = newParameters->GetPhysicsParameter( allNames[nameIndex] );

		//Make bounded or unbounded parameters
		if ( newParameter->GetType() == "Unbounded" || newParameter->GetType() == "GaussianConstrained" )
		{
			minuit->mnparm(nameIndex, allNames[nameIndex], newParameter->GetValue(), STEP_SIZE, 0.0, 0.0, errorFlag);
		}
		else
		{
			minuit->mnparm(nameIndex, allNames[nameIndex], newParameter->GetValue(), STEP_SIZE,
					newParameter->GetMinimum(), newParameter->GetMaximum(), errorFlag);
		}

		//Fix the parameter if required
		if ( newParameter->GetType() == "Fixed" )
		{
			minuit->FixParameter( nameIndex );
		}
		else
		{
			minuit->Release( nameIndex );
		}
	}

	//Set the error analysis
	arguments[0] = NewFunction->UpErrorValue(1);
	minuit->mnexcm("SET ERR", arguments, 1, errorFlag);

	
	arguments[0] = MAXIMUM_MINIMISATION_STEPS;
	arguments[1] = FINAL_GRADIENT_TOLERANCE;
        
	//Now call HESSE to get a rough starting estimate for the error matrix
        //minuit->mnexcm("HESSE", arguments, 1, errorFlag);
	
	//Do the minimisation
	minuit->mnexcm("MIGRAD", arguments, 2, errorFlag);

        //Now call HESSE to properly calculate the error matrix
        minuit->mnexcm("HESSE", arguments, 1, errorFlag);

        //Call MINOS to calculate the non-parabolic errors.
        //MINOS rather than assuming parabolic shape about the minimum, MINOS climbs
        //out of the minimum on either side up until the UP value.
        //minuit->mnexcm("MINOS", arguments, 2, errorFlag);

	//Output time information
	time_t timeNow;
	time(&timeNow);
	cout << "Minuit finished: " << ctime( &timeNow ) << endl;

	//Get the fit status
	double minimumValue = 0.0;
	double fedm = 0.0;
	double errdef = 0.0;
	int variableParameters = 0;
	int parameterNumber = 0;
	int fitStatus = 0;
	minuit->mnstat( minimumValue, fedm, errdef, variableParameters, parameterNumber, fitStatus );

	//Get the fitted parameters
	ResultParameterSet * fittedParameters = new ResultParameterSet( allNames );
	for (int nameIndex = 0; nameIndex < allNames.size(); nameIndex++)
	{
		string parameterName = allNames[nameIndex];
		double parameterValue = 0.0;
		double parameterError = 0.0;
		double xlolim = 0.0;
		double xuplim = 0.0;
		int iuint = 0;
		TString temporaryString;
		minuit->mnpout( nameIndex, temporaryString, parameterValue, parameterError, xlolim, xuplim, iuint );

		PhysicsParameter * oldParameter = newParameters->GetPhysicsParameter( parameterName );
		fittedParameters->SetResultParameter( parameterName, parameterValue, oldParameter->GetOriginalValue(), parameterError,
				-numeric_limits<double>::max(), numeric_limits<double>::max(),
				oldParameter->GetType(), oldParameter->GetUnit() );
	}

	// Get the error matrix and construct a vector containing the correlation coefficients
	int numParams = allNames.size();
	double matrix[numParams][numParams];
	gMinuit->mnemat(&matrix[0][0],numParams);
	vector<double> covarianceMatrix(numParams*(numParams+1)/2);
	for (int row = 0; row < numParams; row++)
	{
		for (int col = 0; col < numParams; col++)
		{
			if(row > col) covarianceMatrix[col+row*(row+1)/2] = matrix[row][col];
			else covarianceMatrix[row+col*(col+1)/2] = matrix[row][col];
		}
	}

	//Make the contour plots
	vector< FunctionContour* > allContours;
	for ( int plotIndex = 0; plotIndex < contours.size(); plotIndex++ )
	{
		//Find the index (in Minuit) of the parameter names requested for the plot
		int xParameterIndex = StringProcessing::VectorContains( &allNames, &( contours[plotIndex].first ) );
		int yParameterIndex = StringProcessing::VectorContains( &allNames, &( contours[plotIndex].second ) );
		if ( xParameterIndex == -1 )
		{
			cerr << "Contour plotting skipped: x parameter \"" << contours[plotIndex].first << "\" not found" << endl;
		}
		else if ( newParameters->GetPhysicsParameter( contours[plotIndex].first )->GetType() == "Fixed" )
		{
			cerr << "Contour plotting skipped: x parameter \"" << contours[plotIndex].first << "\" is fixed" << endl;
		}
		else if ( yParameterIndex == -1 )
		{
			cerr << "Contour plotting skipped: y parameter \"" << contours[plotIndex].second << "\" not found" << endl;
		}
		else if ( newParameters->GetPhysicsParameter( contours[plotIndex].second )->GetType() == "Fixed" )
		{
			cerr << "Contour plotting skipped: y parameter \"" << contours[plotIndex].second << "\" is fixed" << endl;
		}
		else
		{
			//If the parameters have valid indices, ask minuit to plot them
			int numberOfPoints = 40;
			int iErrf;
			double xCoordinates1[numberOfPoints], yCoordinates1[numberOfPoints];
			double xCoordinates2[numberOfPoints], yCoordinates2[numberOfPoints];

			//One sigma contour
			minuit->SetErrorDef( NewFunction->UpErrorValue(1) );
			minuit->mncont( xParameterIndex, yParameterIndex, numberOfPoints, xCoordinates1, yCoordinates1, iErrf );

			//Two sigma contour
			minuit->SetErrorDef( NewFunction->UpErrorValue(2) );
			minuit->mncont( xParameterIndex, yParameterIndex, numberOfPoints, xCoordinates2, yCoordinates2, iErrf );

			//Store the contours
			FunctionContour * newContour = new FunctionContour( contours[plotIndex].first, contours[plotIndex].second, 2 );
			newContour->SetPlot( 1, numberOfPoints, xCoordinates1, yCoordinates1 );
			newContour->SetPlot( 2, numberOfPoints, xCoordinates2, yCoordinates2 );
			allContours.push_back(newContour);
		}
	}

	fitResult = new FitResult( minimumValue, fittedParameters, fitStatus, *( function->GetPhysicsBottle() ), covarianceMatrix, allContours );
}

//The function to pass to Minuit
void Function( int & npar, double * grad, double & fval, double * xval, int iflag )
{
	ParameterSet * temporaryParameters = MinuitWrapper::function->GetParameterSet();
	
	if ( temporaryParameters->SetPhysicsParameters(xval) )
	{
		MinuitWrapper::function->SetParameterSet(temporaryParameters);
		fval = MinuitWrapper::function->Evaluate();
	}
	else
	{
		cerr << "Failed to set physics parameters" << endl;
	}
}

//Return the result of minimisation
FitResult * MinuitWrapper::GetFitResult()
{
	return fitResult;
}

//Request contour plots
void MinuitWrapper::ContourPlots( vector< pair< string, string > > ContourParameters )
{
	contours = ContourParameters;
}
