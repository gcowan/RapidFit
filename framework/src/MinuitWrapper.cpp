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

const double MAXIMUM_MINIMISATION_STEPS = 800.0;
const double FINAL_GRADIENT_TOLERANCE = 0.001;
const double STEP_SIZE = 0.01;
FitFunction * MinuitWrapper::function = 0;
//vector<string> MinuitWrapper::allNames;

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
	//allNames = NewFunction->GetParameterSet()->GetAllNames();
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
		if ( newParameter->GetType() == "Unbounded" )
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
	arguments[0] = NewFunction->UpErrorValue();
	minuit->mnexcm("SET ERR", arguments, 1, errorFlag);

	//Do the minimisation
	arguments[0] = MAXIMUM_MINIMISATION_STEPS;
	arguments[1] = FINAL_GRADIENT_TOLERANCE;
	minuit->mnexcm("MIGRAD", arguments, 2, errorFlag);

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
	// Need to work out a way of pruning the parameters that are not floated in the fit...
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

	// Get the contours
	int numberOfPoints = 40;
	int iErrf;
	double xCoordinates1[numberOfPoints], yCoordinates1[numberOfPoints];
	double xCoordinates2[numberOfPoints], yCoordinates2[numberOfPoints];
	vector< pair< double, double> > oneSigmaContour;
	vector< pair< double, double> > twoSigmaContour;
        vector< vector< pair< double, double> > > contours;
	
	minuit->SetErrorDef(0.5);
	minuit->mncont(0, 1, numberOfPoints, xCoordinates1, yCoordinates1, iErrf);
	minuit->SetErrorDef(2);
	minuit->mncont(0, 1, numberOfPoints, xCoordinates2, yCoordinates2, iErrf);
	
	for ( int point = 0; point < numberOfPoints; point++)
	{
		oneSigmaContour.push_back( pair<double, double> (xCoordinates1[point], yCoordinates1[point]));
		twoSigmaContour.push_back( pair<double, double> (xCoordinates2[point], yCoordinates2[point]));
	}
	contours.push_back(oneSigmaContour);
	contours.push_back(twoSigmaContour);

	fitResult = new FitResult( minimumValue, fittedParameters, fitStatus, *( function->GetPhysicsBottle() ), covarianceMatrix, contours );
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
