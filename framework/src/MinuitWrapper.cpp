/**
  @class MinuitWrapper

  A wrapper to integrate Minuit with RapidFit

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

//	ROOT Headers
#include "Rtypes.h"
//	RapidFit Headers
#include "MinuitWrapper.h"
#include "StringProcessing.h"
#include "FunctionContour.h"
#include "PhysicsBottle.h"
//	System Headers
#include <iostream>
#include <limits>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <float.h>

#define DOUBLE_TOLERANCE DBL_MIN

//const double MAXIMUM_MINIMISATION_STEPS = 100000.0;//800.0;
//const double FINAL_GRADIENT_TOLERANCE = 0.001;//0.001;
//const double STEP_SIZE = 0.001;

//	Required to act as a 'global' pointer to the RapidFit fitfunction member of the MinuitWrapper Class
FitFunction * MinuitWrapper::function = 0;

//Default constructor
MinuitWrapper::MinuitWrapper(): minuit(), fitResult(), contours(), print_verbosity( 0 )
{
	minuit = new TMinuit( 1 );
}

//Constructor with correct argument
MinuitWrapper::MinuitWrapper( int NumberParameters, int output_level ): minuit(), fitResult(), contours(), print_verbosity( output_level )
{
	minuit = new TMinuit( NumberParameters );
}

//Destructor
MinuitWrapper::~MinuitWrapper()
{
	//cout << "hello from MinuitWrapper destructor" << endl;
	delete minuit;
}

void MinuitWrapper::SetSteps( int newSteps )
{
	maxSteps = newSteps;
}

void MinuitWrapper::SetTolerance( double newTolerance )
{
	bestTolerance = newTolerance;
}

void MinuitWrapper::SetOptions( vector<string> newOptions )
{
	Options = newOptions;
}

void MinuitWrapper::SetQuality( int newQuality )
{
	Quality = newQuality;
}

void MinuitWrapper::SetOutputLevel( int output_level )
{
	print_verbosity = output_level;
	Double_t arg=output_level;
	Int_t err;
	if( output_level < 0 )
	{
		arg=-1;
		minuit->mnexcm("SET PRINT", &arg, 1, err);
		minuit->mnexcm("SET NOWarnings", &arg, 1, err );
		err = minuit->Command("SET OUTputfile /dev/zero");
	} else {
		minuit->mnexcm("SET PRINT", &arg, 1, err);
	}
//	minuit->SetPrintLevel( print_verbosity );
}

//Use Migrad to minimise the given function
void MinuitWrapper::Minimise( FitFunction * NewFunction )
{
	function = NewFunction;
	int errorFlag = 0;
	double arguments[2] = {0.0, 0.0};

	//Set the function
	minuit->SetFCN( &MinuitWrapper::Function );

	//Store the parameters
	ParameterSet * newParameters = NewFunction->GetParameterSet();
	vector<string> allNames = newParameters->GetAllNames();
	for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex)
	{
		PhysicsParameter * newParameter = newParameters->GetPhysicsParameter( allNames[nameIndex] );

		double STEP_SIZE= 0.01;
		if( !( fabs( newParameter->GetMaximum() - newParameter->GetMinimum() ) < DOUBLE_TOLERANCE  ) ){
			STEP_SIZE = fabs((newParameter->GetMaximum() - newParameter->GetMinimum()))/10000.0;
		}

		if( newParameter->GetStepSize() > 0 ) { STEP_SIZE = newParameter->GetStepSize(); };

		newParameter->SetStepSize( STEP_SIZE );

		//cout << "STEP SIZE:" << allNames[nameIndex] << "\t"<<nameIndex << "\t"<<STEP_SIZE <<endl;

		//Make bounded or unbounded parameters
		if ( newParameter->GetType() == "Unbounded" || newParameter->GetType() == "GaussianConstrained" )
		{
			minuit->mnparm(nameIndex, allNames[nameIndex], newParameter->GetBlindedValue(), STEP_SIZE, 0.0, 0.0, errorFlag);
		}
		else
		{
			minuit->mnparm(nameIndex, allNames[nameIndex], newParameter->GetBlindedValue(), STEP_SIZE,
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

	//	Syntax for Minuit Commands through the TMinuit Class:
	//
	//	minuit->mnexcm("SOMECOMMAND",&SOMEARGUMENTS,NUMBEROFARGUMENTS,ERRORFLAG);

	//	Set the error analysis
	arguments[0] = NewFunction->UpErrorValue(1);
	minuit->mnexcm("SET ERR", arguments, 1, errorFlag);

	//	Set Migrad Strategy 2
	arguments[0] = Quality;//1;
	minuit->mnexcm("SET STR", arguments, 1, errorFlag);

	minuit->mnexcm("SET NOGradient", arguments, 0, errorFlag);

	string IntOption("Interactive");
	if( StringProcessing::VectorContains( &Options, &IntOption ) != -1 )
	{
		minuit->mnexcm("SET INT", arguments, 0, errorFlag);
	}

	string SeekOption("SeekFirst");
	if( StringProcessing::VectorContains( &Options, &SeekOption ) != -1 )
	{
		arguments[0] = maxSteps;
		arguments[1] = 3;
		minuit->mnexcm("SEEk", arguments, 2, errorFlag);
	}

	string SimplexOption("SimplexFirst");
	if( StringProcessing::VectorContains( &Options, &SimplexOption ) != -1 )
	{
		//	First do Simplex
		minuit->mnexcm("SIMplex",arguments, 2, errorFlag);
	}

	string HesseFirstOption("HesseFirst");
	if( StringProcessing::VectorContains( &Options, &HesseFirstOption ) != -1 )
	{
		minuit->mnexcm("HESSE", arguments, 2, errorFlag);
	}

	arguments[0] = maxSteps;//MAXIMUM_MINIMISATION_STEPS
	arguments[1] = bestTolerance;//FINAL_GRADIENT_TOLERANCE;

	//	Now Do the minimisation
	minuit->mnexcm("MIGRAD", arguments, 2, errorFlag);

	string NoHesse("NoHesse");
	if( StringProcessing::VectorContains( &Options, &NoHesse ) == -1 )
	{
		//	Finally Now call HESSE to properly calculate the error matrix
		minuit->mnexcm("HESSE", arguments, 2, errorFlag);
	}

	string MinosOption("MinosErrors");
	if( StringProcessing::VectorContains( &Options, &MinosOption ) != -1 )
	{
		//	Call MINOS to calculate the non-parabolic errors.
		//	MINOS rather than assuming parabolic shape about the minimum, MINOS climbs
		//	out of the minimum on either side up until the UP value.
		minuit->mnexcm("MINOS", arguments, 2, errorFlag);
	}

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
	for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex)
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
				oldParameter->GetMinimum(), oldParameter->GetMaximum(), oldParameter->GetStepSize(),
				oldParameter->GetType(), oldParameter->GetUnit() );
	}

	// Get the error matrix and construct a vector containing the correlation coefficients
	int numParams = int(allNames.size());

	//	This Causes an outstanding compiler warning due to non standards compliance
	//	I'm taking the approach once used by Pete to solve annoying problems:
	//
	//	A beer to whoever can fix this!		rob.currie@ed.ac.uk

	Double_t matrix[numParams][numParams];

//	Double_t** matrix = new Double_t*[numParams];
//	for( short int i=0; i < numParams; ++i )
//		matrix[i] = new Double_t[numParams];
	gMinuit->mnemat(&matrix[0][0],numParams);
	
	vector<double> covarianceMatrix(unsigned(numParams*(numParams+1)/2));
	for (int row = 0; row < numParams; ++row)
	{
		for (int col = 0; col < numParams; ++col)
		{
			if(row > col) {covarianceMatrix[unsigned(col+row*(row+1)/2)] = matrix[row][col];}
			else {covarianceMatrix[unsigned(row+col*(col+1)/2)] = matrix[row][col];}
		}
	}

//	for( short int i=0; i< numParams; ++i )
//		delete[] matrix[i];
//	delete[] matrix;

	//Make the contour plots
	vector< FunctionContour* > allContours;
	for (unsigned int plotIndex = 0; plotIndex < contours.size(); ++plotIndex )
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
			double* xCoordinates1 = new double[unsigned(numberOfPoints)];
			double* yCoordinates1 = new double[unsigned(numberOfPoints)];
			double* xCoordinates2 = new double[unsigned(numberOfPoints)];
			double* yCoordinates2 = new double[unsigned(numberOfPoints)];

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

	fitResult = new FitResult( minimumValue, fittedParameters, fitStatus, function->GetPhysicsBottle(), covarianceMatrix, allContours );
}

//The function to pass to Minuit
void MinuitWrapper::Function( int & npar, double * grad, double & fval, double * xval, int iflag )
{
	(void) npar;
	(void) grad;
	(void) fval;
	(void) iflag;
	//cout << npar << "\t" << iflag << endl;
	//	CAUTION CAUTION CAUTION
	//	IF YOU EVER USE THIS, IT IS 
	//		!!!UNBLINDED!!!
	//
	//for( int i=0; i< npar; ++i )
	//{
	//	cout << grad[i] << "\t" << xval[i] << endl;
	//}

	ParameterSet * temporaryParameters = function->GetParameterSet();
	
	if ( temporaryParameters->SetPhysicsParameters(xval) )
	{
		function->SetParameterSet(temporaryParameters);
		fval = function->Evaluate();
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
