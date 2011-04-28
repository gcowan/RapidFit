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
#include <cmath>

const double MAXIMUM_MINIMISATION_STEPS = 100000.0;//800.0;
const double FINAL_GRADIENT_TOLERANCE = 0.001;//0.001;
//const double STEP_SIZE = 0.001;
FitFunction * MinuitWrapper::function = 0;

//Default constructor
MinuitWrapper::MinuitWrapper(): print_verbosity( 0 )
{
	minuit = new TMinuit( 1 );
}

//Constructor with correct argument
MinuitWrapper::MinuitWrapper( int NumberParameters, int output_level ): print_verbosity( output_level )
{
	minuit = new TMinuit( NumberParameters );
}

//Destructor
MinuitWrapper::~MinuitWrapper()
{
	delete minuit;
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
	minuit->SetFCN( Function );

	//Store the parameters
	ParameterSet * newParameters = NewFunction->GetParameterSet();
	vector<string> allNames = newParameters->GetAllNames();
	for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex)
	{
		PhysicsParameter * newParameter = newParameters->GetPhysicsParameter( allNames[nameIndex] );

		double STEP_SIZE= 0.01;
		if(!(newParameter->GetMaximum() == newParameter->GetMinimum())){
		STEP_SIZE = fabs((newParameter->GetMaximum() - newParameter->GetMinimum()))/10000.0;
		}
	

		if( allNames[nameIndex] == "gamma" )		STEP_SIZE = 0.01;
		else if( allNames[nameIndex] == "deltaGamma" )	STEP_SIZE = 0.01;
		else if( allNames[nameIndex] == "Aperp_sq" )	STEP_SIZE = 0.01;
		else if( allNames[nameIndex] == "Azero_sq" )	STEP_SIZE = 0.01;
		else if( allNames[nameIndex] == "delta_para" )	STEP_SIZE = 0.1;
		else if( allNames[nameIndex] == "delta_perp" )	STEP_SIZE = 0.1;
		else if( allNames[nameIndex] == "alphaM_pr" )	STEP_SIZE = 0.0001;
		else STEP_SIZE = 0.001;
	
		cout << allNames[nameIndex] << "\t"<<nameIndex << "\t"<<STEP_SIZE <<endl;

		
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

	//	minuit->mnexcm("SOMECOMMAND",&SOMEARGUMENTS,NUMBEROFARGUMENTS,ERRORFLAG);

	//	Set the error analysis
	arguments[0] = NewFunction->UpErrorValue(1);
	minuit->mnexcm("SET ERR", arguments, 1, errorFlag);

	//	Set Migrad Strategy 2
	arguments[0] = 1;
	minuit->mnexcm("SET STR", arguments, 1, errorFlag);

        arguments[0] = MAXIMUM_MINIMISATION_STEPS;
        arguments[1] = FINAL_GRADIENT_TOLERANCE;

	//	First do Simplex
	//minuit->mnexcm("SIMplex",arguments, 2, errorFlag);

	//	Now Do the minimisation
	minuit->mnexcm("MIGRAD", arguments, 2, errorFlag);

	//	Finally Now call HESSE to properly calculate the error matrix
	minuit->mnexcm("HESSE", arguments, 2, errorFlag);

	//	Call MINOS to calculate the non-parabolic errors.
	//	MINOS rather than assuming parabolic shape about the minimum, MINOS climbs
	//	out of the minimum on either side up until the UP value.
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
				-numeric_limits<double>::max(), numeric_limits<double>::max(),
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
	vector<double> covarianceMatrix(numParams*(numParams+1)/2);
	for (int row = 0; row < numParams; ++row)
	{
		for (int col = 0; col < numParams; ++col)
		{
			if(row > col) {covarianceMatrix[col+row*(row+1)/2] = matrix[row][col];}
			else {covarianceMatrix[row+col*(col+1)/2] = matrix[row][col];}
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
			double* xCoordinates1 = new double[numberOfPoints];
			double* yCoordinates1 = new double[numberOfPoints];
			double* xCoordinates2 = new double[numberOfPoints];
			double* yCoordinates2 = new double[numberOfPoints];

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
	int null_i = npar; null_i = 0;
	double null_d = *grad; null_d = 0;
	int null_i2 = iflag; null_i2 = 0;
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
