/*!
 * @class MinuitWrapper
 *
 * A wrapper to integrate Minuit with RapidFit
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

//	ROOT Headers
#include "Rtypes.h"
#include "TVector.h"
#include "TMatrixDSym.h"
#include "TMinuit.h"
#include "TMatrixTUtils.h"
//	RapidFit Headers
#include "MinuitWrapper.h"
#include "StringProcessing.h"
#include "FunctionContour.h"
#include "PhysicsBottle.h"
#include "CorrectedCovariance.h"
//	System Headers
#include <iostream>
#include <limits>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <float.h>
#include <iomanip>
#include <exception>

using namespace::std;

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

//const double MAXIMUM_MINIMISATION_STEPS = 100000.0;//800.0;
//const double FINAL_GRADIENT_TOLERANCE = 0.001;//0.001;
//const double STEP_SIZE = 0.001;

//	Required to act as a 'global' pointer to the RapidFit fitfunction member of the MinuitWrapper Class
IFitFunction * MinuitWrapper::function = 0;

//Default constructor
//MinuitWrapper::MinuitWrapper(): minuit(NULL), fitResult(NULL), contours(), print_verbosity( 0 ), maxSteps(), bestTolerance(), Options(), Quality(), debug(new DebugClass(false) )
//{
//	minuit = new TMinuit( 1 );
//}

//Constructor with correct argument
MinuitWrapper::MinuitWrapper( int NumberParameters, int output_level ) :
	minuit(NULL), fitResult(NULL), contours(), print_verbosity( output_level ), maxSteps(), bestTolerance(), Options(), Quality(), debug(new DebugClass(false) ), nSigma(1)
{
	minuit = new TMinuit( NumberParameters );
}

//Destructor
MinuitWrapper::~MinuitWrapper()
{
	//cout << "hello from MinuitWrapper destructor" << endl;
	delete minuit;
	if( debug != NULL ) delete debug;
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
	double* arg = new double[1];
	arg[0] = output_level;
	Int_t err;
	//cout << "Changing Output Level to: " << output_level << endl;
	if( output_level < 0 )
	{
		arg[0]=-1;
		minuit->mnexcm("SET PRINT", arg, 1, err);
		minuit->mnexcm("SET NOWarnings", arg, 1, err );
		err = minuit->Command("SET OUTputfile /dev/zero");
	} else {
		minuit->mnexcm("SET PRINT", arg, 1, err);
	}
	delete[] arg;
	//	minuit->SetPrintLevel( print_verbosity );
}

void MinuitWrapper::SetupFit( IFitFunction* NewFunction )
{
	function = NewFunction;
	int errorFlag = 0;
	double* arguments = new double[2];// = {0.0, 0.0};

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

		if( newParameter->GetType() != "Fixed" )
		{
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
		}

		//Fix the parameter if required
		//if ( newParameter->GetType() == "Fixed" )
		//{
		//	minuit->FixParameter( nameIndex );
		//}
		//else
		//{
		//	Released on Construction
		//minuit->Release( nameIndex );
		//}
	}

	//      Syntax for Minuit Commands through the TMinuit Class:
	//
	//      minuit->mnexcm("SOMECOMMAND",&SOMEARGUMENTS,NUMBEROFARGUMENTS,ERRORFLAG);

	//      Set the error analysis
	arguments[0] = NewFunction->UpErrorValue( nSigma );
	minuit->mnexcm("SET ERR", arguments, 1, errorFlag);

	//      Set Migrad Strategy 2
	arguments[0] = Quality;//1;
	minuit->mnexcm("SET STR", arguments, 1, errorFlag);

	minuit->mnexcm("SET NOGradient", arguments, 0, errorFlag);

	delete[] arguments;
}

IFitFunction* MinuitWrapper::GetFitFunction()
{
	return function;
}

void MinuitWrapper::FixParameters( vector<double> fix_values, vector<string> ParameterNames )
{
	vector<string> allNames = function->GetParameterSet()->GetAllNames();

	for( unsigned int i=0; i< fix_values.size(); ++i )
	{
		int errorFlag=0;
		int nameIndex = StringProcessing::VectorContains( &allNames, &(ParameterNames[i]) );
		minuit->mnparm(nameIndex, allNames[(unsigned)nameIndex], fix_values[i], 0.001, 0.0, 0.0, errorFlag);
		minuit->FixParameter( nameIndex );
	}
}

void MinuitWrapper::CallHesse()
{
	double* arguments = new double[2];
	int errorFlag=0;
	arguments[0] = maxSteps;//MAXIMUM_MINIMISATION_STEPS
	arguments[1] = bestTolerance;//FINAL_GRADIENT_TOLERANCE;
	//      Now Do the minimisation
	//minuit->mnexcm("MIGRAD", arguments, 2, errorFlag);

	//      Now Do the minimisation
	//              minuit->mnexcm("MIGRAD", arguments, 2, errorFlag);
	arguments[0] = maxSteps;//MAXIMUM_MINIMISATION_STEPS
	arguments[1] = bestTolerance;//FINAL_GRADIENT_TOLERANCE
	minuit->mnexcm("HESSE", arguments, 2, errorFlag);

	delete[] arguments;
}

//Use Migrad to minimise the given function
void MinuitWrapper::Minimise()
{
	currentMinuitInstance = minuit;
	int errorFlag = 0;
	double* arguments = new double[2];// = {0.0, 0.0};
	vector<string> allNames = function->GetParameterSet()->GetAllNames();
	ParameterSet * newParameters = function->GetParameterSet();

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
		minuit->mnhes1();
		//minuit->mnexcm("HESSE", arguments, 2, errorFlag);
	}

	string SetEPS("EPS3");
	if( StringProcessing::VectorContains( &Options, &SetEPS ) != -1 )
	{
		arguments[0] = 1E-5;
		minuit->mnexcm("SET EPS", arguments, 1, errorFlag);
	}

	arguments[0] = maxSteps;//MAXIMUM_MINIMISATION_STEPS
	arguments[1] = bestTolerance;//FINAL_GRADIENT_TOLERANCE;

        time_t timeNow;
        time(&timeNow);

	//	Now Do the minimisation
	string HesseOnly("HesseOnly");
	if( StringProcessing::VectorContains( &Options, &SetEPS ) == -1 )
	{
		minuit->mnexcm("MIGRAD", arguments, 2, errorFlag);
	}
	//Apply an improve step
	//arguments[0] = maxSteps;
	//minuit->mnexcm("IMPROVE", arguments, 2, errorFlag);

        //Get the fit status
        double minimumValue = 0.0;
        double fedm = 0.0;
        double errdef = 0.0;
        int variableParameters = 0;
        int parameterNumber = 0;
        int fitStatus = 0;
        minuit->mnstat( minimumValue, fedm, errdef, variableParameters, parameterNumber, fitStatus );

	//int preHesseStatus = fitStatus;

	time(&timeNow);
	cout << "\nFinal NLL: " << setprecision(15) << minimumValue << "\t\tStatus: " << fitStatus << "\t\t" << ctime(&timeNow) << endl << endl;

	string NoHesse("NoHesse");
	if( StringProcessing::VectorContains( &Options, &NoHesse ) == -1 )
	{
		//	Finally Now call HESSE to properly calculate the error matrix
		this->CallHesse();
	}

	string MinosOption("MinosErrors");
	if( StringProcessing::VectorContains( &Options, &MinosOption ) != -1 )
	{
		vector<string> allFreeNames = function->GetParameterSet()->GetAllFloatNames();
		for( unsigned int i=0; i< Options.size(); ++i )
		{
			string thisOption = Options[i];
			vector<string> thisList = StringProcessing::SplitString( thisOption, ':' );
			string Command = thisList[0];
			if( Command == "MinosFix" )
			{
				for( unsigned int j=1; j< thisList.size(); ++j )
				{
					string thisParam = thisList[j];
					int thisIndex = StringProcessing::VectorContains( &allFreeNames, &thisParam );
					if( thisIndex != -1 )
					{
						cout << "MINOSFIX: FIXING " << thisParam << ". This should only be used for uncorrolated and well defined parameters!" << endl;
						minuit->FixParameter( thisIndex );
					}
					else
					{
						cout << "MINOSFIX: CANNOT FIX PARAMETER " << thisParam << ". It is likely fixed or not in your ParameterSet." << endl;
					}
				}
			}
		}
		//	Call MINOS to calculate the non-parabolic errors.
		//	MINOS rather than assuming parabolic shape about the minimum, MINOS climbs
		//	out of the minimum on either side up until the UP value.
		minuit->mnexcm("MINOS", arguments, 2, errorFlag);
	}

	minuit->mnstat( minimumValue, fedm, errdef, variableParameters, parameterNumber, fitStatus );

	//Output time information
	time(&timeNow);
	cout << "Minuit finished: " << "\t\tStatus: " << fitStatus << "\t\t" << ctime(&timeNow) << endl;

	//Get the final fit status
	minuit->mnstat( minimumValue, fedm, errdef, variableParameters, parameterNumber, fitStatus );

	ResultParameterSet * fittedParameters = this->GetResultParameters( allNames, newParameters );

	RapidFitMatrix* covarianceMatrix = this->GetCovarianceMatrix();

	vector<FunctionContour*> allContours = this->ConstructContours( allNames, newParameters );

	fitResult = new FitResult( minimumValue, fittedParameters, fitStatus, function->GetPhysicsBottle(), covarianceMatrix, allContours );

	string testnewError="TestNewErrors";
	string rooFitError="RooFitErrors";
	if( ( StringProcessing::VectorContains( &Options, &testnewError) != -1 ) || ( StringProcessing::VectorContains( &Options, &rooFitError) != -1 ) )
	{
		//	This will also call ApplyCovarianceMatrix which will change the fitResult
		RapidFitMatrix* newMatrix = CorrectedCovariance::GetCorrectedCovarianceMatrix( this );
		(void) newMatrix;
	}

	delete[] arguments;
}

ResultParameterSet* MinuitWrapper::GetResultParameters( vector<string> allNames, ParameterSet* newParameters )
{
	//Get the fitted parameters
	ResultParameterSet * fittedParameters = new ResultParameterSet( allNames );
	string MinosOption("MinosErrors");
	bool usedMinos = StringProcessing::VectorContains( &Options, &MinosOption ) != -1;
	for( unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex )
	{
		string parameterName = allNames[nameIndex];
		double parameterValue = 0.0;
		double parameterError = 0.0;
		double xlolim = 0.0;
		double xuplim = 0.0;
		int iuint = 0;
		TString temporaryString;
		minuit->mnpout( nameIndex, temporaryString, parameterValue, parameterError, xlolim, xuplim, iuint );

		double err_plus, err_minus, err_para, gcc;

		if( usedMinos )
		{
			minuit->mnerrs( nameIndex, err_plus, err_minus, err_para, gcc );
		}

		PhysicsParameter * oldParameter = newParameters->GetPhysicsParameter( parameterName );
		fittedParameters->SetResultParameter( parameterName, parameterValue, oldParameter->GetOriginalValue(), parameterError,
				xlolim, xuplim, oldParameter->GetType(), oldParameter->GetUnit() );

		if( usedMinos )
		{
			ResultParameter* thisParam = fittedParameters->GetResultParameter( parameterName );
			thisParam->SetAssymErrors( err_plus, fabs(err_minus), err_para );
		}
	}
	return fittedParameters;
}

RapidFitMatrix* MinuitWrapper::GetCovarianceMatrix()
{
	unsigned int numParams = (unsigned)function->GetParameterSet()->GetAllFloatNames().size();
	/*!
	 * This section of code causes some minor warnings against old c++ standards, but the behaviour here, explicitly requires the matrix to be allocated this way
	 *
	 * It's probably possible to mimic the same behaviour with malloc and address allocation, but it's easier to let the compiler 'do it's thing'
	 */
	Double_t matrix[numParams][numParams];
	minuit->mnemat(&matrix[0][0],(int)numParams);

	//cout << "Matrix:" << endl;
	//cout << setprecision(3) << endl;
	TMatrixDSym* covMatrix = new TMatrixDSym( (int)numParams );

	for( unsigned int i=0; i< (unsigned)numParams; ++i )
	{
		for( unsigned int j=0; j< (unsigned)numParams; ++j )
		{
			//cout << "  " << matrix[i][j];
			(*covMatrix)((int)i,(int)j)=matrix[i][j];
		}
		//cout << endl;
	}
	//cout << endl;

	RapidFitMatrix* thisCovMatrix = new RapidFitMatrix();

	thisCovMatrix->thisMatrix = covMatrix;

	thisCovMatrix->theseParameters = function->GetParameterSet()->GetAllFloatNames();

	return thisCovMatrix;
}

vector<double> MinuitWrapper::oldGetCovarianceMatrix( int numParams )
{
	/*!
	 * This section of code causes some minor warnings against old c++ standards, but the behaviour here, explicitly requires the matrix to be allocated this way
	 *
	 * It's probably possible to mimic the same behaviour with malloc and address allocation, but it's easier to let the compiler 'do it's thing'
	 */
	Double_t matrix[numParams][numParams];
	minuit->mnemat(&matrix[0][0],numParams);

	vector<double> covarianceMatrix(unsigned(numParams*(numParams+1)/2));
	for (int row = 0; row < numParams; ++row)
	{
		for (int col = 0; col < numParams; ++col)
		{
			if(row > col) {covarianceMatrix[unsigned(col+row*(row+1)/2)] = matrix[row][col];}
			else {covarianceMatrix[unsigned(row+col*(col+1)/2)] = matrix[row][col];}
		}
	}

	return covarianceMatrix;
}


vector<FunctionContour*> MinuitWrapper::ConstructContours( vector<string> allNames, ParameterSet* newParameters )
{
	//Make the contour plots
	vector<FunctionContour*> allContours;
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
			int numberOfPoints = 20;
			int iErrf;
			double* xCoordinates1 = new double[unsigned(numberOfPoints+1)];
			double* yCoordinates1 = new double[unsigned(numberOfPoints+1)];
			double* xCoordinates2 = new double[unsigned(numberOfPoints+1)];
			double* yCoordinates2 = new double[unsigned(numberOfPoints+1)];
			double* xCoordinates3 = new double[unsigned(numberOfPoints+1)];
			double* yCoordinates3 = new double[unsigned(numberOfPoints+1)];


			//One sigma contour
			minuit->SetErrorDef( 0.5);
			minuit->mncont( xParameterIndex, yParameterIndex, numberOfPoints, xCoordinates1, yCoordinates1, iErrf );
			xCoordinates1[numberOfPoints] = xCoordinates1[0];
			xCoordinates1[numberOfPoints] = yCoordinates1[0];


			//Two sigma contour
			minuit->SetErrorDef( 2.0);
			minuit->mncont( xParameterIndex, yParameterIndex, numberOfPoints, xCoordinates2, yCoordinates2, iErrf );
			xCoordinates2[numberOfPoints] = xCoordinates2[0];
			xCoordinates2[numberOfPoints] = yCoordinates2[0];

			//Three sigma contour
			minuit->SetErrorDef( 4.5);
			minuit->mncont( xParameterIndex, yParameterIndex, numberOfPoints, xCoordinates3, yCoordinates3, iErrf );
			xCoordinates3[numberOfPoints] = xCoordinates3[0];
			xCoordinates3[numberOfPoints] = yCoordinates3[0];

			//Store the contours
			FunctionContour * newContour = new FunctionContour( contours[plotIndex].first, contours[plotIndex].second, 3 );
			newContour->SetPlot( 1, numberOfPoints, xCoordinates1, yCoordinates1 );
			newContour->SetPlot( 2, numberOfPoints, xCoordinates2, yCoordinates2 );
			newContour->SetPlot( 3, numberOfPoints, xCoordinates3, yCoordinates3 );
			allContours.push_back(newContour);

			delete[] xCoordinates1;
			delete[] xCoordinates2;
			delete[] xCoordinates3;
			delete[] yCoordinates1;
			delete[] yCoordinates2;
			delete[] yCoordinates3;
		}
	}
	return allContours;
}

//The function to pass to Minuit
void MinuitWrapper::Function( Int_t & npar, Double_t * grad, Double_t & fval, Double_t * xval, Int_t iflag )
{
	(void) npar;		//	Number of Free Parameters
	(void) grad;		//	Gradient in each Parameter
	(void) iflag;		//	flag for stuff... check this for using derivatives


	//	Very useful code for vry indepth debugging, I advise this be left in to help diagnose FitFunction behaviour!
	//
	//cout << npar << "\t" << iflag << endl;
	//
	//for( int i=0; i< npar; ++i )
	//{
	//	cout << grad[i] << "\t" << xval[i] << endl;
	//}

	//for( int i=0; i< npar; ++i )
	//{
	//	cout << xval[i] << "  " << ((double*)xval)[i] << endl;
	//}

	ParameterSet* test = function->GetParameterSet();

	try
	{
		test->UpdatePhysicsParameters( (double*)xval, npar );
		function->SetParameterSet( test );
		fval = function->Evaluate();
	}
	catch(...)
	{
		cerr << "Failed to set physics parameters or Eval Function" << endl;
		throw(-99999);
	}

	double min, edm, errdef;
	int mnpar, nparx, stat;
	currentMinuitInstance->mnstat( min, edm, errdef, mnpar, nparx, stat );

	cout << "Call: " << left << setw(5) << function->GetCallNum() << " NLL: " << setprecision(10) << fval << " minNLL: " << setprecision(10) << min << " EDM: " << setprecision(3) << setw(5) << edm;
	cout << " Status: " << setw(1) << stat << setw(20) << " " <<  "\r" << flush;
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

void MinuitWrapper::ApplyCovarianceMatrix( RapidFitMatrix* Input )
{
	cout << "Applying from Minuit:  ";
	for( unsigned int i=0; i< (unsigned)Input->theseParameters.size(); ++i )
	{
		cout << Input->theseParameters[i] << "\t";
	}
	cout << endl;
	fitResult->ApplyCovarianceMatrix( Input );
}

void MinuitWrapper::SetDebug( DebugClass* input_debug )
{
	if( debug != NULL ) delete debug;
	debug = new DebugClass( *input_debug );
}

void MinuitWrapper::SetNSigma( int input )
{
	nSigma = input;
}

