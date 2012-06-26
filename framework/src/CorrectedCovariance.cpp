
///	ROOT Headers
#include "TMatrixDSym.h"
#include "TMatrixD.h"
///	RapidFit Headers
#include "IMinimiser.h"
#include "CorrectedCovariance.h"
#include "ParameterSet.h"
#include "StringProcessing.h"
#include "RapidFitMatrix.h"
///	System Headers
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace::std;

void CorrectedCovariance::DumpMatrix( TMatrixD * Input )
{
	cout << setprecision(3) << endl;
	int range = Input->GetNcols();
	for( int i=0; i< range; ++i )
	{
		for( int j=0; j< range; ++j )
		{
			cout << "  " << ((*Input)(i,j));
		}
		cout << endl;
	}
	return;
}

void CorrectedCovariance::DumpMatrix( TMatrixDSym* Input )
{
	cout << setprecision(3) << endl;
	int range = Input->GetNcols();
	for( int i=0; i< range; ++i )
	{
		for( int j=0; j< range; ++j )
		{
			cout << "  " << ((*Input)(i,j));
		}
		cout << endl;
	}
	return;
}

RapidFitMatrix* CorrectedCovariance::GetReducedMatrix( IMinimiser* thisMinimiser, RapidFitMatrix* Covariance_Matrix )
{
	TMatrixDSym* Raw_Covariance_Matrix = Covariance_Matrix->thisMatrix;
	/*!
	 * We now need to correct the Matricies Here to only contain the free parameters from the fit.
	 */
	FitResult* fitResult = thisMinimiser->GetFitResult();

	//      Required Input
	vector<string> RawFloatedParameterNames = fitResult->GetResultParameterSet()->GetAllFloatNames();
	vector<string> ConstrainedParameters = vector<string>();//thisMinimiser->GetFitFunction()->ConstrainedParameter();
	vector<string> FloatedParameterNames;

	for( unsigned int i=0; i< RawFloatedParameterNames.size(); ++i )
	{
		if( StringProcessing::VectorContains( &ConstrainedParameters, &(RawFloatedParameterNames[i]) ) == -1 )
		{
			FloatedParameterNames.push_back( RawFloatedParameterNames[i] );
		}
	}

	unsigned int numFreeParameters = (unsigned)FloatedParameterNames.size();

	//      Fill this vector with ObservableRefs which stores the name and reference number of each Observable
	vector<ObservableRef*> FreeObservableReferences;
	for( unsigned int i=0; i< numFreeParameters; ++i )
	{
		ObservableRef* thisRef = new ObservableRef( FloatedParameterNames[i] );
		fitResult->GetResultParameterSet()->GetResultParameter( *thisRef );
		//cout << string(*thisRef) << "\t" << thisRef->GetIndex() << endl;
		FreeObservableReferences.push_back( thisRef );
	}
	//cout << endl;

	//      Construct and Fill new matricies with the contents of only the free columns in the covariant matrix
	TMatrixDSym* Free_Covariance_Matrix = new TMatrixDSym( (int)FreeObservableReferences.size() );

	for( unsigned int i=0; i< FreeObservableReferences.size(); ++i )
	{
		unsigned int index_i = (unsigned)FreeObservableReferences[i]->GetIndex();

		for( unsigned int j=0; j< FreeObservableReferences.size(); ++j )
		{
			unsigned int index_j = (unsigned)FreeObservableReferences[j]->GetIndex();

			(*Free_Covariance_Matrix)((int)i,(int)j) = (*Raw_Covariance_Matrix)( (int)index_i, (int)index_j );
		}
	}

	while( !FreeObservableReferences.empty() )
	{
		if( FreeObservableReferences.back() != NULL ) delete FreeObservableReferences.back();
		FreeObservableReferences.pop_back();
	}

	RapidFitMatrix* outputMatrix = new RapidFitMatrix();

	outputMatrix->thisMatrix = Free_Covariance_Matrix;

	outputMatrix->theseParameters = FloatedParameterNames;

	return outputMatrix;
}

RapidFitMatrix* CorrectedCovariance::GetCorrectedCovarianceMatrix( IMinimiser* thisMinimiser )
{
	if( thisMinimiser->GetFitResult() == NULL )
	{
		cerr << "This should be Called once the Minimiser has Correctly converged with a sensible Minima" << endl;
		return NULL;
	}
	if( thisMinimiser->GetFitResult()->GetFitStatus() != 3 )
	{
		cerr << "This should be Called once the Minimiser has Correctly converged with a sensible Minima" << endl;
		return NULL;
	}
	if( thisMinimiser->GetFitFunction()->GetWeightsWereUsed() == false )
	{
		//	Silently return the error matrix from the Minimiser as it Cannot be corrected
		return thisMinimiser->GetCovarianceMatrix();
	}


	//      The Matrix currently in the Miniser should be defined by Hesse
	RapidFitMatrix* Raw_Covariance_Matrix = thisMinimiser->GetCovarianceMatrix();

	//	Get the 2nd Error Matrix
	thisMinimiser->GetFitFunction()->SetUseWeightsSquared( true );
	thisMinimiser->CallHesse();
	RapidFitMatrix* Raw_Covariance_Matrix_Weights_Squared = thisMinimiser->GetCovarianceMatrix();
	thisMinimiser->GetFitFunction()->SetUseWeightsSquared( false );

	RapidFitMatrix* Free_Covariance_Matrix = CorrectedCovariance::GetReducedMatrix( thisMinimiser, Raw_Covariance_Matrix );
	RapidFitMatrix* Free_Covariance_Matrix_Weights_Squared = CorrectedCovariance::GetReducedMatrix( thisMinimiser, Raw_Covariance_Matrix_Weights_Squared );


	/*!
	 * The following method to calculate the contents of the corrected matrix is taken from:
	 * http://root.cern.ch/root/html/src/RooAbsPdf.cxx.html#fhMiHE  (09/05/2012)
	 * I take the matrix maths from the section starting at line 1128
	 */

	TMatrixDSym* C = Free_Covariance_Matrix_Weights_Squared->thisMatrix;
	TMatrixDSym* V = Free_Covariance_Matrix->thisMatrix;

	DumpMatrix( V );
	cout << endl;
	DumpMatrix( C );
	cout << endl;

	double det_C = 0.;
	C->Invert( &det_C );

	DumpMatrix( C );

	if( fabs(det_C)<1E-99 )
	{
		cerr << "ERROR!" << endl;
		cerr << "Cannot Invert Original Raw Error Matrix, This is a serious Issue, don't trust the returned Error Matrix!" << endl;
		delete Raw_Covariance_Matrix;
		delete Raw_Covariance_Matrix_Weights_Squared;
		delete Free_Covariance_Matrix;
		delete Free_Covariance_Matrix_Weights_Squared;
		return thisMinimiser->GetCovarianceMatrix();
	}

	TMatrixD* CV = new TMatrixD( *C, TMatrixD::kMult, *V );
	cout << endl;
	DumpMatrix( CV );

	// Calculate corrected covariance matrix = V C-1 V
	TMatrixD* VCV = new TMatrixD( *V, TMatrixD::kMult, TMatrixD( *C, TMatrixD::kMult, *V) );
	cout << endl;
	DumpMatrix(VCV);
	cout << endl;

	int numFreeParameters = (int)Free_Covariance_Matrix->theseParameters.size();

	TMatrixDSym* Correct_Free_Covariance_Matrix = new TMatrixDSym( numFreeParameters );

	for( int i=0 ; i<(int)numFreeParameters; ++i )
	{
		for( int j=i ; j<(int)numFreeParameters; ++j )
		{
			if( i==j )
			{
				(*Correct_Free_Covariance_Matrix)(i,j) = (*VCV)(i,j);
			}
			if( i!=j )
			{
				double deltaRel = ((*VCV)(i,j)-(*VCV)(j,i))/sqrt((*VCV)(i,i)*(*VCV)(j,j));

				if( fabs(deltaRel)>1e-2 )
				{
					cerr << "WARNING: Corrected covariance matrix is not (completely) symmetric: V[" << i << "," << j << "] = ";
					cerr  << (*VCV)(i,j) << " V[" << j << "," << i << "] = " << (*VCV)(j,i) << " explicitly restoring symmetry by inserting average value" << endl;
				}
				(*Correct_Free_Covariance_Matrix)(i,j) = ((*VCV)(i,j)+(*VCV)(j,i)) * 0.5;
			}
		}
	}

	RapidFitMatrix* outputMatrix = new RapidFitMatrix();

	outputMatrix->thisMatrix = Correct_Free_Covariance_Matrix;

	outputMatrix->theseParameters = Free_Covariance_Matrix->theseParameters;

	thisMinimiser->ApplyCovarianceMatrix( outputMatrix );

	delete VCV;
	delete Free_Covariance_Matrix;
	delete Free_Covariance_Matrix_Weights_Squared;
	//delete Correct_Free_Covariance_Matrix;
	delete Raw_Covariance_Matrix;
	delete Raw_Covariance_Matrix_Weights_Squared;
	return outputMatrix;
}

