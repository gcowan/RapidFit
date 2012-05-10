
///	ROOT Headers
#include "TMatrixDSym.h"
#include "TMatrixD.h"
///	RapidFit Headers
#include "IMinimiser.h"
#include "CorrectedCovariance.h"
#include "ParameterSet.h"
#include "StringProcessing.h"
///	System Headers
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace::std;

TMatrixDSym* CorrectedCovariance::GetCorrectedCovarianceMatrix( IMinimiser* thisMinimiser )
{
	ParameterSet* FitParameterSet = thisMinimiser->GetFitFunction()->GetParameterSet();

	//	The Matrix currently in the Miniser should be defined by Hesse
	TMatrixDSym* Raw_Covariance_Matrix = thisMinimiser->GetCovarianceMatrix();

	if( thisMinimiser->GetFitFunction()->GetWeightsWereUsed() == false )
	{
		//	Silently return the error matrix from the Minimiser as it Cannot be corrected
		return Raw_Covariance_Matrix;
	}

	cout << endl << "RAW Covariant Matrix:" << endl;

	for( unsigned int i=0; i< (unsigned)FitParameterSet->GetAllNames().size(); ++i )
	{
		for( unsigned int j=0; j< (unsigned)FitParameterSet->GetAllNames().size(); ++j )
		{
			cout << "  " << (*Raw_Covariance_Matrix)(i,j);
		}
		cout << endl;
	}

	//	Get the 2nd Error Matrix
	thisMinimiser->GetFitFunction()->SetUseWeightsSquared( true );
	thisMinimiser->CallHesse();
	TMatrixDSym* Raw_Covariance_Matrix_Weights_Squared = thisMinimiser->GetCovarianceMatrix();
	thisMinimiser->GetFitFunction()->SetUseWeightsSquared( false );

	/*!
	 * We now need to correct the Matricies Here to only contain the free parameters from the fit.
	 */

	//	Required Input
	vector<string> FloatedParameterNames = FitParameterSet->GetAllFloatNames();
	unsigned int numFreeParameters = (unsigned)FloatedParameterNames.size();

	//	Fill this vector with ObservableRefs which stores the name and reference number of each Observable
	vector<ObservableRef*> FreeObservableReferences;
	cout << endl << "Free Parameters:" << endl;
	for( unsigned int i=0; i< numFreeParameters; ++i )
	{
		ObservableRef* thisRef = new ObservableRef( FloatedParameterNames[i] );
		FitParameterSet->GetPhysicsParameter( *thisRef );
		cout << string(*thisRef) << "\t" << thisRef->GetIndex() << endl;
		FreeObservableReferences.push_back( thisRef );
	}

	//	Construct and Fill new matricies with the contents of only the free columns in the covariant matrix
	TMatrixDSym* Free_Covariance_Matrix = new TMatrixDSym( (int)FreeObservableReferences.size() );
	TMatrixDSym* Free_Covariance_Matrix_Weights_Squared = new TMatrixDSym( (int)FreeObservableReferences.size() );

	cout << endl << "Free Parameter Reduced Covariant Matrix:" << endl;

	for( unsigned int i=0; i< FreeObservableReferences.size(); ++i )
	{
		unsigned int index_i = (unsigned)FreeObservableReferences[i]->GetIndex();

		for( unsigned int j=0; j< FreeObservableReferences.size(); ++j )
		{
			unsigned int index_j = (unsigned)FreeObservableReferences[j]->GetIndex();

			cout << "  " << (*Raw_Covariance_Matrix)( index_i, index_j );

			(*Free_Covariance_Matrix)(i,j) = (*Raw_Covariance_Matrix)( index_i, index_j );
			(*Free_Covariance_Matrix_Weights_Squared)(i,j) = (*Raw_Covariance_Matrix_Weights_Squared)( index_i, index_j );
		}
		cout << endl;
	}



	/*!
	 * The following method to calculate the contents of the corrected matrix is taken from:
	 * http://root.cern.ch/root/html/src/RooAbsPdf.cxx.html#fhMiHE  (09/05/2012)
	 * I take the matrix maths from the section starting at line 1128
	 */

	TMatrixDSym* V = Free_Covariance_Matrix_Weights_Squared;
	TMatrixDSym* C = Free_Covariance_Matrix;

	double det_C = 0.;
	C->Invert( &det_C );

	if( fabs(det_C)<1E-99 )
	{
		cerr << "ERROR!" << endl;
		cerr << "Cannot Invert Original Raw Error Matrix, This is a serious Issue, don't trust the returned Error Matrix!" << endl;
		delete Raw_Covariance_Matrix;
		delete Raw_Covariance_Matrix_Weights_Squared;
		return thisMinimiser->GetCovarianceMatrix();
	}

	// Calculate corrected covariance matrix = V C-1 V
	TMatrixD VCV( *V, TMatrixD::kMult, TMatrixD( *C, TMatrixD::kMult, *V) );

	TMatrixDSym* Correct_Free_Covariance_Matrix = new TMatrixDSym( numFreeParameters );

	for( int i=0 ; i<(int)numFreeParameters; ++i )
	{
		for( int j=i ; j<(int)numFreeParameters; ++j )
		{
			if( i==j )
			{
				(*Correct_Free_Covariance_Matrix)(i,j) = VCV(i,j);
			}
			if( i!=j )
			{
				double deltaRel = (VCV(i,j)-VCV(j,i))/sqrt(VCV(i,i)*VCV(j,j));

				if( fabs(deltaRel)>1e-3 )
				{
					cerr << "WARNING: Corrected covariance matrix is not (completely) symmetric: V[" << i << "," << j << "] = ";
					cerr  << VCV(i,j) << " V[" << j << "," << i << "] = " << VCV(j,i) << " explicitly restoring symmetry by inserting average value" << endl;
				}
				(*Correct_Free_Covariance_Matrix)(i,j) = (VCV(i,j)+VCV(j,i)) * 0.5;
			}
		}
	}

	vector<string> allNames =  FitParameterSet->GetAllNames();
	vector<string> allFreeNames = FitParameterSet->GetAllFloatNames();

	TMatrixDSym* Correct_Covariance_Matrix = new TMatrixDSym( (int)allNames.size() );

	for( unsigned int i=0; i< allNames.size(); ++i )
	{
		bool fixed_i = FitParameterSet->GetPhysicsParameter( allNames[i] )->GetType() == "Fixed";

		for( unsigned int j=0; j< allNames.size(); ++j )
		{
			bool fixed_j = FitParameterSet->GetPhysicsParameter( allNames[j] )->GetType() == "Fixed";

			if( fixed_i || fixed_j )
			{
				(*Correct_Covariance_Matrix)(i,j) = 0;
			}
			else
			{
				int index_i = StringProcessing::VectorContains( &allFreeNames, &(allNames[i]) );
				int index_j = StringProcessing::VectorContains( &allFreeNames, &(allNames[j]) );

				(*Correct_Covariance_Matrix)(i,j) = (*Correct_Free_Covariance_Matrix)( index_i, index_j );
			}
		}
	}

	delete Free_Covariance_Matrix;
	delete Free_Covariance_Matrix_Weights_Squared;
	delete Correct_Free_Covariance_Matrix;
	delete Raw_Covariance_Matrix;
	delete Raw_Covariance_Matrix_Weights_Squared;

	thisMinimiser->ApplyCovarianceMatrix( Correct_Covariance_Matrix );

	return Correct_Covariance_Matrix;
}

