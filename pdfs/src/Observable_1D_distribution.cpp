// $Id: Observable_1D_distribution.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Observable_1D_distribution Observable_1D_distribution.cpp
 *
 *  RapidFit PDF for Bs mass
 *
 *
 *  @author Pete
 *  @date 2011-07-30
 */

#include "Observable_1D_distribution.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"

#include <cmath>
#include <iostream>
#include <string>

PDF_CREATOR( Observable_1D_distribution );

//Constructor
Observable_1D_distribution::Observable_1D_distribution(PDFConfigurator* configurator) :
	// Physics parameters
	BasePDF(),
	wantedObservable( configurator->getObservableToModel().first ),
	givenDataSet( configurator->getObservableToModel().second ),
	wanted_poly( configurator->getFitFunc() ),
	use_function( false )
{
	MakePrototypes();
}

//Make the data point and parameter set
void Observable_1D_distribution::MakePrototypes()
{
	TString Name( wantedObservable ); Name.Append("_Name");
	TString Title( wantedObservable ); Title.Append("_Title");
	TTree* observable_tree = new TTree( Name, Title );
	double new_double = 0.;
	TBranch* new_branch = observable_tree->Branch( "Observable", &new_double, "Observable/D" );
	observable_tree->SetEntries( givenDataSet->GetDataNumber() );
	for( int i=0; i< givenDataSet->GetDataNumber(); ++i )
	{
		new_double = givenDataSet->GetDataPoint( i )->GetObservable( wantedObservable )->GetValue();
		new_branch->Fill();
	}

	observable_tree->SetEstimate( givenDataSet->GetDataNumber() );
	observable_tree->Draw("Observable");
	normalised_histogram = observable_tree->GetHistogram();
	normalised_histogram->SetName("NormalisedHisto_PDF");
	normalised_histogram->SetTitle("NormalisedHisto_PDF");

	normalised_histogram->Scale( 1.0/normalised_histogram->Integral() );

	functional_form = new TF1("Observable_Func", wanted_poly.c_str() );

	if( !wanted_poly.empty() )
	{
		//	Fit a given function to this histogram
		TFitResult* result = new TFitResult( normalised_histogram->Fit( functional_form, "IMLV", "goff" ) );
	
		int result_num = result->Status();
	
		if( result_num == 3 )
		{
			use_function = true;
		}
	}
}

//Destructor
Observable_1D_distribution::~Observable_1D_distribution()
{
}


//Calculate the function value
double Observable_1D_distribution::Evaluate(DataPoint * measurement)
{
	//For this PDF this is effectively an accept reject mechanism to decide if the chosen observable is to be accepted or not

	double input = measurement->GetObservable( wantedObservable )->GetValue();

	double rand = this->GetRandomFunction()->Rndm();

	double accept_reject = -1.;

	if( use_function )
	{
		double lim = functional_form->Eval( input );
		if( rand >= lim ) {	accept_reject = 0.;	}
		else {			accept_reject = 1.;	}
	}
	else
	{
		double lim = normalised_histogram->GetBinContent( normalised_histogram->FindBin( input ) );
		if( rand >= lim ) {	accept_reject = 0.;	}
		else {			accept_reject = 1.;	}
	}

	return accept_reject;
}

// Normalisation
double Observable_1D_distribution::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void)boundary;
	return 1.0;
}

