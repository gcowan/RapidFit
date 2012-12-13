
#include "FlatPDF.h"
#include "Mathematics.h"
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

PDF_CREATOR( FlatPDF );

FlatPDF::FlatPDF( PDFConfigurator* config ) :
	xName( config->getName("x") )
{
}

FlatPDF::~FlatPDF()
{}

bool FlatPDF::SetPhysicsParameters( ParameterSet* input )
{
	(void) input;
	return true;
}

double FlatPDF::Evaluate( DataPoint* input )
{
	PhaseSpaceBoundary* range = input->GetPhaseSpaceBoundary();

	IConstraint* x_const = range->GetConstraint( xName );

	double min = x_const->GetMinimum();
	double max = x_const->GetMaximum();

	double num_range = fabs(max-min);
	return 1./num_range;
}

double FlatPDF::Normalisation( PhaseSpaceBoundary* range )
{
	return 1.;
	IConstraint* x_const = range->GetConstraint( xName );

	double min = x_const->GetMinimum();
	double max = x_const->GetMaximum();

	double num_range = fabs(max-min);
	return num_range;
}

