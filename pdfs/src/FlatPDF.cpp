
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
	(void) input;
	return 1.;
}

double FlatPDF::Normalisation( PhaseSpaceBoundary* range )
{
	IConstraint* x_const = range->GetConstraint( xName );

	double min = x_const->GetMinimum();
	double max = x_const->GetMaximum();

	double num_range = fabs(max-min);

	return (1./num_range);
}

