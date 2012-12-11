
#include "BoxPDF.h"
#include "Mathematics.h"
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

PDF_CREATOR( BoxPDF );

BoxPDF::BoxPDF( PDFConfigurator* config ) :
	xName( config->getName("x") ),
	minName( config->getName("minimum") ),
	maxName( config->getName("maximum") ),
	minimum(0.), maximum(0.), range(0.), norm(0.)
{
	this->makePrototypes();
}

void BoxPDF::makePrototypes()
{
	allObservables.push_back( xName );
	vector<string> parameterNames;
	parameterNames.push_back( minName ); parameterNames.push_back( maxName );
	allParameters = ParameterSet( parameterNames );
}

BoxPDF::~BoxPDF()
{
}

bool BoxPDF::SetPhysicsParameters( ParameterSet* input )
{
	allParameters = *input;

	minimum = input->GetPhysicsParameter( minName )->GetValue();
	maximum = input->GetPhysicsParameter( maxName )->GetValue();

	range = fabs( maximum - minimum );

	norm = 1./range;

	return true;
}

double BoxPDF::Evaluate( DataPoint* input )
{
	double xVal = input->GetObservable( xName )->GetValue();

	double zero = 1E-99;

	if( xVal > minimum )
	{
		if( xVal < maximum )
		{
			return norm;
		}
		else
		{
			return zero;
		}
	}
	else
	{
		return zero;
	}

	return zero;
}

double BoxPDF::Normalisation( PhaseSpaceBoundary* Obsrange )
{
	(void) Obsrange;
	return 1.;
}

