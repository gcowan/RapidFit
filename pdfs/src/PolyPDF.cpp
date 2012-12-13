
#include "PolyPDF.h"
#include "Mathematics.h"
#include <string>
#include <vector>
#include <cmath>

using namespace::std;

PDF_CREATOR( PolyPDF );

PolyPDF::PolyPDF( PDFConfigurator* config ) :
	xName( config->getName("x") ),
	parameterNames(), parameterValues(), order(0)
{
	string orderString = config->getConfigurationValue( "Order" );
	if( !orderString.empty() )
	{
		order = (unsigned)atoi( orderString.c_str() );
	}
	else
	{
		order = 0;
	}

	vector<string> parameterPrimitives;

	for( unsigned int i=0; i<= order; ++i )
	{
		TString polyParam = "Poly_"; polyParam += i;
		parameterNames.push_back( string( polyParam.Data() ) );
		parameterPrimitives.push_back( string( polyParam.Data() ) );
		parameterValues.push_back(0.);
	}

	allParameters = ParameterSet( parameterPrimitives );
}

PolyPDF::~PolyPDF()
{}

bool PolyPDF::SetPhysicsParameters( ParameterSet* input )
{
	allParameters = *input;

	for( unsigned int i=0; i<= order; ++i )
	{
		parameterValues[i] = allParameters.GetPhysicsParameter( parameterNames[i] )->GetValue();
		//cout << string(parameterNames[i]) << ":\t" << parameterValues[i] << endl;
	}

	return true;
}

double PolyPDF::Evaluate( DataPoint* input )
{
	double total = parameterValues[0];

	double x = input->GetObservable( xName )->GetValue();

	for( unsigned int i=1; i<= order; ++i )
	{
		total += parameterValues[i] * pow( x, (double)i );
	}

	if( total <= 1E-99 ) return 1E-99;

	return total;
}

double PolyPDF::Normalisation( PhaseSpaceBoundary* range )
{
	IConstraint* x_const = range->GetConstraint( xName );

	double min = x_const->GetMinimum();
	double max = x_const->GetMaximum();

	double total=0.;

	double const_val = (max - min) * parameterValues[0];

	total+=const_val;

	for( unsigned int i=1; i< order; ++i )
	{
		total += (parameterValues[i] / (double)i) * ( pow( max, (double)(i+1) ) - pow( min, (double)(i+1) ) );
	}

	return total;
}

