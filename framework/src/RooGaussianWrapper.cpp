/**
        @class RooGaussianWrapper

        A wrapper for the RooFit gaussian PDF

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "RooGaussian.h"
#include "RooGaussianWrapper.h"

//Default constructor
RooGaussianWrapper::RooGaussianWrapper()
{
	prototypeDataPoint.push_back( "xValue" );

	prototypeParameterSet.push_back( "mean" );
	prototypeParameterSet.push_back( "sigma" );

	if ( SetUpWrapper( prototypeDataPoint, 1, prototypeParameterSet, 2 ) )
	{
		wrappedPDF = new RooGaussian( "Gaussian", "RooFit gaussian PDF", *observables[0], *variables[0], *variables[1] );
		valid = true;
	}
	else
	{
		valid = false;
	}
}

//Constructor with correct arguments
RooGaussianWrapper::RooGaussianWrapper( vector<string> PDFObservables, vector<string> PDFParameters)
{
	prototypeDataPoint = PDFObservables;
	prototypeParameterSet = PDFParameters;

	if ( SetUpWrapper( PDFObservables, 1, PDFParameters, 2 ) )
	{
		wrappedPDF = new RooGaussian( "Gaussian", "RooFit gaussian PDF", *observables[0], *variables[0], *variables[1] );
		valid = true;
	}
	else
	{
		valid = false;
	}
}

//Destructor
RooGaussianWrapper::~RooGaussianWrapper()
{
}
