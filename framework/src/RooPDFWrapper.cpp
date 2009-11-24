/**
        @class RooPDFWrapper

        The parent class for all wrappers for RooFit PDFs, implementing the basic methods.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "RooPDFWrapper.h"
#include <limits>

//Default constructor
RooPDFWrapper::RooPDFWrapper() : valid( false )
{
}

//Set up the PDF based on the provided name
bool RooPDFWrapper::SetUpWrapper( vector<string> PDFObservables, int NumberObservables, vector<string> PDFParameters, int NumberParameters )
{
	if ( PDFObservables.size() == NumberObservables )
	{
		vector<string>::iterator observable;
		for ( observable = PDFObservables.begin(); observable != PDFObservables.end(); observable++ )
		{
			//Create the RooRealVar
			RooRealVar * newVariable = new RooRealVar( observable->c_str(), observable->c_str(),
					0.0, -numeric_limits<double>::max(), numeric_limits<double>::max() );
			observables.push_back( newVariable );
		}
	}
	else
	{
		cerr << "Incorrect number of observables provided: " << PDFObservables.size() << " not " << NumberObservables << endl;
		return false;
	}

	if ( PDFParameters.size() == NumberParameters )
	{
		pdfParameterSet = new ParameterSet( PDFParameters );

		vector<string>::iterator parameterNames;
		for ( parameterNames = PDFParameters.begin(); parameterNames != PDFParameters.end(); parameterNames++ )
		{
			//Populate the ParameterSet
			pdfParameterSet->SetPhysicsParameter( *parameterNames, new PhysicsParameter() );

			//Create the RooRealVar
			RooRealVar * newVariable = new RooRealVar( parameterNames->c_str(), parameterNames->c_str(),
					0.0, -numeric_limits<double>::max(), numeric_limits<double>::max() );
			variables.push_back( newVariable );
		}
	}
	else
	{
		cerr << "Incorrect number of parameters provided: " << PDFParameters.size() << " not " << NumberParameters << endl;
		return false;
	}

	return true;
}

//Destructor
RooPDFWrapper::~RooPDFWrapper()
{
	delete wrappedPDF;
}

//Indicate whether the function has been set up correctly
bool RooPDFWrapper::IsValid()
{
	return valid;
}

//Set the function parameters
bool RooPDFWrapper::SetPhysicsParameters( ParameterSet * NewParameters )
{
	//Update the internal ParameterSet
	bool worked = pdfParameterSet->SetPhysicsParameters( NewParameters );

	if ( worked )
	{		
		//Pass the values to the wrapped PDF
		for ( int parameterIndex = 0; parameterIndex < prototypeParameterSet.size(); parameterIndex++ )
		{
			PhysicsParameter * newParameter = pdfParameterSet->GetPhysicsParameter( prototypeParameterSet[ parameterIndex ] );
			variables[ parameterIndex ]->setRange( newParameter->GetMinimum(), newParameter->GetMaximum() );
			variables[ parameterIndex ]->setVal( newParameter->GetValue() );
		}
	}
	else
	{
		cerr << "Parameter set does not include all parameters required by PDF" << endl;
	}

	return worked;
}

//Return the integral of the function over the given boundary
double RooPDFWrapper::Integral( DataPoint* NewDataPoint, PhaseSpaceBoundary * NewBoundary )
{
	RooArgSet observableSet;

	//Pass the boundary to the wrapped PDF
	for ( int boundIndex = 0; boundIndex < prototypeDataPoint.size(); boundIndex++ )
	{
		IConstraint * newBound = NewBoundary->GetConstraint( prototypeDataPoint[ boundIndex ] );

		if ( newBound->GetUnit() == "NameNotFoundError" )
		{
			cerr << "Required bound not found in provided phase space boundary: " << prototypeDataPoint[ boundIndex ] << endl;
			return 1.0;
		}
		else
		{
			//Update the observable range
			observables[ boundIndex ]->setRange( newBound->GetMinimum(), newBound->GetMaximum() );
			observableSet.add( *observables[ boundIndex ], false );
		}
	}

	return wrappedPDF->getNorm( observableSet );
}

//Return the function value at the given point
double RooPDFWrapper::Evaluate( DataPoint * NewDataPoint )
{
	//Pass the data point to the wrapped PDF
	for ( int observableIndex = 0; observableIndex < prototypeDataPoint.size(); observableIndex++ )
	{
		Observable * newObservable = NewDataPoint->GetObservable( prototypeDataPoint[ observableIndex ] );

		if ( newObservable->GetUnit() == "NameNotFoundError" )
		{
			cerr << "Required observable not found in provided data point: " << prototypeDataPoint[ observableIndex ] << endl;
			return 0.0;
		}
		else
		{
			observables[ observableIndex ]->setVal( newObservable->GetValue() );
		}
	}

	return wrappedPDF->getVal();
}

//Return a prototype data point
vector<string> RooPDFWrapper::GetPrototypeDataPoint()
{
	return prototypeDataPoint;
}

//Return a prototype set of physics parameters
vector<string> RooPDFWrapper::GetPrototypeParameterSet()
{
	return prototypeParameterSet;
}
