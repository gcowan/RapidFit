/**
        @class JPsiPhiDataGenerator

        A data generator implementing preselection for JPsiPhi compatible events

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "JPsiPhiDataGenerator.h"

//Constructor with correct argument
JPsiPhiDataGenerator::JPsiPhiDataGenerator( PhaseSpaceBoundary * NewBoundary, IPDF * NewPDF ) : AcceptReject( NewBoundary, NewPDF )
{
	moreThanMaximum = 0.11;
}

//Destructor
JPsiPhiDataGenerator::~JPsiPhiDataGenerator()
{
}

//Preselect data points
bool JPsiPhiDataGenerator::Preselection( DataPoint * TestDataPoint, double TestValue )
{
	double time = TestDataPoint->GetObservable( "time" )->GetValue();
	if ( time > 2.0 && TestValue > 0.4 )
	{
		return false;
	}
	else if ( time > 4.0 && TestValue > 0.05 )
	{
		return false;
	}
	else
	{
		return true;
	}
}
