/**
  @class TimeAccRes

  A class for holding a sliced propertime acceptance

  @author Dianne Ferguson
  @data 2013-10-17
  */


#include "TimeAccRes.h"
#include "StringProcessing.h"
#include "AcceptanceSlice.h"
#include "ClassLookUp.h"

#include <stdio.h>
#include <vector>
#include <string>

using namespace::std;

//............................................
// Constructor 
TimeAccRes::TimeAccRes( PDFConfigurator* configurator, bool quiet ) :
	resolutionModel(NULL), timeAcc(NULL), _config( new PDFConfigurator( *configurator ) )
{
	this->ConfigTimeRes( configurator, quiet );
	this->ConfigTimeAcc( configurator, quiet );
}

void TimeAccRes::ConfigTimeRes( PDFConfigurator* configurator, bool quiet )
{
	if( resolutionModel != NULL ) delete resolutionModel;
	string resolutionModelName = configurator->GetResolutionModel();
	resolutionModel = ClassLookUp::LookUpResName( resolutionModelName, configurator, quiet );
}

void TimeAccRes::ConfigTimeAcc( PDFConfigurator* configurator, bool quiet )
{
	if( timeAcc != NULL ) delete timeAcc;

	//.............
	//..............................
	// Configure to use time acceptance machinery
	bool _useTimeAcceptance = configurator->isTrue( "UseTimeAcceptance" ) ;
	if( _useTimeAcceptance )
	{
		if( configurator->hasConfigurationValue( "TimeAcceptanceType", "Upper" ) )
		{
			timeAcc = new SlicedAcceptance( 0., 14.0, 0.00826, quiet) ;
			//timeAcc = new SlicedAcceptance( 0., 14.0, 0.00, quiet) ;
			if( !quiet ) cout << "TimeAccRes:: Constructing timeAcc: Upper time acceptance beta=0.00826 [0 < t < 14] " << endl ;
		}
		else if( configurator->getConfigurationValue( "TimeAcceptanceFile" ) != "" )
		{
			timeAcc = new SlicedAcceptance( "File" , configurator->getConfigurationValue( "TimeAcceptanceFile" ), quiet ) ;
			if( !quiet ) cout << "TimeAccRes:: Constructing timeAcc: using file: " << configurator->getConfigurationValue( "TimeAcceptanceFile" ) << endl ;
		}
	}

	if( timeAcc == NULL )
	{
		timeAcc = new SlicedAcceptance( 0., 20., quiet ) ;
		if( !quiet ) cout << "TimeAccRes:: Constructing timeAcc: DEFAULT FLAT [0 < t < 20]  " << endl ;
	}
}

TimeAccRes::~TimeAccRes()
{
	if( timeAcc != NULL ) delete timeAcc;
	if( resolutionModel != NULL ) delete resolutionModel;
}

TimeAccRes::TimeAccRes( const TimeAccRes& input ) : resolutionModel(NULL), timeAcc(NULL), _config(NULL)
{
	if( input._config != NULL )
	{
		_config = new PDFConfigurator( *(input._config) );
		this->ConfigTimeRes( _config, true );
		this->ConfigTimeAcc( _config, true );
	}
}

//..........................
//This method allows the instance to add the parameters it needs to the list
void TimeAccRes::addParameters( vector<string> & parameterNames )
{
	resolutionModel->addParameters( parameterNames );
	return;
}

//..........................
//To take the current value of a parameter into the instance
void TimeAccRes::setParameters( ParameterSet & parameters )
{
	resolutionModel->setParameters( parameters );
	return;
}

bool TimeAccRes::CacheValid() const
{
	return resolutionModel->CacheValid();
}

//..........................
//This method allows the instance to add the specific observables it needs to the list
void TimeAccRes::addObservables( vector<string> & observableNames )
{
	resolutionModel->addObservables( observableNames );
	return;
}

//..........................
//To take the current value of an obserable into the instance
void TimeAccRes::setObservables( DataPoint * measurement )
{
	resolutionModel->setObservables( measurement );
	return;
}

//..........................
//To take the current value of an obserable into the instance
bool TimeAccRes::isPerEvent( ) {  return resolutionModel->isPerEvent(); }

double TimeAccRes::Exp( double time, double gamma )
{
	return resolutionModel->Exp( time, gamma ) * timeAcc->getValue( time );
}

double TimeAccRes::ExpInt( double tlow, double thigh, double gamma )
{
	double returnable_ExpInt=0.;

	double tlo_boundary = tlow;
	double thi_boundary = thigh;

	double thisExpInt = 0.;

	double tlo=0.;
	double thi=0.;

	for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
	{
		AcceptanceSlice* thisSlice = timeAcc->getSlice(islice);

		const double slice_lo = thisSlice->tlow();
		const double slice_hi = thisSlice->thigh();

		tlo = tlo_boundary > slice_lo ? tlo_boundary : slice_lo;
		thi = thi_boundary < slice_hi ? thi_boundary : slice_hi;
		if( thi > tlo )
		{
			thisExpInt = resolutionModel->ExpInt( tlo, thi, gamma );
			returnable_ExpInt += thisExpInt * thisSlice->height();
		}
	}

	return returnable_ExpInt;
}

double TimeAccRes::ExpSin( double time, double gamma, double dms )
{
	return resolutionModel->ExpSin( time, gamma, dms ) * timeAcc->getValue( time );
}

double TimeAccRes::ExpSinInt( double tlow, double thigh, double gamma, double dms )
{
	double returnable_ExpSinInt=0.;

	double tlo_boundary = tlow;
	double thi_boundary = thigh;

	double thisExpSinInt = 0.;

	double tlo=0.;
	double thi=0.;

	for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
	{
		AcceptanceSlice* thisSlice = timeAcc->getSlice(islice);

		const double slice_lo = thisSlice->tlow();
		const double slice_hi = thisSlice->thigh();

		tlo = tlo_boundary > slice_lo ? tlo_boundary : slice_lo;
		thi = thi_boundary < slice_hi ? thi_boundary : slice_hi;
		if( thi > tlo )
		{
			thisExpSinInt = resolutionModel->ExpSinInt( tlo, thi, gamma, dms );
			returnable_ExpSinInt += thisExpSinInt * thisSlice->height();
		}
	}

	return returnable_ExpSinInt;
}

double TimeAccRes::ExpCos( double time, double gamma, double dms )
{
	return resolutionModel->ExpCos( time, gamma, dms ) * timeAcc->getValue( time );
}

double TimeAccRes::ExpCosInt( double tlow, double thigh, double gamma, double dms )
{
	double returnable_ExpCosInt=0.;

	double tlo_boundary = tlow;
	double thi_boundary = thigh;

	double thisExpCosInt = 0.;

	double tlo=0.;
	double thi=0.;

	for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
	{
		AcceptanceSlice* thisSlice = timeAcc->getSlice(islice);

		const double slice_lo = thisSlice->tlow();
		const double slice_hi = thisSlice->thigh();

		tlo = tlo_boundary > slice_lo ? tlo_boundary : slice_lo;
		thi = thi_boundary < slice_hi ? thi_boundary : slice_hi;
		if( thi > tlo )
		{
			thisExpCosInt = resolutionModel->ExpSinInt( tlo, thi, gamma, dms );
			returnable_ExpCosInt += thisExpCosInt * thisSlice->height();
		}
	}

	return returnable_ExpCosInt;
}


