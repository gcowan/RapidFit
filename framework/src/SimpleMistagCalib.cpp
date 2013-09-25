
#include "ParameterSet.h"
#include "DataPoint.h"
#include "PDFConfigurator.h"
#include "DataPoint.h"
#include "ParameterSet.h"
#include "IMistagCalib.h"
#include "SimpleMistagCalib.h"

#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace::std;

SimpleMistagCalib::SimpleMistagCalib( PDFConfigurator* configurator ) : IMistagCalib(),
	_tag(), _mistag(), _mistagP0(), _mistagP1(), _mistagSetPoint(),
	_mistagDeltaP1(), _mistagDeltaP0(), _mistagDeltaSetPoint(),
	tagName( configurator->getName("tag") ),
	mistagName( configurator->getName("mistag") ),
	mistagP1Name( configurator->getName("mistagP1") ),
	mistagP0Name( configurator->getName("mistagP0") ),
	mistagSetPointName( configurator->getName("mistagSetPoint") ),
	mistagDeltaP1Name( configurator->getName("mistagDeltaP1") ),
	mistagDeltaP0Name( configurator->getName("mistagDeltaP0") ),
	mistagDeltaSetPointName( configurator->getName("mistagDeltaSetPoint") )
{
}

SimpleMistagCalib::~SimpleMistagCalib()
{
}

void SimpleMistagCalib::addParameters( vector<string>& parameterNames ) const
{
	parameterNames.push_back( mistagP1Name );
	parameterNames.push_back( mistagP0Name );
	parameterNames.push_back( mistagSetPointName );
	parameterNames.push_back( mistagDeltaP1Name );
	parameterNames.push_back( mistagDeltaP0Name );
	parameterNames.push_back( mistagDeltaSetPointName );
}

void SimpleMistagCalib::setParameters( const ParameterSet& parameters )
{
	_mistagP1 = parameters.GetPhysicsParameter( mistagP1Name )->GetValue();
	_mistagP0 = parameters.GetPhysicsParameter( mistagP0Name )->GetValue();
	_mistagSetPoint = parameters.GetPhysicsParameter( mistagSetPointName )->GetValue();
	_mistagDeltaP1 = parameters.GetPhysicsParameter( mistagDeltaP1Name )->GetValue();
	_mistagDeltaP0 = parameters.GetPhysicsParameter( mistagDeltaP0Name )->GetValue();
	_mistagDeltaSetPoint = parameters.GetPhysicsParameter( mistagDeltaSetPointName )->GetValue();
}


void SimpleMistagCalib::addObservables( vector<string>& observableNames ) const
{
	observableNames.push_back( tagName );
	observableNames.push_back( mistagName );
}

void SimpleMistagCalib::setObservables( const DataPoint* measurement )
{
	_tag = measurement->GetObservable( tagName )->GetValue();
	_mistag = measurement->GetObservable( mistagName )->GetValue();
}

double SimpleMistagCalib::q() const
{
	return _tag;
}

/*
double SimpleMistagCalib::mistag() const
{
	double returnValue = -1000.;

	if( (fabs(q()) < 0.5) || (fabs(q()) > 1.) )
	{
		returnValue = 0.5;
	}
	else if( (_mistag>=0.0) && (_mistag <= 0.5) )
	{
		//Normal case
		returnValue =  _mistagP0 + _mistagP1*(_mistag - _mistagSetPoint ) ;
		if( returnValue < 0 )  returnValue = 0;
		if( returnValue > 0.5) returnValue = 0.5;
	}
	else if( _mistag < 0.0 )
	{
		//PDF_THREAD_LOCK
		cout << "Bs2JpsiPhi_Signal_v6::mistag() : _mistag < 0 so set to 0 " << endl;
		//PDF_THREAD_UNLOCK
		returnValue = 0;
	}
	else if( _mistag > 0.5 )
	{
		//PDF_THREAD_LOCK
		cout << "Bs2JpsiPhi_Signal_v6::mistag() : _mistag > 0.5 so set to 0.5 "  << endl;
		//PDF_THREAD_UNLOCK
		returnValue = 0.5;
	}
	else
	{
		//PDF_THREAD_LOCK
		cout << "Bs2JpsiPhi_Signal_v6::mistag() : WARNING ******If you got here you dont know what you are doing  "  << endl;
		//PDF_THREAD_UNLOCK
		exit(1);
	}
	return returnValue;
}
*/

double SimpleMistagCalib::mistagBbar() const
{
	double returnValue = -1000.;

	if( fabs(q()) < 0.5 ) {
		returnValue = 0.5 ;
	}
	else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
		//Normal case
		returnValue =  _mistagP0-(_mistagDeltaP0*0.5) + (_mistagP1-(_mistagDeltaP1*0.5))*(_mistag - (_mistagSetPoint-(_mistagDeltaSetPoint*0.5)) );
		//if( true ) returnValue =   _mistagDeltaP0 + (_mistagDeltaP1)*(_mistag - (_mistagDeltaSetPoint) ) ;// to mock up independent P1/P0 for each tag
		if( returnValue < 0 )  returnValue = 0;
		if( returnValue > 0.5) returnValue = 0.5;
	}
	else if( _mistag < 0.0 ) {
		//PDF_THREAD_LOCK
		cout << "Bs2JpsiPhi_Signal_v6::mistagBbar() : _mistag < 0 so deltaMistag set to 0 also " << endl;
		//PDF_THREAD_UNLOCK
		returnValue = 0;
	}
	else if( _mistag > 0.5 ) {
		//PDF_THREAD_LOCK
		cout << "Bs2JpsiPhi_Signal_v6::mistagBbar() : _mistag > 0.5 so so deltaMistag set to 0.5 also "  << endl;
		//PDF_THREAD_UNLOCK
		returnValue = 0.5;
	}
	else {
		//PDF_THREAD_LOCK
		cout << "Bs2JpsiPhi_Signal_v6::mistagBbar() : WARNING ******If you got here you dont know what you are doing  "  << endl;
		//PDF_THREAD_UNLOCK
		exit(1);
	}
	return returnValue;
}

double SimpleMistagCalib::mistagB() const
{
	double returnValue = -1000.;

	if( (fabs(q()) < 0.5) || (fabs(q()) > 1.) ) {
		returnValue = 0.5;
	}
	else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
		//Normal case
		returnValue =  _mistagP0+(_mistagDeltaP0*0.5) + (_mistagP1+(_mistagDeltaP1*0.5))*(_mistag - (_mistagSetPoint+(_mistagDeltaSetPoint*0.5)) );
		//if( true ) returnValue =  _mistagP0 + (_mistagP1)*(_mistag - (_mistagSetPoint) ) ;  // to mock up independent P1/P0 for each tag
		if( returnValue < 0 )  returnValue = 0;
		if( returnValue > 0.5) returnValue = 0.5;
	}
	else if( _mistag < 0.0 ) {
		//PDF_THREAD_LOCK
		cout << "Bs2JpsiPhi_Signal_v6::mistagB() : _mistag < 0 so deltaMistag set to 0 also " << endl;
		//PDF_THREAD_UNLOCK
		returnValue = 0;
	}
	else if( _mistag > 0.5 ) {
		//PDF_THREAD_LOCK
		cout << "Bs2JpsiPhi_Signal_v6::mistagB() : _mistag > 0.5 so so deltaMistag set to 0.5 also "  << endl;
		//PDF_THREAD_UNLOCK
		returnValue = 0.5;
	}
	else {
		//PDF_THREAD_LOCK
		cout << "Bs2JpsiPhi_Signal_v6::mistagB() : WARNING ******If you got here you dont know what you are doing  "  << endl;
		//PDF_THREAD_UNLOCK
		exit(1);
	}
	return returnValue;
}

double SimpleMistagCalib::D1() const
{
	return 1.0 - this->q()*( this->mistagB() - this->mistagBbar() );
}

double SimpleMistagCalib::D2() const
{
	return this->q()*( 1.0 - this->mistagB() - this->mistagBbar() );
}

void SimpleMistagCalib::Print() const
{
	cout << endl;
	cout << "_mistagP0            " << _mistagP0 << "\t" << string(mistagP0Name) << endl;
	cout << "_mistagDeltaP0       " << _mistagDeltaP0 << "\t" << string(mistagDeltaP0Name) << endl;
	cout << "_mistagP1            " << _mistagP1 << "\t" << string(mistagP1Name) << endl;
	cout << "_mistagDeltaP1       " << _mistagDeltaP1 << "\t" << string(mistagDeltaP1Name) << endl;
	cout << "_mistag              " << _mistag << "\t" << string(mistagName) << endl;
	cout << "_mistagSetPoint      " << _mistagSetPoint << "\t" << string(mistagSetPointName) << endl;
	cout << "_mistagDeltaSetPoint " << _mistagDeltaSetPoint << "\t" << string(mistagDeltaSetPointName) << endl;
	cout << "this->mistagB()      " << this->mistagB() << endl;
	cout << "this->mistagBbar()   " << this->mistagBbar() << endl;
}


