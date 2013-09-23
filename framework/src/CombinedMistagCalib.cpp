
#include "ParameterSet.h"
#include "DataPoint.h"
#include "PDFConfigurator.h"
#include "DataPoint.h"
#include "ParameterSet.h"
#include "IMistagCalib.h"
#include "CombinedMistagCalib.h"
#include "ObservableRef.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>

using namespace::std;

CombinedMistagCalib::CombinedMistagCalib( PDFConfigurator* configurator ) : IMistagCalib(),
	_tagOS(), _tagSS(), _mistagOS(), _mistagSS(),
	_mistagP0_OS(), _mistagP1_OS(), _mistagSetPoint_OS(), _mistagDeltaP1_OS(), _mistagDeltaP0_OS(), _mistagDeltaSetPoint_OS(),
	_mistagP0_SS(), _mistagP1_SS(), _mistagSetPoint_SS(), _mistagDeltaP1_SS(), _mistagDeltaP0_SS(), _mistagDeltaSetPoint_SS(),
	_mistagP0_OSSS(), _mistagP1_OSSS(), _mistagSetPoint_OSSS(), _mistagDeltaP1_OSSS(), _mistagDeltaP0_OSSS(), _mistagDeltaSetPoint_OSSS(),
	_eta(), _combinedtag(),
	//	Observable Names
	tagOSName( configurator->getName("tagdecision_os") ),
	tagSSName( configurator->getName("tagdecision_ss") ),
	tagCombName( configurator->getName("tagdecision") ),
	mistagOSName( configurator->getName("tagomega_os") ),
	mistagSSName( configurator->getName("tagomega_ss") ),
	//	Physics Parameter Names
	mistagP1Name_OS( configurator->getName("mistagP1_OS") ),
	mistagP0Name_OS( configurator->getName("mistagP0_OS") ),
	mistagSetPointName_OS( configurator->getName("mistagSetPoint_OS") ),
	mistagDeltaP1Name_OS( configurator->getName("mistagDeltaP1_OS") ),
	mistagDeltaP0Name_OS( configurator->getName("mistagDeltaP0_OS") ),
	mistagDeltaSetPointName_OS( configurator->getName("mistagDeltaSetPoint_OS") ),
	mistagP1Name_SS( configurator->getName("mistagP1_SS") ),
	mistagP0Name_SS( configurator->getName("mistagP0_SS") ),
	mistagSetPointName_SS( configurator->getName("mistagSetPoint_SS") ),
	mistagDeltaP1Name_SS( configurator->getName("mistagDeltaP1_SS") ),
	mistagDeltaP0Name_SS( configurator->getName("mistagDeltaP0_SS") ),
	mistagDeltaSetPointName_SS( configurator->getName("mistagDeltaSetPoint_SS") ),
	mistagP1Name_OSSS( configurator->getName("mistagP1_OSSS") ),
	mistagP0Name_OSSS( configurator->getName("mistagP0_OSSS") ),
	mistagSetPointName_OSSS( configurator->getName("mistagSetPoint_OSSS") ),
	mistagDeltaP1Name_OSSS( configurator->getName("mistagDeltaP1_OSSS") ),
	mistagDeltaP0Name_OSSS( configurator->getName("mistagDeltaP0_OSSS") ),
	mistagDeltaSetPointName_OSSS( configurator->getName("mistagDeltaSetPoint_OSSS") )
{
}

CombinedMistagCalib::~CombinedMistagCalib()
{
}

void CombinedMistagCalib::addParameters( vector<string>& parameterNames ) const
{
	parameterNames.push_back( mistagP1Name_OS );
	parameterNames.push_back( mistagP0Name_OS );
	parameterNames.push_back( mistagSetPointName_OS );
	parameterNames.push_back( mistagDeltaP1Name_OS );
	parameterNames.push_back( mistagDeltaP0Name_OS );
	parameterNames.push_back( mistagDeltaSetPointName_OS );
	parameterNames.push_back( mistagP1Name_SS );
	parameterNames.push_back( mistagP0Name_SS );
	parameterNames.push_back( mistagSetPointName_SS );
	parameterNames.push_back( mistagDeltaP1Name_SS );
	parameterNames.push_back( mistagDeltaP0Name_SS );
	parameterNames.push_back( mistagDeltaSetPointName_SS );
	parameterNames.push_back( mistagP1Name_OSSS );
	parameterNames.push_back( mistagP0Name_OSSS );
	parameterNames.push_back( mistagSetPointName_OSSS );
	parameterNames.push_back( mistagDeltaP1Name_OSSS );
	parameterNames.push_back( mistagDeltaP0Name_OSSS );
	parameterNames.push_back( mistagDeltaSetPointName_OSSS );
}

void CombinedMistagCalib::setParameters( const ParameterSet& parameters )
{
	_mistagP1_OS = parameters.GetPhysicsParameter( mistagP1Name_OS )->GetValue();
	_mistagP0_OS = parameters.GetPhysicsParameter( mistagP0Name_OS )->GetValue();
	_mistagSetPoint_OS = parameters.GetPhysicsParameter( mistagSetPointName_OS )->GetValue();
	_mistagDeltaP1_OS = parameters.GetPhysicsParameter( mistagDeltaP1Name_OS )->GetValue();
	_mistagDeltaP0_OS = parameters.GetPhysicsParameter( mistagDeltaP0Name_OS )->GetValue();
	_mistagDeltaSetPoint_OS = parameters.GetPhysicsParameter( mistagDeltaSetPointName_OS )->GetValue();

	_mistagP1_SS = parameters.GetPhysicsParameter( mistagP1Name_SS )->GetValue();
	_mistagP0_SS = parameters.GetPhysicsParameter( mistagP0Name_SS )->GetValue();
	_mistagSetPoint_SS = parameters.GetPhysicsParameter( mistagSetPointName_SS )->GetValue();
	_mistagDeltaP1_SS = parameters.GetPhysicsParameter( mistagDeltaP1Name_SS )->GetValue();
	_mistagDeltaP0_SS = parameters.GetPhysicsParameter( mistagDeltaP0Name_SS )->GetValue();
	_mistagDeltaSetPoint_SS = parameters.GetPhysicsParameter( mistagDeltaSetPointName_SS )->GetValue();

	_mistagP1_OSSS = parameters.GetPhysicsParameter( mistagP1Name_OSSS )->GetValue();
	_mistagP0_OSSS = parameters.GetPhysicsParameter( mistagP0Name_OSSS )->GetValue();
	_mistagSetPoint_OSSS = parameters.GetPhysicsParameter( mistagSetPointName_OSSS )->GetValue();
	_mistagDeltaP1_OSSS = parameters.GetPhysicsParameter( mistagDeltaP1Name_OSSS )->GetValue();
	_mistagDeltaP0_OSSS = parameters.GetPhysicsParameter( mistagDeltaP0Name_OSSS )->GetValue();
	_mistagDeltaSetPoint_OSSS = parameters.GetPhysicsParameter( mistagDeltaSetPointName_OSSS )->GetValue();
}


void CombinedMistagCalib::addObservables( vector<string>& observableNames ) const
{
	observableNames.push_back( tagOSName );
	observableNames.push_back( mistagOSName );
	observableNames.push_back( tagSSName );
	observableNames.push_back( mistagSSName );
	observableNames.push_back( tagCombName );
}

void CombinedMistagCalib::setObservables( const DataPoint* measurement )
{
	_tagOS = (int) measurement->GetObservable( tagOSName )->GetValue();
	_mistagOS = measurement->GetObservable( mistagOSName )->GetValue();
	_tagSS = (int) measurement->GetObservable( tagSSName )->GetValue();
	_mistagSS = measurement->GetObservable( mistagSSName )->GetValue();

	_combinedtag = (int)measurement->GetObservable( tagCombName )->GetValue();

	if( !_tagOS && !_tagSS ) _eta = this->getEta();
}

double CombinedMistagCalib::getEta() const
{
	double thisEta = 0.;
	if( _tagOS == _tagSS )
	{
		if( _tagOS == -1 )
		{
			thisEta = ( this->mistagOSBbar() *this->mistagSSB() ) / ( this->mistagOSBbar()*this->mistagSSBbar() + (1.-this->mistagOSBbar())*(1.-this->mistagSSBbar()) );
		}
		else
		{
			thisEta = ( this->mistagOSB() * this->mistagSSB() ) / ( this->mistagOSB()*this->mistagSSB() + (1.-this->mistagOSB())*(1.-this->mistagSSB()) );
		}
	}
	else
	{
		if( _combinedtag == -1 )
		{
			if( _mistagSS > _mistagOS )
			{
				thisEta = ( this->mistagOSBbar() * ( 1. - this->mistagSSBbar() ) ) / ( this->mistagOSBbar() * ( 1. - this->mistagSSBbar() ) + ( 1. - this->mistagOSBbar() ) * this->mistagSSBbar() );
			}
			else
			{
				thisEta = ( ( 1. - this->mistagOSBbar() ) * this->mistagSSBbar() ) / ( ( 1. - this->mistagOSBbar() ) * this->mistagSSBbar() + this->mistagOSBbar() * ( 1. - this->mistagSSBbar() ) );
			}
		}
		else
		{
			if( _mistagSS > _mistagOS )
			{
				thisEta = ( this->mistagOSB() * ( 1. - this->mistagSSB() ) ) / ( this->mistagOSB() * ( 1. - this->mistagSSB() ) + ( 1. - this->mistagOSB() ) * this->mistagSSB() );
			}
			else
			{
				thisEta = ( ( 1. - this->mistagOSB() ) * this->mistagSSB() ) / ( ( 1. - this->mistagOSB() ) * this->mistagSSB() + this->mistagOSB() * ( 1. - this->mistagSSB() ) );
			}
		}
	}

	return thisEta;
}

double CombinedMistagCalib::q() const
{
	if( _tagOS && !_tagSS )
	{
		return (double)_tagOS;
	}
	else if( !_tagSS && _tagOS )
	{
		return (double)_tagSS;
	}
	else
	{
		return (double)_combinedtag;
	}
}

double CombinedMistagCalib::mistag() const
{
	/*
	//Normal case
	returnValue =  _mistagP0 + _mistagP1*(_mistag - _mistagSetPoint ) ;
	if( returnValue < 0 )  returnValue = 0;
	if( returnValue > 0.5) returnValue = 0.5;
	return returnValue; */
	return 0.5;
}

double CombinedMistagCalib::mistagBbar() const
{
	double returnValue=0.;
	if( !_tagOS && !_tagSS )
	{
		returnValue = 0.5;
	}
	else if( _tagOS && !_tagSS )
	{
		returnValue = this->mistagOSBbar();
	}
	else if( !_tagOS && _tagSS )
	{
		returnValue = this->mistagSSBbar();
	}
	else
	{
		returnValue = _mistagP0_OSSS-(_mistagDeltaP0_OSSS*0.5) + (_mistagP1_OSSS-(_mistagDeltaP1_OSSS*0.5))*( _eta - (_mistagSetPoint_OSSS-(_mistagDeltaSetPoint_OSSS*0.5)) );
	}

	return returnValue;
}

double CombinedMistagCalib::mistagOSBbar() const
{
	return _mistagP0_OS-(_mistagDeltaP0_OS*0.5) + (_mistagP1_OS-(_mistagDeltaP1_OS*0.5))*( _mistagOS - (_mistagSetPoint_OS-(_mistagDeltaSetPoint_OS*0.5)) );
}

double CombinedMistagCalib::mistagSSBbar() const
{
	return _mistagP0_SS-(_mistagDeltaP0_SS*0.5) + (_mistagP1_SS-(_mistagDeltaP1_SS*0.5))*( _mistagSS - (_mistagSetPoint_SS-(_mistagDeltaSetPoint_SS*0.5)) );
}

double CombinedMistagCalib::mistagOSB() const
{
	return _mistagP0_OS+(_mistagDeltaP0_OS*0.5) + (_mistagP1_OS+(_mistagDeltaP1_OS*0.5))*( _mistagOS - (_mistagSetPoint_OS+(_mistagDeltaSetPoint_OS*0.5)) );
}

double CombinedMistagCalib::mistagSSB() const
{
	return _mistagP0_SS+(_mistagDeltaP0_SS*0.5) + (_mistagP1_SS+(_mistagDeltaP1_SS*0.5))*( _mistagSS - (_mistagSetPoint_SS+(_mistagDeltaSetPoint_SS*0.5)) );
}

double CombinedMistagCalib::mistagB() const
{
	double returnValue=0.;
	if( !_tagOS && !_tagSS )
	{
		returnValue = 0.5;
	}
	else if( _tagOS && !_tagSS )
	{
		returnValue = this->mistagOSB();
	}
	else if( !_tagOS && _tagSS )
	{
		returnValue = this->mistagSSB();
	}
	else
	{
		returnValue = _mistagP0_OSSS+(_mistagDeltaP0_OSSS*0.5) + (_mistagP1_OSSS+(_mistagDeltaP1_OSSS*0.5))*( _eta - (_mistagSetPoint_OSSS+(_mistagDeltaSetPoint_OSSS*0.5)) );
	}

	return returnValue;
}

double CombinedMistagCalib::D1() const
{
	return 1.0 - this->q()*( this->mistagB() - this->mistagBbar() );
}

double CombinedMistagCalib::D2() const
{
	return this->q()*( 1.0 - this->mistagB() - this->mistagBbar() );
}



