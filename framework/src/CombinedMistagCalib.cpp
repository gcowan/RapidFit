
#include "ParameterSet.h"
#include "DataPoint.h"
#include "PDFConfigurator.h"
#include "DataPoint.h"
#include "ParameterSet.h"
#include "IMistagCalib.h"
#include "CombinedMistagCalib.h"
#include "SimpleMistagCalib.h"
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
	_mistagOSSS(), _combinedtag(), _OSTagged(false), _SSTagged(false), _OSSSTagged(false), _storedD1(), _storedD2(),
	//	Observable Names
	tagOSName( configurator->getName("tagdecision_os") ),
	tagSSName( configurator->getName("tagdecision_ss") ),
	tagCombName( configurator->getName("tagdecision") ),
	mistagOSName( configurator->getName("tagomega_os") ),
	mistagSSName( configurator->getName("tagomega_ss") ),
	mistagOSSSName( configurator->getName("tagomega") ),
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
	mistagDeltaSetPointName_OSSS( configurator->getName("mistagDeltaSetPoint_OSSS") ),
	_IMistagCalib_OS(NULL), _IMistagCalib_SS(NULL), _IMistagCalib_OSSS(NULL),
	_debugMistag(false), _onTuple(false), _useFixedEta(false)
{
	_debugMistag = configurator->isTrue( "DebugMistagModel" );
	_onTuple = configurator->isTrue( "MistagUseTuple" );

	_useFixedEta = configurator->isTrue( "useFixedCombination" );

	if( _debugMistag )
	{
		PDFConfigurator* OSCalib = new PDFConfigurator( *configurator );
		OSCalib->appendParameterNames( "mistagP0:mistagP1:mistagSetPoint:mistagDeltaP0:mistagDeltaP1:mistagDeltaSetPoint:_OS" );
		OSCalib->addParameterSubstitution( "tag:tagdecision" );
		OSCalib->addParameterSubstitution( "mistag:tagomega" );
		PDFConfigurator* SSCalib = new PDFConfigurator( *configurator );
		SSCalib->appendParameterNames( "mistagP0:mistagP1:mistagSetPoint:mistagDeltaP0:mistagDeltaP1:mistagDeltaSetPoint:_SS" );
		SSCalib->addParameterSubstitution( "tag:tagdecision" );
		SSCalib->addParameterSubstitution( "mistag:tagomega" );
		PDFConfigurator* OSSSCalib = new PDFConfigurator( *configurator );
		OSSSCalib->appendParameterNames( "mistagP0:mistagP1:mistagSetPoint:mistagDeltaP0:mistagDeltaP1:mistagDeltaSetPoint:_OSSS" );
		OSSSCalib->addParameterSubstitution( "tag:tagdecision" );
		OSSSCalib->addParameterSubstitution( "mistag:tagomega" );

		_IMistagCalib_OS = new SimpleMistagCalib( OSCalib );
		_IMistagCalib_SS = new SimpleMistagCalib( SSCalib );
		_IMistagCalib_OSSS = new SimpleMistagCalib( OSSSCalib );

		delete OSCalib;
		delete SSCalib;
		delete OSSSCalib;
	}
}

CombinedMistagCalib::~CombinedMistagCalib()
{
	if( _IMistagCalib_OS != NULL ) delete _IMistagCalib_OS;
	if( _IMistagCalib_SS != NULL ) delete _IMistagCalib_SS;
	if( _IMistagCalib_OSSS != NULL ) delete _IMistagCalib_OSSS;
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

	if( _debugMistag )
	{
		_IMistagCalib_OS->setParameters( parameters );
		_IMistagCalib_SS->setParameters( parameters );
		_IMistagCalib_OSSS->setParameters( parameters );
	}
}


void CombinedMistagCalib::addObservables( vector<string>& observableNames ) const
{
	observableNames.push_back( tagOSName );
	observableNames.push_back( mistagOSName );
	observableNames.push_back( tagSSName );
	observableNames.push_back( mistagSSName );
	observableNames.push_back( tagCombName );
	if( _debugMistag || _onTuple ) observableNames.push_back( mistagOSSSName );
}

void CombinedMistagCalib::setObservables( const DataPoint* measurement )
{
	_tagOS = (int) measurement->GetObservable( tagOSName )->GetValue();
	_mistagOS = measurement->GetObservable( mistagOSName )->GetValue();
	_tagSS = (int) measurement->GetObservable( tagSSName )->GetValue();
	_mistagSS = measurement->GetObservable( mistagSSName )->GetValue();

	_combinedtag = (int)measurement->GetObservable( tagCombName )->GetValue();

	if( _tagOS != 0 && _tagSS == 0 ) _OSTagged = true;
	else _OSTagged = false;

	if( _tagOS == 0 && _tagSS != 0 ) _SSTagged = true;
	else _SSTagged = false;

	if( _tagOS != 0 && _tagSS != 0 && _combinedtag != 0 ) _OSSSTagged = true;
	else _OSSSTagged = false;

	if( _onTuple )
	{
		_mistagOS = measurement->GetObservable( mistagOSSSName )->GetValue();
		_mistagSS = measurement->GetObservable( mistagOSSSName )->GetValue();
		_mistagOSSS = measurement->GetObservable( mistagOSSSName )->GetValue();
	}
	else
	{
		//_OSSSTagged = _OSTagged && _SSTagged;
		if( _OSSSTagged )
		{
			if( _useFixedEta ) _mistagOSSS = this->getFixedEta();
			else _mistagOSSS = this->getFloatedEta();
			_mistagOS = 0.5;
			_mistagSS = 0.5;
		}
	}

	if( _OSSSTagged )
	{
		_OSTagged = false;
		_SSTagged = false;
	}
	else
	{
		_mistagOSSS = 0.5;
	}

	if( _OSTagged ) _mistagSS = 0.5;
	if( _SSTagged ) _mistagOS = 0.5;

	_storedD1 = this->RealD1();
	_storedD2 = this->RealD2();


	if( _debugMistag )
	{

		if( _OSTagged )
		{
			double thisEta = measurement->GetObservable( mistagOSSSName )->GetValue();
			if( fabs( thisEta - _mistagOS ) > 1E-6 )
			{
				cout << thisEta << "\t\t" << _mistagOS << endl << endl;
				exit(0);
			}
		}

		if( _SSTagged )
		{
			double thisEta = measurement->GetObservable( mistagOSSSName )->GetValue();
			if( fabs( thisEta - _mistagSS ) > 1E-6 )
			{
				cout << thisEta << "\t\t" << _mistagSS << endl << endl;
				exit(0);
			}
		}


		int thisTag = 0;

		if( _OSTagged ) thisTag = _tagOS;
		if( _SSTagged ) thisTag = _tagSS;
		if( _OSSSTagged ) thisTag = _combinedtag;

		bool diff = false;
		if( thisTag != _combinedtag )
		{
			cout << "this" << "\t\t" << "combined" << endl;
			cout << thisTag << "\t\t" << _combinedtag << endl;
			diff = true;
		}

		if( this->OSTagged() )
		{
			_IMistagCalib_OS->setObservables( measurement );
			double OS_D1 = _IMistagCalib_OS->D1();
			double OS_D2 = _IMistagCalib_OS->D2();
			double OS_q = _IMistagCalib_OS->q();
			double OS_mB = _IMistagCalib_OS->mistagB();
			double OS_mBb = _IMistagCalib_OS->mistagBbar();
			if( fabs( OS_D1 - _storedD1 ) > 1E-6 )
			{
				cout << " Diff OST D1 !!" << endl;
				cout << _storedD1 << "\t" << OS_D1 << endl;
				diff = true;
			}
			if( fabs( OS_D2 - _storedD2 ) > 1E-6 )
			{
				cout << " Diff OST D2 !!" << endl;
				cout << _storedD2 << "\t" << OS_D2 << endl;
				diff = true;
			}
			if( fabs( OS_q - _tagOS ) > 1E-6 )
			{
				cout << " Diff in OS Tag!!" << endl;
				cout << _tagOS << "\t" << OS_q << endl;
				diff = true;
			}
			if( fabs( OS_mB - this->mistagB() ) > 1E-6 )
			{
				cout << " Diff in OS mistagB!!" << endl;
				cout << OS_mB << "\t" << this->mistagB() << endl;
				diff = true;
			}
			if( fabs( OS_mBb - this->mistagBbar() ) > 1E-6 )
			{
				cout << " Diff in OS mistagBar!!" << endl;
				cout << OS_mBb << "\t" << this->mistagBbar() << endl;
				diff = true;
			}
			if( diff )
			{
				if( _onTuple ) cout << "TRUE" << endl;
				else cout << "FALSE" << endl;
				_IMistagCalib_OS->Print();
				this->Print();
			}
		}
		if( this->SSTagged() )
		{
			_IMistagCalib_SS->setObservables( measurement );
			double OS_D1 = _IMistagCalib_SS->D1();
			double OS_D2 = _IMistagCalib_SS->D2();
			double OS_q = _IMistagCalib_SS->q();
			double OS_mB = _IMistagCalib_SS->mistagB();
			double OS_mBb = _IMistagCalib_SS->mistagBbar();
			if( fabs( OS_D1 - _storedD1 ) > 1E-6 )
			{
				cout << " Diff SST D1 !!" << endl;
				diff = true;
			}
			if( fabs( OS_D2 - _storedD2 ) > 1E-6 )
			{
				cout << " Diff SST D2 !!" << endl;
				diff = true;
			}
			if( fabs( OS_q - _tagSS ) > 1E-6 )
			{
				cout << " Diff in OS Tag!!" << endl;
				cout << _tagSS << "\t" << OS_q << endl;
				diff = true;
			}
			if( fabs( OS_mB - this->mistagB() ) > 1E-6 )
			{
				cout << " Diff in SS mistagB!!" << endl;
				cout << OS_mB << "\t" << this->mistagB() << endl;
				diff = true;
			}
			if( fabs( OS_mBb - this->mistagBbar() ) > 1E-6 )
			{
				cout << " Diff in SS mistagBar!!" << endl;
				cout << OS_mBb << "\t" << this->mistagBbar() << endl;
				diff = true;
			}
			if( diff )
			{
				if( _onTuple ) cout << "TRUE" << endl;
				else cout << "FALSE" << endl;
				_IMistagCalib_SS->Print();
				this->Print();
			}
		}
		if( this->OSSSTagged() )
		{
			_IMistagCalib_OSSS->setObservables( measurement );
			double OS_D1 = _IMistagCalib_OSSS->D1();
			double OS_D2 = _IMistagCalib_OSSS->D2();
			double OS_q = _IMistagCalib_OSSS->q();
			double OS_mB = _IMistagCalib_OSSS->mistagB();
			double OS_mBb = _IMistagCalib_OSSS->mistagBbar();
			if( fabs( OS_D1 - _storedD1 ) > 1E-6 )
			{
				cout << " Diff OSSS D1 !!" << endl;
				cout << _storedD1 << "\t" << OS_D1 << endl;
				diff = true;
			}
			if( fabs( OS_D2 - _storedD2 ) > 1E-6 )
			{
				cout << " Diff OSSS D2 !!" << endl;
				cout << _storedD2 << "\t" << OS_D2 << endl;
				diff = true;
			}
			if( fabs( OS_q - _combinedtag ) > 1E-6 )
			{
				cout << " Diff in OS Tag!!" << endl;
				cout << _combinedtag << "\t" << OS_q << endl;
				diff = true;
			}
			if( fabs( OS_mB - this->mistagB() ) > 1E-6 )
			{
				cout << " Diff in OSSS mistagB!!" << endl;
				cout << OS_mB << "\t" << this->mistagB() << endl;
				diff = true;
			}
			if( fabs( OS_mBb - this->mistagBbar() ) > 1E-6 )
			{
				cout << " Diff in OSSS mistagBar!!" << endl;
				cout << OS_mBb << "\t" << this->mistagBbar() << endl;
				diff = true;
			}
			if( diff )
			{
				if( _onTuple ) cout << "TRUE" << endl;
				else cout << "FALSE" << endl;
				_IMistagCalib_OSSS->Print();
				this->Print();
			}
		}

		if( diff )
		{
			measurement->Print();
			exit(0);
		}
	}
}

bool CombinedMistagCalib::OSTagged() const
{
	return _OSTagged;
}

bool CombinedMistagCalib::SSTagged() const
{
	return _SSTagged;
}

bool CombinedMistagCalib::OSSSTagged() const
{
	return _OSSSTagged;
}

double CombinedMistagCalib::getFixedEta() const
{
	double thisEta = 0.;
	if( _tagOS == _tagSS )
	{
		thisEta = ( _mistagOS *_mistagSS ) / ( _mistagOS*_mistagSS + (1.-_mistagOS)*(1.-_mistagSS) );
	}
	else
	{
		if( _mistagSS > _mistagOS )
		{
			thisEta = ( _mistagOS * ( 1. - _mistagSS ) ) / ( _mistagOS * ( 1. - _mistagSS ) + ( 1. - _mistagOS ) * _mistagSS );
		}
		else
		{
			thisEta = ( ( 1. - _mistagOS ) * _mistagSS ) / ( ( 1. - _mistagOS ) * _mistagSS + _mistagOS * ( 1. - _mistagSS ) );
		}
	}

	return thisEta;
}


double CombinedMistagCalib::getFloatedEta() const
{
	double thisEta = 0.;
	if( _tagOS == _tagSS )
	{
		if( _tagOS == -1 )
		{
			thisEta = ( this->mistagOSBbar() *this->mistagSSBbar() ) / ( this->mistagOSBbar()*this->mistagSSBbar() + (1.-this->mistagOSBbar())*(1.-this->mistagSSBbar()) );
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
	double returnValue=0.;

	if( this->OSTagged() )
	{
		returnValue = (double)_tagOS;
	}
	else if( this->SSTagged() )
	{
		returnValue = (double)_tagSS;
	}
	else if( this->OSSSTagged() )
	{
		returnValue = (double)_combinedtag;
	}
	else
	{
		returnValue = 0.;
	}

	return returnValue;
}

double CombinedMistagCalib::mistagBbar() const
{
	double returnValue=0.;

	if( this->OSTagged() )
	{
		returnValue = this->mistagOSBbar();
	}
	else if( this->SSTagged() )
	{
		returnValue = this->mistagSSBbar();
	}
	else if( this->OSSSTagged() )
	{
		returnValue = this->mistagOSSSBbar();
	}
	else
	{
		returnValue = 0.5;
	}

	if( returnValue > 0.5 )
	{
		returnValue = 0.5;
	}
	else if( returnValue < 0. )
	{
		returnValue = 0.;
	}

	return returnValue;
}

double CombinedMistagCalib::mistagOSBbar() const
{
	if( fabs(this->q()) < 0.5 ) return 0.5;

	double returnValue = 0.;
	if( _mistagOS > 0.5 ) returnValue = 0.5;
	else if( _mistagOS < 0. ) returnValue = 0.;
	else returnValue = _mistagP0_OS-(_mistagDeltaP0_OS*0.5) + (_mistagP1_OS-(_mistagDeltaP1_OS*0.5))*( _mistagOS - (_mistagSetPoint_OS-(_mistagDeltaSetPoint_OS*0.5)) );
	return returnValue;
}

double CombinedMistagCalib::mistagSSBbar() const
{
	if( fabs(this->q()) < 0.5 ) return 0.5;

	double returnValue = 0.;
	if( _mistagSS > 0.5 ) returnValue = 0.5;
	else if( _mistagSS < 0. ) returnValue = 0.;
	else returnValue = _mistagP0_SS-(_mistagDeltaP0_SS*0.5) + (_mistagP1_SS-(_mistagDeltaP1_SS*0.5))*( _mistagSS - (_mistagSetPoint_SS-(_mistagDeltaSetPoint_SS*0.5)) );
	return returnValue;
}

double CombinedMistagCalib::mistagOSSSBbar() const
{
	if( fabs(this->q()) < 0.5 ) return 0.5;

	double returnValue = 0.;
	if( _mistagOSSS > 0.5 ) returnValue = 0.5;
	else if( _mistagOSSS < 0. ) returnValue = 0.;
	else returnValue = _mistagP0_OSSS-(_mistagDeltaP0_OSSS*0.5) + (_mistagP1_OSSS-(_mistagDeltaP1_OSSS*0.5))*( _mistagOSSS - (_mistagSetPoint_OSSS-(_mistagDeltaSetPoint_OSSS*0.5)) );
	return returnValue;
}

double CombinedMistagCalib::mistagOSB() const
{
	if( fabs(this->q()) < 0.5 ) return 0.5;

	double returnValue = 0.;
	if( _mistagOS > 0.5 ) returnValue = 0.5;
	else if( _mistagOS < 0. ) returnValue = 0.;
	else returnValue = _mistagP0_OS+(_mistagDeltaP0_OS*0.5) + (_mistagP1_OS+(_mistagDeltaP1_OS*0.5))*( _mistagOS - (_mistagSetPoint_OS+(_mistagDeltaSetPoint_OS*0.5)) );
	return returnValue;
}

double CombinedMistagCalib::mistagSSB() const
{
	if( fabs(this->q()) < 0.5 ) return 0.5;

	double returnValue = 0.;
	if( _mistagSS > 0.5 ) returnValue = 0.5;
	else if( _mistagSS < 0. ) returnValue = 0.;
	else returnValue = _mistagP0_SS+(_mistagDeltaP0_SS*0.5) + (_mistagP1_SS+(_mistagDeltaP1_SS*0.5))*( _mistagSS - (_mistagSetPoint_SS+(_mistagDeltaSetPoint_SS*0.5)) );
	return returnValue;
}

double CombinedMistagCalib::mistagOSSSB() const
{
	if( fabs(this->q()) < 0.5 ) return 0.5;

	double returnValue = 0.;
	if( _mistagOSSS > 0.5 ) returnValue = 0.5;
	else if( _mistagOSSS < 0. ) returnValue = 0.;
	else returnValue = _mistagP0_OSSS+(_mistagDeltaP0_OSSS*0.5) + (_mistagP1_OSSS+(_mistagDeltaP1_OSSS*0.5))*( _mistagOSSS - (_mistagSetPoint_OSSS+(_mistagDeltaSetPoint_OSSS*0.5)) ); 
	return returnValue;
}

double CombinedMistagCalib::mistagB() const
{
	double returnValue=0.;

	if( this->OSTagged() )
	{
		returnValue = this->mistagOSB();
	}
	else if( this->SSTagged() )
	{
		returnValue = this->mistagSSB();
	}
	else if( this->OSSSTagged() )
	{
		returnValue = this->mistagOSSSB();
	}
	else
	{
		returnValue = 0.5;
	}

	if( returnValue > 0.5 )
	{
		returnValue = 0.5;
	}
	else if( returnValue < 0. )
	{
		returnValue = 0.;
	}

	return returnValue;
}

double CombinedMistagCalib::RealD1() const
{
	return 1.0 - this->q()*( this->mistagB() - this->mistagBbar() );
}

double CombinedMistagCalib::RealD2() const
{
	return this->q()*( 1.0 - this->mistagB() - this->mistagBbar() );
}

double CombinedMistagCalib::D1() const
{
	return _storedD1;
}

double CombinedMistagCalib::D2() const
{
	return _storedD2;
}

void CombinedMistagCalib::Print() const
{
	cout << endl;
	cout << "_mistagP0_OS              " << _mistagP0_OS << endl;
	cout << "_mistagP0_SS              " << _mistagP0_SS << endl;
	cout << "_mistagP0_OSSS            " << _mistagP0_OSSS << endl;
	cout << "_mistagDeltaP0_OS         " << _mistagDeltaP0_OS << endl;
	cout << "_mistagDeltaP0_SS         " << _mistagDeltaP0_SS << endl;
	cout << "_mistagDeltaP0_OSSS       " << _mistagDeltaP0_OSSS << endl;
	cout << "_mistagP1_OS              " << _mistagP1_OS << endl;
	cout << "_mistagP1_SS              " << _mistagP1_SS << endl;
	cout << "_mistagP1_OSSS            " << _mistagP1_OSSS << endl;
	cout << "_mistagDeltaP1_OS         " << _mistagDeltaP1_OS << endl;
	cout << "_mistagDeltaP1_SS         " << _mistagDeltaP1_SS << endl;
	cout << "_mistagDeltaP1_OSSS       " << _mistagDeltaP1_OSSS << endl;
	cout << "_mistagOS                 " << _mistagOS << endl;
	cout << "_mistagSS                 " << _mistagSS << endl;
	cout << "_mistagOSSS               " << _mistagOSSS << endl;
	cout << "_mistagSetPoint_OS        " << _mistagSetPoint_OS << endl;
	cout << "_mistagSetPoint_OSSS      " << _mistagSetPoint_OSSS << endl;
	cout << "_mistagDeltaSetPoint_OS   " << _mistagDeltaSetPoint_OS << endl;
	cout << "_mistagDeltaSetPoint_SS   " << _mistagDeltaSetPoint_SS << endl;
	cout << "_mistagDeltaSetPoint_OSSS " << _mistagDeltaSetPoint_OS << endl;
	cout << "this->mistagB()           " << this->mistagB() << endl;
	cout << "this->mistagBbar()        " << this->mistagBbar() << endl;
}

