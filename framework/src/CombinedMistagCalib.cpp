
#include "ParameterSet.h"
#include "DataPoint.h"
#include "PDFConfigurator.h"
#include "DataPoint.h"
#include "ParameterSet.h"
#include "IMistagCalib.h"
#include "CombinedMistagCalib.h"
#include "SimpleMistagCalib.h"
#include "ObservableRef.h"
#include "DebugClass.h"

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
	_combinedtag(), _OSTagged(false), _SSTagged(false), _OSSSTagged(false), _storedD1(), _storedD2(),
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
	_debugMistag(false), _onTuple(false), _floatCalib(false), _untagged(false)
{
	_debugMistag = configurator->isTrue( "DebugMistagModel" );
	_onTuple = ! configurator->isTrue( "Mistag3fbModel" );
	_floatCalib = configurator->isTrue( "FloatCombinedCalib" );

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

	if( _onTuple || _debugMistag )
	{
		parameterNames.push_back( mistagP1Name_OSSS );
		parameterNames.push_back( mistagP0Name_OSSS );
		parameterNames.push_back( mistagSetPointName_OSSS );
		parameterNames.push_back( mistagDeltaP1Name_OSSS );
		parameterNames.push_back( mistagDeltaP0Name_OSSS );
		parameterNames.push_back( mistagDeltaSetPointName_OSSS );
	}
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

	if( _onTuple || _debugMistag )
	{
		_mistagP1_OSSS = parameters.GetPhysicsParameter( mistagP1Name_OSSS )->GetValue();
		_mistagP0_OSSS = parameters.GetPhysicsParameter( mistagP0Name_OSSS )->GetValue();
		_mistagSetPoint_OSSS = parameters.GetPhysicsParameter( mistagSetPointName_OSSS )->GetValue();
		_mistagDeltaP1_OSSS = parameters.GetPhysicsParameter( mistagDeltaP1Name_OSSS )->GetValue();
		_mistagDeltaP0_OSSS = parameters.GetPhysicsParameter( mistagDeltaP0Name_OSSS )->GetValue();
		_mistagDeltaSetPoint_OSSS = parameters.GetPhysicsParameter( mistagDeltaSetPointName_OSSS )->GetValue();
	}
}


void CombinedMistagCalib::addObservables( vector<string>& observableNames ) const
{
	observableNames.push_back( tagOSName );
	observableNames.push_back( mistagOSName );
	observableNames.push_back( tagSSName );
	observableNames.push_back( mistagSSName );
	if( _onTuple || _debugMistag ) observableNames.push_back( tagCombName );
	if( _onTuple || _debugMistag ) observableNames.push_back( mistagOSSSName );
}

int CombinedMistagCalib::GetCombinedTag() const
{
	if( _tagSS == _tagOS ) return _tagSS;

	int decision = 0;

	double p_B_OS = this->GetProbBOS();
	double p_B_SS = this->GetProbBSS();
	double p_B = p_B_OS * p_B_SS;

	double p_Bb_OS = this->GetProbBbarOS();
	double p_Bb_SS = this->GetProbBbarSS();
	double p_Bb = p_Bb_OS * p_Bb_SS;

	if( p_B  > p_Bb ) decision = +1;
	else if( p_B < p_Bb ) decision = -1;
	else decision = 0;

	return decision;
}

double CombinedMistagCalib::GetProbBOS() const
{
	if( !_floatCalib )
	{
		return (1.-(double)_tagOS)*0.5 + ((double)_tagOS)*(1.-_mistagOS);
	}
	else
	{
		return (1.-(double)_tagOS)*0.5 + ((double)_tagOS)*(1.-this->mistagOSB());
	}
}

double CombinedMistagCalib::GetProbBbarOS() const
{
	if( !_floatCalib )
	{
		return (1.+(double)_tagOS)*0.5 - ((double)_tagOS)*(1.-_mistagOS);
	}
	else
	{
		return (1.+(double)_tagOS)*0.5 - ((double)_tagOS)*(1.-this->mistagOSBbar());
	}
}

double CombinedMistagCalib::GetProbBSS() const
{
	if( !_floatCalib )
	{
		return (1.-(double)_tagSS)*0.5 + ((double)_tagSS)*(1.-_mistagSS);
	}
	else
	{
		return (1.-(double)_tagSS)*0.5 + ((double)_tagSS)*(1.-this->mistagSSB());
	}
}

double CombinedMistagCalib::GetProbBbarSS() const
{
	if( !_floatCalib )
	{
		return (1.+(double)_tagSS)*0.5 - ((double)_tagSS)*(1.-_mistagSS);
	}
	else
	{
		return (1.+(double)_tagSS)*0.5 - ((double)_tagSS)*(1.-this->mistagSSBbar());
	}
}

void CombinedMistagCalib::setObservables( DataPoint* measurement )
{
	double readTagOS = measurement->GetObservable( tagOSName )->GetValue();
	_tagOS = (readTagOS>=0.)?(int)ceil(readTagOS):(int)floor(readTagOS);
	_mistagOS = measurement->GetObservable( mistagOSName )->GetValue();
	double readTagSS = (int) measurement->GetObservable( tagSSName )->GetValue();
	_tagSS = (readTagSS>=0.)?(int)ceil(readTagSS):(int)floor(readTagSS);
	_mistagSS = measurement->GetObservable( mistagSSName )->GetValue();

	if( _tagOS != 0 && _tagSS == 0 ) _OSTagged = true;
	else _OSTagged = false;

	if( _tagOS == 0 && _tagSS != 0 ) _SSTagged = true;
	else _SSTagged = false;

	if( _tagOS != 0 && _tagSS != 0 ) _OSSSTagged = true;
	else _OSSSTagged = false;

	if( _OSSSTagged )
	{
		_OSTagged = false;
		_SSTagged = false;
		_combinedtag = this->GetCombinedTag();
	}

	if( (_tagOS == _tagSS) && (_tagOS == 0) ) _untagged = true;
	else _untagged = false;

	_storedD1 = this->RealD1();
	_storedD2 = this->RealD2();

	if( _debugMistag )
	{

		double readTagOSSS = measurement->GetObservable( tagCombName )->GetValue();
		int _tagOSSS = (readTagOSSS>=0.)?(int)ceil(readTagOSSS):(int)floor(readTagOSSS);

		double this_q = this->q();
		int _calculated = (this_q>=0.)?(int)ceil(this_q):(int)floor(this_q);

		vector<double> allData;

		vector<int> tagdecisions;
		tagdecisions.push_back( _tagOS );	allData.push_back( _tagOS );
		tagdecisions.push_back( _tagSS );	allData.push_back( _tagSS );
		tagdecisions.push_back( _tagOSSS );	allData.push_back( _tagOSSS );
		tagdecisions.push_back( _calculated ); 	allData.push_back( _calculated );

		DebugClass::AppendToFile( "decisions.txt", tagdecisions );

		_mistagOSSS = measurement->GetObservable( mistagOSSSName )->GetValue();

		double calculated_mistag_2011 = this->getFixedEta();
		double calculated_mistag_2013 = this->getFloatedMistag();

		double oldval1, oldval2, oldval3, oldval4, oldval5, oldval6;

		oldval1 = _mistagDeltaSetPoint_OS;
		oldval2 = _mistagDeltaP1_OS;
		oldval3 = _mistagDeltaP0_OS;

		oldval4 = _mistagDeltaSetPoint_SS;
		oldval5 = _mistagDeltaP1_SS;
		oldval6 = _mistagDeltaP0_SS;

		_mistagDeltaSetPoint_OS = 0.;
		_mistagDeltaP1_OS = 0.;
		_mistagDeltaP0_OS = 0.;
		_mistagDeltaSetPoint_SS = 0.;
		_mistagDeltaP1_SS = 0.;
		_mistagDeltaP0_SS = 0.;

		double calculated_mistag_2013_woAssym = this->getFloatedMistag();

		_mistagDeltaSetPoint_OS = oldval1;
		_mistagDeltaP1_OS = oldval2;
		_mistagDeltaP0_OS = oldval3;
		_mistagDeltaSetPoint_SS = oldval4;
		_mistagDeltaP1_SS = oldval5;
		_mistagDeltaP0_SS = oldval6;

		vector<double> mistagvalues;
		mistagvalues.push_back( _mistagOS );				allData.push_back( _mistagOS );
		mistagvalues.push_back( _mistagSS );				allData.push_back( _mistagSS );
		mistagvalues.push_back( _mistagOSSS );				allData.push_back( _mistagOSSS );
		mistagvalues.push_back( calculated_mistag_2011 );		allData.push_back( calculated_mistag_2011 );
		mistagvalues.push_back( calculated_mistag_2013 );		allData.push_back( calculated_mistag_2013 );
		mistagvalues.push_back( calculated_mistag_2013_woAssym );	allData.push_back( calculated_mistag_2013_woAssym );

		DebugClass::AppendToFile( "mistags.txt", tagdecisions );

		DebugClass::AppendToFile( "AllData.txt", allData );
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


double CombinedMistagCalib::getFloatedMistag() const
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
		if( _mistagSS > _mistagOS )
		{
			if( _tagOS == -1 )
			{
				if( _tagSS == -1 )
				{
					thisEta = ( this->mistagOSBbar() * ( 1. - this->mistagSSBbar() ) ) / ( this->mistagOSBbar() * ( 1. - this->mistagSSBbar() ) + ( 1. - this->mistagOSBbar() ) * this->mistagSSBbar() );
				}
				else
				{
					thisEta = ( this->mistagOSBbar() * ( 1. - this->mistagSSB() ) ) / ( this->mistagOSBbar() * ( 1. - this->mistagSSB() ) + ( 1. - this->mistagOSBbar() ) * this->mistagSSB() );
				}
			}
			else
			{
				if( _tagSS == -1 )
				{
					thisEta = ( this->mistagOSB() * ( 1. - this->mistagSSBbar() ) ) / ( this->mistagOSB() * ( 1. - this->mistagSSBbar() ) + ( 1. - this->mistagOSB() ) * this->mistagSSBbar() );
				}
				else
				{
					thisEta = ( this->mistagOSB() * ( 1. - this->mistagSSB() ) ) / ( this->mistagOSB() * ( 1. - this->mistagSSB() ) + ( 1. - this->mistagOSB() ) * this->mistagSSB() );
				}
			}
		}
		else
		{
			if( _tagOS == -1 )
			{
				if( _tagSS == -1 )
				{
					thisEta = ( ( 1. - this->mistagOSBbar() ) * this->mistagSSBbar() ) / ( ( 1. - this->mistagOSBbar() ) * this->mistagSSBbar() + this->mistagOSBbar() * ( 1. - this->mistagSSBbar() ) );
				}
				else
				{
					thisEta = ( ( 1. - this->mistagOSBbar() ) * this->mistagSSB() ) / ( ( 1. - this->mistagOSBbar() ) * this->mistagSSB() + this->mistagOSBbar() * ( 1. - this->mistagSSB() ) );
				}
			}
			else
			{
				if( _tagSS == -1 )
				{
					thisEta = ( ( 1. - this->mistagOSB() ) * this->mistagSSBbar() ) / ( ( 1. - this->mistagOSB() ) * this->mistagSSBbar() + this->mistagOSB() * ( 1. - this->mistagSSBbar() ) );
				}
				else
				{
					thisEta = ( ( 1. - this->mistagOSB() ) * this->mistagSSB() ) / ( ( 1. - this->mistagOSB() ) * this->mistagSSB() + this->mistagOSB() * ( 1. - this->mistagSSB() ) );
				}
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
		if( _onTuple ) returnValue = (double)_combinedtag;
		else returnValue = (double)this->GetCombinedTag();
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
	//if( fabs(this->q()) < 0.5 ) return 0.5;

	double returnValue = 0.;
	if( _mistagOS > 0.5 ) returnValue = 0.5;
	else if( _mistagOS < 0. ) returnValue = 0.;
	else returnValue = _mistagP0_OS-(_mistagDeltaP0_OS*0.5) + (_mistagP1_OS-(_mistagDeltaP1_OS*0.5))*( _mistagOS - (_mistagSetPoint_OS-(_mistagDeltaSetPoint_OS*0.5)) );

	if( returnValue > 0.5 ) returnValue = 0.5;
	else if( returnValue < 0. ) returnValue = 0.;
	return returnValue;
}

double CombinedMistagCalib::mistagSSBbar() const
{
	//if( fabs(this->q()) < 0.5 ) return 0.5;

	double returnValue = 0.;
	if( _mistagSS > 0.5 ) returnValue = 0.5;
	else if( _mistagSS < 0. ) returnValue = 0.;
	else returnValue = _mistagP0_SS-(_mistagDeltaP0_SS*0.5) + (_mistagP1_SS-(_mistagDeltaP1_SS*0.5))*( _mistagSS - (_mistagSetPoint_SS-(_mistagDeltaSetPoint_SS*0.5)) );

	if( returnValue > 0.5 ) returnValue = 0.5;
	else if( returnValue < 0. ) returnValue = 0.;
	return returnValue;
}

double CombinedMistagCalib::mistagOSSSBbar() const
{
	//if( fabs(this->q()) < 0.5 ) return 0.5;

	double returnValue = 0.;
	if( _onTuple )
	{
		//      Mistag calculated from nTuple and so needs calibrating
		returnValue = _mistagP0_OSSS-(_mistagDeltaP0_OSSS*0.5) + (_mistagP1_OSSS-(_mistagDeltaP1_OSSS*0.5))*( this->getFixedEta() - (_mistagSetPoint_OSSS-(_mistagDeltaSetPoint_OSSS*0.5)) );
	}
	else
	{
		//      Mistag calculated using calibrated OS and SS so doesn't need calibrating
		returnValue = this->getFloatedMistag();
	}
	if( returnValue > 0.5 ) returnValue = 0.5;
	else if( returnValue < 0. ) returnValue = 0.;
	return returnValue;
}

double CombinedMistagCalib::mistagOSB() const
{
	//if( fabs(this->q()) < 0.5 ) return 0.5;

	double returnValue = 0.;
	if( _mistagOS > 0.5 ) returnValue = 0.5;
	else if( _mistagOS < 0. ) returnValue = 0.;
	else returnValue = _mistagP0_OS+(_mistagDeltaP0_OS*0.5) + (_mistagP1_OS+(_mistagDeltaP1_OS*0.5))*( _mistagOS - (_mistagSetPoint_OS+(_mistagDeltaSetPoint_OS*0.5)) );

	if( returnValue > 0.5 ) returnValue = 0.5;
	else if( returnValue < 0. ) returnValue = 0.;
	return returnValue;
}

double CombinedMistagCalib::mistagSSB() const
{
	//if( fabs(this->q()) < 0.5 ) return 0.5;

	double returnValue = 0.;
	if( _mistagSS > 0.5 ) returnValue = 0.5;
	else if( _mistagSS < 0. ) returnValue = 0.;
	else returnValue = _mistagP0_SS+(_mistagDeltaP0_SS*0.5) + (_mistagP1_SS+(_mistagDeltaP1_SS*0.5))*( _mistagSS - (_mistagSetPoint_SS+(_mistagDeltaSetPoint_SS*0.5)) );

	if( returnValue > 0.5 ) returnValue = 0.5;
	else if( returnValue < 0. ) returnValue = 0.;
	return returnValue;
}

double CombinedMistagCalib::mistagOSSSB() const
{
	//if( fabs(this->q()) < 0.5 ) return 0.5;

	double returnValue = 0.;
	if( _onTuple )
	{
		//	Mistag calculated from nTuple and so needs calibrating
		returnValue = _mistagP0_OSSS+(_mistagDeltaP0_OSSS*0.5) + (_mistagP1_OSSS+(_mistagDeltaP1_OSSS*0.5))*( this->getFixedEta() - (_mistagSetPoint_OSSS+(_mistagDeltaSetPoint_OSSS*0.5)) ); 
	}
	else
	{
		//	Mistag calculated using calibrated OS and SS so doesn't need calibrating
		returnValue = this->getFloatedMistag();
	}

	if( returnValue > 0.5 ) returnValue = 0.5;
	else if( returnValue < 0. ) returnValue = 0.;
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
	if( !_untagged )
	{
		return 1.0 - this->q()*( this->mistagB() - this->mistagBbar() );
	}
	else
	{
		return 1.;
	}
}

double CombinedMistagCalib::RealD2() const
{
	if( !_untagged )
	{
		return this->q()*( 1.0 - this->mistagBbar() - this->mistagB() );
	}
	else
	{
		return 0.;
	}
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
	cout << "_mistagSetPoint_OS        " << _mistagSetPoint_OS << endl;
	cout << "_mistagSetPoint_OSSS      " << _mistagSetPoint_OSSS << endl;
	cout << "_mistagDeltaSetPoint_OS   " << _mistagDeltaSetPoint_OS << endl;
	cout << "_mistagDeltaSetPoint_SS   " << _mistagDeltaSetPoint_SS << endl;
	cout << "_mistagDeltaSetPoint_OSSS " << _mistagDeltaSetPoint_OS << endl;
	cout << "this->mistagB()           " << this->mistagB() << endl;
	cout << "this->mistagBbar()        " << this->mistagBbar() << endl;
}

bool CombinedMistagCalib::eventIsTagged() const
{
	return !_untagged;
}

