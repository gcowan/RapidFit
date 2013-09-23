
#ifndef _COMBINED_MISTAGCALIB_H_
#define _COMBINED_MISTAGCALIB_H_

#include "IMistagCalib.h"
#include "PDFConfigurator.h"
#include "ParameterSet.h"
#include "DataPoint.h"

#include <string>

using namespace::std;

class CombinedMistagCalib : public IMistagCalib
{
	public:
		CombinedMistagCalib( PDFConfigurator* configurator );
		~CombinedMistagCalib();

		void addParameters( vector<string>& parameterNames ) const;
		void setParameters( const ParameterSet& parameters );

		void addObservables( vector<string>& observableNames ) const;
		void setObservables( const DataPoint* measurement );

		double q() const;

		double mistag() const;

		double mistagBbar() const;

		double mistagB() const;

		double D1() const;

		double D2() const;

	protected:

		double mistagSSB() const;
		double mistagOSB() const;
		double mistagSSBbar() const;
		double mistagOSBbar() const;

		double getEta() const;

		int _tagOS, _tagSS, _combinedtag;
		double _mistagOS, _mistagSS, _eta;
		double _mistagP0_OS, _mistagP1_OS, _mistagSetPoint_OS, _mistagDeltaP1_OS, _mistagDeltaP0_OS, _mistagDeltaSetPoint_OS;
		double _mistagP0_SS, _mistagP1_SS, _mistagSetPoint_SS, _mistagDeltaP1_SS, _mistagDeltaP0_SS, _mistagDeltaSetPoint_SS;
		double _mistagP0_OSSS, _mistagP1_OSSS, _mistagSetPoint_OSSS, _mistagDeltaP1_OSSS, _mistagDeltaP0_OSSS, _mistagDeltaSetPoint_OSSS;

		ObservableRef tagOSName, tagSSName, mistagOSName, mistagSSName, tagCombName;
		ObservableRef mistagP1Name_OS, mistagP0Name_OS, mistagSetPointName_OS, mistagDeltaP1Name_OS, mistagDeltaP0Name_OS, mistagDeltaSetPointName_OS;
		ObservableRef mistagP1Name_SS, mistagP0Name_SS, mistagSetPointName_SS, mistagDeltaP1Name_SS, mistagDeltaP0Name_SS, mistagDeltaSetPointName_SS;
		ObservableRef mistagP1Name_OSSS, mistagP0Name_OSSS, mistagSetPointName_OSSS, mistagDeltaP1Name_OSSS, mistagDeltaP0Name_OSSS, mistagDeltaSetPointName_OSSS;

};

#endif

