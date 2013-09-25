
#ifndef _SIMPLE_MISTAGCALIB_H_
#define _SIMPLE_MISTAGCALIB_H_

#include "IMistagCalib.h"
#include "PDFConfigurator.h"
#include "ParameterSet.h"
#include "DataPoint.h"

using namespace::std;

class SimpleMistagCalib : public IMistagCalib
{
	public:
		SimpleMistagCalib( PDFConfigurator* configurator );
		~SimpleMistagCalib();

		void addParameters( vector<string>& parameterNames ) const;
		void setParameters( const ParameterSet& parameters );

		void addObservables( vector<string>& observableNames ) const;
		void setObservables( const DataPoint* measurement );

		double D1() const;

		double D2() const;

		void Print() const;

		double q() const;

		double mistagBbar() const;

		double mistagB() const;

	private:

		double _tag, _mistag, _mistagP0, _mistagP1, _mistagSetPoint, _mistagDeltaP1, _mistagDeltaP0, _mistagDeltaSetPoint;
		ObservableRef tagName, mistagName, mistagP1Name, mistagP0Name, mistagSetPointName;
		ObservableRef mistagDeltaP1Name, mistagDeltaP0Name, mistagDeltaSetPointName;

		//double mistag() const;

};

#endif

