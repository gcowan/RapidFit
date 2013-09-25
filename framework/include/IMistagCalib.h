
#ifndef _I_MISTAGCALIB_H_
#define _I_MISTAGCALIB_H_

#include "ParameterSet.h"
#include "DataPoint.h"

using namespace::std;

class IMistagCalib
{
	public:
		virtual ~IMistagCalib() {};

		virtual void addParameters( vector<string> & parameterNames ) const = 0;
		virtual void setParameters( const ParameterSet & parameters ) = 0;

		virtual void addObservables( vector<string> & observableNames ) const = 0;
		virtual void setObservables( const DataPoint * measurement ) = 0;

		virtual double D1() const = 0;

		virtual double D2() const = 0;

		virtual void Print() const = 0;

		virtual double q() const = 0;

		virtual double mistagBbar() const = 0;

		virtual double mistagB() const = 0;

	protected:

		IMistagCalib() {};

};

#endif

