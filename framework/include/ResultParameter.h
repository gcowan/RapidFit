/**
        @class ResultParameter

        A physics parameter after it has been fitted

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef RESULT_PARAMETER_H
#define RESULT_PARAMETER_H

//	RapidFit Headers
#include "PhysicsParameter.h"
//	System Headers
#include <string>

using namespace::std;

class ResultParameter
{
	public:
		ResultParameter( const PhysicsParameter* );
		ResultParameter( string, double, double, double, double, double, string, string );
		~ResultParameter();

		bool GetAssym() const;
		void SetAssymErrors( const double err_plus, const double err_minus, const double err_para );
		double GetErrLow() const;
		double GetErrHi() const;

		string GetName() const;
		double GetValue() const;
		double GetOriginalValue() const;
		double GetError() const;
		void SetError( double );
		double GetPull() const;
		double GetMinimum() const;
		double GetMaximum() const;
		string GetType() const;
		void ForceType( string );
		string GetUnit() const;
		void ForceOriginalValue( double );
		void ForcePullValue( double );
		PhysicsParameter* GetDummyPhysicsParameter() const;
		double GetStepSize() const;
		void SetScanStatus( bool );
		bool GetScanStatus() const;

		void Print() const;

		string FitXML() const;
		string ToyXML() const;
		string ConstraintXML() const;

	private:
		string name;
		double value;
		double originalValue;
		double error;
		double minimum;
		double maximum;
		double stepSize;
		string type;
		string unit;
		bool ScanStatus;

		string XML( const bool=true ) const;

		bool assym_err;
		double err_hi, err_low;
};

#endif

