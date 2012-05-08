/**
        @class PhysicsParameter

        A parameter to be adjusted by the fitter, with a starting value and limits

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef PHYSICS_PARAMETER_H
#define PHYSICS_PARAMETER_H

//	System Headers
#include <string>

using namespace::std;

class PhysicsParameter
{
	public:
		PhysicsParameter( string );
		PhysicsParameter( string, double, double, double, double, string, string );
		PhysicsParameter( string, double, double, string, string );
		~PhysicsParameter();

		void SetName( string );
		string GetName() const;

		double GetValue() const;
		void SetValue(double);

		double GetBlindedValue() const;
		void SetBlindedValue(double);

		double GetTrueValue() const;
		void SetTrueValue(double);

		double GetMinimum() const;
		void SetMinimum(double);

		double GetMaximum() const;
		void SetMaximum(double);

		void SetLimits(double, double);
	
		void SetBlindOffset( double ) ;
		void SetBlinding( bool ) ;

		string GetType() const;
		void SetType(string);

		double GetOriginalValue() const;
		void ForceOriginalValue( double );

		string GetUnit() const;

		double GetStepSize() const;
		void SetStepSize( double );

		void Print() const;

		string XML() const;

		void SetBlindingInfo( string, double );

	private:
		string name;
		double value;
		double originalValue;
		double minimum;
		double maximum;
		double stepSize;
		string type;
		string unit;
	
		bool toBeBlinded;
		double blindOffset;

		string blindString;
		double blindScale;
};

#endif

