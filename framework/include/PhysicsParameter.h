/**
        @class PhysicsParameter

        A parameter to be adjusted by the fitter, with a starting value and limits

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef PHYSICS_PARAMETER_H
#define PHYSICS_PARAMETER_H

//	System Headers
#include <string>

using namespace std;

class PhysicsParameter
{
	public:
		PhysicsParameter();
		PhysicsParameter( string, double, double, double, string, string );
		PhysicsParameter( string, double, string, string );
		~PhysicsParameter();

		double GetValue();
		void SetValue(double);
		double GetBlindedValue();
		void SetBlindedValue(double);
		double GetTrueValue();
		void SetTrueValue(double);

		double GetMinimum();
		void SetMinimum(double);

		double GetMaximum();
		void SetMaximum(double);

		void SetLimits(double, double);
	
		void SetBlindOffset( double ) ;
		void SetBlinding( bool ) ;

		string GetType();
		void SetType(string);

		double GetOriginalValue();
		void ForceOriginalValue( double );
		string GetUnit();
	
		void print() ;

	private:
		double value;
		double originalValue;
		double minimum;
		double maximum;
		string type;
		string unit;
	
		bool toBeBlinded ;
		double blindOffset ;
};

#endif
