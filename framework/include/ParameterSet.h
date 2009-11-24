/**
        @class ParameterSet

        A collection of physics parameters

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef PARAMETER_SET_H
#define PARAMETER_SET_H

#include <vector>
#include <string>
#include "PhysicsParameter.h"

using namespace std;

class ParameterSet
{
	public:
		ParameterSet();
		ParameterSet( vector<string> );
		~ParameterSet();

		vector<string> GetAllNames();
		PhysicsParameter * GetPhysicsParameter(string);
		bool SetPhysicsParameter( string, PhysicsParameter* );
		bool SetPhysicsParameter( string, double, double, double, string, string );
		bool SetPhysicsParameters( ParameterSet* );

		//Not very nice in OO programming terms, and unsafe. Much faster though
		bool SetPhysicsParameters( double* );
		bool SetPhysicsParameters( vector<double> );

	private:
		vector<PhysicsParameter> allParameters;
		vector<string> allNames;
};

#endif
