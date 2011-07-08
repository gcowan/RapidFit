/**
        @class ParameterSet

        A collection of physics parameters

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef PARAMETER_SET_H
#define PARAMETER_SET_H

//	ROOT Headers
#include "TString.h"
//	RapidFit Headers
#include "ObservableRef.h"
#include "PhysicsParameter.h"
//	System Headers
#include <vector>
#include <string>

using namespace std;

class ParameterSet
{
	public:
		ParameterSet();
		ParameterSet( vector<string> );
		~ParameterSet();

		void SET_ID( TString );
		void SET_ID( string );
		string GET_ID();
		vector<string> GetAllNames();
		vector<string> GetAllFloatNames();
		vector<string> GetAllFixedNames();
		PhysicsParameter * GetPhysicsParameter( pair<string,int>* );
		PhysicsParameter * GetPhysicsParameter( string );
		PhysicsParameter * GetPhysicsParameter( ObservableRef& );
		bool SetPhysicsParameter( string, PhysicsParameter* );
		bool SetPhysicsParameter( string, double, double, double, double, string, string );
		bool SetPhysicsParameters( ParameterSet* );

		//Not very nice in OO programming terms, and unsafe. Much faster though
		bool SetPhysicsParameters( double* );
		bool SetPhysicsParameters( vector<double> );
	
		void print() ;

	private:
		string stored_id;
		vector<PhysicsParameter> allParameters;
		vector<string> allNames;
};

#endif

