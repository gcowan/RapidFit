/**
        @class ResultParameterSet

        A set of physics parameters after fitting

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef RESULT_PARAMETER_SET_H
#define RESULT_PARAMETER_SET_H

//	RapidFit Headers
#include "ResultParameter.h"
#include "ParameterSet.h"
#include "RapidFitMatrix.h"
//	System Headers
#include <vector>
#include <string>

using namespace::std;

class ResultParameterSet
{
	public:
		//Constructor
		ResultParameterSet( vector<string> );
		ResultParameterSet( const ResultParameterSet& );

		//Destructor
		~ResultParameterSet();

		//Get the list of parameters within
		vector<string> GetAllNames() const;

		vector<string> GetAllFloatNames() const;

		vector<string> GetAllFixedNames() const;

		//Get the result parameter according to this index, or this name
		ResultParameter * GetResultParameter( int ) const;
		ResultParameter * GetResultParameter( string ) const;
		ResultParameter * GetResultParameter( const ObservableRef& ) const;

		//Set the result parameter to be equal to:
		bool SetResultParameter( string, ResultParameter* );
		bool SetResultParameter( string, double, double, double, double, double, string, string );

		//Forcibly change or add a new parameter to this object
		bool ForceNewResultParameter( ResultParameter* );
		bool ForceNewResultParameter( string, double, double, double, double, double, string, string );

		//Get a ParameterSet constructed from the values of the parameters in this class
		ParameterSet* GetDummyParameterSet() const;

		void Print() const;

		string FitXML() const;
		string ToyXML() const;

		void ApplyCovarianceMatrix( RapidFitMatrix* Input );
	private:
		ResultParameterSet operator=( const ResultParameterSet& input );

		vector<ResultParameter*> allParameters;
		vector<string> allNames;

		string XML( const bool=true ) const;
};

#endif

