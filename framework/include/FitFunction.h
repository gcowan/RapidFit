/**
        @class FitFunction

        Parent class for the function to minimise
        Overload the evaluate methods and UP value for Chi2, NLL, etc.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef FIT_FUNCTION_H
#define FIT_FUNCTION_H

//	ROOT Headers
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
//	RapidFit Headers
#include "PhysicsBottle.h"
#include "RapidFitIntegrator.h"
#include "ObservableRef.h"

class FitFunction
{
	public:
		FitFunction();
		~FitFunction();

		void SetupTrace( TString FileName );
		void SetPhysicsBottle( PhysicsBottle* );
		PhysicsBottle * GetPhysicsBottle();
		bool SetParameterSet( ParameterSet* );
		ParameterSet * GetParameterSet();
		double Evaluate();
		void Finalise();
		void UseEventWeights(string);

		//Overload this function in child classes
		virtual double UpErrorValue(int);

	protected:
		//	Uncopyable!
		FitFunction ( const FitFunction& );
		FitFunction& operator = ( const FitFunction& );
		//Overload these functions in child classes
		virtual double EvaluateDataSet( IPDF*, IDataSet*, RapidFitIntegrator* );

		PhysicsBottle * allData;
		vector< RapidFitIntegrator* > allIntegrators;
		double testDouble;
		bool useWeights;
		string weightObservableName;

		//	This is fr traing Minuit
		//	(this could give some VERY oool graphs in ResultSpace :D )
		TFile* Fit_File;
		TTree* Fit_Tree;
		vector<Double_t> branch_objects;
		vector<ObservableRef> branch_names;
		Double_t fit_calls;
};

#endif
