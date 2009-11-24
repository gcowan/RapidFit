/**
        @class InvalidObject

        A generic return type indicating the requested class could not be instantiated

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef INVALID_OBJECT_H
#define INVALID_OBJECT_H

#include "IPDF.h"
#include "FitFunction.h"
#include "IMinimiser.h"
#include "IDataSet.h"
#include "IDataGenerator.h"

class InvalidObject : public IPDF, public FitFunction, public IMinimiser, public IDataSet, public IDataGenerator
{
	public:
		InvalidObject();
		InvalidObject(string);
		~InvalidObject();

		//IPDF
		virtual bool IsValid();
		virtual bool SetPhysicsParameters( ParameterSet* );
		virtual double Integral( DataPoint*, PhaseSpaceBoundary* );
		virtual double Evaluate( DataPoint* );
		virtual vector<string> GetPrototypeDataPoint();
		virtual vector<string> GetPrototypeParameterSet();
		virtual vector<string> GetDoNotIntegrateList();

		//IMinimiser
                virtual void Minimise( FitFunction* );
                virtual FitResult * GetFitResult();

		//IDataSet
		virtual DataPoint * GetDataPoint(int);
		virtual void AddDataPoint( DataPoint* );
		virtual int GetDataNumber();
		virtual PhaseSpaceBoundary * GetBoundary();

		//IDataGenerator
                virtual int GenerateData(int);
		virtual IDataSet * GetDataSet();

                //FitFunction
		virtual double UpErrorValue();

	protected:
                //FitFunction
		virtual double EvaluateDataSet( IPDF*, IDataSet* );
		virtual double EvaluateParameterSet( ParameterSet* );

	private:
		string errorMessage;
};

#endif
