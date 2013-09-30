//	Class designed to contain common structs/functions required for multi-threading the fits in RapidFit

#pragma once
#ifndef RAPIDFIT_THREADING_H
#define RAPIDFIT_THREADING_H

#include "DataPoint.h"
#include "IDataSet.h"
#include "RapidFitIntegrator.h"
#include "ComponentRef.h"

#include <vector>
#include <string>

using namespace::std;

class IPDF;
class IDataSet;

//      Threading Struct which contains all of the objects required for running multiple concurrent fits to data subsets
//	This object is useful as multiple bits of information need to be provided to the running thread
struct Fitting_Thread{
	explicit Fitting_Thread() :
		dataSubSet(), fittingPDF(NULL), useWeights(false), dataPoint_Result(), FitBoundary(NULL),
		stored_integral(0.), weightsSquared(false), dataSet(NULL), thisComponent(NULL)
	{}

	vector<DataPoint*> dataSubSet;		/*!	DataPoints to be evaluated by this thread		*/
	IDataSet* dataSet;			/*!	DataSet containtaining the DataPoints			*/
	IPDF* fittingPDF;			/*!	Pointer to the PDF instance to be used by this thread	*/
	bool useWeights;			/*!	Are we performing a weighted fit?			*/
	vector<double> dataPoint_Result;	/*!	Result for evaluating each datapoint			*/
	PhaseSpaceBoundary* FitBoundary;	/*!	PhaseSpaceBoundary containing all data			*/
	double stored_integral;			/*!	Stored Integral for Numerical Integral fits		*/
	bool weightsSquared;			/*!	Are we using Weight Squared?				*/

	ComponentRef* thisComponent;

	private:
		Fitting_Thread(const Fitting_Thread&);
		Fitting_Thread& operator=(const Fitting_Thread&);
};

class Threading
{
	public:
		//	Number of cores on machine this is compiled for
		static int numCores();

		//	Split the data into subset(s) with a safe default
		static vector<vector<DataPoint*> > divideData( IDataSet*, int=1 );

		static vector<IDataSet*> divideDataSet( IDataSet* input, unsigned int subsets=1 );

		//	Function to divide the data values used in the threaded GSL Norm function
		static vector<vector<double*> > divideDataNormalise( vector<double*> input, int subsets=1 );

	private:

		//	Cannot Construct this class, it's simply a collection of static methods
		Threading();
		~Threading();
};

#endif

