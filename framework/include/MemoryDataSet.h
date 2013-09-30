/**
  @class MemoryDataSet

  A data set which simply stores a vector of pointers to datapoint objects

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

#pragma once
#ifndef MEMORY_DATA_SET_H
#define MEMORY_DATA_SET_H

//	RapidFit Headers
#include "IDataSet.h"
#include "DataPoint.h"
//	System Headers
#include <vector>

using namespace::std;

class MemoryDataSet : public IDataSet
{
	public:
		MemoryDataSet( PhaseSpaceBoundary*, vector<DataPoint*> );
		MemoryDataSet( PhaseSpaceBoundary* );
		~MemoryDataSet();

		//Interface functions
		virtual DataPoint * GetDataPoint(int) const;
		virtual bool AddDataPoint( DataPoint* );
		virtual int GetDataNumber( DataPoint* templateDataPoint =NULL ) const;
		virtual PhaseSpaceBoundary * GetBoundary() const;
		virtual void SetBoundary( const PhaseSpaceBoundary* );

		virtual void SortBy( string );

		virtual IDataSet* GetDiscreteDataSet( const vector<ObservableRef> discreteParam, const vector<double> discreteVal ) const;

		virtual vector<DataPoint*> GetDiscreteSubSet( const vector<ObservableRef> discreteParam, const vector<double> discreteVal ) const;
		virtual vector<DataPoint*> GetDiscreteSubSet( const vector<string> discreteParam, const vector<double> discreteVal ) const;

		void Clear();

		/*!                      
		 * @brief Returns an estimate of the total Yield
		 */                      
		double Yield();
		double YieldError();

		string GetWeightName() const;
		bool GetWeightsWereUsed() const;
		void UseEventWeights( const string Name );

		double GetSumWeights() const;
		double GetSumWeightsSq() const;
		void ApplyAlpha( const double, const double );
		double GetAlpha() const;

		void ApplyExternalAlpha( const string alphaName );

		void NormaliseWeights();

		virtual void Print() const;

		void ClearAllPseudoObservables();

		void SafeAddDataPoint( DataPoint* NewDataPoint );

	private:
		//	Uncopyable!
		MemoryDataSet ( const MemoryDataSet& );
		MemoryDataSet& operator = ( const MemoryDataSet& );
		vector<DataPoint*> allData;
		PhaseSpaceBoundary * dataBoundary;
		mutable vector<int> allSubSets;

		bool useWeights;
		string WeightName;
		double alpha;
		string alphaName;

		bool canDelete;
};

bool compare_datapoints ( pair<DataPoint,ObservableRef> first, pair<DataPoint,ObservableRef> second );

#endif

