/**
        @class RootFileDataSet

        A collection of data points stored in a Root NTuple

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef ROOT_FILE_DATA_SET_H
#define ROOT_FILE_DATA_SET_H

#include "TFile.h"
#include <string>
#include <iostream>
#include "TNtuple.h"
#include "IDataSet.h"
#include "TBranch.h"

using namespace std;

class RootFileDataSet : public IDataSet
{
	public:
		RootFileDataSet();
		RootFileDataSet( string, PhaseSpaceBoundary* );
		RootFileDataSet( string, string, PhaseSpaceBoundary* );
		~RootFileDataSet();

		//Interface functions
		virtual DataPoint * GetDataPoint(int);
		virtual bool AddDataPoint( DataPoint* );
		virtual int GetDataNumber();
		virtual PhaseSpaceBoundary * GetBoundary();
	
		bool CheckTNtupleWithBoundary( TNtuple*, PhaseSpaceBoundary* );
		bool SetBranches();
	
	private:
		TFile * inputFile;
		TNtuple * rootNTuple;
		PhaseSpaceBoundary * dataBoundary;
		vector<TBranch *> branches;
		vector<Float_t> observableValues;
		DataPoint * outputDataPoint;
};

#endif
