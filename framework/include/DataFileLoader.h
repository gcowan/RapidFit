/**
        @class DataFileLoader

        Holds methods for loading given file formats into an IDataSet

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef DATA_FILE_LOADER_H
#define DATA_FILE_LOADER_H

#include "IDataSet.h"
#include "PhaseSpaceBoundary.h"
#include "TNtuple.h"
#include <string>
#include <iostream>

using namespace std;

class DataFileLoader
{
	public:
		DataFileLoader();
		DataFileLoader( string, PhaseSpaceBoundary*, long );
		~DataFileLoader();

		IDataSet * GetDataSet();
	private:
		void LoadAsciiFileIntoMemory( string, long );
		void LoadRootFileIntoMemory( string, string, long );
		bool CheckTNtupleWithBoundary( TNtuple *, PhaseSpaceBoundary * );
		IDataSet * data;
};

#endif
