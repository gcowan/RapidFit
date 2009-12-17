/**
  @class DataSetConfiguration

  A class for holding the data to create a data set

  @author Benjamin Wynne
  @data 2009-12-16
  */

#ifndef DATA_SET_CONFIGURATION_H
#define DATA_SET_CONFIGURATION_H

#include "TNtuple.h"
#include "IPDF.h"
#include "IDataSet.h"
#include <string>
#include <vector>

class DataSetConfiguration
{
	public:
		DataSetConfiguration();
		DataSetConfiguration( string, long, vector<string>, vector<string> );
		DataSetConfiguration( string, long, vector<string>, vector<string>, IPDF* );
		~DataSetConfiguration();

		bool SetPhysicsParameters( ParameterSet* );
		IDataSet * MakeDataSet( PhaseSpaceBoundary*, IPDF* );

	private:
		IDataSet * LoadDataFile( vector<string>, vector<string>, PhaseSpaceBoundary*, long );
		IDataSet * LoadAsciiFileIntoMemory( string, long, PhaseSpaceBoundary* );
		IDataSet * LoadRootFileIntoMemory( string, string, long, PhaseSpaceBoundary* );
		bool CheckTNtupleWithBoundary( TNtuple *, PhaseSpaceBoundary * );

		string source;
		long numberEvents;
		vector<string> arguments, argumentNames;
		IPDF * generatePDF;
		bool separateGeneratePDF, parametersAreSet;
};

#endif
