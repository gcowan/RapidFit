/**
  @class DataSetConfiguration

  A class for holding the data to create a data set

  @author Benjamin Wynne
  @data 2009-12-16
  */

#ifndef DATA_SET_CONFIGURATION_H
#define DATA_SET_CONFIGURATION_H

//	ROOT Headers
#include "TNtuple.h"
//	RapidFit Headers
#include "IPDF.h"
#include "IDataSet.h"
//	System Headers
#include <string>
#include <vector>

class DataSetConfiguration
{
	public:
		DataSetConfiguration();
		DataSetConfiguration( string, long, string, vector<string>, vector<string> );
		DataSetConfiguration( string, long, string, vector<string>, vector<string>, IPDF* );
		~DataSetConfiguration();

		bool SetPhysicsParameters( vector<ParameterSet*> );
		bool SetSource( string );
		IDataSet * MakeDataSet( PhaseSpaceBoundary*, IPDF*, int=-1 );
		IPDF* GetGenerationPDF();

	private:
		//	Uncopyable!
		DataSetConfiguration ( const DataSetConfiguration& );
		DataSetConfiguration& operator = ( const DataSetConfiguration& );
		IDataSet * LoadDataFile( vector<string>, vector<string>, PhaseSpaceBoundary*, long );
		IDataSet * LoadAsciiFileIntoMemory( string, long, PhaseSpaceBoundary* );
		IDataSet * LoadRootFileIntoMemory( string, string, long, PhaseSpaceBoundary* );
		bool CheckTNtupleWithBoundary( TNtuple *, PhaseSpaceBoundary * );

		string source;
		string cutString;
		long numberEvents;
		vector<string> arguments, argumentNames;
		IPDF * generatePDF;
		bool separateGeneratePDF, parametersAreSet;
};

#endif
