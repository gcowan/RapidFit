//	RapidFit Headers
#include "XMLTag.h"
#include "ParameterSet.h"
#include "PDFWithData.h"
#include "FitFunctionConfiguration.h"
#include "MinimiserConfiguration.h"
#include "OutputConfiguration.h"
#include "DataSetConfiguration.h"
#include "ComponentPlotter.h"
#include "ConstraintFunction.h"
#include "ScanParam.h"
#include "PrecalculatorConfig.h"
#include "DebugClass.h"
#include "MultiXMLConfigReader.h"
#include "XMLConfigReader.h"
#include "I_XMLConfigReader.h"
//	System Headers
#include <vector>
#include <string>

using namespace::std;

MultiXMLConfigReader::MultiXMLConfigReader( vector<string> fileNames, DebugClass* thisDebug ) : I_XMLConfigReader(), XMLReaders(), storedSeed(-1), storedRepeats(-1)
{
	for( unsigned int i=0; i< fileNames.size(); ++i )
	{
		XMLReaders.push_back( new XMLConfigReader( fileNames[i], NULL, thisDebug ) );
	}
	debug = new DebugClass( thisDebug );
}

MultiXMLConfigReader::~MultiXMLConfigReader()
{
	if( debug != NULL ) delete debug;
	while( !XMLReaders.empty() )
	{
		if( XMLReaders.back() != NULL ) delete XMLReaders.back();
		XMLReaders.pop_back();
	}
}

bool MultiXMLConfigReader::IsValid() const
{
	bool test = true;
	for( unsigned int i=0; i< XMLReaders.size(); ++i )	test=test&&XMLReaders[i]->IsValid();
	return test;
}

vector<string> MultiXMLConfigReader::GetXML() const
{
	return vector<string>();
}

ParameterSet* MultiXMLConfigReader::GetFitParameters( vector<string> CommandLineParam )
{
	vector<ParameterSet*> allSets;
	for( unsigned int i=0; i< XMLReaders.size(); ++i )	allSets.push_back( XMLReaders[i]->GetFitParameters( CommandLineParam ) );

	ParameterSet* thisSet = new ParameterSet( allSets, true );

	return thisSet;
}

MinimiserConfiguration* MultiXMLConfigReader::GetMinimiserConfiguration()
{
	if( XMLReaders.empty() ) return NULL;
	else return XMLReaders[0]->GetMinimiserConfiguration();
}

FitFunctionConfiguration* MultiXMLConfigReader::GetFitFunctionConfiguration()
{
	if( XMLReaders.empty() ) return NULL;
	else return XMLReaders[0]->GetFitFunctionConfiguration();
}

OutputConfiguration* MultiXMLConfigReader::GetOutputConfiguration()
{
	if( XMLReaders.empty() ) return NULL;
	else return XMLReaders[0]->GetOutputConfiguration();
}

vector<PDFWithData*> MultiXMLConfigReader::GetPDFsAndData( vector<int> StartingValues )
{
	(void) StartingValues;
	vector<PDFWithData*> fullSet;
	for( unsigned int i=0; i< XMLReaders.size(); ++i )
	{
		vector<PDFWithData*> thisSet = XMLReaders[i]->GetPDFsAndData();
		for( unsigned int j=0; j< thisSet.size(); ++j )
		{
			fullSet.push_back( thisSet[j] );
		}
	}
	return fullSet;
}

vector<ConstraintFunction* > MultiXMLConfigReader::GetConstraints()
{
	vector<ConstraintFunction*> fullSet;
	for( unsigned int i=0; i< XMLReaders.size(); ++i )
	{
		vector<ConstraintFunction*> thisSet = XMLReaders[i]->GetConstraints();
		for( unsigned int j=0; j< thisSet.size(); ++j )
		{
			fullSet.push_back( thisSet[j] );
		}
	}
	return fullSet;
}

vector<PhaseSpaceBoundary*> MultiXMLConfigReader::GetPhaseSpaceBoundaries()
{
	vector<PhaseSpaceBoundary*> fullSet;
	for( unsigned int i=0; i< XMLReaders.size(); ++i )
	{
		vector<PhaseSpaceBoundary*> thisSet = XMLReaders[i]->GetPhaseSpaceBoundaries();
		for( unsigned int j=0; j< thisSet.size(); ++j )
		{
			fullSet.push_back( thisSet[j] );
		}
	}
	return fullSet;
}

PrecalculatorConfig* MultiXMLConfigReader::GetPrecalculatorConfig()
{
	if( XMLReaders.empty() ) return NULL;
	else return XMLReaders[0]->GetPrecalculatorConfig();
}

vector<int> MultiXMLConfigReader::GetAllDataSetSizes()
{
	cout << "MultiXMLConfigReader DOES NOT SUPPORT THIS FUNCTION YET!!!" << endl;
	return vector<int>();
}

vector<int> MultiXMLConfigReader::GetAllStartEntries()
{
	cout << "MultiXMLConfigReader DOES NOT SUPPORT THIS FUNCTION YET!!!" << endl;
	return vector<int>();
}

int MultiXMLConfigReader::GetNumberRepeats()
{
	if( storedRepeats > 0 ) return storedRepeats;
	else return XMLReaders[0]->GetNumberRepeats();
}

unsigned int MultiXMLConfigReader::GetSeed()
{
	if( storedSeed >= 0 ) return (unsigned)storedSeed;
	else return XMLReaders[0]->GetSeed();
}

void MultiXMLConfigReader::SetSeed( unsigned int new_seed )
{
	storedSeed = new_seed;
}

void MultiXMLConfigReader::SetDebug( DebugClass* input_debug )
{
	for( unsigned int i=0; i< XMLReaders.size(); ++i )
	{
		XMLReaders[i]->SetDebug( input_debug );
	}
        if( debug != NULL ) delete debug;
        debug = new DebugClass( *input_debug );
}

