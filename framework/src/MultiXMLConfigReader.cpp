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
#include "StringProcessing.h"
//	System Headers
#include <vector>
#include <string>
#include <sstream>

using namespace::std;

MultiXMLConfigReader::MultiXMLConfigReader( vector<string> fileNames, vector<pair<string, string> >* OverrideXML ) : I_XMLConfigReader(), XMLReaders(), storedSeed(-1), storedRepeats(-1), _fileNames( fileNames )
{
	for( unsigned int i=0; i< fileNames.size(); ++i )
	{
		TString thisNum; thisNum+=i+1;
		XMLReaders.push_back( new XMLConfigReader( fileNames[i], OverrideXML, string(thisNum.Data()) ) );
	}
	if( DebugClass::DebugThisClass( "XMLConfigReader_TAG" ) ) if( !DebugClass::GetClassNames().empty() ) exit(0);
}

MultiXMLConfigReader::~MultiXMLConfigReader()
{
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
	vector<string> total_XML;

	for( unsigned int i=0; i< XMLReaders.size(); ++i )
	{
		stringstream thisStream;
		thisStream << "\n##########   " << _fileNames[i] << "   ##########\n\n\n";
		total_XML.push_back( thisStream.str() );

		vector<string> thisFile = XMLReaders[i]->GetXML();

		for( unsigned int j=0; j< thisFile.size(); ++j )
		{
			total_XML.push_back( thisFile[j] );
		}

		total_XML.push_back( "\n\n\n##########   EOF   ##########\n\n\n" );
	}

	return total_XML;
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
		vector<PDFWithData*> thisSet;
		try
		{
			thisSet = XMLReaders[i]->GetPDFsAndData();
		}
		catch(...)
		{
			thisSet = vector<PDFWithData*>();
		}

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
		vector<ConstraintFunction*> thisSet;
		try
		{
			thisSet = XMLReaders[i]->GetConstraints();
		}
		catch(...)
		{
			thisSet = vector<ConstraintFunction*>();
		}

		for( unsigned int j=0; j< thisSet.size(); ++j )
		{
			vector<string> knownNames;
			for( unsigned int k=0; k< fullSet.size(); ++k )
			{
				vector<string> allNames_Full = fullSet[k]->GetConstraintNames();
				for( unsigned int a=0; a< allNames_Full.size(); ++a )
				{
					knownNames.push_back( allNames_Full[a] );
					//cout << "known: " << allNames_Full[a] << endl;
				}
			}

			vector<string> theseNames = thisSet[j]->GetConstraintNames();
			bool good = true;
			for( unsigned int k=0; k< theseNames.size(); ++k )
			{
				if( StringProcessing::VectorContains( &knownNames, &(theseNames[k]) ) != -1 )
				{
					cout << "Rejecting multiple Constraint: " << theseNames[k] << endl;
					good = false;
					//break;
				}
			}
			if( good )
			{
				fullSet.push_back( thisSet[j] );
			}
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
	for( unsigned int i=0; i< XMLReaders.size(); ++i )
	{
		XMLReaders[i]->SetSeed( new_seed );
	}
}

unsigned int MultiXMLConfigReader::GetOriginalSeed() const
{
	return XMLReaders[0]->GetOriginalSeed();
}

