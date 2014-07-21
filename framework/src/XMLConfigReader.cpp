/**
  @class XMLConfigReader

  Opens an xml config file and uses it to create RapidFit data objects

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

//	RapidFit Headers
#include "I_XMLConfigReader.h"
#include "XMLConfigReader.h"
#include "ClassLookUp.h"
#include "SumPDF.h"
#include "NormalisedSumPDF.h"
#include "ProdPDF.h"
#include "StringProcessing.h"
#include "AcceptReject.h"
#include "IConstraint.h"
#include "ObservableContinuousConstraint.h"
#include "ObservableDiscreteConstraint.h"
#include "Blinder.h"
#include "ScanParam.h"
#include "PDFConfigurator.h"
#include "RapidFitIntegrator.h"
#include "RapidRun.h"
#include "XMLObjectGenerator.h"
//	System Headers
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <float.h>
#include <vector>

using namespace::std;

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6


void XMLConfigReader::PrintTag( const XMLTag* thisTag )
{
	cout << "PATH:\t" << thisTag->GetPath() << endl;
	cout << "NAME:\t" << thisTag->GetName() << endl;
	cout << "VALUE:\t" << endl;

	vector<string> values = thisTag->GetRAWValue();

	for( unsigned int i=0; i< values.size(); ++i )
	{
		cout << values[i] << endl;
	}
	cout << endl;

	vector<XMLTag*> theseChildren = thisTag->GetChildren();
	for( unsigned int i=0; i< theseChildren.size(); ++i )
	{
		cout << "CHILD:\t" << i << "\t" << endl;
		XMLConfigReader::PrintTag( theseChildren[i] );
		cout << endl;
	}
}

//Constructor with file name argument
XMLConfigReader::XMLConfigReader( string FileName, vector<pair<string, string> >* OverrideXML, string GlobalPrepend ) : I_XMLConfigReader(),
	fileName( FileName ), fileTags(), wholeFile(), All_XML_Tags(NULL), children(), seed(-1), XMLValid(false), _original_seed(-999)
{
	All_XML_Tags = new XMLTag(OverrideXML, GlobalPrepend);

	//Open the config file
	ifstream configFile( FileName.c_str() );
	if ( !configFile.is_open() )
	{
		cerr << "Failed to open config file \"" << FileName << "\"" << endl;
		exit(1);
	}

	//Read the whole file into a vector
	while( configFile.good() )
	{
		string newLine;
		getline( configFile, newLine );
		wholeFile.push_back( newLine );
	}

	StringProcessing::RemoveWhiteSpace(wholeFile);
	StringProcessing::RemoveLeadingWhiteSpace(wholeFile);

	vector<string> value;
	vector<XMLTag*> File_Tags = All_XML_Tags->FindTagsInContent( wholeFile, value );
	children = File_Tags[0]->GetChildren();
	fileTags = File_Tags;

	XMLValid = this->TestXML();

	unsigned int PDFNum_Counter=0;
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "ToFit" )
		{
			for( unsigned int i=0; i< children[childIndex]->GetChildren().size(); ++i )
			{
				if( children[childIndex]->GetChildren()[i]->GetName() == "PDF" || children[childIndex]->GetChildren()[i]->GetName() == "CommonPDF" )
				{
					TString ThisNum=""; ThisNum+=PDFNum_Counter;
					children[childIndex]->AppendPath( ThisNum.Data() );
					++PDFNum_Counter;
					break;
				}
			}
		}
	}

	cout << endl << endl << endl << "Finished Updating PDF Paths" << endl << endl << endl;

	if( DebugClass::DebugThisClass( "XMLConfigReader_TAG" ) )
	{
		for( unsigned int childIndex=0; childIndex<children.size(); ++childIndex )
		{
			XMLConfigReader::PrintTag( children[childIndex] );
		}
		if( !DebugClass::GetClassNames().empty() && GlobalPrepend.empty() ) exit(0);
	}
}

bool XMLConfigReader::IsValid() const
{
	return XMLValid;
}

bool XMLConfigReader::TestXML()
{
	if( fileTags.size() != 1 )
	{
		cout << "XMLConfigReader: Possible Error! your XML should contain only 1 top level <RapidFit> tag!" << endl;
		cout << "XMLConfigReader: Everything must be included within these tags." << endl;
		return false;
	}
	else
	{
		unsigned int numParamSets=0;
		unsigned int numToFits=0;
		unsigned int numMinimisers=0;
		unsigned int numFitFunctions=0;
		unsigned int numOutputs=0;
		unsigned int numPrecalculators=0;
		unsigned int numSeeds=0;
		unsigned int numRepeats=0;
		unsigned int numCommonPDFs=0;
		unsigned int numCommonPhaseSpaces=0;
		for( unsigned int i=0; i< children.size(); ++i )
		{
			string name = children[i]->GetName();
			if( name == "ParameterSet" ) ++numParamSets;
			else if( name == "Minimiser" ) ++numMinimisers;
			else if( name == "FitFunction" ) ++numFitFunctions;
			else if( name == "ToFit" ) ++numToFits;
			else if( name == "Output" ) ++numOutputs;
			else if( name == "Precalculator" ) ++numPrecalculators;
			else if( name == "Seed" ) ++numSeeds;
			else if( name == "NumberRepeats" ) ++numRepeats;
			else if( name == "CommonPDF" ) ++numCommonPDFs;
			else if( name == "CommonPhaseSpace" ) ++numCommonPhaseSpaces;
			else
			{
				cout << endl << endl;
				cout << "!!!WARNING!!!" << endl;
				cout << "XMLConfigReader: Possible Error, Tag: <" << name << "> unknown. Ignoring!" << endl;
				cout << "!!!WARNING!!!" << endl;
				cout << endl << endl;
			}
		}
		bool returnable=true;
		if( numParamSets != 1 )
		{
			if( numParamSets == 0 )
			{
				cout << "XMLConfigReader: You need to include a <ParameterSet> in your XML!" << endl;
			}
			if( numParamSets > 1 )
			{
				cout << "XMLConfigReader: You need to have just 1 <ParameterSet> in your XML" << endl;
			}
			returnable=false;
		}
		if( numToFits == 0 )
		{
			cout << "XMLConfigReader: You Must include at least 1 <ToFit> segment in your XML" << endl;
			returnable=false;
		}
		if( numMinimisers != 1 )
		{
			if( numMinimisers == 0 )
			{
				cout << "XMLConfigReader: You need to include a <Minimiser> in your XML!" << endl;
			}
			if( numMinimisers > 1 )
			{
				cout << "XMLConfigReader: You need to have just 1 <Minimiser> in your XML" << endl;
			}
			returnable=false;
		}
		if( numFitFunctions != 1 )
		{
			if( numFitFunctions == 0 )
			{
				cout << "XMLConfigReader: You need to include a <FitFunction> in your XML!" << endl;
			}
			if( numFitFunctions > 1 )
			{
				cout << "XMLConfigReader: You need to have just 1 <FitFunction> in your XML" << endl;
			}
			returnable=false;
		}
		if( numOutputs > 1 )
		{
			cout << "XMLConfigReader: You need to have at most 1 <Output> in your XML" << endl;
			returnable=false;
		}
		if( numPrecalculators > 1 )
		{
			cout << "XMLConfigReader: You need to have at most 1 <Precalculator> in your XML" << endl;
			returnable=false;
		}
		if( numSeeds > 1 )
		{
			cout << "XMLConfigReader: You need to have at most 1 <Seed> in your XML" << endl;
			returnable=false;
		}
		if( numRepeats > 1 )
		{
			cout << "XMLConfigReader: You need to have at most 1 <NumberRepeat> in your XML" << endl;
			returnable=false;
		}
		if( numCommonPDFs > 1 )
		{
			cout << "XMLConfigReader: You need to have at most 1 <CommonPDF> in your XML" << endl;
		}
		if( numCommonPhaseSpaces > 1 )
		{
			cout << "XMLConfigReader: You need to have at most 1 <CommonPhaseSpace> in your XML" << endl;
		}
		return returnable;
	}
	return true;
}

//Destructor
XMLConfigReader::~XMLConfigReader()
{
	if( All_XML_Tags != NULL ) delete All_XML_Tags;
	//cout << "Hello from XMLConfigReader destructor" << endl;
	/*	while( !fileTags.empty() )
		{
		if( fileTags.back() != NULL ) delete fileTags.back();
		fileTags.pop_back();
		}
		while( !children.empty() )
		{
		if( children.back() != NULL ) delete children.back();
		children.pop_back();
		}*/
}

vector<string> XMLConfigReader::GetXML() const
{
	return wholeFile;
}

ParameterSet* XMLConfigReader::GetFitParameters( vector<string> CommandLineParam )
{
	ParameterSet* RawParameters = this->GetRawFitParameters();

	if( CommandLineParam.empty() ) return RawParameters;

	for( unsigned int i=0; i < CommandLineParam.size() ; ++i )
	{
		vector<string> input_args = StringProcessing::SplitString( CommandLineParam[i], ',' );
		if( input_args.size() != 6 )
		{
			cerr << "Cannot understand Physics Parameter info you passed at runtime" << endl;
			cerr << "They should be defined as:\t\tgamma,value,min,max,stepsize,type" << endl;
			exit(-598);
		}
		bool found_parameter=false;

		bool local_find = false;
		PhysicsParameter* try_get_param=NULL;
		try{
			try_get_param = RawParameters->GetPhysicsParameter( input_args[0] );
			local_find = true;
		}
		catch( int e )
		{
			local_find = false;
		}

		if( local_find )
		{
			string Unit = try_get_param->GetUnit();

			double step = strtod(input_args[4].c_str(),NULL);
			double max = strtod(input_args[3].c_str(),NULL);
			double min = strtod(input_args[2].c_str(),NULL);
			double val = strtod(input_args[1].c_str(),NULL);
			string type = input_args[5];

			if( ( ( fabs( min - max ) < 	1E-6 ) && ( max < 1E-6 ) ) || (type=="Fixed") )
			{
				//		Unbounded				name	value	step	type	unit
				PhysicsParameter* new_param = new PhysicsParameter( input_args[0], val , step, type, Unit);
				RawParameters->SetPhysicsParameter( input_args[0], new_param );
			}
			else
			{
				//					name	val	min	max	step	type	unit
				RawParameters->SetPhysicsParameter( input_args[0], val, min, max, step, type, Unit );
			}
			found_parameter = local_find;
		}

		if( !found_parameter )
		{
			cerr << "Couldn't find parameter: " << input_args[0] << " exiting!" << endl;
			exit(-345);
		}
	}


	return RawParameters;
}

ParameterSet* XMLConfigReader::GetRawFitParameters()
{
	vector<ParameterSet*> All_Parameters;
	//Find the ParameterSet tag
	for( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "ParameterSet" )
		{
			All_Parameters.push_back( XMLObjectGenerator::GetParameterSet( children[childIndex] ) );
		}
	}
	if( All_Parameters.empty() ) All_Parameters.push_back( new ParameterSet( vector<string>()) );

	ParameterSet* output = new ParameterSet( All_Parameters );
	for( vector<ParameterSet*>::iterator param_i = All_Parameters.begin(); param_i != All_Parameters.end(); ++param_i )
	{
		delete *param_i;
	}
	return output;
}

//Return the minimiser for the fit
MinimiserConfiguration * XMLConfigReader::GetMinimiserConfiguration()
{
	//Find the Minimiser tag
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "Minimiser" )
		{
			return XMLObjectGenerator::MakeMinimiser( children[childIndex], GetOutputConfiguration() );
		}
	}

	//If no such tag is found, fail
	cerr << "Minimiser tag not found in config file" << endl;
	exit(1);
}


//Return the output configuration for the fit
OutputConfiguration * XMLConfigReader::GetOutputConfiguration()
{
	//Find the Output tag
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "Output" )
		{
			return XMLObjectGenerator::MakeOutputConfiguration( children[childIndex] );
		}
	}

	//If no such tag is found, make default
	//cout << "Output tag not found in config file - using default" << endl;
	return new OutputConfiguration( vector<pair<string,string> >(), string("None"), vector<ScanParam*>(), vector<pair<ScanParam*,ScanParam*> >(), vector<CompPlotter_config*>()  );
}

//Return the number of repeats for the fit
int XMLConfigReader::GetNumberRepeats()
{
	//Find the NumberRepeats tag
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "NumberRepeats" )
		{
			return XMLTag::GetIntegerValue( children[childIndex] );
		}
	}

	//If no such tag is found, fail
	cerr << "NumberRepeats tag not found in config file" << endl;
	return 1;
}

//Return the function to minimise
FitFunctionConfiguration * XMLConfigReader::GetFitFunctionConfiguration()
{
	//Find the FitFunction tag
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "FitFunction" )
		{
			return XMLObjectGenerator::MakeFitFunction( children[childIndex] );
		}
	}

	//If no such tag is found, fail
	cerr << "FitFunction tag not found in config file" << endl;
	exit(1);
}

//Organise all the PDFs and DataSets
vector< PDFWithData* > XMLConfigReader::GetPDFsAndData( vector<int> Starting_Value )
{
	//Collect all ToFit elements
	vector< XMLTag* > toFits;
	vector< PDFWithData* > pdfsAndData;
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "ToFit" )
		{
			toFits.push_back( children[childIndex] );
		}
	}

	//Go through the collected ToFit elements
	if ( toFits.size() == 0 )
	{
		cerr << "No ToFit tags found in config file" << endl;
		throw(1);
	}
	else
	{
		unsigned int PDF_NUM=0;
		//Loop over all ToFits
		for ( unsigned int fitIndex = 0; fitIndex < toFits.size(); ++fitIndex )
		{
			XMLTag * pdfTag=NULL;
			XMLTag * dataTag=NULL;
			XMLTag* overloadConfigurator=NULL;
			bool foundPDF = false;
			bool foundData = false;
			bool foundConstraint = false;

			//Find the PDF and data set configuration
			vector< XMLTag* > fitComponents = toFits[fitIndex]->GetChildren();
			for ( unsigned int componentIndex = 0; componentIndex < fitComponents.size(); ++componentIndex )
			{
				string name = fitComponents[componentIndex]->GetName();
				if ( name == "PDF" || name == "SumPDF" || name == "NormalisedSumPDF" || name == "ProdPDF" )
				{
					pdfTag = fitComponents[componentIndex];
					foundPDF = true;
				}
				else if ( name == "DataSet" )
				{
					dataTag = fitComponents[componentIndex];
					foundData = true;
				}
				else if ( name == "ConstraintFunction" )
				{
					//Just so can check for  bad ToFit tags
					foundConstraint = true;
				}
				else if ( name == "CommonPDF" )
				{
					string value = XMLTag::GetStringValue( fitComponents[componentIndex] );
					if( value == "True" )
					{
						pdfTag = this->FindCommonPDFXML();
						if( pdfTag == NULL ) foundPDF=false;
						else foundPDF=true;
					}
				}
				else if( name == "PDFConfigurator" )
				{
					overloadConfigurator = fitComponents[componentIndex];
				}
				else
				{
					cerr << "Unrecognised fit component: " << name << endl;
					exit(1);
				}
			}


			//Make the data set, and populate the physics bottle
			if( foundData && foundPDF )
			{
				if( !Starting_Value.empty() )
				{
					int s_val_index = (int) pdfsAndData.size(); ++s_val_index;
					XMLTag* common = this->FindCommonPDFXML();
					ParameterSet* thisSet = this->GetFitParameters();
					PhaseSpaceBoundary* thisPhaseSpaceBoundary = this->FindCommonPhaseSpace();
					pdfsAndData.push_back( XMLObjectGenerator::GetPDFWithData( dataTag, pdfTag, Starting_Value[(unsigned)s_val_index], overloadConfigurator, common, thisSet, thisPhaseSpaceBoundary ) );
					delete thisPhaseSpaceBoundary;
					delete thisSet;
				} else {
					XMLTag* common = this->FindCommonPDFXML();
					ParameterSet* thisSet = this->GetFitParameters();
					PhaseSpaceBoundary* thisPhaseSpaceBoundary = this->FindCommonPhaseSpace();
					pdfsAndData.push_back( XMLObjectGenerator::GetPDFWithData( dataTag, pdfTag, -1, overloadConfigurator, common, thisSet, thisPhaseSpaceBoundary ) );
					delete thisPhaseSpaceBoundary;
					delete thisSet;
				}
				TString ThisLabel = pdfsAndData.back()->GetPDF()->GetLabel();
				ThisLabel.Append( "_PDF" ); ThisLabel+=PDF_NUM; ++PDF_NUM;
				pdfsAndData.back()->GetPDF()->SetLabel( ThisLabel.Data() );
			}
			else if( !foundConstraint )
			{
				cerr << "A ToFit xml tag is incomplete:";
				if ( !foundData && foundPDF )
				{
					cerr << " Data set configuration missing";
				}
				else if ( !foundPDF && foundData )
				{
					cerr << " PDF configuration missing";
				}
				else
				{
					cerr << " No configuration found";
				}
				cerr << endl;
				exit(1);
			}
		}
	}

	return pdfsAndData;
}

//Organise all the PDFs and DataSets
vector< ConstraintFunction* > XMLConfigReader::GetConstraints()
{
	//Collect all ToFit elements
	vector< XMLTag* > toFits;
	vector< ConstraintFunction* > constraints;
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "ToFit" )
		{
			toFits.push_back( children[childIndex] );
		}
	}

	//Go through the collected ToFit elements
	if ( toFits.size() == 0 )
	{
		cerr << "No ToFit tags found in config file(2)" << endl;
		throw(1);
	}
	else
	{
		//Loop over all ToFits
		for ( unsigned int fitIndex = 0; fitIndex < toFits.size(); ++fitIndex )
		{
			XMLTag* constraintTag=NULL;
			bool foundConstraint = false;

			//Find the PDF and data set configuration
			vector< XMLTag* > fitComponents = toFits[fitIndex]->GetChildren();
			for ( unsigned int componentIndex = 0; componentIndex < fitComponents.size(); ++componentIndex )
			{
				if ( fitComponents[componentIndex]->GetName() == "ConstraintFunction" )
				{
					constraintTag = fitComponents[componentIndex];
					foundConstraint = true;
				}
			}

			if (foundConstraint)
			{
				constraints.push_back( XMLObjectGenerator::GetConstraintFunction(constraintTag) );
			}
		}
	}

	return constraints;
}

//  Return the PhaseSpaceBoundaries defined in the XML
//  I don't want to add bloat and store the PhaseSpaceBoudaries as the file is loaded,
//  this just leads to storing lots of not wanted data within this class
//  I try to make use of the structure correctly
vector<PhaseSpaceBoundary*> XMLConfigReader::GetPhaseSpaceBoundaries()
{
	vector<PhaseSpaceBoundary*> PhaseSpaceBoundary_vec;
	//Find the ParameterSets
	vector< XMLTag* > toFits;
	//  Children Defined globally on initialization!
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "ToFit" )
		{
			toFits.push_back( children[childIndex] );
		}
	}
	//  Loop over all of the fits
	for ( unsigned short int fitnum=0; fitnum < toFits.size(); ++fitnum )
	{
		vector<XMLTag*> allChildren = toFits[fitnum]->GetChildren();
		for( unsigned short int Childnum=0; Childnum < allChildren.size(); ++Childnum )
		{
			if( allChildren[Childnum]->GetName() == "DataSet" )
			{
				vector<XMLTag*> allFitDataSets = allChildren[Childnum]->GetChildren();
				for( unsigned short int DataSetNum=0; DataSetNum < allFitDataSets.size(); ++DataSetNum )
				{
					if( allFitDataSets[DataSetNum]->GetName() == "PhaseSpaceBoundary" )
					{
						PhaseSpaceBoundary_vec.push_back( XMLObjectGenerator::GetPhaseSpaceBoundary( allFitDataSets[DataSetNum] ) );
					}
				}
			}

		}
	}
	return PhaseSpaceBoundary_vec;
}

PhaseSpaceBoundary* XMLConfigReader::FindCommonPhaseSpace( XMLTag* PhaseTag )
{
	vector<XMLTag*> phaseChildren;
	if( PhaseTag != NULL ) phaseChildren = PhaseTag->GetChildren();

	PhaseSpaceBoundary* returnablePhaseSpace = NULL;

	for( unsigned int childIndex=0; childIndex< children.size(); ++childIndex )
	{
		if( children[childIndex]->GetName() == "CommonPhaseSpace" )
		{
			returnablePhaseSpace = XMLObjectGenerator::GetPhaseSpaceBoundary( children[childIndex]->GetChildren()[0] );
		}
	}

	if( returnablePhaseSpace == NULL ) returnablePhaseSpace = new PhaseSpaceBoundary( vector<string>() );

	for( unsigned int phaseChild=0; phaseChild< phaseChildren.size(); ++phaseChild )
	{
		if( phaseChildren[phaseChild]->GetName() == "Observable" )
		{
			string temp; (void) temp;
			IConstraint* new_observable = XMLObjectGenerator::GetConstraint( phaseChildren[phaseChild], temp );
			//new_observable->Print();
			//	Add Constraint if non exists and Overwrite if it already exists
			returnablePhaseSpace->AddConstraint( new_observable->GetName(), new_observable, true );
			delete new_observable;
		}
	}

	return returnablePhaseSpace;
}

/*
   IPDF* XMLConfigReader::FindCommonPDF( PhaseSpaceBoundary* override )
   {
   for( unsigned int childIndex=0; childIndex< children.size(); ++childIndex )
   {
   if( children[childIndex]->GetName() == "CommonPDF" )
   {
   if( override == NULL )
   {
   return XMLObjectGenerator::GetPDF( children[childIndex]->GetChildren()[0], this->FindCommonPhaseSpace(), NULL );
   }
   else
   {
   return XMLObjectGenerator::GetPDF( children[childIndex]->GetChildren()[0], override, NULL );
   }
   }
   }
   return NULL;
   }
   */

XMLTag* XMLConfigReader::FindCommonPDFXML()
{
	for( unsigned int childIndex=0; childIndex< children.size(); ++childIndex )
	{
		if( children[childIndex]->GetName() == "CommonPDF" )
		{
			unsigned int subParamNum=0;
			unsigned int appendParamNum=0;
			unsigned int configParamNum=0;
			for( unsigned int i=0; i< children[childIndex]->GetChildren().size(); ++i )
			{
				AppendCommonNames( children[childIndex]->GetChildren()[i], subParamNum, appendParamNum, configParamNum );
			}

			return children[childIndex]->GetChildren()[0];
		}
	}
	return NULL;
}

void XMLConfigReader::AppendCommonNames( XMLTag* input, unsigned int &subParamNum, unsigned int &appendParamNum, unsigned int &configParamNum )
{
	for( unsigned int i=0; i< input->GetChildren().size(); ++i )
	{
		if( input->GetChildren()[i]->GetName() == "ParameterSubstitution" )
		{
			TString ThisNum("Common"); ThisNum+=subParamNum;
			input->GetChildren()[i]->AppendPath( ThisNum.Data() );
			++subParamNum;
		}
		else if ( input->GetChildren()[i]->GetName() == "AppendParameterNames" )
		{
			TString ThisNum("Common"); ThisNum+=appendParamNum;
			input->GetChildren()[i]->AppendPath( ThisNum.Data() );
			++appendParamNum;
		}
		else if ( input->GetChildren()[i]->GetName() == "ConfigurationParameter" )
		{
			TString ThisNum("Common"); ThisNum+=configParamNum;
			input->GetChildren()[i]->AppendPath( ThisNum.Data() );
			++configParamNum;
		}

		for( unsigned int j=0; j< input->GetChildren()[i]->GetChildren().size(); ++j )
		{
			AppendCommonNames( input->GetChildren()[i]->GetChildren()[j], subParamNum, appendParamNum, configParamNum );
		}
	}
}

XMLTag* XMLConfigReader::FindCommonPDFConfiguratorXML()
{
	for( unsigned int childIndex=0; childIndex< children.size(); ++childIndex )
	{
		if( children[childIndex]->GetName() == "CommonPDFConfigurator" )
		{
			return children[childIndex]->GetChildren()[0];
		}
	}
	return NULL;
}

//Make a precalculator for a data set
PrecalculatorConfig* XMLConfigReader::GetPrecalculatorConfig( )
{
	PrecalculatorConfig* preconfig = NULL;

	for( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		//Check the correct tag has been provided
		if ( children[childIndex]->GetName() == "Precalculator" )
		{
			preconfig = new PrecalculatorConfig();
			//Parse the tag components
			string precalculatorName = "Unspecified";
			string weightName = "Unspecified";
			string filename = "WeightedFile";
			bool useAlpha=false;
			int config=1;
			vector< XMLTag* > newchildren = children[childIndex]->GetChildren();
			for ( unsigned int newchildIndex = 0; newchildIndex < newchildren.size(); ++newchildIndex )
			{
				string name = newchildren[newchildIndex]->GetName();
				if ( name == "Name" )
				{
					precalculatorName = XMLTag::GetStringValue( newchildren[newchildIndex] );
				}
				else if ( name == "WeightName" )
				{
					weightName = XMLTag::GetStringValue( newchildren[newchildIndex] );
				}
				else if ( name == "Config" )
				{
					config = XMLTag::GetIntegerValue( newchildren[newchildIndex] );
				}
				else if ( name == "OutputFile" )
				{
					filename = XMLTag::GetStringValue( newchildren[newchildIndex] );
				}
				else if ( name == "UseAlpha" )
				{
					useAlpha = XMLTag::GetBooleanValue( newchildren[newchildIndex] );
				}
				else
				{
					cerr << "Unrecognised Precalculator component: " << name << endl;
					exit(1);
				}
			}
			preconfig->SetCalculatorName( precalculatorName );
			preconfig->SetWeightName( weightName );
			preconfig->SetConfig( (unsigned)config );
			preconfig->SetFileName( filename );
			preconfig->SetAlpha( useAlpha );
			return preconfig;
		}
	}

	return NULL;
}


//	Return the Integer Seed used in RapidFit XML tag <Seed>SomeInt</Seed>
unsigned int XMLConfigReader::GetSeed()
{
	unsigned int this_seed=0;
	if( seed < 0 )
	{
		//Find the NumberRepeats tag
		for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
		{
			if ( children[childIndex]->GetName() == "Seed" )
			{
				seed = XMLTag::GetIntegerValue( children[childIndex] );
				cout << "Using seed: " << seed << " from input file." << endl;
				return unsigned(seed);
			}
		}
		seed = 0 ;
		//If no such tag is found, report
		cout << "Seed tag not found in config file, defaulting to TRandom3(0)." << endl;
	} else  this_seed = unsigned(seed);

	if( _original_seed < 0 ) _original_seed = this_seed;
	return this_seed;
}

//	Set a new TRandom seed that is returned by the XMLFile
void XMLConfigReader::SetSeed( unsigned int new_seed )
{
	_original_seed = (int)new_seed;
	seed = (int)new_seed;
}

int XMLConfigReader::GetTotalDataSetSize()
{
	int DataSetSize=0;
	//Find the ParameterSets
	vector< XMLTag* > toFits;
	//  Children Defined globally on initialization!
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "ToFit" )
		{
			toFits.push_back( children[childIndex] );
		}
	}
	//  Loop over all of the fits
	for ( unsigned short int fitnum=0; fitnum < toFits.size(); ++fitnum )
	{
		vector<XMLTag*> allChildren = toFits[fitnum]->GetChildren();
		for( unsigned short int Childnum=0; Childnum < allChildren.size(); ++Childnum )
		{
			if( allChildren[Childnum]->GetName() == "DataSet" )
			{
				vector< XMLTag* > dataComponents = allChildren[Childnum]->GetChildren();
				for ( unsigned int dataIndex = 0; dataIndex < dataComponents.size(); ++dataIndex )
				{
					if ( dataComponents[dataIndex]->GetName() == "NumberEvents" )
					{
						DataSetSize += XMLTag::GetIntegerValue( dataComponents[dataIndex] );
					}
				}
			}
		}
	}

	return DataSetSize;
}

vector<int> XMLConfigReader::GetAllDataSetSizes()
{
	vector<int> DataSetSizes;
	//Find the ParameterSets
	vector< XMLTag* > toFits;
	//  Children Defined globally on initialization!
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "ToFit" )
		{
			toFits.push_back( children[childIndex] );
		}
	}
	//  Loop over all of the fits
	for ( unsigned short int fitnum=0; fitnum < toFits.size(); ++fitnum )
	{
		vector<XMLTag*> allChildren = toFits[fitnum]->GetChildren();
		for( unsigned short int Childnum=0; Childnum < allChildren.size(); ++Childnum )
		{
			if( allChildren[Childnum]->GetName() == "DataSet" )
			{
				vector< XMLTag* > dataComponents = allChildren[Childnum]->GetChildren();
				for ( unsigned int dataIndex = 0; dataIndex < dataComponents.size(); ++dataIndex )
				{
					if ( dataComponents[dataIndex]->GetName() == "NumberEvents" )
					{
						DataSetSizes.push_back( XMLTag::GetIntegerValue( dataComponents[dataIndex] ) );
					}
				}
			}
		}
	}
	return DataSetSizes;
}

vector<int> XMLConfigReader::GetAllStartEntries()
{
	vector<int> StartEntries;
	//Find the ParameterSets
	vector< XMLTag* > toFits;
	//  Children Defined globally on initialization!
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "ToFit" )
		{
			toFits.push_back( children[childIndex] );
		}
	}
	//  Loop over all of the fits
	for ( unsigned short int fitnum=0; fitnum < toFits.size(); ++fitnum )
	{
		vector<XMLTag*> allChildren = toFits[fitnum]->GetChildren();
		for( unsigned short int Childnum=0; Childnum < allChildren.size(); ++Childnum )
		{
			if( allChildren[Childnum]->GetName() == "DataSet" )
			{
				vector< XMLTag* > dataComponents = allChildren[Childnum]->GetChildren();
				bool found_one = false;
				for ( unsigned int dataIndex = 0; dataIndex < dataComponents.size(); ++dataIndex )
				{
					if ( dataComponents[dataIndex]->GetName() == "NumberEvents" )
					{
						found_one = true;
						StartEntries.push_back( XMLTag::GetIntegerValue( dataComponents[dataIndex] ) );
					}
				}
				if( !found_one )
				{
					StartEntries.push_back( 0 );
				}
			}
		}
	}
	return StartEntries;
}

unsigned int XMLConfigReader::GetOriginalSeed() const
{
	return (unsigned)_original_seed;
}

