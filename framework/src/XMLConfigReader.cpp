/**
  @class XMLConfigReader

  Opens an xml config file and uses it to create RapidFit data objects

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

#include <stdlib.h>
#include "XMLConfigReader.h"
#include "ClassLookUp.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "SumPDF.h"
#include "NormalisedSumPDF.h"
#include "ProdPDF.h"
#include "StringProcessing.h"
#include "AcceptReject.h"
#include "InvalidObject.h"
#include "DataFileLoader.h"
#include "ObservableContinuousConstraint.h"
#include "ObservableDiscreteConstraint.h"

//Default constructor
XMLConfigReader::XMLConfigReader() : isLoaded(false)
{
}

//Constructor with file name argument
XMLConfigReader::XMLConfigReader( string FileName ) : isLoaded(false)
{
	//Open the config file
	ifstream configFile( FileName.c_str() );
	if ( !configFile.is_open() )
	{
		cerr << "Failed to open config file \"" << FileName << "\"" << endl;
		return;
	}

	//Read the whole file into a vector
	vector<string> wholeFile;
	while ( configFile.good() )
	{
		string newLine;
		getline( configFile, newLine );
		wholeFile.push_back( newLine );
	}
	StringProcessing::RemoveWhiteSpace(wholeFile);

	//Check the root tag is correct
	if ( wholeFile[0] == "<RapidFit>" && wholeFile[ wholeFile.size() - 1 ] == "</RapidFit>" )
	{
		//Parse the first level tags
		vector<string> middleOfFile;
		for ( int lineIndex = 1; lineIndex < wholeFile.size() - 1; lineIndex++ )
		{
			middleOfFile.push_back( wholeFile[lineIndex] );
		}
		vector<string> value;
		children = XMLTag::FindTagsInContent( middleOfFile, value );

		//Anything in "value" is considered debug data
		for ( int valueIndex = 0; valueIndex < value.size(); valueIndex++ )
		{
			cout << value[valueIndex] << endl;
		}
	}
	else
	{
		cerr << "Config file in wrong format: root tag is \"" << wholeFile[0] << "\" not \"<RapidFit>\"" << endl;
	}

	isLoaded = true;
}

//Destructor
XMLConfigReader::~XMLConfigReader()
{
}

//Return whether file is loaded
bool XMLConfigReader::IsLoaded()
{
	return isLoaded;
}

//Return the parameter set
ParameterSet * XMLConfigReader::GetFitParameters()
{
	//Find the ParameterSet tag
	for ( int childIndex = 0; childIndex < children.size(); childIndex++ )
	{
		if ( children[childIndex]->GetName() == "ParameterSet" )
		{
			return GetParameterSet( children[childIndex] );
		}
	}

	//If no such tag is found, error
	cerr << "Config file does not specify parameters" << endl;
	return new ParameterSet();
}

//Return the minimiser for the fit
string XMLConfigReader::GetMinimiserName()
{
	//Get the name of the minimiser
	string minimiserName = "Uninitialised";
	for ( int childIndex = 0; childIndex < children.size(); childIndex++ )
	{
		if ( children[childIndex]->GetName() == "Minimiser" )
		{
			minimiserName = children[childIndex]->GetValue()[0];
		}
	}

	return minimiserName;
}

//Return the number of repeats for the fit
int XMLConfigReader::GetNumberRepeats()
{
	//Get the name of the minimiser
	string repeatString = "0";
	for ( int childIndex = 0; childIndex < children.size(); childIndex++ )
	{
		if ( children[childIndex]->GetName() == "NumberRepeats" )
		{
			repeatString = children[childIndex]->GetValue()[0];
		}
	}

	return atoi( repeatString.c_str() );
}

//Return the function to minimise
FitFunction * XMLConfigReader::GetFitFunction()
{
	//Find the FitFunction tag
	for ( int childIndex = 0; childIndex < children.size(); childIndex++ )
	{
		if ( children[childIndex]->GetName() == "FitFunction" )
		{
			return MakeFitFunction( children[childIndex] );
		}
	}

	//Tag not found, fail
	cerr << "FitFunction tag not found" << endl;
	exit(1);
}

//Make the FitFunction
FitFunction * XMLConfigReader::MakeFitFunction( XMLTag * FunctionTag )
{
	if ( FunctionTag->GetName() == "FitFunction" )
	{
		string functionName = "Uninitialised";
		string weightName = "Uninitialised";
		bool hasWeight = false;
		vector< XMLTag* > functionInfo = FunctionTag->GetChildren();
		if ( functionInfo.size() == 0 )
		{
			//Old style - just specifies the function name
			functionName = FunctionTag->GetValue()[0];
		}
		else
		{
			//New style - can have weights
			for ( int childIndex = 0; childIndex < functionInfo.size(); childIndex++ )
			{
				if ( functionInfo[childIndex]->GetName() == "FunctionName" )
				{
					functionName = functionInfo[childIndex]->GetValue()[0];
				}
				else if ( functionInfo[childIndex]->GetName() == "WeightName" )
				{
					hasWeight = true;
					weightName = functionInfo[childIndex]->GetValue()[0];
				}
				else
				{
					cerr << "Unrecognised FitFunction component: \"" << functionInfo[childIndex]->GetName() << endl;
					exit(1);
				}
			}
		}

		//Make the function
		if (hasWeight)
		{
			return ClassLookUp::LookUpFitFunctionName( functionName, weightName );
		}
		else
		{
			return ClassLookUp::LookUpFitFunctionName( functionName, "" );
		}
	}
	else
	{
		cerr << "Tag name incorrect: \"" << FunctionTag->GetName() << " not \"FitFunction\"" << endl;
		exit(1);
	}
}

//Organise all the PDFs and DataSets
vector< PDFWithData* > XMLConfigReader::GetPDFsAndData()
{
	//Collect all ToFit elements
	vector< XMLTag* > toFits;
	vector< PDFWithData* > pdfsAndData;
	for ( int childIndex = 0; childIndex < children.size(); childIndex++ )
	{
		if ( children[childIndex]->GetName() == "ToFit" )
		{
			toFits.push_back( children[childIndex] );
		}
	}

	if ( toFits.size() == 0 )
	{
		cerr << "No fits configured" << endl;
	}
	else
	{
		//Loop over all ToFits
		for ( int fitIndex = 0; fitIndex < toFits.size(); fitIndex++ )
		{
			XMLTag * pdfTag;
			XMLTag * dataTag;
			bool foundPDF = false;
			bool foundData = false;

			//Find the PDF and data set configuration
			vector< XMLTag* > fitComponents = toFits[fitIndex]->GetChildren();
			for ( int componentIndex = 0; componentIndex < fitComponents.size(); componentIndex++ )
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
				else
				{
					cerr << "Unrecognised fit component: " << name << endl;
				}
			}

			//Make the data set, and populate the physics bottle
			if ( foundData && foundPDF )
			{
				pdfsAndData.push_back( GetPDFWithData( dataTag, pdfTag ) );
			}
			else
			{
				cerr << "A ToFit xml tag is incomplete:";
				if ( !foundData )
				{
					cerr << " Data set configuration missing.";
				}
				if ( !foundPDF )
				{
					cerr << " PDF configuration missing.";
				}
				cerr << endl;
			}
		}
	}

	return pdfsAndData;
}

//Create a ParameterSet from the appropriate xml tag
ParameterSet * XMLConfigReader::GetParameterSet( XMLTag * InputTag )
{
	//Load data from the xml file
	if ( InputTag->GetName() == "ParameterSet" )
	{
		vector< PhysicsParameter* > physicsParameters;
		vector<string> names;
		string name = "";

		//Create each physics parameter
		vector< XMLTag* > parameters = InputTag->GetChildren();
		for ( int parameterIndex = 0; parameterIndex < parameters.size(); parameterIndex++ )
		{
			PhysicsParameter * newParameter = GetPhysicsParameter( parameters[parameterIndex], name );
			if ( newParameter->GetType() != "Invalid" )
			{
				physicsParameters.push_back( newParameter );
				names.push_back(name);
			}
		}

		//Create the parameter set
		ParameterSet * newParameters = new ParameterSet(names);
		for ( int nameIndex = 0; nameIndex < names.size(); nameIndex++ )
		{
			newParameters->SetPhysicsParameter( names[nameIndex], physicsParameters[nameIndex] );

		}
		return newParameters;
	}
	else
	{
		cout << "Incorrect xml tag provided: \"" << InputTag->GetName() << "\" not \"ParameterSet\"" << endl;
		return new ParameterSet();
	}
}

//Make a physics parameter from an appropriate XML tag
PhysicsParameter * XMLConfigReader::GetPhysicsParameter( XMLTag * InputTag, string & ParameterName )
{
	//Check the tag is actually a physics parameter
	if ( InputTag->GetName() == "PhysicsParameter" )
	{
		//Create some default values;
		ParameterName = "Uninitialised";
		string type = "Uninitialised";
		string unit = "Uninitialised";
		double value = 0.0;
		double minimum = 0.0;
		double maximum = 0.0;

		//Loop over the tag children, which correspond to the parameter elements
		vector< XMLTag* > elements = InputTag->GetChildren();
		for ( int elementIndex = 0; elementIndex < elements.size(); elementIndex++ )
		{
			string name = elements[elementIndex]->GetName();
			if ( name == "Name" )
			{
				ParameterName = elements[elementIndex]->GetValue()[0];
			}
			else if ( name == "Value" )
			{
				value = strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL );
			}
			else if ( name == "Minimum" )
			{
				minimum = strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL );
			}
			else if ( name == "Maximum" )
			{
				maximum = strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL );
			}
			else if ( name == "Type" )
			{
				type = elements[elementIndex]->GetValue()[0];
			}
			else if ( name == "Unit" )
			{
				unit = elements[elementIndex]->GetValue()[0];
			}
			else
			{
				cerr << "Unrecognised physics parameter configuration: " << name << endl;
			}
		}

		//Now construct the physics parameter
		if ( ( maximum == 0.0 && minimum == 0.0 ) || type == "Unbounded" )
		{
			//Unbounded parameter
			return new PhysicsParameter( ParameterName, value, type, unit );
		}
		else
		{
			//Bounded parameter
			return new PhysicsParameter( ParameterName, value, minimum, maximum, type, unit );
		}
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << InputTag->GetName() << "\" not \"PhysicsParameter\"" << endl;
		return new PhysicsParameter( "Invalid", 0.0, 0.0, 0.0, "Invalid", "Invalid" );
	}
}

//Collect the information needed to make a data set
PDFWithData * XMLConfigReader::GetPDFWithData( XMLTag * DataTag, XMLTag * FitPDFTag )
{
	//Check the tag actually is a data set
	if ( DataTag->GetName() == "DataSet" )
	{
		string dataSource = "Uninitialised";
		long numberEvents = 0;
		PhaseSpaceBoundary * dataBoundary;
		vector<string> dataArguments;
		bool boundaryFound = false;
		XMLTag * generatePDFTag;
		bool generatePDFFlag = false;

		//Retrieve the data set config
		vector< XMLTag* > dataComponents = DataTag->GetChildren();
		for ( int dataIndex = 0; dataIndex < dataComponents.size(); dataIndex++ )
		{
			string name = dataComponents[dataIndex]->GetName();
			if ( name == "Source" )
			{
				dataSource = dataComponents[dataIndex]->GetValue()[0];
			}
			else if ( name == "FileName" )
			{
				dataArguments.push_back( dataComponents[dataIndex]->GetValue()[0] );
			}
			else if ( name == "NumberEvents" )
			{
				numberEvents = strtol( dataComponents[dataIndex]->GetValue()[0].c_str(), NULL, 10 );
			}
			else if ( name == "PhaseSpaceBoundary" )
			{
				boundaryFound = true;
				dataBoundary = GetPhaseSpaceBoundary( dataComponents[dataIndex] );
			}
			else if ( name == "PDF" || name == "SumPDF" || name == "NormalisedSumPDF" || name == "ProdPDF" )
			{
				generatePDFTag = dataComponents[dataIndex];
				generatePDFFlag = true;
			}
			else
			{
				cerr << "Unrecognised data set component: " << name << endl;
			}
		}

		//Return the collection of configuration information - data generation will happen later
		if (!boundaryFound)
		{
			cerr << "DataSet defined without PhaseSpaceBoundary" << endl;
			return new PDFWithData();
		}

		//Make the fit PDF
		IPDF * fitPDF = GetPDF( FitPDFTag, dataBoundary );

		if (generatePDFFlag)
		{
			//Make a separate data generation PDF
			IPDF * generatePDF = GetPDF( generatePDFTag, dataBoundary );

			return new PDFWithData( generatePDF, fitPDF, dataSource, numberEvents, dataArguments, dataBoundary );
		}
		else
		{
			return new PDFWithData( fitPDF, dataSource, numberEvents, dataArguments, dataBoundary );
		}
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << DataTag->GetName() << "\" not \"DataSet\"" << endl;
		return new PDFWithData();
	}
}

//Make a PhaseSpaceBoundary from the appropriate xml tag
PhaseSpaceBoundary * XMLConfigReader::GetPhaseSpaceBoundary( XMLTag * InputTag )
{
	//Check the tag is actually a PhaseSpaceBoundary
	if ( InputTag->GetName() == "PhaseSpaceBoundary" )
	{
		vector< IConstraint* > constraints;
		vector<string> names;
		string name;

		//Create each single bound
		vector< XMLTag* > constraintTags = InputTag->GetChildren();
		for ( int boundIndex = 0; boundIndex < constraintTags.size(); boundIndex++ )
		{
			IConstraint * newConstraint = GetConstraint( constraintTags[boundIndex], name );
			if ( newConstraint->GetUnit() != "Invalid" )
			{
				constraints.push_back(newConstraint);
				names.push_back(name);
			}
		}

		//Create the parameter set
		PhaseSpaceBoundary * newBoundary = new PhaseSpaceBoundary(names);
		for ( int nameIndex = 0; nameIndex < names.size(); nameIndex++ )
		{
			newBoundary->SetConstraint( names[nameIndex], constraints[nameIndex] );

		}
		return newBoundary;
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << InputTag->GetName() << "\" not \"PhaseSpaceBoundary\"" << endl;
		return new PhaseSpaceBoundary();
	}
}

//Make an IConstraint from the appropriate xml tag
IConstraint * XMLConfigReader::GetConstraint( XMLTag * InputTag, string & Name )
{
	//Check the tag is actually a single bound
	if ( InputTag->GetName() == "Observable" )
	{
		//Create some default values;
		Name = "Uninitialised";
		string unit = "Uninitialised";
		double minimum = 0.0;
		double maximum = 0.0;
		vector<double> allValues;

		//Loop over the tag children, which correspond to the parameter elements
		vector< XMLTag* > elements = InputTag->GetChildren();
		for ( int elementIndex = 0; elementIndex < elements.size(); elementIndex++ )
		{
			if ( elements[elementIndex]->GetName() == "Name" )
			{
				Name = elements[elementIndex]->GetValue()[0];
			}
			else if ( elements[elementIndex]->GetName() == "Minimum" )
			{
				minimum = strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL );
			}
			else if ( elements[elementIndex]->GetName() == "Maximum" )
			{
				maximum = strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL );
			}
			else if ( elements[elementIndex]->GetName() == "Value" )
			{
				allValues.push_back( strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL ) );
			}
			else if ( elements[elementIndex]->GetName() == "Unit" )
			{
				unit = elements[elementIndex]->GetValue()[0];
			}
			else
			{
				cerr << "Unrecognised constraint configuration: " <<  elements[elementIndex]->GetName() << endl;
			}
		}

		//If there are discrete values, make a discrete constraint
		if ( allValues.size() > 0 )
		{
			return new ObservableDiscreteConstraint( Name, allValues, unit );
		}
		else
		{
			return new ObservableContinuousConstraint( Name, minimum, maximum, unit );
		}
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << InputTag->GetName() << "\" not \"Observable\"" << endl;
		return new ObservableContinuousConstraint( "Invalid", 0.0, 0.0, "Invalid" );
	}
}

//Create a PDF from an appropriate xml tag
IPDF * XMLConfigReader::GetNamedPDF( XMLTag * InputTag )
{
	//Check the tag actually is a PDF
	if ( InputTag->GetName() == "PDF" )
	{
		vector< XMLTag* > pdfConfig = InputTag->GetChildren();
		string name;
		vector<string> observableNames, parameterNames;

		//Load the PDF configuration
		for ( int configIndex = 0; configIndex < pdfConfig.size(); configIndex++ )
		{
			if ( pdfConfig[configIndex]->GetName() == "Name" )
			{
				name = pdfConfig[configIndex]->GetValue()[0];
			}
			else if ( pdfConfig[configIndex]->GetName() == "ObservableName" )
			{
				observableNames.push_back( pdfConfig[configIndex]->GetValue()[0] );
			}
			else if ( pdfConfig[configIndex]->GetName() == "ParameterName" )
			{
				parameterNames.push_back( pdfConfig[configIndex]->GetValue()[0] );
			}
			else
			{
				cerr << "Unrecognised PDF configuration: " << pdfConfig[configIndex]->GetName() << endl;
			}
		}

		//Check if the name is recognised as a PDF
		return ClassLookUp::LookUpPDFName( name, observableNames, parameterNames );
	}
	else
	{
		cerr << "Incorrect tag provided: \"" << InputTag->GetName() << "\" not \"PDF\"" << endl;
		return new InvalidObject( "Incorrect tag provided: \"" + InputTag->GetName() + "\" not \"PDF\"" );
	}
}

//Create a SumPDF from an appropriate xml tag
IPDF * XMLConfigReader::GetSumPDF( XMLTag * InputTag, PhaseSpaceBoundary * InputBoundary )
{
	//Check the tag actually is a PDF
	if ( InputTag->GetName() == "SumPDF" )
	{
		vector< XMLTag* > pdfConfig = InputTag->GetChildren();
		string fractionName = "unspecified";
		vector< IPDF* > componentPDFs;

		//Load the PDF configuration
		for ( int configIndex = 0; configIndex < pdfConfig.size(); configIndex++ )
		{
			if ( pdfConfig[configIndex]->GetName() == "FractionName" )
			{
				fractionName = pdfConfig[configIndex]->GetValue()[0];
			}
			else
			{
				componentPDFs.push_back( GetPDF( pdfConfig[configIndex], InputBoundary ) );
			}
		}

		//Check there are two component PDFs to sum
		if ( componentPDFs.size() == 2 )
		{
			if ( fractionName == "unspecified" )
			{
				return new SumPDF( componentPDFs[0], componentPDFs[1], InputBoundary );
			}
			else
			{
				return new SumPDF( componentPDFs[0], componentPDFs[1], InputBoundary, fractionName );
			}
		}
		else
		{
			cerr << "Incorrect number of PDFs to sum: " << componentPDFs.size() << " not 2" << endl;
			return new InvalidObject( "SumPDF can only sum over 2 PDFs" );
		}
	}
	else
	{
		cerr << "Incorrect tag provided: \"" << InputTag->GetName() << "\" not \"SumPDF\"" << endl;
		return new InvalidObject( "Incorrect tag provided: \"" + InputTag->GetName() + "\" not \"SumPDF\"" );
	}
}

//Create a NormalisedSumPDF from an appropriate xml tag
IPDF * XMLConfigReader::GetNormalisedSumPDF( XMLTag * InputTag, PhaseSpaceBoundary * InputBoundary )
{
	//Check the tag actually is a PDF
	if ( InputTag->GetName() == "NormalisedSumPDF" )
	{
		vector< XMLTag* > pdfConfig = InputTag->GetChildren();
		string fractionName = "unspecified";
		vector< IPDF* > componentPDFs;

		//Load the PDF configuration
		for ( int configIndex = 0; configIndex < pdfConfig.size(); configIndex++ )
		{
			if ( pdfConfig[configIndex]->GetName() == "FractionName" )
			{
				fractionName = pdfConfig[configIndex]->GetValue()[0];
			}
			else
			{
				componentPDFs.push_back( GetPDF( pdfConfig[configIndex], InputBoundary ) );
			}
		}

		//Check there are two component PDFs to sum
		if ( componentPDFs.size() == 2 )
		{
			if ( fractionName == "unspecified" )
			{
				return new NormalisedSumPDF( componentPDFs[0], componentPDFs[1], InputBoundary );
			}
			else
			{
				return new NormalisedSumPDF( componentPDFs[0], componentPDFs[1], InputBoundary, fractionName );
			}
		}
		else
		{
			cerr << "Incorrect number of PDFs to sum: " << componentPDFs.size() << " not 2" << endl;
			return new InvalidObject( "NormalisedSumPDF can only sum over 2 PDFs" );
		}
	}
	else
	{
		cerr << "Incorrect tag provided: \"" << InputTag->GetName() << "\" not \"NormalisedSumPDF\"" << endl;
		return new InvalidObject( "Incorrect tag provided: \"" + InputTag->GetName() + "\" not \"NormalisedSumPDF\"" );
	}
}

//Create a ProdPDF from an appropriate xml tag
IPDF * XMLConfigReader::GetProdPDF( XMLTag * InputTag, PhaseSpaceBoundary * InputBoundary )
{
	//Check the tag actually is a PDF
	if ( InputTag->GetName() == "ProdPDF" )
	{
		vector< XMLTag* > pdfConfig = InputTag->GetChildren();
		vector< IPDF* > componentPDFs;

		//Load the PDF configuration
		for ( int configIndex = 0; configIndex < pdfConfig.size(); configIndex++ )
		{
			componentPDFs.push_back( GetPDF( pdfConfig[configIndex], InputBoundary ) );
		}

		//Check there are two component PDFs to sum
		if ( componentPDFs.size() == 2 )
		{
			return new ProdPDF( componentPDFs[0], componentPDFs[1] );
		}
		else
		{
			cerr << "Incorrect number of PDFs to multiply: " << componentPDFs.size() << " not 2" << endl;
			return new InvalidObject( "ProdPDF can only multiply 2 PDFs" );
		}
	}
	else
	{
		cerr << "Incorrect tag provided: \"" << InputTag->GetName() << "\" not \"ProdPDF\"" << endl;
		return new InvalidObject( "Incorrect tag provided: \"" + InputTag->GetName() + "\" not \"ProdPDF\"" );
	}
}

//Choose one of the PDF instantiation methods
IPDF * XMLConfigReader::GetPDF( XMLTag * InputTag, PhaseSpaceBoundary * InputBoundary )
{
	if ( InputTag->GetName() == "PDF" )
	{
		return GetNamedPDF(InputTag);
	}
	else if ( InputTag->GetName() == "SumPDF" )
	{
		return GetSumPDF( InputTag, InputBoundary );
	}
	else if ( InputTag->GetName() == "NormalisedSumPDF" )
	{
		return GetNormalisedSumPDF( InputTag, InputBoundary );        
	}
	else if ( InputTag->GetName() == "ProdPDF" )
	{
		return GetProdPDF( InputTag, InputBoundary );
	}
	else
	{
		cerr << "Unrecognised PDF configuration: " << InputTag->GetName() << endl;
	}
}
