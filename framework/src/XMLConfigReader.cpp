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
		exit(1);
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
		exit(1);
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

	//If no such tag is found, fail
	cerr << "ParameterSet tag not found in config file" << endl;
	exit(1);
}

//Return the minimiser for the fit
MinimiserConfiguration * XMLConfigReader::GetMinimiserConfiguration()
{
	//Find the Minimiser tag
	for ( int childIndex = 0; childIndex < children.size(); childIndex++ )
	{
		if ( children[childIndex]->GetName() == "Minimiser" )
		{
			return MakeMinimiser( children[childIndex] );
		}
	}

	//If no such tag is found, fail
	cerr << "Minimiser tag not found in config file" << endl;
	exit(1);
}

//Make a minimiser configuration object
MinimiserConfiguration * XMLConfigReader::MakeMinimiser( XMLTag * MinimiserTag )
{
	if ( MinimiserTag->GetName() == "Minimiser" )
	{
		//Examine all minimiser components
		string minimiserName = "Uninitialised";
		vector< XMLTag* > minimiserComponents = MinimiserTag->GetChildren();
		vector<string> valueLines = MinimiserTag->GetValue();
		if ( valueLines.size() == 1 )
		{
			minimiserName = MinimiserTag->GetValue()[0];
		}
		else
		{
			cerr << "Minimiser tag contains " << valueLines.size() << " lines, not 1" << endl;
			exit(1);
		}
		return new MinimiserConfiguration( minimiserName, GetOutputConfiguration() );
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << MinimiserTag->GetName() << "\" not \"Minimiser\"" << endl;
		exit(1);
	}
}

//Return the output configuration for the fit
OutputConfiguration * XMLConfigReader::GetOutputConfiguration()
{
	//Find the Output tag
	for ( int childIndex = 0; childIndex < children.size(); childIndex++ )
	{
		if ( children[childIndex]->GetName() == "Output" )
		{
			return MakeOutputConfiguration( children[childIndex] );
		}
	}

	//If no such tag is found, make default
	cout << "Output tag not found in config file - using default" << endl;
	return new OutputConfiguration();
}

//Make an output configuration object
OutputConfiguration * XMLConfigReader::MakeOutputConfiguration( XMLTag * OutputTag )
{
	if ( OutputTag->GetName() == "Output" )
	{
		vector< pair< string, string > > contourPlots;
		vector<string> projections;
		string pullType = "None";
		vector< XMLTag* > outputComponents = OutputTag->GetChildren();
		for ( int childIndex = 0; childIndex < outputComponents.size(); childIndex++ )
		{
			if ( outputComponents[childIndex]->GetName() == "ContourPlot" )
			{
				contourPlots.push_back( MakeContourPlot( outputComponents[childIndex] ) );
			}
			else if ( outputComponents[childIndex]->GetName() == "Projection" )
			{
				projections.push_back( outputComponents[childIndex]->GetValue()[0] );
			}
			else if ( outputComponents[childIndex]->GetName() == "DoPullPlots" )
			{
				pullType = outputComponents[childIndex]->GetValue()[0];
			}
			else
			{
				cerr << "Unrecognised output component: " << outputComponents[childIndex]->GetName() << endl;
				exit(1);
			}
		}

		return new OutputConfiguration( contourPlots, projections, pullType );
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << OutputTag->GetName() << "\" not \"Output\"" << endl;
		exit(1);
	}
}

//Return the pair of observables to plot the function contours for
pair< string, string > XMLConfigReader::MakeContourPlot( XMLTag * PlotTag )
{
	if ( PlotTag->GetName() == "ContourPlot" )
	{
		vector< XMLTag* > plotComponents = PlotTag->GetChildren();
		if ( plotComponents.size() == 2 )
		{
			//Retrieve the names of the parameters to plot
			string xName, yName;
			bool hasX = false;
			bool hasY = false;
			for ( int childIndex = 0; childIndex < plotComponents.size(); childIndex++ )
			{
				if ( plotComponents[childIndex]->GetName() == "XParameter" )
				{
					xName = plotComponents[childIndex]->GetValue()[0];
					hasX = true;
				}
				else if ( plotComponents[childIndex]->GetName() == "YParameter" )
				{
					yName = plotComponents[childIndex]->GetValue()[0];
					hasY = true;
				}
				else
				{
					cerr << "Unrecognised ContourPlot component: " << plotComponents[childIndex]->GetName() << endl;
					exit(1);
				}
			}

			//Check both parameters are specified
			if ( hasX && hasY )
			{
				return make_pair( xName, yName );
			}
			else
			{
				cerr << "ContourPlot tag is missing parameter name for ";
				if ( !hasX )
				{
					cerr << "x";
				}
				else
				{
					cerr << "y";
				}
				cerr << " axis" << endl;
				exit(1);
			}
		}
		else
		{
			cerr << "ContourPlot tag should only contain the two parameters to plot contours for" << endl;
			exit(1);
		}
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << PlotTag << "\" not \"ContourPlot\"" << endl;
		exit(1);
	}
}

//Return the number of repeats for the fit
int XMLConfigReader::GetNumberRepeats()
{
	//Find the NumberRepeats tag
	for ( int childIndex = 0; childIndex < children.size(); childIndex++ )
	{
		if ( children[childIndex]->GetName() == "NumberRepeats" )
		{
			return atoi( children[childIndex]->GetValue()[0].c_str() );
		}
	}

	//If no such tag is found, fail
	cerr << "NumberRepeats tag not found in config file" << endl;
	exit(1);
}

//Return the function to minimise
FitFunctionConfiguration * XMLConfigReader::GetFitFunctionConfiguration()
{
	//Find the FitFunction tag
	for ( int childIndex = 0; childIndex < children.size(); childIndex++ )
	{
		if ( children[childIndex]->GetName() == "FitFunction" )
		{
			return MakeFitFunction( children[childIndex] );
		}
	}

	//If no such tag is found, fail
	cerr << "FitFunction tag not found in config file" << endl;
	exit(1);
}

//Make FitFunction configuration object
FitFunctionConfiguration * XMLConfigReader::MakeFitFunction( XMLTag * FunctionTag )
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
			return new FitFunctionConfiguration( functionName, weightName );
		}
		else
		{
			return new FitFunctionConfiguration( functionName );
		}
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << FunctionTag->GetName() << " not \"FitFunction\"" << endl;
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
		cerr << "No ToFit tags found in config file" << endl;
		exit(1);
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
					exit(1);
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
				exit(1);
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
		cerr << "Incorrect xml tag provided: \"" << InputTag->GetName() << "\" not \"ParameterSet\"" << endl;
		exit(1);
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
		double mean = 0.0;
		double sigma = 0.0;
		bool hasValue = false;
		bool hasMaximum = false;
		bool hasMinimum = false;
		bool hasMean = false;
		bool hasSigma = false;

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
				hasValue = true;
				value = strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL );
			}
			else if ( name == "Minimum" )
			{
				hasMinimum = true;
				minimum = strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL );
			}
			else if ( name == "Maximum" )
			{
				hasMaximum = true;
				maximum = strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL );
			}
			else if ( name == "Mean" )
			{
				hasMean = true;
				mean = strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL );
			}
			else if ( name == "Sigma" )
			{
				hasSigma = true;
				sigma = strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL );
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
				exit(1);
			}
		}

		//Now construct the physics parameter
		if (hasValue)
		{
			if ( hasMaximum && hasMinimum )
			{
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
				//Check for ambiguity
				if ( hasMaximum || hasMinimum )
				{
					cerr << "Ambiguous parameter definition: " << ParameterName << " has value and ";
					if ( hasMaximum )
					{
						cerr << "maximum, but not minimum";
					}
					if ( hasMinimum )
					{
						cerr << "minimum, but not maximum";
					}
					cerr << " defined" << endl;
					exit(1);
				}
				else
				{
					//Unbounded parameter
					return new PhysicsParameter( ParameterName, value, type, unit );
				}
			}
		}
		else if ( hasMean && hasSigma )
		{
			//Gaussian constrained parameter
			return new PhysicsParameter( ParameterName, mean, mean - sigma, mean + sigma, type, unit );
		}
		else
		{
			cerr << "Ambiguous definition for parameter " << ParameterName << endl;
			exit(1);
		}
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << InputTag->GetName() << "\" not \"PhysicsParameter\"" << endl;
		exit(1);
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
		vector<string> dataArguments, argumentNames;
		bool boundaryFound = false;
		vector< IPrecalculator* > dataProcessors;
		vector< DataSetConfiguration* > dataSetMakers;
		bool generatePDFFlag = false;
		IPDF * generatePDF;

		//Retrieve the data set config
		vector< XMLTag* > dataComponents = DataTag->GetChildren();
		for ( int dataIndex = 0; dataIndex < dataComponents.size(); dataIndex++ )
		{
			string name = dataComponents[dataIndex]->GetName();
			if ( name == "Source" )
			{
				dataSource = dataComponents[dataIndex]->GetValue()[0];
			}
			else if ( name == "Subset" )
			{
				dataSetMakers.push_back( MakeDataSetConfiguration( dataComponents[dataIndex], dataBoundary ) );
			}
			else if ( name == "FileName" || name == "NTuplePath" )
			{
				argumentNames.push_back(name);
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
				generatePDF = GetPDF( dataComponents[dataIndex], dataBoundary );
				generatePDFFlag = true;
			}
			else if ( name == "Precalculator" )
			{
				dataProcessors.push_back( MakePrecalculator( dataComponents[dataIndex], dataBoundary ) );
			}
			else
			{
				cerr << "Unrecognised data set component: " << name << endl;
				exit(1);
			}
		}

		//Return the collection of configuration information - data generation will happen later
		if (!boundaryFound)
		{
			cerr << "DataSet defined without PhaseSpaceBoundary" << endl;
			exit(1);
		}

		//If there are no separate data sources, go for backwards compatibility
		if ( dataSetMakers.size() < 1 )
		{
			DataSetConfiguration * oldStyleConfig;

			if (generatePDFFlag)
			{
				oldStyleConfig = new DataSetConfiguration( dataSource, numberEvents, dataArguments, argumentNames, generatePDF );
			}
			else
			{
				oldStyleConfig = new DataSetConfiguration( dataSource, numberEvents, dataArguments, argumentNames );
			}

			dataSetMakers.push_back(oldStyleConfig);
		}

		//Make the objects
		IPDF * fitPDF = GetPDF( FitPDFTag, dataBoundary );
		return new PDFWithData( fitPDF, dataBoundary, dataSetMakers, dataProcessors);
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << DataTag->GetName() << "\" not \"DataSet\"" << endl;
		exit(1);
	}
}

//Collect the information needed to make a data set
DataSetConfiguration * XMLConfigReader::MakeDataSetConfiguration( XMLTag * DataTag, PhaseSpaceBoundary * DataBoundary )
{
	//Check the tag actually is a data set
	if ( DataTag->GetName() == "Subset" )
	{
		string dataSource = "Uninitialised";
		long numberEvents = 0;
		vector<string> dataArguments, argumentNames;
		bool generatePDFFlag = false;
		IPDF * generatePDF;

		//Retrieve the data set config
		vector< XMLTag* > dataComponents = DataTag->GetChildren();
		for ( int dataIndex = 0; dataIndex < dataComponents.size(); dataIndex++ )
		{
			string name = dataComponents[dataIndex]->GetName();
			if ( name == "Source" )
			{
				dataSource = dataComponents[dataIndex]->GetValue()[0];
			}
			else if ( name == "FileName" || name == "NTuplePath" )
			{
				argumentNames.push_back(name);
				dataArguments.push_back( dataComponents[dataIndex]->GetValue()[0] );
			}
			else if ( name == "NumberEvents" )
			{
				numberEvents = strtol( dataComponents[dataIndex]->GetValue()[0].c_str(), NULL, 10 );
			}
			else if ( name == "PDF" || name == "SumPDF" || name == "NormalisedSumPDF" || name == "ProdPDF" )
			{
				generatePDF = GetPDF( dataComponents[dataIndex], DataBoundary );
				generatePDFFlag = true;
			}
			else
			{
				cerr << "Unrecognised data set component: " << name << endl;
				exit(1);
			}
		}

		if (generatePDFFlag)
		{
			return new DataSetConfiguration( dataSource, numberEvents, dataArguments, argumentNames, generatePDF );
		}
		else
		{
			return new DataSetConfiguration( dataSource, numberEvents, dataArguments, argumentNames );
		}
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << DataTag->GetName() << "\" not \"Subset\"" << endl;
		exit(1);
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
		exit(1);
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
				exit(1);
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
		exit(1);
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
				exit(1);
			}
		}

		//Check if the name is recognised as a PDF
		return ClassLookUp::LookUpPDFName( name, observableNames, parameterNames );
	}
	else
	{
		cerr << "Incorrect tag provided: \"" << InputTag->GetName() << "\" not \"PDF\"" << endl;
		exit(1);
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
			exit(1);
		}
	}
	else
	{
		cerr << "Incorrect tag provided: \"" << InputTag->GetName() << "\" not \"SumPDF\"" << endl;
		exit(1);
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
			exit(1);
		}
	}
	else
	{
		cerr << "Incorrect tag provided: \"" << InputTag->GetName() << "\" not \"NormalisedSumPDF\"" << endl;
		exit(1);
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
			exit(1);
		}
	}
	else
	{
		cerr << "Incorrect tag provided: \"" << InputTag->GetName() << "\" not \"ProdPDF\"" << endl;
		exit(1);
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
		exit(1);
	}
}

//Make a precalculator for a data set
IPrecalculator * XMLConfigReader::MakePrecalculator( XMLTag * InputTag, PhaseSpaceBoundary * InputBoundary )
{
	//Check the correct tag has been provided
	if ( InputTag->GetName() == "Precalculator" )
	{
		//Parse the tag components
		vector< IPDF* > componentPDFs;
		string precalculatorName = "Unspecified";
		string weightName = "Unspecified";
		vector< XMLTag* > children = InputTag->GetChildren();
		for ( int childIndex = 0; childIndex < children.size(); childIndex++ )
		{
			string name = children[childIndex]->GetName();
			if ( name == "Name" )
			{
				precalculatorName = children[childIndex]->GetValue()[0];
			}
			else if ( name == "WeightName" )
			{
				weightName = children[childIndex]->GetValue()[0];
			}
			else if ( name == "PDF" || name == "SumPDF" || name == "NormalisedSumPDF" || name == "ProdPDF" )
			{
				componentPDFs.push_back( GetPDF( children[childIndex], InputBoundary ) );
			}
			else
			{
				cerr << "Unrecognised Precalculator component: " << name << endl;
				exit(1);
			}
		}

		//Check a signal and background PDF have been provided
		if ( componentPDFs.size() == 2 )
		{
			return ClassLookUp::LookUpPrecalculator( precalculatorName, componentPDFs[0], componentPDFs[1], GetFitParameters(), weightName );
		}
		else
		{
			cerr << "Precalculator expecting 2 PDFs, not " << componentPDFs.size() << endl;
			exit(1);
		}
	}
	else
	{
		cerr << "Incorrect tag provided: \"" << InputTag->GetName() << "\" not \"Precalculator\"" << endl;
		exit(1);
	}
}
