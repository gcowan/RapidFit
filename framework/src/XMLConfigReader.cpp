/**
  @class XMLConfigReader

  Opens an xml config file and uses it to create RapidFit data objects

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

//	RapidFit Headers
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
//	System Headers
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <float.h>

using namespace::std;

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

//Constructor with file name argument
XMLConfigReader::XMLConfigReader( string FileName, vector<pair<string, string> >* OverrideXML ) : wholeFile(), All_XML_Tags(new XMLTag(OverrideXML)), children(), seed(-1), debug(new DebugClass(false) )
{
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

	vector<string> value;
	vector<XMLTag*> File_Tags = All_XML_Tags->FindTagsInContent( wholeFile, value );
	children = File_Tags[0]->GetChildren();

	if( ( File_Tags.size() != 1) || ( children.size() == 0 ) )
	{
		cerr << "Error processing XMLFile" << endl;
		exit( -234 );
	}
}

//Destructor
XMLConfigReader::~XMLConfigReader()
{
	delete All_XML_Tags;
	//cout << "Hello from XMLConfigReader destructor" << endl;
	if( debug != NULL ) delete debug;
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

			//		min			-		max			< 	1E-5

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
			All_Parameters.push_back( GetParameterSet( children[childIndex] ) );
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
		//vector<string> valueLines = MinimiserTag->GetValue();
		int MAXIMUM_MINIMISATION_STEPS = 1000000;
		double FINAL_GRADIENT_TOLERANCE = 0.001;
		MinimiserConfiguration* returnableConfig = NULL;
		vector<string> minimiserOptions;
		int Quality = 1;
		bool MultiMini=false;
		if ( minimiserComponents.size() == 0 )
		{
			minimiserName = MinimiserTag->GetValue()[0];
		}
		else
		{
			//New style - can have weights
			for ( unsigned int childIndex = 0; childIndex < minimiserComponents.size(); ++childIndex )
			{
				if ( minimiserComponents[childIndex]->GetName() == "MinimiserName" )
				{
					minimiserName = minimiserComponents[childIndex]->GetValue()[0];
				}
				else if( minimiserComponents[childIndex]->GetName() == "MaxSteps" )
				{
					MAXIMUM_MINIMISATION_STEPS = atoi( minimiserComponents[childIndex]->GetValue()[0].c_str() );
				}
				else if( minimiserComponents[childIndex]->GetName() == "GradTolerance" )
				{
					FINAL_GRADIENT_TOLERANCE = strtod( minimiserComponents[childIndex]->GetValue()[0].c_str(), NULL );
				}
				else if( minimiserComponents[childIndex]->GetName() == "ConfigureMinimiser" )
				{
					minimiserOptions.push_back( minimiserComponents[childIndex]->GetValue()[0] );
				}
				else if( minimiserComponents[childIndex]->GetName() == "Quality" )
				{
					Quality = atoi( minimiserComponents[childIndex]->GetValue()[0].c_str() );
				}
				else if( minimiserComponents[childIndex]->GetName() == "MultiMini" )
				{
					if( minimiserComponents[childIndex]->GetValue()[0] == "true" )
					{
						MultiMini = true;
					}
				}
				else
				{
					cerr << "Minimiser not properly configured" << endl;
					exit(9234);
				}
			}
		}
		if( debug->GetStatus() ) cout << "XMLConfigReader:\tMinimiser: " << minimiserName << " requested in XML, creating" << endl;
		returnableConfig = new MinimiserConfiguration( minimiserName, GetOutputConfiguration() );
		returnableConfig->SetSteps( MAXIMUM_MINIMISATION_STEPS );
		returnableConfig->SetTolerance( FINAL_GRADIENT_TOLERANCE );
		returnableConfig->SetOptions( minimiserOptions );
		returnableConfig->SetQuality( Quality );
		returnableConfig->SetMultiMini( MultiMini );
		if( debug->GetStatus() )
		{
			cout << "XMLConfigReader:\tMinimiser Configured With:\t" << endl;
			cout << "XMLConfigReader:\tSteps:\t" << MAXIMUM_MINIMISATION_STEPS << endl;
			cout << "XMLConfigReader:\tTolerance:\t" << FINAL_GRADIENT_TOLERANCE << endl;
			cout << "XMLConfigReader:\tOptions:" << endl;
			for( unsigned int i=0; i< minimiserOptions.size(); ++i )
			{
				cout << "\t\t\t" << minimiserOptions[i] << endl;
			}
			cout << "XMLConfigReader:\tQuality:\t" << Quality << endl;
			cout << "XMLConfigReader:\tMultiMini:\t" << string((MultiMini==true)?"True":"False") << endl;
		}
		return returnableConfig;
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
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "Output" )
		{
			return MakeOutputConfiguration( children[childIndex] );
		}
	}

	//If no such tag is found, make default
	//cout << "Output tag not found in config file - using default" << endl;
	return new OutputConfiguration( vector<pair<string,string> >(), string("None"), vector<ScanParam*>(), vector<pair<ScanParam*,ScanParam*> >(), vector<CompPlotter_config*>()  );
}

//Make an output configuration object
OutputConfiguration * XMLConfigReader::MakeOutputConfiguration( XMLTag * OutputTag )
{
	if ( OutputTag->GetName() == "Output" )
	{
		vector< pair< string, string > > contourPlots;
		vector<string> componentprojections;
		string pullType = "None";
		vector<ScanParam*> ScanParameters;
		vector<pair<ScanParam*, ScanParam*> > _2DScanParameters;

		vector<CompPlotter_config*> compPlotter_vec;

		vector< XMLTag* > outputComponents = OutputTag->GetChildren();
		for ( unsigned int childIndex = 0; childIndex < outputComponents.size(); ++childIndex )
		{
			if ( outputComponents[childIndex]->GetName() == "ContourPlot" )
			{
				contourPlots.push_back( MakeContourPlot( outputComponents[childIndex] ) );
			}
			else if ( outputComponents[childIndex]->GetName() == "Projection" )
			{
				compPlotter_vec.push_back( getCompPlotterConfigs( outputComponents[childIndex] ) );
			}
			else if ( outputComponents[childIndex]->GetName() == "ComponentProjection" )
			{
				compPlotter_vec.push_back( getCompPlotterConfigs( outputComponents[childIndex] ) );
			}
			else if ( outputComponents[childIndex]->GetName() == "DoPullPlots" )
			{
				pullType = outputComponents[childIndex]->GetValue()[0];
			}
			else if ( outputComponents[childIndex]->GetName() == "Scan" )
			{
				ScanParam* temp_SParam = GetScanParam( outputComponents[childIndex] );
				ScanParameters.push_back( temp_SParam );
			}
			else if ( outputComponents[childIndex]->GetName() == "TwoDScan" )
			{
				pair<ScanParam*, ScanParam*> temp_2DScan = Get2DScanParam( outputComponents[childIndex] );
				_2DScanParameters.push_back( temp_2DScan );
			}
			else if ( outputComponents[childIndex]->GetName() == "2DScan" )
			{
				cerr << "PLEASE MOVE YOUR XML TO THE NEW SYNTAX \'<TwoDScan>\' TO BE STANDARDS COMPLIANT!" << endl;
				pair<ScanParam*, ScanParam*> temp_2DScan = Get2DScanParam( outputComponents[childIndex] );
				_2DScanParameters.push_back( temp_2DScan );
			}
			else
			{
				cerr << "Unrecognised output component: " << outputComponents[childIndex]->GetName() << endl;
				exit(1);
			}
		}

		if( debug->GetStatus() )
		{
			cout << "XMLConfigReader:\tOutput Configuration Requested, Contructing." << endl;
		}

		return new OutputConfiguration( contourPlots, pullType, ScanParameters, _2DScanParameters, compPlotter_vec );
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << OutputTag->GetName() << "\" not \"Output\"" << endl;
		exit(1);
	}
}

CompPlotter_config* XMLConfigReader::getCompPlotterConfigs( XMLTag* CompTag )
{
	vector< XMLTag* > projComps = CompTag->GetChildren();

	CompPlotter_config* returnable_config = new CompPlotter_config();

	returnable_config->observableName = "Undefined";
	returnable_config->data_bins = 100;
	returnable_config->PDF_points = 128;
	returnable_config->logY = false;

	returnable_config->xmin = -999;
	returnable_config->xmax = -999;
	returnable_config->ymin = -999;
	returnable_config->ymax = -999;

	returnable_config->CalcChi2 = false;
	returnable_config->Chi2Value = -999;

	returnable_config->OnlyZero = false;

	//	If we have no extra config just take the name and use defaults
	if( projComps.size() == 0 )
	{
		returnable_config->observableName = CompTag->GetValue()[0];
	}

	if( CompTag->GetName() == "Projection" ) returnable_config->OnlyZero = true;
	else returnable_config->OnlyZero = false;

	for( unsigned int childIndex = 0; childIndex < projComps.size(); ++childIndex )
	{
		if( projComps[childIndex]->GetName() == "DataBins" )
		{
			returnable_config->data_bins = atoi( projComps[childIndex]->GetValue()[0].c_str() );
		}
		else if( projComps[childIndex]->GetName() == "PDFpoints" )
		{
			returnable_config->PDF_points = atoi( projComps[childIndex]->GetValue()[0].c_str() );
		}
		else if( projComps[childIndex]->GetName() == "LogY" )
		{
			string value =  projComps[childIndex]->GetValue()[0];
			if( value == "True" )
			{
				returnable_config->logY = true;
			}
		}
		else if( projComps[childIndex]->GetName() == "Name" )
		{
			returnable_config->observableName = projComps[childIndex]->GetValue()[0];
		}
		else if( projComps[childIndex]->GetName() == "WidthKey" )
		{
			vector<string> widths = StringProcessing::SplitString( projComps[childIndex]->GetValue()[0], ':' );
			if( widths.empty() ) returnable_config->width_key.push_back( atoi(projComps[childIndex]->GetValue()[0].c_str()) );
			for( vector<string>::iterator width_i = widths.begin(); width_i != widths.end(); ++width_i ) returnable_config->width_key.push_back( atoi( width_i->c_str() ) );
		}
		else if( projComps[childIndex]->GetName() == "ColorKey" )
		{
			vector<string> colors = StringProcessing::SplitString( projComps[childIndex]->GetValue()[0], ':' );
			if( colors.empty() ) returnable_config->style_key.push_back( atoi(projComps[childIndex]->GetValue()[0].c_str()) );
			for( vector<string>::iterator color_i = colors.begin(); color_i != colors.end(); ++color_i ) returnable_config->color_key.push_back( atoi( color_i->c_str() ) );
		}
		else if( projComps[childIndex]->GetName() == "StyleKey" )
		{
			vector<string> styles = StringProcessing::SplitString( projComps[childIndex]->GetValue()[0], ':' );
			if( styles.empty() ) returnable_config->style_key.push_back( atoi(projComps[childIndex]->GetValue()[0].c_str()) );
			for( vector<string>::iterator style_i = styles.begin(); style_i != styles.end(); ++style_i ) returnable_config->style_key.push_back( atoi( style_i->c_str() ) );
		}
		else if( projComps[childIndex]->GetName() == "Title" )
		{
			returnable_config->PlotTitle = projComps[childIndex]->GetValue()[0];
		}
		else if( projComps[childIndex]->GetName() == "CompNames" )
		{
			returnable_config->component_names = StringProcessing::SplitString( projComps[childIndex]->GetValue()[0], ':' );
			if( returnable_config->component_names.empty() ) returnable_config->component_names.push_back( projComps[childIndex]->GetValue()[0] );
		}
		else if( projComps[childIndex]->GetName() == "Xmax" )
		{
			returnable_config->xmax = strtod( projComps[childIndex]->GetValue()[0].c_str(), NULL );
		}
		else if( projComps[childIndex]->GetName() == "Xmin" )
		{
			returnable_config->xmin = strtod( projComps[childIndex]->GetValue()[0].c_str(), NULL );
		}
		else if( projComps[childIndex]->GetName() == "Ymax" )
		{
			returnable_config->ymax = strtod( projComps[childIndex]->GetValue()[0].c_str(), NULL );
		}
		else if( projComps[childIndex]->GetName() == "Ymin" )
		{
			returnable_config->ymin = strtod( projComps[childIndex]->GetValue()[0].c_str(), NULL );
		}
		else if( projComps[childIndex]->GetName() == "XTitle" )
		{
			returnable_config->xtitle = TString(projComps[childIndex]->GetValue()[0]);
		}
		else if( projComps[childIndex]->GetName() == "YTitle" )
		{
			returnable_config->ytitle = TString(projComps[childIndex]->GetValue()[0]);
		}
		else if( projComps[childIndex]->GetName() == "TrustNumerical" )
		{
			string value = projComps[childIndex]->GetValue()[0];
			if( value == "True" )
			{
				returnable_config->ScaleNumerical = false;
			}
		}
		else if( projComps[childIndex]->GetName() == "CalcChi2" )
		{
			string value = projComps[childIndex]->GetValue()[0];
			if( value == "True" )
			{
				returnable_config->CalcChi2 = true;
			}
		}
		else if( projComps[childIndex]->GetName() == "DrawPull" )
		{
			string value = projComps[childIndex]->GetValue()[0];
			if( value == "True" )
			{
				returnable_config->DrawPull = true;
			}
		}
		else if( projComps[childIndex]->GetName() == "AddLHCb" )
		{
			string value = projComps[childIndex]->GetValue()[0];
			if( value == "True" )
			{
				returnable_config->addLHCb = true;
			}
		}
		else
		{
			cerr << "Sorry Don't understand: " << projComps[childIndex]->GetName() << " ignoring!" << endl;
		}
	}


	return returnable_config;
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
			for ( unsigned int childIndex = 0; childIndex < plotComponents.size(); ++childIndex )
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
	for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
	{
		if ( children[childIndex]->GetName() == "NumberRepeats" )
		{
			return atoi( children[childIndex]->GetValue()[0].c_str() );
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
		bool want_Trace = false;
		bool change_style = false;
		string Trace_FileName;
		string Strategy;
		int Threads = -1;
		bool integratorTest = true;
		bool NormaliseWeights=false;
		vector< XMLTag* > functionInfo = FunctionTag->GetChildren();
		if ( functionInfo.size() == 0 )
		{
			//Old style - just specifies the function name
			functionName = FunctionTag->GetValue()[0];
		}
		else
		{
			//New style - can have weights
			for ( unsigned int childIndex = 0; childIndex < functionInfo.size(); ++childIndex )
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
				else if ( functionInfo[childIndex]->GetName() == "Trace" )
				{
					want_Trace = true;
					Trace_FileName = functionInfo[childIndex]->GetValue()[0];
				}
				else if ( functionInfo[childIndex]->GetName() == "Strategy" )
				{
					change_style = true;
					Strategy = functionInfo[childIndex]->GetValue()[0];
				}
				else if ( functionInfo[childIndex]->GetName() == "Threads" )
				{
					Threads = atoi( functionInfo[childIndex]->GetValue()[0].c_str() );
				}
				else if ( functionInfo[childIndex]->GetName() == "SetIntegratorTest" )
				{
					if( functionInfo[childIndex]->GetValue()[0] == "True" )
					{
						integratorTest = true;
					}
					else
					{
						integratorTest = false;
					}
				}
				else if ( functionInfo[childIndex]->GetName() == "NormaliseWeights" )
				{
					if( functionInfo[childIndex]->GetValue()[0] == "True" )
					{
						NormaliseWeights = true;
					}
					else
					{
						NormaliseWeights = false;
					}
				}
				else
				{
					cerr << "Unrecognised FitFunction component: " << functionInfo[childIndex]->GetName() << endl;
					exit(1);
				}
			}
		}

		FitFunctionConfiguration* returnable_function = NULL;

		//Make the function
		if (hasWeight)
		{
			cerr <<"Weighted events have been asked for in XML using:\t\t" << weightName << endl ;
			returnable_function = new FitFunctionConfiguration( functionName, weightName );
		}
		else
		{
			returnable_function = new FitFunctionConfiguration( functionName );
		}

		if( want_Trace )
		{
			returnable_function->SetupTrace( Trace_FileName );
		}

		if( change_style )
		{
			returnable_function->SetStrategy( Strategy );
		}

		returnable_function->SetThreads( Threads );
		returnable_function->SetNormaliseWeights( NormaliseWeights );
		returnable_function->SetIntegratorTest( integratorTest );

		return returnable_function;
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << FunctionTag->GetName() << " not \"FitFunction\"" << endl;
		exit(1);
	}
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
		exit(1);
	}
	else
	{
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
					string value = fitComponents[componentIndex]->GetValue()[0];
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
					pdfsAndData.push_back( GetPDFWithData( dataTag, pdfTag, Starting_Value[(unsigned)s_val_index], overloadConfigurator ) );
				} else {
					pdfsAndData.push_back( GetPDFWithData( dataTag, pdfTag, -1, overloadConfigurator ) );
				}
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
		cerr << "No ToFit tags found in config file" << endl;
		exit(1);
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
				constraints.push_back( GetConstraintFunction(constraintTag) );
			}
		}
	}

	return constraints;
}

//Create a ConstraintFunction for the appropriate xml tag
ConstraintFunction * XMLConfigReader::GetConstraintFunction( XMLTag * InputTag )
{
	if ( InputTag->GetName() == "ConstraintFunction" )
	{
		vector< ExternalConstraint* > constraints;
		vector< XMLTag* > functionComponents = InputTag->GetChildren();
		for ( unsigned int componentIndex = 0; componentIndex < functionComponents.size(); ++componentIndex )
		{
			if ( functionComponents[componentIndex]->GetName() == "ExternalConstraint" )
			{
				constraints.push_back( GetExternalConstraint( functionComponents[componentIndex] ) );
			}
		}

		ConstraintFunction* returnable = new ConstraintFunction(constraints);
		while( !constraints.empty() )
		{
			if( constraints.back() != NULL ) delete constraints.back();
			constraints.pop_back();
		}
		return returnable;
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << InputTag->GetName() << "\" not \"ConstraintFunction\"" << endl;
		exit(1);
	}
}

//Create an ExternalConstraint for the appropriate xml tag
ExternalConstraint * XMLConfigReader::GetExternalConstraint( XMLTag * InputTag )
{
	if ( InputTag->GetName() == "ExternalConstraint" )
	{
		string name;
		double value=0.;
		double error=0.;
		vector< XMLTag* > externalComponents = InputTag->GetChildren();
		for ( unsigned int componentIndex = 0; componentIndex < externalComponents.size(); ++componentIndex )
		{
			if ( externalComponents[componentIndex]->GetName() == "Name" )
			{
				name = externalComponents[componentIndex]->GetValue()[0];
			}
			else if ( externalComponents[componentIndex]->GetName() == "Value" )
			{
				value = strtod( externalComponents[componentIndex]->GetValue()[0].c_str(), NULL );
			}
			else if ( externalComponents[componentIndex]->GetName() == "Error" )
			{
				error = strtod( externalComponents[componentIndex]->GetValue()[0].c_str(), NULL );
			}
			else
			{
				cerr << "Unrecognised constraint component: " << externalComponents[componentIndex]->GetName() << endl;
				exit(1);
			}
		}

		return new ExternalConstraint( name, value, error );
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << InputTag->GetName() << "\" not \"ExternalConstraint\"" << endl;
		exit(1);
	}
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
		for ( unsigned int parameterIndex = 0; parameterIndex < parameters.size(); ++parameterIndex )
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
		for ( unsigned int nameIndex = 0; nameIndex < names.size(); ++nameIndex )
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
		string blindString = "Uninitialised";
		double value = 0.0;
		double minimum = 0.0;
		double maximum = 0.0;
		double blindScale = 0.0 ;
		double blindOffset = 0.0 ;
		double stepSize = -1.;
		bool hasValue = false;
		bool hasMaximum = false;
		bool hasMinimum = false;
		bool hasBlindString = false ;
		bool hasBlindScale = false ;

		//Loop over the tag children, which correspond to the parameter elements
		vector< XMLTag* > elements = InputTag->GetChildren();
		for ( unsigned int elementIndex = 0; elementIndex < elements.size(); ++elementIndex )
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
			else if ( name == "Type" )
			{
				type = elements[elementIndex]->GetValue()[0];
			}
			else if ( name == "Unit" )
			{
				unit = elements[elementIndex]->GetValue()[0];
			}
			else if ( name == "BlindString" )
			{
				hasBlindString = true ;
				blindString = elements[elementIndex]->GetValue()[0];
			}
			else if ( name == "BlindScale" )
			{
				hasBlindScale =  true ;
				blindScale = strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL );
			}
			else if ( name == "StepSize" )
			{
				stepSize = strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL );
			}
			else
			{
				cerr << "Unrecognised physics parameter configuration: " << name << endl;
				exit(1);
			}
		}

		//See of blinding has been specified, and if so construct the blinding offset
		if( hasBlindString && hasBlindScale )
		{
			blindOffset = Blinder::getBlindOffset( blindString.c_str(), blindScale ) ;
		}
		else if( (hasBlindString && ! hasBlindScale) || (! hasBlindString && hasBlindScale) )
		{
			cerr << "Blinding information incomplete for parameter: " << ParameterName << " Ignoring blinding for this parameter" << endl;
		}

		//Now construct the physics parameter
		if (hasValue)
		{
			if ( hasMaximum && hasMinimum )
			{
				if ( ( ( fabs(maximum - 0.0) < DOUBLE_TOLERANCE ) && ( ( fabs(minimum - 0.0) < DOUBLE_TOLERANCE ) ) ) || type == "Unbounded" )
				{
					//Unbounded parameter
					PhysicsParameter * p = new PhysicsParameter( ParameterName, value, stepSize, type, unit );
					if( hasBlindString && hasBlindScale )
					{
						p->SetBlindOffset( blindOffset );
						p->SetBlindingInfo( blindString.c_str(), blindScale );
					}
					return p;
				}
				else
				{
					//Bounded parameter
					PhysicsParameter * p =  new PhysicsParameter( ParameterName, value, minimum, maximum, stepSize, type, unit );
					if( hasBlindString && hasBlindScale )
					{
						p->SetBlindOffset( blindOffset ) ;
						p->SetBlindingInfo( blindString.c_str(), blindScale );
					}
					return p;
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
					PhysicsParameter * p =  new PhysicsParameter( ParameterName, value, stepSize, type, unit );
					if( hasBlindString && hasBlindScale )
					{
						p->SetBlindOffset( blindOffset );
						p->SetBlindingInfo( blindString.c_str(), blindScale );
					}
					return p ;
				}
			}
		}
		else
		{
			cerr << "Ambiguous definition for parameter " << ParameterName << endl;
			cerr << "\t\tTry Adding A <Value> tag to the Parameter definition." << endl;
			exit(1);
		}
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << InputTag->GetName() << "\" not \"PhysicsParameter\"" << endl;
		exit(1);
	}
	exit(-99);
}

//Collect the information needed to make a data set
PDFWithData * XMLConfigReader::GetPDFWithData( XMLTag * DataTag, XMLTag * FitPDFTag, int Starting_Value, XMLTag* overloadConfigurator )
{
	//Check the tag actually is a data set
	if ( DataTag->GetName() == "DataSet" )
	{
		string dataSource = "Uninitialised";
		long numberEvents = 0;
		string cutString = "";
		PhaseSpaceBoundary * dataBoundary=NULL;
		vector<string> dataArguments, argumentNames;
		bool boundaryFound = false;
		vector< DataSetConfiguration* > dataSetMakers;
		bool generatePDFFlag = false;
		IPDF * generatePDF=NULL;
		XMLTag* generatePDFXML=NULL;
		XMLTag* pdfOptionXML=NULL;

		//Retrieve the data set config
		vector< XMLTag* > dataComponents = DataTag->GetChildren();
		for ( unsigned int dataIndex = 0; dataIndex < dataComponents.size(); ++dataIndex )
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
			else if ( name == "CutString" )
			{
				cutString = dataComponents[dataIndex]->GetValue()[0];
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
				if( dataBoundary == NULL ) dataBoundary = GetPhaseSpaceBoundary( dataComponents[dataIndex] );
				else
				{
					cerr << "CANNOT DO MULTIPLE PHASESPACES YET" << endl;
					exit(8972);
				}
			}
			else if ( name == "CommonPhaseSpace" )
			{
				boundaryFound = true;
				if( dataBoundary == NULL ) dataBoundary = this->FindCommonPhaseSpace( dataComponents[dataIndex] );
				else
				{
					cerr << "CANNOT DO MULTIPLE PHASESPACES YET" << endl;
					exit(8972);
				}
			}
			else if ( name == "PDFConfiguration" )
			{
				pdfOptionXML = dataComponents[dataIndex];
			}
			else if ( name == "PDF" || name == "SumPDF" || name == "NormalisedSumPDF" || name == "ProdPDF" )
			{
				generatePDFXML = dataComponents[dataIndex];
				generatePDFFlag = true;
			}
			else if ( name == "StartingEntry" )
			{
				if( Starting_Value < 0 )
				{
					Starting_Value = atoi( dataComponents[dataIndex]->GetValue()[0].c_str() );
				}
			}
			else
			{
				cerr << "Unrecognised data set component: " << name << endl;
				exit(1);
			}
		}

		//dataBoundary->Print();

		if( generatePDFFlag == true )
		{
			generatePDF = GetPDF( generatePDFXML, dataBoundary, pdfOptionXML );
		}

		//Return the collection of configuration information - data generation will happen later
		if(!boundaryFound)
		{
			cerr << "DataSet defined without PhaseSpaceBoundary" << endl;
			exit(1);
		}

		//If there are no separate data sources, go for backwards compatibility
		if( dataSetMakers.size() < 1 )
		{
			DataSetConfiguration * oldStyleConfig;

			if(generatePDFFlag)
			{
				oldStyleConfig = new DataSetConfiguration( dataSource, numberEvents, cutString, dataArguments, argumentNames, generatePDF, dataBoundary );
			}
			else
			{
				if( Starting_Value > 0 )
				{
					oldStyleConfig = new DataSetConfiguration( dataSource, numberEvents, cutString, dataArguments, argumentNames, Starting_Value, dataBoundary );
				}
				else
				{
					oldStyleConfig = new DataSetConfiguration( dataSource, numberEvents, cutString, dataArguments, argumentNames, 0 , dataBoundary );
				}
			}

			dataSetMakers.push_back(oldStyleConfig);
		}

		//Make the objects
		IPDF * fitPDF = GetPDF( FitPDFTag, dataBoundary, overloadConfigurator );
		if( generatePDF != NULL ) delete generatePDF;
		PDFWithData* returnable = new PDFWithData( fitPDF, dataBoundary, dataSetMakers );
		if( fitPDF != NULL ) delete fitPDF;
		while( !dataSetMakers.empty() )
		{
			if( dataSetMakers.back() != NULL ) delete dataSetMakers.back();
			dataSetMakers.pop_back();
		}
		if( dataBoundary != NULL ) delete dataBoundary;
		return returnable;
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
		string cutString = "";
		long numberEvents = 0;
		vector<string> dataArguments, argumentNames;
		bool generatePDFFlag = false;
		IPDF * generatePDF=NULL;
		XMLTag* PDFXML=NULL;
		XMLTag* pdfOptions=NULL;

		//Retrieve the data set config
		vector< XMLTag* > dataComponents = DataTag->GetChildren();
		for ( unsigned int dataIndex = 0; dataIndex < dataComponents.size(); ++dataIndex )
		{
			string name = dataComponents[dataIndex]->GetName();
			if ( name == "Source" )
			{
				dataSource = dataComponents[dataIndex]->GetValue()[0];
			}
			else if ( name == "CutString" )
			{
				cutString = dataComponents[dataIndex]->GetValue()[0];
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
				PDFXML = dataComponents[dataIndex];
				generatePDFFlag = true;
			}
			else if( name == "CommonPDF" )
			{
				PDFXML = this->FindCommonPDFXML();
				generatePDFFlag = true;
			}
			else if( name == "PDFConfiguration" )
			{
				pdfOptions = dataComponents[dataIndex];
			}
			else
			{
				cerr << "Unrecognised data set component: " << name << endl;
				exit(1);
			}
		}

		if( generatePDFFlag )	generatePDF = GetPDF( PDFXML, DataBoundary, pdfOptions );

		if (generatePDFFlag)
		{
			return new DataSetConfiguration( dataSource, numberEvents, cutString, dataArguments, argumentNames, generatePDF );
		}
		else
		{
			return new DataSetConfiguration( dataSource, numberEvents, cutString, dataArguments, argumentNames );
		}
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << DataTag->GetName() << "\" not \"Subset\"" << endl;
		exit(1);
	}
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
						PhaseSpaceBoundary_vec.push_back( GetPhaseSpaceBoundary( allFitDataSets[DataSetNum] ) );
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
			returnablePhaseSpace = this->GetPhaseSpaceBoundary( children[childIndex]->GetChildren()[0] );
		}
	}

	if( returnablePhaseSpace == NULL ) returnablePhaseSpace = new PhaseSpaceBoundary( vector<string>() );

	for( unsigned int phaseChild=0; phaseChild< phaseChildren.size(); ++phaseChild )
	{
		if( phaseChildren[phaseChild]->GetName() == "Observable" )
		{
			string temp; (void) temp;
			IConstraint* new_observable = this->GetConstraint( phaseChildren[phaseChild], temp );
			//new_observable->Print();
			//	Add Constraint if non exists and Overwrite if it already exists
			returnablePhaseSpace->AddConstraint( new_observable->GetName(), new_observable, true );
			delete new_observable;
		}
	}

	return returnablePhaseSpace;
}

IPDF* XMLConfigReader::FindCommonPDF( PhaseSpaceBoundary* override )
{
	for( unsigned int childIndex=0; childIndex< children.size(); ++childIndex )
	{
		if( children[childIndex]->GetName() == "CommonPDF" )
		{
			if( override == NULL )
			{
				return this->GetPDF( children[childIndex]->GetChildren()[0], this->FindCommonPhaseSpace(), NULL );
			}
			else
			{
				return this->GetPDF( children[childIndex]->GetChildren()[0], override, NULL );
			}
		}
	}
	return NULL;
}

XMLTag* XMLConfigReader::FindCommonPDFXML()
{
	for( unsigned int childIndex=0; childIndex< children.size(); ++childIndex )
	{
		if( children[childIndex]->GetName() == "CommonPDF" )
		{
			return children[childIndex]->GetChildren()[0];
		}
	}
	return NULL;
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

//Make a PhaseSpaceBoundary from the appropriate xml tag
PhaseSpaceBoundary * XMLConfigReader::GetPhaseSpaceBoundary( XMLTag* InputTag )
{
	//Check the tag is actually a PhaseSpaceBoundary
	if ( InputTag->GetName() == "PhaseSpaceBoundary" )
	{
		vector< IConstraint* > constraints;
		vector<string> names;
		string name;

		//Create each single bound
		vector< XMLTag* > constraintTags = InputTag->GetChildren();
		for ( unsigned int boundIndex = 0; boundIndex < constraintTags.size(); ++boundIndex )
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
		for( unsigned int nameIndex = 0; nameIndex < names.size(); ++nameIndex )
		{
			newBoundary->AddConstraint( names[nameIndex], constraints[nameIndex] );
		}
		while(!constraints.empty())
		{
			if( constraints.back() != NULL ) delete constraints.back();
			constraints.pop_back();
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
		string tf1;
		double minimum = -9999.0;
		double maximum = -9999.0;
		vector<double> allValues;

		//Loop over the tag children, which correspond to the parameter elements
		vector< XMLTag* > elements = InputTag->GetChildren();
		for ( unsigned int elementIndex = 0; elementIndex < elements.size(); ++elementIndex )
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
			else if ( elements[elementIndex]->GetName() == "TF1" )
			{
				tf1 = elements[elementIndex]->GetValue()[0];
			}
			else
			{
				cerr << "Unrecognised constraint configuration: " <<  elements[elementIndex]->GetName() << endl;
				exit(1);
			}
		}

		if( tf1.empty() ) tf1 = Name;

		IConstraint* returnable_const=NULL;

		//If there are discrete values, make a discrete constraint
		if ( allValues.size() > 0 )
		{
			returnable_const = new ObservableDiscreteConstraint( Name, allValues, unit, tf1 );
		}
		else
		{
			returnable_const = new ObservableContinuousConstraint( Name, minimum, maximum, unit, tf1 );
		}

		return returnable_const;
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << InputTag->GetName() << "\" not \"Observable\"" << endl;
		exit(18);
	}
}

//Create a PDF from an appropriate xml tag
IPDF * XMLConfigReader::GetNamedPDF( XMLTag * InputTag, XMLTag* overloadConfigurator, bool print )
{
	IPDF* returnable_NamedPDF;
	//Check the tag actually is a PDF
	if( InputTag->GetName() == "PDF" )
	{
		vector< XMLTag* > pdfConfig = InputTag->GetChildren();
		XMLTag* common = this->FindCommonPDFConfiguratorXML();
		if( common != NULL )
		{
			vector<XMLTag*> optional = common->GetChildren();
			for( unsigned int i=0; i< optional.size(); ++i )
			{
				pdfConfig.push_back( optional[i] );
			}
		}
		string name;
		vector<string> observableNames, parameterNames;
		PDFConfigurator* configurator = new PDFConfigurator();

		//Load the PDF configuration
		for ( unsigned int configIndex = 0; configIndex < pdfConfig.size(); ++configIndex )
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
			else if ( pdfConfig[configIndex]->GetName() == "ParameterSubstitution" )
			{
				configurator->addParameterSubstitution( pdfConfig[configIndex]->GetValue()[0] );
			}
			else if ( pdfConfig[configIndex]->GetName() == "AppendParameterNames" )
			{
				configurator->appendParameterNames( pdfConfig[configIndex]->GetValue()[0] );
			}
			else if ( pdfConfig[configIndex]->GetName() == "ConfigurationParameter" )
			{
				configurator->addConfigurationParameter( pdfConfig[configIndex]->GetValue()[0] );
			}
			else if ( pdfConfig[configIndex]->GetName() == "FractionName" )
			{
				//	Silently ignore this option
			}
			else
			{
				cerr << "(1)Unrecognised PDF configuration: " << pdfConfig[configIndex]->GetName() << endl;
				exit(1);
			}
		}

		if( !configurator->empty() )
		{
			cout << "PDFConfigurator: ";
			if( overloadConfigurator != NULL ) cout << "Additional ";
			cout << "Options Passed to this PDF" << endl;
			configurator->Print();
		}

                if( overloadConfigurator != NULL )
                {
                        pdfConfig = overloadConfigurator->GetChildren();
                }

		//Load the PDF configuration
		for ( unsigned int configIndex = 0; configIndex < pdfConfig.size(); ++configIndex )
		{
			if ( pdfConfig[configIndex]->GetName() == "ParameterSubstitution" )
			{
				configurator->addParameterSubstitution( pdfConfig[configIndex]->GetValue()[0] );
			}
			else if ( pdfConfig[configIndex]->GetName() == "AppendParameterNames" )
			{
				configurator->appendParameterNames( pdfConfig[configIndex]->GetValue()[0] );
			}
			else if ( pdfConfig[configIndex]->GetName() == "ConfigurationParameter" )
			{
				configurator->addConfigurationParameter( pdfConfig[configIndex]->GetValue()[0] );
			}
			else if ( pdfConfig[configIndex]->GetName() == "FractionName" )
			{
				//      Silently ignore this option
			}
			else
			{
				cerr << "(1)Unrecognised PDF configuration: " << pdfConfig[configIndex]->GetName() << endl;
				exit(1);
			}
		}

		//Check if the name is recognised as a PDF
		returnable_NamedPDF = ClassLookUp::LookUpPDFName( name, configurator );
	}
	else
	{
		cerr << "Incorrect tag provided: \"" << InputTag->GetName() << "\" not \"PDF\"" << endl;
		exit(1);
	}

	//	returnable_NamedPDF->SetRandomFunction( GetSeed() );
	return returnable_NamedPDF;
}

//Create a SumPDF from an appropriate xml tag
IPDF * XMLConfigReader::GetSumPDF( XMLTag * InputTag, PhaseSpaceBoundary * InputBoundary, XMLTag* overloadConfigurator, bool print )
{
	IPDF* returnable_SUMPDF;
	vector< IPDF* > componentPDFs;
	//Check the tag actually is a PDF
	if ( InputTag->GetName() == "SumPDF" )
	{
		vector< XMLTag* > pdfConfig = InputTag->GetChildren();
		string fractionName = "unspecified";

		if( overloadConfigurator != NULL )
		{
			vector<XMLTag*> optional = overloadConfigurator->GetChildren();
			vector<string> names;
			string fracstr="FractionName";
			for( unsigned int i=0; i< optional.size(); ++i ) names.push_back( optional[i]->GetName() );
			for( int i=int(optional.size()-1); i != -1; --i )
			{
				if( optional[i]->GetName() == "FractionName" )
				{
					cout << optional[i]->GetValue()[0] << endl;
					fractionName = optional[(unsigned)i]->GetValue()[0];
					overloadConfigurator->RemoveChild( i );
					break;
				}
				cout << optional[i]->GetName() << endl;
			}
		}

		unsigned int i=1;
		//Load the PDF configuration
		for ( unsigned int configIndex = 0; configIndex < pdfConfig.size(); ++configIndex )
		{
			if( pdfConfig[configIndex]->GetName() == "FractionName" && fractionName == "unspecified" )
			{
				fractionName = pdfConfig[configIndex]->GetValue()[0];
			}
			else
			{
				//cout << "--------------------------------------------------------------------------------------------------------------------" << endl;
				//cout << "Constructing SumPDF Component: " << i << endl;
				//cout << endl;
				componentPDFs.push_back( GetPDF( pdfConfig[configIndex], InputBoundary, overloadConfigurator, print ) );
				++i;
			}
		}

		//Check there are two component PDFs to sum
		if( componentPDFs.size() == 2 )
		{
			if ( fractionName == "unspecified" )
			{
				cerr << "Please provide fractionName" << endl;
				exit(-87632);
			}
			else
			{
				returnable_SUMPDF = new SumPDF( componentPDFs[0], componentPDFs[1], InputBoundary, fractionName );
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

	returnable_SUMPDF->SetName( "Sum" );
	returnable_SUMPDF->SetLabel( "Sum_("+componentPDFs[0]->GetLabel()+")+("+componentPDFs[1]->GetLabel()+")" );

	delete componentPDFs[0]; delete componentPDFs[1];
	//	returnable_SUMPDF->SetRandomFunction( GetSeed() );
	return returnable_SUMPDF;
}

//Create a NormalisedSumPDF from an appropriate xml tag
IPDF * XMLConfigReader::GetNormalisedSumPDF( XMLTag * InputTag, PhaseSpaceBoundary * InputBoundary, XMLTag* overloadConfigurator, bool print )
{
	IPDF* returnable_NormPDF = NULL;
	vector< IPDF* > componentPDFs;
	//Check the tag actually is a PDF
	if ( InputTag->GetName() == "NormalisedSumPDF" )
	{
		vector< XMLTag* > pdfConfig = InputTag->GetChildren();
		string fractionName = "unspecified";

		if( overloadConfigurator != NULL )
		{
			vector<XMLTag*> optional = overloadConfigurator->GetChildren();
			vector<string> names;
			string fracstr="FractionName";
			for( unsigned int i=0; i< optional.size(); ++i ) names.push_back( optional[i]->GetName() );
			for( int i=int(optional.size()-1); i > -1; --i )
			{
				if( optional[i]->GetName() == "FractionName" )
				{
					//cout << optional[i]->GetValue()[0] << endl;
					fractionName = optional[(unsigned)i]->GetValue()[0];
					overloadConfigurator->RemoveChild( i );
					break;
				}
			}
		}

		unsigned int i=1;
		//Load the PDF configuration
		for ( unsigned int configIndex = 0; configIndex < pdfConfig.size(); ++configIndex )
		{
			if( pdfConfig[configIndex]->GetName() == "FractionName" && fractionName == "unspecified" )
			{
				fractionName = pdfConfig[configIndex]->GetValue()[0];
			}
			else
			{
				//cout << "--------------------------------------------------------------------------------------------------------------------" << endl;
				//cout << "Constructing NormalisedSum Component: " << i << endl;
				//cout << endl;
				componentPDFs.push_back( GetPDF( pdfConfig[configIndex], InputBoundary, overloadConfigurator, print ) );
				++i;
			}
		}

		//Check there are two component PDFs to sum
		if ( componentPDFs.size() == 2 )
		{
			if ( fractionName == "unspecified" )
			{
				cerr << "Error Must Provide a FractionName for this PDF!\tNormalisedSum_(" << componentPDFs[0]->GetLabel()<<")+("<<componentPDFs[1]->GetLabel()<< ")" << endl;
				exit(1);
			}
			else
			{
				returnable_NormPDF =  new NormalisedSumPDF( componentPDFs[0], componentPDFs[1], InputBoundary, fractionName );
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

	//cout << returnable_NormPDF->GetName() << "\t\t" << returnable_NormPDF->GetLabel() << endl;

	delete componentPDFs[0]; delete componentPDFs[1];
	//	returnable_NormPDF->SetRandomFunction( GetSeed() );
	return returnable_NormPDF;
}

//Create a ProdPDF from an appropriate xml tag
IPDF * XMLConfigReader::GetProdPDF( XMLTag * InputTag, PhaseSpaceBoundary * InputBoundary, XMLTag* overloadConfigurator, bool print )
{
	IPDF* returnable_ProdPDF;
	vector< IPDF* > componentPDFs;
	//Check the tag actually is a PDF
	if ( InputTag->GetName() == "ProdPDF" )
	{
		unsigned int i=1;
		vector< XMLTag* > pdfConfig = InputTag->GetChildren();
		//Load the PDF configuration
		for ( unsigned int configIndex = 0; configIndex < pdfConfig.size(); ++configIndex )
		{
			//cout << "--------------------------------------------------------------------------------------------------------------------" << endl;
			//cout << "Constructing ProdPDF Component: " << i << endl;
			//cout << endl;
			componentPDFs.push_back( GetPDF( pdfConfig[configIndex], InputBoundary, overloadConfigurator, print ) );
			++i;
		}

		//Check there are two component PDFs to sum
		if ( componentPDFs.size() == 2 )
		{
			returnable_ProdPDF = new ProdPDF( componentPDFs[0], componentPDFs[1] );
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

	delete componentPDFs[0]; delete componentPDFs[1];
	//	returnable_ProdPDF->SetRandomFunction( GetSeed() );
	return returnable_ProdPDF;
}

//Choose one of the PDF instantiation methods
IPDF * XMLConfigReader::GetPDF( XMLTag * InputTag, PhaseSpaceBoundary * InputBoundary, XMLTag* overloadConfigurator, bool print )
{
	if( overloadConfigurator != NULL && print )
	{
		vector<XMLTag*> pdfConfig = overloadConfigurator->GetChildren();

		PDFConfigurator* configurator = new PDFConfigurator();

		//Load the PDF configuration
		for( unsigned int configIndex = 0; configIndex < pdfConfig.size(); ++configIndex )
		{
			if ( pdfConfig[configIndex]->GetName() == "ParameterSubstitution" )
			{
				configurator->addParameterSubstitution( pdfConfig[configIndex]->GetValue()[0] );
			}
			else if ( pdfConfig[configIndex]->GetName() == "AppendParameterNames" )
			{
				configurator->appendParameterNames( pdfConfig[configIndex]->GetValue()[0] );
			}
			else if ( pdfConfig[configIndex]->GetName() == "ConfigurationParameter" )
			{
				configurator->addConfigurationParameter( pdfConfig[configIndex]->GetValue()[0] );
			}
		}

		if( ! configurator->empty() )
		{
			cout << "====================================================================================================================" << endl;
			cout << "PDFConfigurator: Config Parameters Passed to CommonPDF" << endl;
			configurator->Print();
			cout << endl;
		}
		delete configurator;
	}

	IPDF* returnable_pdf;
	if ( InputTag->GetName() == "PDF" )
	{
		returnable_pdf = GetNamedPDF( InputTag, overloadConfigurator, false );
	}
	else if ( InputTag->GetName() == "SumPDF" )
	{
		returnable_pdf = GetSumPDF( InputTag, InputBoundary, overloadConfigurator, false );
	}
	else if ( InputTag->GetName() == "NormalisedSumPDF" )
	{
		returnable_pdf = GetNormalisedSumPDF( InputTag, InputBoundary, overloadConfigurator, false );
	}
	else if ( InputTag->GetName() == "ProdPDF" )
	{
		returnable_pdf = GetProdPDF( InputTag, InputBoundary, overloadConfigurator, false );
	}
	else
	{
		cerr << "(2)Unrecognised PDF configuration: " << InputTag->GetName() << endl;
		exit(1);
	}

	returnable_pdf->SetRandomFunction( int(GetSeed()) );
	returnable_pdf->SetMCCacheStatus( false );
	return returnable_pdf;
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
			int config=1;
			vector< XMLTag* > newchildren = children[childIndex]->GetChildren();
			for ( unsigned int newchildIndex = 0; newchildIndex < newchildren.size(); ++newchildIndex )
			{
				string name = newchildren[newchildIndex]->GetName();
				if ( name == "Name" )
				{
					precalculatorName = newchildren[newchildIndex]->GetValue()[0];
				}
				else if ( name == "WeightName" )
				{
					weightName = newchildren[newchildIndex]->GetValue()[0];
				}
				else if ( name == "Config" )
				{
					config = atoi( newchildren[newchildIndex]->GetValue()[0].c_str() );
				}
				else if ( name == "OutputFile" )
				{
					filename = newchildren[newchildIndex]->GetValue()[0];
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
			return preconfig;
		}
	}

	return NULL;
}


//	Return the Integer Seed used in RapidFit XML tag <Seed>SomeInt</Seed>
unsigned int XMLConfigReader::GetSeed()
{
	if( seed < 0 )
	{
		//Find the NumberRepeats tag
		for ( unsigned int childIndex = 0; childIndex < children.size(); ++childIndex )
		{
			if ( children[childIndex]->GetName() == "Seed" )
			{
				seed = abs( atoi( children[childIndex]->GetValue( )[0].c_str() ) );
				cout << "Using seed: " << seed << " from input file." << endl;
				return unsigned(seed);
			}
		}
		seed = 0 ;
		//If no such tag is found, report
		cout << "Seed tag not found in config file, defaulting to TRandom3(0)." << endl;
	} else  return unsigned(seed);
	return 0;
}

//	Set a new TRandom seed that is returned by the XMLFile
void XMLConfigReader::SetSeed( unsigned int new_seed )
{
	seed = (int)new_seed;
}

pair<ScanParam*, ScanParam*> XMLConfigReader::Get2DScanParam( XMLTag * InputTag )
{
	pair<ScanParam*, ScanParam*> Returnable_Pair;
	ScanParam* Param_1=NULL;
	ScanParam* Param_2=NULL;
	if ( ( InputTag->GetName() == "2DScan" ) || ( InputTag->GetName() == "TwoDScan" ) )
	{
		bool param1 = false;
		bool param2 = false;
		//Loop over the tag children, which correspond to the parameter elements
		vector< XMLTag* > elements = InputTag->GetChildren();
		for ( unsigned int elementIndex = 0; elementIndex < elements.size(); ++elementIndex )
		{
			string name = elements[elementIndex]->GetName();
			if ( name == "X_Param" )  {
				Param_1 = XMLConfigReader::GetScanParam( elements[elementIndex] );
				param1 = true;
			} else if ( name == "Y_Param" )  {
				Param_2 = XMLConfigReader::GetScanParam( elements[elementIndex] );
				param2 = true;
			} else {
				cerr << "Bad Configuration for 2DScan, XMLTag " << name << " is bad." << endl;
				exit(-1);
			}
		}

		if( !param1 || !param2 )
		{
			cerr << "Sorry 2DScan not correctly defined" <<endl;
			exit(-3);
		}

		Returnable_Pair.first = Param_1;
		Returnable_Pair.second = Param_2;
		return Returnable_Pair;
	}
	return Returnable_Pair;
}

//Make a physics parameter from an appropriate XML tag
ScanParam * XMLConfigReader::GetScanParam( XMLTag * InputTag )
{
	ScanParam * ScanParamLocal;
	//Check the tag is actually a physics parameter
	if ( (InputTag->GetName() == "Scan") || (InputTag->GetName() == "X_Param") || (InputTag->GetName() == "Y_Param") )
	{
		//Create some default values;
		vector<string> name_tag, type;
		vector<double> maximum, minimum;
		vector<int> points, sigma;

		//Loop over the tag children, which correspond to the parameter elements
		vector< XMLTag* > elements = InputTag->GetChildren();
		for ( unsigned int elementIndex = 0; elementIndex < elements.size(); ++elementIndex )
		{
			string name = elements[elementIndex]->GetName();
			if ( name == "Name" )
			{
				name_tag.push_back( elements[elementIndex]->GetValue()[0].c_str() );
			}
			else if ( name == "Minimum" )
			{
				minimum.push_back( strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL ) );
			}
			else if ( name == "Maximum" )
			{
				maximum.push_back( strtod( elements[elementIndex]->GetValue()[0].c_str(), NULL ) );
			}
			else if ( name == "Points" )
			{
				points.push_back( atoi( elements[elementIndex]->GetValue()[0].c_str() ) );
			}
			else if ( name == "Sigma" )
			{
				sigma.push_back( atoi( elements[elementIndex]->GetValue()[0].c_str() ) );
			}
			else
			{
				cerr << "Unrecognised Scan configuration: " << name << endl;
				exit(1);
			}
		}

		//Check for ambiguity
		if ( ( maximum.empty() || minimum.empty() ) && ( sigma.empty() ) )
		{
			cerr << "Ambiguous Scan definition: " << name_tag[0] << " has value and ";
			if ( maximum.empty() )
			{
				cerr << " maximum, but not minimum";
			}
			if ( minimum.empty() )
			{
				cerr << " minimum, but not maximum";
			}
			cerr << " defined" << endl;
			exit(1);
		}

		if( ( !maximum.empty() || !minimum.empty() ) && !sigma.empty() )
		{
			cerr << "Ambiguous Scan definition:\t" << name_tag[0] <<"\t" << type[0] << endl;
			cerr << "has both numerical and sigma based limits, defaulting to numerical!" << endl;
			while( !sigma.empty() ) { sigma.pop_back(); }
		}

		ScanParamLocal = new ScanParam( name_tag, maximum, minimum, sigma, points );
		return ScanParamLocal;

	} else	{ cerr << "Unreconised Scan Config: " << InputTag->GetName() << endl; exit(1); }

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
						DataSetSize += (int) atoi( dataComponents[dataIndex]->GetValue()[0].c_str() );
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
						DataSetSizes.push_back( (int) atoi( dataComponents[dataIndex]->GetValue()[0].c_str() ) );
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
						StartEntries.push_back( (int) atoi( dataComponents[dataIndex]->GetValue()[0].c_str() ) );
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

void XMLConfigReader::SetDebug( DebugClass* input_debug )
{
	if( input_debug != NULL )
	{
		if( debug != NULL ) delete debug;
		debug = new DebugClass( *input_debug );

		if( debug->DebugThisClass("XMLConfigReader") )
		{
			debug->SetStatus(true);
			cout << "XMLConfigReader: Debugging Enabled!" << endl;
		}
		else
		{
			debug->SetStatus(false);
		}
	}
	else
	{
		if( debug != NULL ) delete debug;
		debug = new DebugClass(false);
	}
}

