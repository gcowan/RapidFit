
#include "XMLTag.h"
#include "Blinder.h"
#include "IConstraint.h"
#include "ObservableDiscreteConstraint.h"
#include "ObservableContinuousConstraint.h"
#include "PhysicsParameter.h"
#include "PhaseSpaceBoundary.h"
#include "XMLObjectGenerator.h"
#include "ParameterSet.h"
#include "FitFunctionConfiguration.h"
#include "RapidRun.h"
#include "IPDF.h"
#include "ClassLookUp.h"
#include "StringProcessing.h"
#include "DataSetConfiguration.h"
#include "MinimiserConfiguration.h"
#include "OutputConfiguration.h"
#include "ExternalConstraint.h"
#include "ExternalConstMatrix.h"
#include "ConstraintFunction.h"
#include "ComponentPlotter_config.h"

#include <string>
#include <vector>
#include <stdio.h>
#include <iostream>

#define DOUBLE_TOLERANCE 1E-6

using namespace::std;

//Make a physics parameter from an appropriate XML tag
PhysicsParameter * XMLObjectGenerator::GetPhysicsParameter( XMLTag * InputTag, string & ParameterName )
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
				ParameterName = XMLTag::GetStringValue( elements[elementIndex] );
			}
			else if ( name == "Value" )
			{
				hasValue = true;
				value = XMLTag::GetDoubleValue( elements[elementIndex] );
			}
			else if ( name == "Minimum" )
			{
				hasMinimum = true;
				minimum = XMLTag::GetDoubleValue( elements[elementIndex] );
			}
			else if ( name == "Maximum" )
			{
				hasMaximum = true;
				maximum = XMLTag::GetDoubleValue( elements[elementIndex] );
			}
			else if ( name == "Type" )
			{
				type = XMLTag::GetStringValue( elements[elementIndex] );
			}
			else if ( name == "Unit" )
			{
				unit = XMLTag::GetStringValue( elements[elementIndex] );
			}
			else if ( name == "BlindString" )
			{
				hasBlindString = true ;
				blindString = XMLTag::GetStringValue( elements[elementIndex] );
			}
			else if ( name == "BlindScale" )
			{
				hasBlindScale =  true ;
				blindScale = XMLTag::GetDoubleValue( elements[elementIndex] );
			}
			else if ( name == "StepSize" )
			{
				stepSize = XMLTag::GetDoubleValue( elements[elementIndex] );
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

//Make a PhaseSpaceBoundary from the appropriate xml tag
PhaseSpaceBoundary* XMLObjectGenerator::GetPhaseSpaceBoundary( XMLTag* InputTag )
{
	//Check the tag is actually a PhaseSpaceBoundary
	if ( InputTag->GetName() == "PhaseSpaceBoundary" || InputTag->GetName() == "CommonPhaseSpace" )
	{
		vector< IConstraint* > constraints;
		vector<string> names;
		string name;

		//Create each single bound
		vector< XMLTag* > constraintTags = InputTag->GetChildren();
		for ( unsigned int boundIndex = 0; boundIndex < constraintTags.size(); ++boundIndex )
		{
			IConstraint * newConstraint = XMLObjectGenerator::GetConstraint( constraintTags[boundIndex], name );
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
IConstraint * XMLObjectGenerator::GetConstraint( XMLTag * InputTag, string & Name )
{
	//Check the tag is actually a single bound
	if ( InputTag->GetName() == "Observable" )
	{
		//Create some default values;
		Name = "Uninitialised";
		string unit = "Uninitialised";
		string tf1;
		double minimum = -99999.0;
		double maximum = -99999.0;
		vector<double> allValues;
		bool gotValue=false;
		bool gotMin=false;
		bool gotMax=false;
		//Loop over the tag children, which correspond to the parameter elements
		vector< XMLTag* > elements = InputTag->GetChildren();
		for ( unsigned int elementIndex = 0; elementIndex < elements.size(); ++elementIndex )
		{
			if ( elements[elementIndex]->GetName() == "Name" )
			{
				Name = XMLTag::GetStringValue( elements[elementIndex] );
			}
			else if ( elements[elementIndex]->GetName() == "Minimum" )
			{
				gotMin=true;
				minimum = XMLTag::GetDoubleValue( elements[elementIndex] );
			}
			else if ( elements[elementIndex]->GetName() == "Maximum" )
			{
				gotMax=true;
				maximum = XMLTag::GetDoubleValue( elements[elementIndex] );
			}
			else if ( elements[elementIndex]->GetName() == "Value" )
			{
				gotValue=true;
				allValues.push_back( XMLTag::GetDoubleValue( elements[elementIndex] ) );
			}
			else if ( elements[elementIndex]->GetName() == "Unit" )
			{
				unit = XMLTag::GetStringValue( elements[elementIndex] );
			}
			else if ( elements[elementIndex]->GetName() == "TF1" )
			{
				tf1 = XMLTag::GetStringValue( elements[elementIndex] );
			}
			else
			{
				cerr << "Unrecognised constraint configuration: " <<  elements[elementIndex]->GetName() << endl;
				exit(1);
			}
		}

		if( tf1.empty() ) tf1 = Name;

		IConstraint* returnable_const=NULL;
		if(!gotValue)
		{
			if(gotMax&&gotMin){}
			else{
				cerr << "Observable needs a range! " << Name << endl;
				exit(1);
			}

		}
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

//Create a ParameterSet from the appropriate xml tag
ParameterSet * XMLObjectGenerator::GetParameterSet( XMLTag * InputTag )
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
			PhysicsParameter * newParameter = XMLObjectGenerator::GetPhysicsParameter( parameters[parameterIndex], name );
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

//Make FitFunction configuration object
FitFunctionConfiguration * XMLObjectGenerator::MakeFitFunction( XMLTag * FunctionTag )
{
	if ( FunctionTag->GetName() == "FitFunction" )
	{
		string functionName = "Uninitialised";
		string weightName = "Uninitialised";
		string alphaName = "undefined";
		bool hasWeight = false;
		bool hasAlpha = false;
		bool want_Trace = false;
		bool change_style = false;
		string Trace_FileName;
		string Strategy;
		int Threads = -1;
		bool integratorTest = true;
		bool NormaliseWeights = false;
		bool SingleNormaliseWeights = false;
		bool OffSetNLL = false;//true;
		vector<string> ParameterSortList;
		vector< XMLTag* > functionInfo = FunctionTag->GetChildren();
		RapidFitIntegratorConfig* thisConfig = new RapidFitIntegratorConfig();
		if ( functionInfo.size() == 0 )
		{
			//Old style - just specifies the function name
			functionName = XMLTag::GetStringValue( FunctionTag );
		}
		else
		{
			//New style - can have weights
			for ( unsigned int childIndex = 0; childIndex < functionInfo.size(); ++childIndex )
			{
				if ( functionInfo[childIndex]->GetName() == "FunctionName" )
				{
					functionName = XMLTag::GetStringValue( functionInfo[childIndex] );
				}
				else if ( functionInfo[childIndex]->GetName() == "UseGSLNumericalIntegration" )
				{
					thisConfig->useGSLIntegrator = XMLTag::GetBooleanValue( functionInfo[childIndex] );
				}
				else if ( functionInfo[childIndex]->GetName() == "FixedIntegrationPoints" )
				{
					thisConfig->FixedIntegrationPoints = (unsigned)XMLTag::GetIntegerValue( functionInfo[childIndex] );
				}
				else if ( functionInfo[childIndex]->GetName() == "WeightName" )
				{
					hasWeight = true;
					weightName = XMLTag::GetStringValue( functionInfo[childIndex] );
				}
				else if ( functionInfo[childIndex]->GetName() == "AlphaName" )
				{
					hasAlpha = true;
					alphaName = XMLTag::GetStringValue( functionInfo[childIndex] );
				}
				else if ( functionInfo[childIndex]->GetName() == "Trace" )
				{
					want_Trace = true;
					Trace_FileName = XMLTag::GetStringValue( functionInfo[childIndex] );
				}
				else if ( functionInfo[childIndex]->GetName() == "Strategy" )
				{
					change_style = true;
					Strategy = XMLTag::GetStringValue( functionInfo[childIndex] );
				}
				else if ( functionInfo[childIndex]->GetName() == "Threads" )
				{
					Threads = XMLTag::GetIntegerValue( functionInfo[childIndex] );
				}
				else if ( functionInfo[childIndex]->GetName() == "SetIntegratorTest" )
				{
					integratorTest = XMLTag::GetBooleanValue( functionInfo[childIndex] );
				}
				else if ( functionInfo[childIndex]->GetName() == "NormaliseWeights" )
				{
					NormaliseWeights = XMLTag::GetBooleanValue( functionInfo[childIndex] );
				}
				else if ( functionInfo[childIndex]->GetName() == "SingleNormaliseWeights" )
				{
					SingleNormaliseWeights = XMLTag::GetBooleanValue( functionInfo[childIndex] );
				}
				else if ( functionInfo[childIndex]->GetName() == "OffSetNLL" )
				{
					OffSetNLL = XMLTag::GetBooleanValue( functionInfo[childIndex] );
				}
				else if ( functionInfo[childIndex]->GetName() == "RequiredParameterOrder" )
				{
					string thisOrder = XMLTag::GetStringValue( functionInfo[childIndex] );
					ParameterSortList = StringProcessing::SplitString( thisOrder, ':' );
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
			cout <<"Weighted events have been asked for in XML using:\t\t" << weightName << endl << endl;
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

		if( hasAlpha )
		{
			returnable_function->SetAlphaName( alphaName );
		}

		if( !RapidRun::isGridified() ) returnable_function->SetThreads( Threads );
		else returnable_function->SetThreads( 1 );
		returnable_function->SetNormaliseWeights( NormaliseWeights );
		returnable_function->SetSingleNormaliseWeights( SingleNormaliseWeights );
		returnable_function->SetIntegratorTest( integratorTest );
		returnable_function->SetIntegratorConfig( thisConfig );
		returnable_function->SetOffSetNLL( OffSetNLL );
		returnable_function->SetFloatedParameterList( ParameterSortList );

		delete thisConfig;
		return returnable_function;
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << FunctionTag->GetName() << " not \"FitFunction\"" << endl;
		exit(1);
	}
}

//Create a PDF from an appropriate xml tag
IPDF * XMLObjectGenerator::GetNamedPDF( XMLTag * InputTag, PhaseSpaceBoundary* InputBoundary, XMLTag* overloadConfigurator, ParameterSet* thisParameterSet, bool print )
{
	(void) print;
	IPDF* returnable_NamedPDF=NULL;

	vector< XMLTag* > pdfConfig = InputTag->GetChildren();

	if( overloadConfigurator != NULL )
	{
		vector<XMLTag*> optional = overloadConfigurator->GetChildren();
		for( unsigned int i=0; i< optional.size(); ++i )
		{
			pdfConfig.push_back( optional[i] );
		}
	}

	string name = InputTag->GetName();

	if( name != "PDF" && name != "SumPDF" && name != "ProdPDF" && name != "NormalisedSumPDF" )
	{
		cerr << "Unrecognised PDF Tag: " << name << endl << endl;
		exit(-1234567);
	}

	vector<string> observableNames, parameterNames;
	PDFConfigurator* configurator = new PDFConfigurator();

	unsigned int appendParamNum=0;
	unsigned int configParamNum=0;
	unsigned int subParamNum=0;

	//Load the PDF configuration
	for ( unsigned int configIndex = 0; configIndex < pdfConfig.size(); ++configIndex )
	{
		if ( pdfConfig[configIndex]->GetName() == "Name" )
		{
			name = XMLTag::GetStringValue( pdfConfig[configIndex] );
		}
		else if ( pdfConfig[configIndex]->GetName() == "Label" )
		{
			configurator->SetPDFLabel( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
		}
		else if ( pdfConfig[configIndex]->GetName() == "ObservableName" )
		{
			observableNames.push_back( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
		}
		else if ( pdfConfig[configIndex]->GetName() == "ParameterName" )
		{
			parameterNames.push_back( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
		}
		else if ( pdfConfig[configIndex]->GetName() == "ParameterSubstitution" )
		{
			TString ThisNum; ThisNum+=subParamNum;
			pdfConfig[configIndex]->AppendPath( ThisNum.Data() );
			configurator->addParameterSubstitution( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
			++subParamNum;
		}
		else if ( pdfConfig[configIndex]->GetName() == "AppendParameterNames" )
		{
			TString ThisNum; ThisNum+=appendParamNum;
			pdfConfig[configIndex]->AppendPath( ThisNum.Data() );
			configurator->appendParameterNames( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
			++appendParamNum;
		}
		else if ( pdfConfig[configIndex]->GetName() == "ConfigurationParameter" )
		{
			TString ThisNum; ThisNum+=configParamNum;
			pdfConfig[configIndex]->AppendPath( ThisNum.Data() );
			configurator->addConfigurationParameter( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
			++configParamNum;
		}
		else if ( pdfConfig[configIndex]->GetName() == "AppendAllOthers" )
		{
			cerr << "AppendAllOthers has not been implemented yet!" << endl;
			exit(-5623);
			/*
			   TString ThisNum; ThisNum+=appendParamNum;
			   pdfConfig[configIndex]->AppendPath( ThisNum.Data() );
			   configurator->appendOtherParameterNames( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
			   ++appendParamNum;
			 */
		}
		else if ( pdfConfig[configIndex]->GetName() == "FractionName" )
		{
			TString ThisNum; ThisNum+=configParamNum;
			pdfConfig[configIndex]->AppendPath( ThisNum.Data() );
			configurator->AddFractionName( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
			++configParamNum;
		}
		else if( pdfConfig[configIndex]->GetName() == "PDF" || pdfConfig[configIndex]->GetName() == "SumPDF" || pdfConfig[configIndex]->GetName() == "NormalisedSumPDF" || pdfConfig[configIndex]->GetName() == "ProdPDF" )
		{
			IPDF* thisPDF = XMLObjectGenerator::GetPDF( pdfConfig[configIndex], InputBoundary, overloadConfigurator, thisParameterSet, false );
			configurator->AddDaughterPDF( thisPDF );
			delete thisPDF;
		}
		else
		{
			cerr << "(1)Unrecognised PDF configuration: " << pdfConfig[configIndex]->GetName() << endl;
			exit(1);
		}
	}

	cout << endl << endl << endl;

	if( !configurator->empty() )
	{
		cout << "PDFConfigurator: ";
		if( overloadConfigurator != NULL ) cout << "Additional ";
		cout << "Options Passed to this PDF" << endl;
		configurator->Print();
	}


	//      XMLTags from the CommonPDF XML object

	if( overloadConfigurator != NULL )
	{
		pdfConfig = overloadConfigurator->GetChildren();
	}
	else
	{
		pdfConfig = vector<XMLTag*>();
	}

	//Load the PDF configuration
	for ( unsigned int configIndex = 0; configIndex < pdfConfig.size(); ++configIndex )
	{
		if ( pdfConfig[configIndex]->GetName() == "ParameterSubstitution" )
		{
			configurator->addParameterSubstitution( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
		}
		else if ( pdfConfig[configIndex]->GetName() == "AppendParameterNames" )
		{
			configurator->appendParameterNames( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
		}
		else if ( pdfConfig[configIndex]->GetName() == "ConfigurationParameter" )
		{
			configurator->addConfigurationParameter( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
		}
		else if ( pdfConfig[configIndex]->GetName() == "Label" )
		{
			configurator->SetPDFLabel( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
		}
		else
		{
			cerr << "(1b)Unrecognised PDF configuration: " << pdfConfig[configIndex]->GetName() << endl;
			exit(1);
		}
	}

	configurator->SetPhaseSpaceBoundary( InputBoundary );

	//Check if the name is recognised as a PDF
	returnable_NamedPDF = ClassLookUp::LookUpPDFName( name, configurator );

	return returnable_NamedPDF;
}

//Choose one of the PDF instantiation methods
IPDF * XMLObjectGenerator::GetPDF( XMLTag * InputTag, PhaseSpaceBoundary * InputBoundary, XMLTag* overloadConfigurator, ParameterSet* thisParameterSet, bool print )
{
	if( overloadConfigurator != NULL && print )
	{
		vector<XMLTag*> pdfConfig = overloadConfigurator->GetChildren();

		PDFConfigurator* configurator = new PDFConfigurator();

		unsigned int appendParamNum=0;
		unsigned int configParamNum=0;
		unsigned int subParamNum=0;

		//Load the PDF configuration
		for( unsigned int configIndex = 0; configIndex < pdfConfig.size(); ++configIndex )
		{
			if ( pdfConfig[configIndex]->GetName() == "ParameterSubstitution" )
			{
				TString ThisNum; ThisNum+=subParamNum;
				pdfConfig[configIndex]->AppendPath( ThisNum.Data() );
				configurator->addParameterSubstitution( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
				++subParamNum;
			}
			else if ( pdfConfig[configIndex]->GetName() == "AppendParameterNames" )
			{
				TString ThisNum; ThisNum+=appendParamNum;
				pdfConfig[configIndex]->AppendPath( ThisNum.Data() );
				configurator->appendParameterNames( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
				++appendParamNum;
			}
			else if ( pdfConfig[configIndex]->GetName() == "AppendAllOthers" )
			{
				cerr << "AppendAllOthers has not been implemented yet!" << endl;
				exit(-2367);
				/*
				   TString ThisNum; ThisNum+=appendParamNum;
				   pdfConfig[configIndex]->AppendPath( ThisNum.Data() );
				   configurator->appendOtherParameterNames( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
				   ++appendParamNum;
				 */
			}
			else if ( pdfConfig[configIndex]->GetName() == "ConfigurationParameter" )
			{
				TString ThisNum; ThisNum+=configParamNum;
				pdfConfig[configIndex]->AppendPath( ThisNum.Data() );
				configurator->addConfigurationParameter( XMLTag::GetStringValue( pdfConfig[configIndex] ) );
				++configParamNum;
			}
		}

		if( !configurator->empty() )
		{
			cout << "====================================================================================================================" << endl;
			cout << "PDFConfigurator: Config Parameters Passed to CommonPDF" << endl;
			configurator->Print();
			cout << endl;
		}
		delete configurator;
	}

	IPDF* returnable_pdf = XMLObjectGenerator::GetNamedPDF( InputTag, InputBoundary, overloadConfigurator, thisParameterSet );
	cout << "XMLConfigReader:: Constructed " << returnable_pdf->GetLabel() << " PDF" << endl;

	vector<string> haveNames = thisParameterSet->GetAllNames();
	vector<string> wantedNames = returnable_pdf->GetPrototypeParameterSet();
	bool can_be_defined=true;
	for( unsigned int i=0; i< wantedNames.size(); ++i )
	{
		if( StringProcessing::VectorContains( &haveNames, &(wantedNames[i]) ) == -1 )
		{
			can_be_defined = false;
			break;
		}
	}

	if( can_be_defined ) returnable_pdf->UpdatePhysicsParameters( thisParameterSet );

	return returnable_pdf;
}
//Collect the information needed to make a data set
DataSetConfiguration * XMLObjectGenerator::MakeDataSetConfiguration( XMLTag * DataTag, PhaseSpaceBoundary * DataBoundary, XMLTag* common, ParameterSet* thisParameterSet )
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
				dataSource = XMLTag::GetStringValue( dataComponents[dataIndex] );
			}
			else if ( name == "CutString" )
			{
				cutString = XMLTag::GetStringValue( dataComponents[dataIndex] );
			}
			else if ( name == "FileName" || name == "NTuplePath" )
			{
				argumentNames.push_back(name);
				dataArguments.push_back( XMLTag::GetStringValue( dataComponents[dataIndex] ) );
			}
			else if ( name == "NumberEvents" )
			{
				numberEvents = XMLTag::GetIntegerValue( dataComponents[dataIndex] );
			}
			else if ( name == "PDF" || name == "SumPDF" || name == "NormalisedSumPDF" || name == "ProdPDF" )
			{
				PDFXML = dataComponents[dataIndex];
				generatePDFFlag = true;
			}
			else if( name == "CommonPDF" )
			{
				PDFXML = common;
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

		if( generatePDFFlag )   generatePDF = XMLObjectGenerator::GetPDF( PDFXML, DataBoundary, pdfOptions, thisParameterSet );
		generatePDF->SetMCCacheStatus( false );
		delete thisParameterSet;

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

//Make a minimiser configuration object
MinimiserConfiguration * XMLObjectGenerator::MakeMinimiser( XMLTag * MinimiserTag, OutputConfiguration* thisOutput )
{
	if ( MinimiserTag->GetName() == "Minimiser" )
	{
		//Examine all minimiser components
		string minimiserName = "Uninitialised";
		vector< XMLTag* > minimiserComponents = MinimiserTag->GetChildren();
		//vector<string> valueLines = XMLTag::GetStringValue( MinimiserTag );
		int MAXIMUM_MINIMISATION_STEPS = 1000000;
		double FINAL_GRADIENT_TOLERANCE = 0.001;
		MinimiserConfiguration* returnableConfig = NULL;
		vector<string> minimiserOptions;
		int Quality = 1;
		int nSigma = 1;
		bool MultiMini=false;
		if ( minimiserComponents.size() == 0 )
		{
			minimiserName = XMLTag::GetStringValue( MinimiserTag );
		}
		else
		{
			//New style - can have weights
			for ( unsigned int childIndex = 0; childIndex < minimiserComponents.size(); ++childIndex )
			{
				if ( minimiserComponents[childIndex]->GetName() == "MinimiserName" )
				{
					minimiserName = XMLTag::GetStringValue( minimiserComponents[childIndex] );
				}
				else if ( minimiserComponents[childIndex]->GetName() == "NSigma" )
				{
					nSigma = XMLTag::GetIntegerValue( minimiserComponents[childIndex] );
				}
				else if( minimiserComponents[childIndex]->GetName() == "MaxSteps" )
				{
					MAXIMUM_MINIMISATION_STEPS = XMLTag::GetIntegerValue( minimiserComponents[childIndex] );
				}
				else if( minimiserComponents[childIndex]->GetName() == "GradTolerance" )
				{
					FINAL_GRADIENT_TOLERANCE = XMLTag::GetDoubleValue( minimiserComponents[childIndex] );
				}
				else if( minimiserComponents[childIndex]->GetName() == "ConfigureMinimiser" )
				{
					minimiserOptions.push_back( XMLTag::GetStringValue( minimiserComponents[childIndex] ) );
				}
				else if( minimiserComponents[childIndex]->GetName() == "Quality" )
				{
					Quality = XMLTag::GetIntegerValue( minimiserComponents[childIndex] );
				}
				else if( minimiserComponents[childIndex]->GetName() == "MultiMini" )
				{
					MultiMini = XMLTag::GetBooleanValue( minimiserComponents[childIndex] );
				}
				else
				{
					cerr << "Minimiser not properly configured" << endl;
					exit(9234);
				}
			}
		}
		if( DebugClass::DebugThisClass( "XMLConfigReader" ) ) cout << "XMLConfigReader:\tMinimiser: " << minimiserName << " requested in XML, creating" << endl;
		returnableConfig = new MinimiserConfiguration( minimiserName, thisOutput );
		returnableConfig->SetSteps( MAXIMUM_MINIMISATION_STEPS );
		returnableConfig->SetTolerance( FINAL_GRADIENT_TOLERANCE );
		returnableConfig->SetOptions( minimiserOptions );
		returnableConfig->SetQuality( Quality );
		returnableConfig->SetMultiMini( MultiMini );
		returnableConfig->SetNSigma( nSigma );
		if( DebugClass::DebugThisClass( "XMLConfigReader" ) )
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

//Create an ExternalConstraint for the appropriate xml tag
ExternalConstraint * XMLObjectGenerator::GetExternalConstraint( XMLTag * InputTag )
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
				name = XMLTag::GetStringValue( externalComponents[componentIndex] );
			}
			else if ( externalComponents[componentIndex]->GetName() == "Value" )
			{
				value = XMLTag::GetDoubleValue( externalComponents[componentIndex] );
			}
			else if ( externalComponents[componentIndex]->GetName() == "Error" )
			{
				error = XMLTag::GetDoubleValue( externalComponents[componentIndex] );
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

//Create an ExternalConstMatrix for the appropriate xml tag
ExternalConstMatrix * XMLObjectGenerator::GetExternalConstMatrix( XMLTag * InputTag )
{
	if ( InputTag->GetName() == "ExternalConstMatrix" )
	{
		string names, values, errors, correlations;
		bool hasNames=false, hasValues=false, hasErrors=false, hasCorrelations=false;
		vector< XMLTag* > externalComponents = InputTag->GetChildren();
		for ( unsigned int componentIndex = 0; componentIndex < externalComponents.size(); ++componentIndex )
		{
			if ( externalComponents[componentIndex]->GetName() == "Names" )
			{
				names = XMLTag::GetStringValue( externalComponents[componentIndex] );
				hasNames=true;
			}
			else if ( externalComponents[componentIndex]->GetName() == "Values" )
			{
				values = XMLTag::GetStringValue( externalComponents[componentIndex] );
				hasValues=true;
			}
			else if ( externalComponents[componentIndex]->GetName() == "Errors" )
			{
				errors = XMLTag::GetStringValue( externalComponents[componentIndex] );
				hasErrors=true;
			}
			else if ( externalComponents[componentIndex]->GetName() == "Correlations" )
			{
				correlations = XMLTag::GetStringValue( externalComponents[componentIndex] );
				hasCorrelations=true;
			}
			else
			{
				cerr << "Unrecognised constraint component: " << externalComponents[componentIndex]->GetName() << endl;
				exit(1);
			}
		}
		if( !hasNames || !hasValues || !hasErrors || !hasCorrelations )
		{
			cerr << "Incorrect xml tags provided to assemble \"ExternalConstMatrix\"" << endl;
			exit(-87623);
		}
		return new ExternalConstMatrix( names, values, errors, correlations );
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << InputTag->GetName() << "\" not \"ExternalConstMatrix\"" << endl;
		exit(1);
	}
}

//Create a ConstraintFunction for the appropriate xml tag
ConstraintFunction * XMLObjectGenerator::GetConstraintFunction( XMLTag * InputTag )
{
	if ( InputTag->GetName() == "ConstraintFunction" )
	{
		vector< IConstraintFunction* > constraints;
		vector< XMLTag* > functionComponents = InputTag->GetChildren();
		for ( unsigned int componentIndex = 0; componentIndex < functionComponents.size(); ++componentIndex )
		{
			if ( functionComponents[componentIndex]->GetName() == "ExternalConstraint" )
			{
				constraints.push_back( (IConstraintFunction*)GetExternalConstraint( functionComponents[componentIndex] ) );
			}
			else if( functionComponents[componentIndex]->GetName() == "ExternalConstMatrix" )
			{
				constraints.push_back( (IConstraintFunction*)GetExternalConstMatrix( functionComponents[componentIndex] ) );
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

//Make a physics parameter from an appropriate XML tag
ScanParam * XMLObjectGenerator::GetScanParam( XMLTag * InputTag )
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
				name_tag.push_back( XMLTag::GetStringValue( elements[elementIndex] ) );
			}
			else if ( name == "Minimum" )
			{
				minimum.push_back( XMLTag::GetDoubleValue( elements[elementIndex] ) );
			}
			else if ( name == "Maximum" )
			{
				maximum.push_back( XMLTag::GetDoubleValue( elements[elementIndex] ) );
			}
			else if ( name == "Points" )
			{
				points.push_back( XMLTag::GetIntegerValue( elements[elementIndex] ) );
			}
			else if ( name == "Sigma" )
			{
				sigma.push_back( XMLTag::GetIntegerValue( elements[elementIndex] ) );
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

	} else  { cerr << "Unreconised Scan Config: " << InputTag->GetName() << endl; exit(1); }

}

pair<ScanParam*, ScanParam*> XMLObjectGenerator::Get2DScanParam( XMLTag * InputTag )
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
				Param_1 = XMLObjectGenerator::GetScanParam( elements[elementIndex] );
				param1 = true;
			} else if ( name == "Y_Param" )  {
				Param_2 = XMLObjectGenerator::GetScanParam( elements[elementIndex] );
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

CompPlotter_config* XMLObjectGenerator::getCompPlotterConfigs( XMLTag* CompTag )
{
	vector< XMLTag* > projComps = CompTag->GetChildren();

	CompPlotter_config* returnable_config = new CompPlotter_config();

	returnable_config->observableName = "Undefined";
	returnable_config->data_bins = 100;
	returnable_config->PDF_points = 128;
	returnable_config->logY = false;

	returnable_config->xmin = -99999;
	returnable_config->xmax = -99999;
	returnable_config->ymin = -99999;
	returnable_config->ymax = -99999;

	returnable_config->CalcChi2 = false;
	returnable_config->Chi2Value = -99999;

	returnable_config->OnlyZero = false;

	RapidFitIntegratorConfig* projectionIntegratorConfig = returnable_config->integratorConfig;

#ifdef __RAPIDFIT_USE_GSL
	//      When GSL is available we prefer it for Projections
	//      This Integrator is highly user configurable and is multi-thread safe
	projectionIntegratorConfig->useGSLIntegrator = true;
#endif

	//      If we have no extra config just take the name and use defaults
	if( projComps.size() == 0 )
	{
		returnable_config->observableName = XMLTag::GetStringValue( CompTag );
	}

	if( CompTag->GetName() == "Projection" ) returnable_config->OnlyZero = true;
	else returnable_config->OnlyZero = false;

	for( unsigned int childIndex = 0; childIndex < projComps.size(); ++childIndex )
	{
		if( projComps[childIndex]->GetName() == "DataBins" )
		{
			returnable_config->data_bins = XMLTag::GetIntegerValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "PDFpoints" )
		{
			returnable_config->PDF_points = XMLTag::GetIntegerValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "LogY" )
		{
			returnable_config->logY = XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "LogX" )
		{
			returnable_config->logX = XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "Name" )
		{
			returnable_config->observableName = XMLTag::GetStringValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "WidthKey" )
		{
			vector<string> widths = StringProcessing::SplitString(
					XMLTag::GetStringValue( projComps[childIndex] ), ':' );

			if( widths.empty() )
			{
				returnable_config->width_key.push_back( XMLTag::GetIntegerValue( projComps[childIndex] ) );
			}

			for( vector<string>::iterator width_i = widths.begin(); width_i != widths.end(); ++width_i )
				returnable_config->width_key.push_back( atoi( width_i->c_str() ) );

		}
		else if( projComps[childIndex]->GetName() == "ColorKey" )
		{
			vector<string> colors = StringProcessing::SplitString(
					XMLTag::GetStringValue( projComps[childIndex] ), ':' );

			if( colors.empty() )
			{
				returnable_config->style_key.push_back( XMLTag::GetIntegerValue( projComps[childIndex] ) );
			}

			for( vector<string>::iterator color_i = colors.begin(); color_i != colors.end(); ++color_i )
				returnable_config->color_key.push_back( atoi( color_i->c_str() ) );

		}
		else if( projComps[childIndex]->GetName() == "StyleKey" )
		{
			vector<string> styles = StringProcessing::SplitString(
					XMLTag::GetStringValue( projComps[childIndex] ), ':' );
			if( styles.empty() )
			{
				returnable_config->style_key.push_back( XMLTag::GetIntegerValue( projComps[childIndex] ) );
			}
			for( vector<string>::iterator style_i = styles.begin(); style_i != styles.end(); ++style_i )
				returnable_config->style_key.push_back( atoi( style_i->c_str() ) );

		}
		else if( projComps[childIndex]->GetName() == "Title" )
		{
			returnable_config->PlotTitle = XMLTag::GetStringValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "CompNames" || projComps[childIndex]->GetName() == "ComponentNames" )
		{
			returnable_config->component_names = StringProcessing::SplitString(
					XMLTag::GetStringValue( projComps[childIndex] ), ':' );

			if( returnable_config->component_names.empty() )
				returnable_config->component_names.push_back( XMLTag::GetStringValue( projComps[childIndex] ) );
		}
		else if( projComps[childIndex]->GetName() == "CombinationNames" )
		{
			returnable_config->combination_names = StringProcessing::SplitString(
					XMLTag::GetStringValue( projComps[childIndex] ), ':' );

			if( returnable_config->combination_names.empty() )
				returnable_config->combination_names.push_back( XMLTag::GetStringValue( projComps[childIndex] ) );
		}
		else if( projComps[childIndex]->GetName() == "Xmax" )
		{
			returnable_config->xmax = XMLTag::GetDoubleValue( projComps[childIndex] );;
		}
		else if( projComps[childIndex]->GetName() == "Xmin" )
		{
			returnable_config->xmin = XMLTag::GetDoubleValue( projComps[childIndex] );;
		}
		else if( projComps[childIndex]->GetName() == "Ymax" )
		{
			returnable_config->ymax = XMLTag::GetDoubleValue( projComps[childIndex] );;
		}
		else if( projComps[childIndex]->GetName() == "Ymin" )
		{
			returnable_config->ymin = XMLTag::GetDoubleValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "XTitle" )
		{
			returnable_config->xtitle = TString( XMLTag::GetStringValue( projComps[childIndex] ) );
		}
		else if( projComps[childIndex]->GetName() == "YTitle" )
		{
			returnable_config->ytitle = TString( XMLTag::GetStringValue( projComps[childIndex] ) );
		}
		else if( projComps[childIndex]->GetName() == "TrustNumerical" )
		{
			returnable_config->ScaleNumerical = XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "CalcChi2" )
		{
			returnable_config->CalcChi2 = XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "DrawPull" || projComps[childIndex]->GetName() == "DrawPulls" )
		{
			returnable_config->DrawPull = XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "LimitPulls" )
		{
			returnable_config->LimitPulls = XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "AddLHCb" )
		{
			returnable_config->addLHCb = XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "AddRightLHCb" )
		{
			returnable_config->addRightLHCb = XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "LegendTextSize" )
		{
			returnable_config->LegendTextSize =  XMLTag::GetDoubleValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "TopRightLegend" )
		{
			returnable_config->TopRightLegend =  XMLTag::GetBooleanValue( projComps[childIndex] );
			returnable_config->TopLeftLegend=false;
			returnable_config->BottomRightLegend=false;
			returnable_config->BottomLeftLegend=false;
		}
		else if( projComps[childIndex]->GetName() == "TopLeftLegend" )
		{
			returnable_config->TopRightLegend=false;
			returnable_config->TopLeftLegend =  XMLTag::GetBooleanValue( projComps[childIndex] );
			returnable_config->BottomRightLegend=false;
			returnable_config->BottomLeftLegend=false;
		}
		else if( projComps[childIndex]->GetName() == "BottomRightLegend" )
		{
			returnable_config->TopRightLegend=false;
			returnable_config->TopLeftLegend=false;
			returnable_config->BottomRightLegend =  XMLTag::GetBooleanValue( projComps[childIndex] );
			returnable_config->BottomLeftLegend=false;
		}
		else if( projComps[childIndex]->GetName() == "BottomLeftLegend" )
		{
			returnable_config->TopRightLegend=false;
			returnable_config->TopLeftLegend=false;
			returnable_config->BottomRightLegend=false;
			returnable_config->BottomLeftLegend =  XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "NoLegend" )
		{
			returnable_config->useLegend =  !XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "UseSpline" )
		{
			returnable_config->useSpline =  XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "Threads" )
		{
			if( !RapidRun::isGridified() ) projectionIntegratorConfig->numThreads =  (unsigned)XMLTag::GetIntegerValue( projComps[childIndex] );
			else projectionIntegratorConfig->numThreads = 1;
		}
		else if( projComps[childIndex]->GetName() == "FixedIntegrationPoints" )
		{
			projectionIntegratorConfig->FixedIntegrationPoints = (unsigned)XMLTag::GetIntegerValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "UseGSLNumericalIntegration" )
		{
			projectionIntegratorConfig->useGSLIntegrator =  XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "PlotAllCombinatons" )
		{
			returnable_config->plotAllCombinations = XMLTag::GetBooleanValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "DefaultDiscreteValue" )
		{
			returnable_config->defaultCombinationValue = XMLTag::GetDoubleValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "XaxisTitleScale" )
		{
			returnable_config->XaxisTitleScale = XMLTag::GetDoubleValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "XaxisLabelScale" )
		{
			returnable_config->XaxisLabelScale = XMLTag::GetDoubleValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "YaxisTitleScale" )
		{
			returnable_config->YaxisTitleScale = XMLTag::GetDoubleValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "YaxisLabelScale" )
		{
			returnable_config->YaxisLabelScale = XMLTag::GetDoubleValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "OnlyUseCombination" )
		{
			returnable_config->ForceCombinationNumber = XMLTag::GetIntegerValue( projComps[childIndex] );
		}
		else if( projComps[childIndex]->GetName() == "OnlyUseComponent" )
		{
			returnable_config->ForceComponentNumber = XMLTag::GetIntegerValue( projComps[childIndex] );
		}
		else
		{
			cerr << "XMLConfigurationReader: Sorry Don't understand XMLTag: " << projComps[childIndex]->GetName() << " ignoring!" << endl;
		}
	}

	return returnable_config;
}

//Return the pair of observables to plot the function contours for
pair< string, string > XMLObjectGenerator::MakeContourPlot( XMLTag * PlotTag )
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
					xName = XMLTag::GetStringValue( plotComponents[childIndex] );
					hasX = true;
				}
				else if ( plotComponents[childIndex]->GetName() == "YParameter" )
				{
					yName = XMLTag::GetStringValue( plotComponents[childIndex] );
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

//Make an output configuration object
OutputConfiguration * XMLObjectGenerator::MakeOutputConfiguration( XMLTag * OutputTag )
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
				compPlotter_vec.push_back( XMLObjectGenerator::getCompPlotterConfigs( outputComponents[childIndex] ) );
			}
			else if ( outputComponents[childIndex]->GetName() == "ComponentProjection" )
			{
				compPlotter_vec.push_back( XMLObjectGenerator::getCompPlotterConfigs( outputComponents[childIndex] ) );
			}
			else if ( outputComponents[childIndex]->GetName() == "DoPullPlots" )
			{
				pullType = XMLTag::GetStringValue( outputComponents[childIndex] );
			}
			else if ( outputComponents[childIndex]->GetName() == "Scan" )
			{
				ScanParam* temp_SParam = XMLObjectGenerator::GetScanParam( outputComponents[childIndex] );
				ScanParameters.push_back( temp_SParam );
			}
			else if ( outputComponents[childIndex]->GetName() == "TwoDScan" )
			{
				pair<ScanParam*, ScanParam*> temp_2DScan = XMLObjectGenerator::Get2DScanParam( outputComponents[childIndex] );
				_2DScanParameters.push_back( temp_2DScan );
			}
			else if ( outputComponents[childIndex]->GetName() == "2DScan" )
			{
				cerr << "PLEASE MOVE YOUR XML TO THE NEW SYNTAX \'<TwoDScan>\' TO BE STANDARDS COMPLIANT!" << endl;
				pair<ScanParam*, ScanParam*> temp_2DScan = XMLObjectGenerator::Get2DScanParam( outputComponents[childIndex] );
				_2DScanParameters.push_back( temp_2DScan );
			}
			else
			{
				cerr << "Unrecognised output component: " << outputComponents[childIndex]->GetName() << endl;
				exit(1);
			}
		}

		if( DebugClass::DebugThisClass( "XMLConfigReader" ) )
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

//Collect the information needed to make a data set
PDFWithData * XMLObjectGenerator::GetPDFWithData( XMLTag * DataTag, XMLTag * FitPDFTag, int Starting_Value, XMLTag* overloadConfigurator, XMLTag* common, ParameterSet* thisParameterSet, PhaseSpaceBoundary* thisPhaseSpaceBoundary )
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
		DataSetConfiguration* dataSetMaker = NULL;
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
				dataSource = XMLTag::GetStringValue( dataComponents[dataIndex] );
			}
			else if ( name == "Subset" )
			{
				dataSetMaker = XMLObjectGenerator::MakeDataSetConfiguration( dataComponents[dataIndex], dataBoundary, common, thisParameterSet );
			}
			else if ( name == "CutString" )
			{
				cutString = XMLTag::GetStringValue( dataComponents[dataIndex] );
			}
			else if ( name == "FileName" || name == "NTuplePath" )
			{
				argumentNames.push_back(name);
				dataArguments.push_back( XMLTag::GetStringValue( dataComponents[dataIndex] ) );
			}
			else if ( name == "NumberEvents" )
			{
				numberEvents = XMLTag::GetIntegerValue( dataComponents[dataIndex] );
			}
			else if ( name == "PhaseSpaceBoundary" )
			{
				boundaryFound = true;
				if( dataBoundary == NULL ) dataBoundary = XMLObjectGenerator::GetPhaseSpaceBoundary( dataComponents[dataIndex] );
				else
				{
					cerr << "CANNOT DO MULTIPLE PHASESPACES YET" << endl;
					exit(8972);
				}
			}
			else if ( name == "CommonPhaseSpace" )
			{
				boundaryFound = true;
				if( dataBoundary == NULL )
				{
					PhaseSpaceBoundary* tempBound = XMLObjectGenerator::GetPhaseSpaceBoundary( dataComponents[dataIndex] );
					dataBoundary = new PhaseSpaceBoundary(*thisPhaseSpaceBoundary);
					vector<string> additionalNames = tempBound->GetAllNames();
					for( unsigned int i=0; i< additionalNames.size(); ++i )
					{
						dataBoundary->AddConstraint( additionalNames[i], tempBound->GetConstraint( additionalNames[i] ), true );
					}
				}
				else
				{
					cerr << "CANNOT DO MULTIPLE PHASESPACES YET" << endl;
					exit(8972);
				}
			}
			else if ( name == "PDFConfigurator" )
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
					Starting_Value = XMLTag::GetIntegerValue( dataComponents[dataIndex] );
				}
			}
			else if ( name == "CommonPDF" )
			{
				string value = XMLTag::GetStringValue( dataComponents[dataIndex] );
				if( value == "True" )
				{
					generatePDFXML = common;
					generatePDFFlag = true;
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
			generatePDF = XMLObjectGenerator::GetPDF( generatePDFXML, dataBoundary, pdfOptionXML, thisParameterSet );
			generatePDF->SetMCCacheStatus( false );
		}

		//Return the collection of configuration information - data generation will happen later
		if(!boundaryFound)
		{
			cerr << "DataSet defined without PhaseSpaceBoundary" << endl;
			exit(1);
		}

		//If there are no separate data sources, go for backwards compatibility
		if( dataSetMaker == NULL )
		{
			DataSetConfiguration * oldStyleConfig = NULL;

			if( generatePDFFlag )
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

			dataSetMaker = oldStyleConfig;
		}
		//Make the objects
		IPDF * fitPDF = XMLObjectGenerator::GetPDF( FitPDFTag, dataBoundary, overloadConfigurator, thisParameterSet );
		fitPDF->SetMCCacheStatus( false );

		if( generatePDF != NULL ) delete generatePDF;
		PDFWithData* returnable = new PDFWithData( fitPDF, dataBoundary, dataSetMaker );
		delete dataSetMaker;
		if( fitPDF != NULL ) delete fitPDF;
		if( dataBoundary != NULL ) delete dataBoundary;
		return returnable;
	}
	else
	{
		cerr << "Incorrect xml tag provided: \"" << DataTag->GetName() << "\" not \"DataSet\"" << endl;
		exit(1);
	}
}

