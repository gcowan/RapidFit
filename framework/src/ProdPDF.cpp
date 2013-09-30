
///	RapidFit Headers
#include "BasePDF.h"
#include "ProdPDF.h"
#include "StringProcessing.h"
#include "ClassLookUp.h"
///	System Headers
#include <iostream>
#include <cstdlib>
#include <sstream>

using namespace::std;

PDF_CREATOR( ProdPDF );

//Constructor not specifying fraction parameter name
//ProdPDF::ProdPDF( IPDF * FirstPDF, IPDF * SecondPDF ) : BasePDF(), prototypeDataPoint(), prototypeParameterSet(), doNotIntegrateList(), firstPDF( ClassLookUp::CopyPDF(FirstPDF) ), secondPDF( ClassLookUp::CopyPDF(SecondPDF) )
ProdPDF::ProdPDF( PDFConfigurator* config ) : BasePDF(), prototypeDataPoint(), prototypeParameterSet(), doNotIntegrateList(), firstPDF( NULL ), secondPDF( NULL )
{
	if( config->GetDaughterPDFs().size() != 2 )
	{
		cerr << "ProdPDF requires ONLY 2 daughter PDFs" << endl;
		exit(-54262);
	}
	else
	{
		firstPDF = ClassLookUp::CopyPDF( config->GetDaughterPDFs()[0] );
		secondPDF = ClassLookUp::CopyPDF( config->GetDaughterPDFs()[1] );
	}

	this->SetName( "ProdPDF" );
	this->SetLabel( "ProdPDF["+firstPDF->GetLabel()+"]x["+secondPDF->GetLabel()+"]" );

	cout << endl;
	cout << "Constructing ProdPDF " << this->GetLabel() << endl;
	cout << endl;

	MakePrototypes();

	firstPDF->SetDebugMutex( this->DebugMutex(), false );
	secondPDF->SetDebugMutex( this->DebugMutex(), false );

	//      Only one PDF has components beyond here
	//      Add an index to allow for the ProdPDF to return the correct component in future

	//      By definition the observable you project over is only in 1 of the pdfs

	this->TurnThisCachingOff();

	//This by design will create a ParameterSet with the same structure as the prototypeParameterSet list
	allParameters.AddPhysicsParameters( firstPDF->GetPhysicsParameters(), false );
	allParameters.AddPhysicsParameters( secondPDF->GetPhysicsParameters(), false );

	//this->SetComponentStatus( firstPDF->GetComponentStatus() || secondPDF->GetComponentStatus() );
	//secondPDF->SetComponentStatus( firstPDF->GetComponentStatus() || secondPDF->GetComponentStatus() );
	//firstPDF->SetComponentStatus( firstPDF->GetComponentStatus() || secondPDF->GetComponentStatus() );
}

void ProdPDF::SetComponentStatus( const bool input )
{
	firstPDF->SetComponentStatus( input );
	secondPDF->SetComponentStatus( input );
	this->ReallySetComponentStatus( input );
}

bool ProdPDF::GetComponentStatus() const
{
	return this->ReallyGetComponentStatus();
}

vector<string> ProdPDF::PDFComponents()
{
	vector<string> first_components;
	for( unsigned int i=0; i< firstPDF->PDFComponents().size(); ++i )
	{
		first_components.push_back( firstPDF->PDFComponents()[i] );
	}
	vector<string> second_components;
	for( unsigned int j=0; j< secondPDF->PDFComponents().size(); ++j )
	{
		second_components.push_back( secondPDF->PDFComponents()[j] );
	}

	string zero="0";

	if( first_components.empty() || (StringProcessing::VectorContains( &first_components, &zero ) == -1 ) ) first_components.push_back("0");
	if( second_components.empty() || (StringProcessing::VectorContains( &second_components, &zero ) == -1 ) ) second_components.push_back("0");

	first_components = StringProcessing::MoveElementToStart( first_components, "0" );
	second_components = StringProcessing::MoveElementToStart( second_components, "0" );

	component_list = vector<string>();

	if( (first_components.size() != 1) && (second_components.size() == 1) )
	{
		//component_list.pop_back();
		for( unsigned int i=0; i< first_components.size(); ++i )
		{
			component_list.push_back( StringProcessing::AddNumberToLeft( first_components[i], 1 ) );
		}
		string total="10";
		if( StringProcessing::VectorContains( &component_list, &total ) == -1 )
		{
			component_list.push_back( "0" );
		}
		else
		{
			component_list[ (unsigned)StringProcessing::VectorContains( &component_list, &total ) ] = "0";
		}
	}
	else if( (second_components.size() != 1) && (first_components.size() == 1) )
	{
		//component_list.pop_back();
		for( unsigned int j=0; j< second_components.size(); ++j )
		{
			component_list.push_back( StringProcessing::AddNumberToLeft( second_components[j], 2 ) );
		}
		string total="20";
		if( StringProcessing::VectorContains( &component_list, &total ) == -1 )
		{
			component_list.push_back( "0" );
		}
		else
		{
			component_list[ (unsigned)StringProcessing::VectorContains( &component_list ,&total ) ] = "0";
		}
	}
	else
	{
		if( (first_components.size() != 1) && (second_components.size() != 1) )
		{
			cerr << "CANNOT UNDERSTAND MULTIPLE PDFS WITH MULTIPLE COMPONENTS INSIDE THE SAME PDF, HENCE I WILL ONLY CLAIM 1 COMPONENT" << endl;
		}
		component_list.push_back( "0" );
	}

	component_list = StringProcessing::MoveElementToStart( component_list, "0" );

	return component_list;
}

ProdPDF::ProdPDF( const ProdPDF& input ) : BasePDF( (BasePDF) input ),
	prototypeDataPoint( input.prototypeDataPoint ),
	prototypeParameterSet( input.prototypeParameterSet ),
	doNotIntegrateList( input.doNotIntegrateList ),
	firstPDF( ClassLookUp::CopyPDF( input.firstPDF ) ), secondPDF( ClassLookUp::CopyPDF( input.secondPDF ) )
{
	firstPDF->SetDebugMutex( this->DebugMutex(), false );
	secondPDF->SetDebugMutex( this->DebugMutex(), false );
}

void ProdPDF::SetUpIntegrator( const RapidFitIntegratorConfig* input )
{
	firstPDF->SetUpIntegrator( input );
	secondPDF->SetUpIntegrator( input );
}

void ProdPDF::TurnThisCachingOff()
{
	this->ReallyTurnCachingOff();
}

void ProdPDF::TurnCachingOff()
{
	this->TurnThisCachingOff();
	firstPDF->TurnCachingOff();
	secondPDF->TurnCachingOff();
}

//Assemble the vectors of parameter/observable names needed
void ProdPDF::MakePrototypes()
{
	prototypeParameterSet = StringProcessing::CombineUniques( firstPDF->GetPrototypeParameterSet(), secondPDF->GetPrototypeParameterSet() );
	prototypeDataPoint = StringProcessing::CombineUniques( firstPDF->GetPrototypeDataPoint(), secondPDF->GetPrototypeDataPoint() );
	doNotIntegrateList = StringProcessing::CombineUniques( firstPDF->GetDoNotIntegrateList(), secondPDF->GetDoNotIntegrateList() );
}

//Destructor
ProdPDF::~ProdPDF()
{
	//cout << "Hello from Product destructor" << endl;
	if( firstPDF != NULL ) delete firstPDF;
	if( secondPDF != NULL ) delete secondPDF;
}

//Set the function parameters
bool ProdPDF::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	firstPDF->UpdatePhysicsParameters( NewParameterSet );
	secondPDF->UpdatePhysicsParameters( NewParameterSet );
	bool output = allParameters.SetPhysicsParameters( NewParameterSet );
	return output;
}

//Return the integral of the function over the given boundary
double ProdPDF::Normalisation( DataPoint* NewDataPoint, PhaseSpaceBoundary * NewBoundary )
{
	//Note that this is almost certainly wrong. However, I don't know a good analytical solution.
	//In cases that the formula is incorrect, it will be caught by the numerical integration check.
	double termOne = firstPDF->Integral( NewDataPoint, NewBoundary );
	double termTwo = secondPDF->Integral( NewDataPoint, NewBoundary );

	double prod = termOne * termTwo;

	if( std::isnan(prod) || prod <= 0. )
	{
		PDF_THREAD_LOCK

			cout << "ProdPDF::Normalisation\t\t" << termOne << "\tx\t" << termTwo << endl;

		cout << firstPDF->GetLabel() << "\t\t\t\t\t" << secondPDF->GetLabel() << endl << endl;

		PDF_THREAD_UNLOCK
			throw(-653102);
	}

	return prod;
}

//Return the function value at the given point
double ProdPDF::Evaluate( DataPoint * NewDataPoint )
{
	double termOne = firstPDF->Evaluate( NewDataPoint );
	double termTwo = secondPDF->Evaluate( NewDataPoint );

	double prod = termOne * termTwo;

	if( std::isnan(prod) || prod <= 0. )
	{
		PDF_THREAD_LOCK

			cout << "ProdPDF::Evaluate\t\t" << termOne << "\tx\t" << termTwo << endl;

		cout << firstPDF->GetLabel() << "\t\t\t\t\t" << secondPDF->GetLabel() << endl << endl;

		PDF_THREAD_UNLOCK
			throw(-653102);
	}

	return prod;
}

//Return the function value at the given point for numerical integration
double ProdPDF::EvaluateForNumericIntegral( DataPoint * NewDataPoint )
{
	double termOne = firstPDF->EvaluateForNumericIntegral( NewDataPoint );
	double termTwo = secondPDF->EvaluateForNumericIntegral( NewDataPoint );

	double prod = termOne * termTwo;

	if( std::isnan(prod) || prod <= 0. )
	{
		PDF_THREAD_LOCK

			cout << "ProdPDF::EvaluateForNumerical\t\t" << termOne << "\tx\t" << termTwo << endl;

		cout << firstPDF->GetLabel() << "\t\t\t\t\t" << secondPDF->GetLabel() << endl << endl;

		PDF_THREAD_UNLOCK
			throw(-653102);
	}

	return prod;
}

//Return a prototype data point
vector<string> ProdPDF::GetPrototypeDataPoint()
{
	return prototypeDataPoint;
}

//Return a prototype set of physics parameters
vector<string> ProdPDF::GetPrototypeParameterSet()
{
	return prototypeParameterSet;
}

//Return a list of parameters not to be integrated
vector<string> ProdPDF::GetDoNotIntegrateList()
{
	return doNotIntegrateList;
}

bool ProdPDF::GetNumericalNormalisation() const
{
	return firstPDF->GetNumericalNormalisation() || secondPDF->GetNumericalNormalisation();
}

//Return the function value at the given point
double ProdPDF::EvaluateComponent( DataPoint* NewDataPoint, ComponentRef* componentIndexObj )
{
	string componentIndex = componentIndexObj->getComponentName();
	int component_num = componentIndexObj->getComponentNumber();
	if( component_num == -1 )
	{
		if( componentIndex=="" || componentIndex=="0" ) component_num = 0;
		else if( componentIndex=="10" || componentIndex=="20" ) component_num = 0;
		else component_num = StringProcessing::GetNumberOnLeft( componentIndex );

		componentIndexObj->setComponentNumber( component_num );
		componentIndexObj->addSubComponent( StringProcessing::RemoveFirstNumber( componentIndex ) );
	}

	//cout << "Called:" << component_num << "\t" << StringProcessing::RemoveFirstNumber( componentIndex ) << endl;

	double returnable_value = -1;

	double termOne = -1;
	double termTwo = -1;

	//      By definition only 1 of the PDFs should have components to plot
	switch( component_num )
	{
		case 0:
			termOne = firstPDF->EvaluateForNumericIntegral( NewDataPoint );
			termTwo = secondPDF->EvaluateForNumericIntegral( NewDataPoint );
			break;

		case 1:
			termOne = firstPDF->EvaluateComponent( NewDataPoint, componentIndexObj->getSubComponent() );
			termTwo = secondPDF->EvaluateForNumericIntegral( NewDataPoint );
			break;

		case 2:
			termOne = firstPDF->EvaluateForNumericIntegral( NewDataPoint );
			termTwo = secondPDF->EvaluateComponent( NewDataPoint, componentIndexObj->getSubComponent() );
			break;
		default:
			return 0;
	}

	//cout << "ProdPDF: " << termOne * termTwo << "\t\t" << termOne << " , " << termTwo << "\t" << this->GetLabel() << endl;
	returnable_value = termOne * termTwo;

	//if( returnable_value < 1E-10 ) NewDataPoint->Print();

	return returnable_value;
}

IPDF* ProdPDF::GetFirstPDF() const
{
	return firstPDF;
}

IPDF* ProdPDF::GetSecondPDF() const
{
	return secondPDF;
}

bool ProdPDF::GetCachingEnabled() const
{
	return firstPDF->GetCachingEnabled() && secondPDF->GetCachingEnabled();
}

void ProdPDF::SetCachingEnabled( bool input )
{
	firstPDF->SetCachingEnabled( input );
	secondPDF->SetCachingEnabled( input );
}

string ProdPDF::XML() const
{
	stringstream xml;

	xml << "<ProdPDF>" << endl;
	xml << firstPDF->XML();
	xml << secondPDF->XML();
	xml << "</ProdPDF>" << endl;

	return xml.str();
}

void ProdPDF::SetDebugMutex( pthread_mutex_t* Input, bool can_remove )
{
	can_remove_mutex = can_remove;
	if( debug_mutex != NULL && can_remove_mutex ) delete debug_mutex;
	firstPDF->SetDebugMutex( Input, false );
	secondPDF->SetDebugMutex( Input, false );
	debug_mutex = Input;
}

void ProdPDF::SetDebug( DebugClass* input_debug )
{
	BasePDF::SetDebug( input_debug );
	firstPDF->SetDebug( input_debug );
	secondPDF->SetDebug( input_debug );
}

string ProdPDF::GetComponentName( ComponentRef* componentIndexObj )
{
	if( componentIndexObj == NULL ) return "Unknown";
	else
	{
		string componentIndex = componentIndexObj->getComponentName();
		int component_num = componentIndexObj->getComponentNumber();

		if( component_num == -1 || componentIndexObj->getSubComponent() == NULL )
		{
			if( componentIndex=="" || componentIndex=="0" ) component_num = 0;
			else if( componentIndex=="10" || componentIndex=="20" ) component_num = 0;
			else component_num = StringProcessing::GetNumberOnLeft( componentIndex );

			componentIndexObj->setComponentNumber( component_num );
			componentIndexObj->addSubComponent( StringProcessing::RemoveFirstNumber( componentIndex ) );

		}

		//      All components are:
		//                              0       :       total of the whole NormalisedSumPDF
		//                              1xx     :       all of the components in the first PDF
		//                              2yy     :       all of the components in the second PDF

		//      !!!!!   NORMALISE THE VALUE     !!!!!
		switch( component_num )
		{
			case 1:
				//      We want a component from the first PDF
				return firstPDF->GetComponentName( componentIndexObj->getSubComponent() );

			case 2:
				//      We want a component from the second PDF
				return secondPDF->GetComponentName( componentIndexObj->getSubComponent() );

			case 0:
				//      We wanr this PDF
				return this->GetLabel();

			default:
				return "Unknown-Prod";
		}

	}
	return "Invalid-Prob";
}

