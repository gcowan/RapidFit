/*!
 * @class SumPDF
 *
 * An implementation of IPDF for adding the values of two other IPDFs
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

///	RapidFit Headers
#include "SumPDF.h"
#include "StringProcessing.h"
#include "ClassLookUp.h"
#include "ComponentRef.h"
///	System Headers
#include <iostream>
#include <sstream>
#include <float.h>

using namespace::std;

PDF_CREATOR( SumPDF );

SumPDF::SumPDF( const SumPDF& input ) : BasePDF( (BasePDF) input ), prototypeDataPoint(input.prototypeDataPoint), prototypeParameterSet(input.prototypeParameterSet), doNotIntegrateList(input.doNotIntegrateList),
	firstPDF(ClassLookUp::CopyPDF(input.firstPDF) ), secondPDF( ClassLookUp::CopyPDF(input.secondPDF) ), firstFraction(input.firstFraction), firstIntegralCorrection(input.firstIntegralCorrection),
	secondIntegralCorrection(input.secondIntegralCorrection), fractionName(input.fractionName)
{
	firstPDF->SetDebugMutex( this->DebugMutex(), false );
	secondPDF->SetDebugMutex( this->DebugMutex(), false );
}

//Constructor specifying fraction parameter name
//SumPDF::SumPDF( IPDF * FirstPDF, IPDF * SecondPDF, PhaseSpaceBoundary * InputBoundary, string FractionName ) : prototypeDataPoint(), prototypeParameterSet(), doNotIntegrateList(), firstPDF( ClassLookUp::CopyPDF(FirstPDF) ), secondPDF( ClassLookUp::CopyPDF(SecondPDF) ), firstFraction(0.5), firstIntegralCorrection(), secondIntegralCorrection(), fractionName(FractionName)

SumPDF::SumPDF( PDFConfigurator* config ) : BasePDF(), prototypeDataPoint(), prototypeParameterSet(), doNotIntegrateList(), firstPDF(NULL), secondPDF(NULL), firstFraction(0.5),
	firstIntegralCorrection(), secondIntegralCorrection(), fractionName(), integrationBoundary(NULL)
{
	if( config->GetDaughterPDFs().size() != 2 )
	{
		cerr << "SumPDF requires ONLY 2 daughter PDFs" << endl;
		exit(-54261);
	}
	else
	{
		firstPDF = ClassLookUp::CopyPDF( config->GetDaughterPDFs()[0] );
		secondPDF = ClassLookUp::CopyPDF( config->GetDaughterPDFs()[1] );
		fractionName = config->getName( config->GetFractionName() );
		integrationBoundary = new PhaseSpaceBoundary( *(config->GetPhaseSpaceBoundary()) );
	}

	cout << endl;
	cout << "Constructing SumPDF" << endl;
	cout << endl;
	this->SetName( "SumPDF" );
	this->SetLabel( "SumPDF_("+firstPDF->GetLabel()+")+("+firstPDF->GetLabel()+")" );
	MakePrototypes(integrationBoundary);

	firstPDF->SetDebugMutex( this->DebugMutex(), false );
	secondPDF->SetDebugMutex( this->DebugMutex(), false );

	//this->SetComponentStatus( firstPDF->GetComponentStatus() || secondPDF->GetComponentStatus() );
	//secondPDF->SetComponentStatus( firstPDF->GetComponentStatus() || secondPDF->GetComponentStatus() );
	//firstPDF->SetComponentStatus( firstPDF->GetComponentStatus() || secondPDF->GetComponentStatus() );
}

void SumPDF::SetComponentStatus( const bool input )
{
	firstPDF->SetComponentStatus( input );
	secondPDF->SetComponentStatus( input );
	this->ReallySetComponentStatus( input );
}

bool SumPDF::GetComponentStatus() const
{
	return this->ReallyGetComponentStatus();
}

vector<string> SumPDF::PDFComponents()
{
	vector<string> firstpdf_components;
	for( unsigned int i=0; i< firstPDF->PDFComponents().size(); ++i )
	{
		firstpdf_components.push_back( string( firstPDF->PDFComponents()[i] ) );
	}
	vector<string> secondpdf_components;
	for( unsigned int j=0; j< secondPDF->PDFComponents().size(); ++j )
	{
		secondpdf_components.push_back( string( secondPDF->PDFComponents()[j] ) );
	}

	string zero="0";

	if( firstPDF->GetName() != "NormalisedSum" || firstPDF->GetName() != "Sum" )
	{
		if( StringProcessing::VectorContains( &firstpdf_components, &zero ) == -1 ) firstpdf_components.push_back( "0" );
	}
	if( secondPDF->GetName() != "NormalisedSum" || secondPDF->GetName() != "Sum" )
	{
		if( StringProcessing::VectorContains( &secondpdf_components, &zero ) == -1 ) secondpdf_components.push_back( "0" );
	}

	firstpdf_components = StringProcessing::MoveElementToStart( firstpdf_components, "0" );
	secondpdf_components = StringProcessing::MoveElementToStart( secondpdf_components, "0" );

	//      All components are:
	//                              0       :       total of the whole NormalisedSumPDF
	//                              1xx     :       all of the components in the first PDF
	//                              2yy     :       all of the components in the second PDF

	vector<string> empty;
	component_list.swap(empty);

	for( unsigned int i=0; i< firstpdf_components.size(); ++i )
	{
		component_list.push_back( StringProcessing::AddNumberToLeft( firstpdf_components[i], 1 ) );
	}
	for( unsigned int j=0; j< secondpdf_components.size(); ++j )
	{
		component_list.push_back( StringProcessing::AddNumberToLeft( secondpdf_components[j], 2 ) );
	}

	component_list = StringProcessing::MoveElementToStart( component_list, "0" );

	return component_list;
}

void SumPDF::SetUpIntegrator( const RapidFitIntegratorConfig* input )
{
	firstPDF->SetUpIntegrator( input );
	secondPDF->SetUpIntegrator( input );
}

void SumPDF::TurnThisCachingOff()
{
	this->ReallyTurnCachingOff();
}

void SumPDF::TurnCachingOff()
{
	this->TurnThisCachingOff();
	firstPDF->TurnCachingOff();
	secondPDF->TurnCachingOff();
}

//Assemble the vectors of parameter/observable names needed
void SumPDF::MakePrototypes( PhaseSpaceBoundary * InputBoundary )
{
	//Make sure the ratio of the two PDFs is included
	vector<string> secondParameterSet = secondPDF->GetPrototypeParameterSet();
	secondParameterSet.push_back(fractionName);

	//Make the prototype parameter set
	prototypeParameterSet = StringProcessing::CombineUniques( firstPDF->GetPrototypeParameterSet(), secondParameterSet );

	//Make the prototype data point
	vector<string> firstObservables = firstPDF->GetPrototypeDataPoint();
	vector<string> secondObservables = secondPDF->GetPrototypeDataPoint();
	prototypeDataPoint = StringProcessing::CombineUniques( firstObservables, secondObservables );

	//Make the do not integrate list
	doNotIntegrateList = StringProcessing::CombineUniques( firstPDF->GetDoNotIntegrateList(), secondPDF->GetDoNotIntegrateList() );

	//Make the correctionss to the integrals for observables unused by only one PDF
	vector<string>::iterator observableIterator;
	IConstraint * inputConstraint;
	firstIntegralCorrection = 1.0;
	secondIntegralCorrection = 1.0;
	bool doIntegrate=false;
	for ( observableIterator = firstObservables.begin(); observableIterator != firstObservables.end(); ++observableIterator )
	{
		if ( StringProcessing::VectorContains( &secondObservables, &(*observableIterator) ) == -1 )
		{
			//The first PDF uses this observable, the second doesn't
			inputConstraint = InputBoundary->GetConstraint( *observableIterator );
			doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrateList, &(*observableIterator) ) == -1 );

			//Update this integral correction
			if ( !inputConstraint->IsDiscrete() && doIntegrate )
			{
				secondIntegralCorrection *= ( inputConstraint->GetMaximum() - inputConstraint->GetMinimum() );
			}
		}
	}
	for ( observableIterator = secondObservables.begin(); observableIterator != secondObservables.end(); ++observableIterator )
	{
		if ( StringProcessing::VectorContains( &firstObservables, &(*observableIterator) ) == -1 )
		{
			//The second PDF uses this observable, the first doesn't
			inputConstraint = InputBoundary->GetConstraint( *observableIterator );
			doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrateList, &(*observableIterator) ) == -1 );

			//Update this integral correction
			if ( !inputConstraint->IsDiscrete() && doIntegrate )
			{
				firstIntegralCorrection *= ( inputConstraint->GetMaximum() - inputConstraint->GetMinimum() );
			}
		}
	}
}

//Destructor
SumPDF::~SumPDF()
{
	if( firstPDF != NULL ) delete firstPDF;
	if( secondPDF != NULL ) delete secondPDF;
}

bool SumPDF::GetNumericalNormalisation() const
{
	return firstPDF->GetNumericalNormalisation() || secondPDF->GetNumericalNormalisation();
}

//Set the function parameters
bool SumPDF::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	PhysicsParameter * newFraction = NewParameterSet->GetPhysicsParameter(fractionName);
	if ( newFraction->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Parameter \"" << fractionName << "\" expected but not found" << endl;
		return false;
	}
	else
	{
		double newFractionValue = newFraction->GetValue();

		//Stupidity check
		if ( newFractionValue > 1.0 || newFractionValue < 0.0 )
		{
			cerr << "Requested impossible fraction: " << newFractionValue << endl;
			return false;
		}
		else
		{
			firstFraction = newFractionValue;
			firstPDF->UpdatePhysicsParameters( NewParameterSet );
			secondPDF->UpdatePhysicsParameters( NewParameterSet );
			return allParameters.SetPhysicsParameters( NewParameterSet );
		}
	}
}

//Return the integral of the function over the given boundary
double SumPDF::Normalisation( DataPoint* NewDataPoint, PhaseSpaceBoundary * NewBoundary )
{
	//Get the PDFs' integrals, weighted by firstFraction and corrected for unused observables
	double termOne = firstPDF->Integral( NewDataPoint, NewBoundary ) * firstFraction * firstIntegralCorrection;
	double termTwo = secondPDF->Integral( NewDataPoint, NewBoundary ) * ( 1 - firstFraction ) * secondIntegralCorrection;
	return termOne + termTwo;
}

//Return the function value at the given point
double SumPDF::Evaluate( DataPoint * NewDataPoint )
{
	if( firstFraction > 1.0 || firstFraction < 0.0 )                  
	{
		cerr << "Requested impossible fraction: " << firstFraction << endl;
		return DBL_MAX;
	}
	//Get the PDFs' values, weighted by firstFraction
	double termOne = firstPDF->Evaluate( NewDataPoint ) * firstFraction;
	double termTwo = secondPDF->Evaluate( NewDataPoint ) * ( 1 - firstFraction );
	return termOne + termTwo;
}


//Return a prototype data point
vector<string> SumPDF::GetPrototypeDataPoint()
{
	return prototypeDataPoint;
}

//Return a prototype set of physics parameters
vector<string> SumPDF::GetPrototypeParameterSet()
{
	return prototypeParameterSet;
}

//Return a list of parameters not to be integrated
vector<string> SumPDF::GetDoNotIntegrateList()
{
	return doNotIntegrateList;
}

//Return the components of the function value at the given point
double SumPDF::EvaluateComponent( DataPoint* NewDataPoint, ComponentRef* componentIndexObj )
{
	string componentIndex = componentIndexObj->getComponentName();
	int component_num = componentIndexObj->getComponentNumber();
	if( component_num == -1 )
	{
		component_num = StringProcessing::GetNumberOnLeft( componentIndex );
		componentIndexObj->setComponentNumber( component_num );
		componentIndexObj->addSubComponent( StringProcessing::RemoveFirstNumber( componentIndex ) );
	}

	double returnable_value = -1;

	double termOne = -1;
	double termTwo = -1;

	//      All components are:
	//                              0       :       total of the whole NormalisedSumPDF
	//                              1xx     :       all of the components in the first PDF
	//                              2yy     :       all of the components in the second PDF

	switch( component_num )
	{
		case 1:
			//      We want a component from the first PDF, integrate over the second PDF
			termOne = firstPDF->EvaluateComponent( NewDataPoint, componentIndexObj->getSubComponent() );
			termTwo = 0;
			break;

		case 2:
			//      We want a component from the second PDF, integrate over the first PDF
			termOne = 0;
			termTwo = secondPDF->EvaluateComponent( NewDataPoint, componentIndexObj->getSubComponent() );
			break;

		case 0:
			termOne = firstPDF->EvaluateComponent( NewDataPoint, componentIndexObj->getSubComponent() );
			termTwo = secondPDF->EvaluateComponent( NewDataPoint, componentIndexObj->getSubComponent() );
			break;

		default:
			//cerr << "\n\n\t\tSumPDF can't understand component number: " << component_num << endl<<endl;
			return 0.;
	}

	returnable_value = termOne * firstFraction + termTwo * (1 - firstFraction);

	return returnable_value;
}

IPDF* SumPDF::GetFirstPDF() const
{
	return firstPDF;
}

IPDF* SumPDF::GetSecondPDF() const
{
	return secondPDF;
}

void SumPDF::SetCachingEnabled( bool input )
{
	firstPDF->SetCachingEnabled( input );
	secondPDF->SetCachingEnabled( input );
}

bool SumPDF::GetCachingEnabled() const
{
	return firstPDF->GetCachingEnabled() && secondPDF->GetCachingEnabled();
}

string SumPDF::XML() const
{
	stringstream xml;

	xml << "<SumPDF>" << endl;
	xml << firstPDF->XML() << endl;
	xml << secondPDF->XML() << endl;
	xml << "</SumPDF>" << endl;

	return xml.str();
}

void SumPDF::SetDebugMutex( pthread_mutex_t* Input, bool can_remove )
{
	can_remove_mutex = can_remove;
	if( debug_mutex != NULL && can_remove_mutex ) delete debug_mutex;
	firstPDF->SetDebugMutex( Input, false );
	secondPDF->SetDebugMutex( Input, false );
	debug_mutex = Input;
}

void SumPDF::SetDebug( DebugClass* input_debug )
{
	BasePDF::SetDebug( input_debug );
	firstPDF->SetDebug( input_debug );
	secondPDF->SetDebug( input_debug );
}

string SumPDF::GetComponentName( ComponentRef* componentIndexObj )
{
	if( componentIndexObj == NULL ) return "Unknown";
	else
	{
		string componentIndex = componentIndexObj->getComponentName();
		int component_num = componentIndexObj->getComponentNumber();

		if( component_num == -1 || componentIndexObj->getSubComponent() == NULL )
		{
			component_num = StringProcessing::GetNumberOnLeft( componentIndex );
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
				return "Unknown-Sum";
		}

	}
	return "Invalid-Sum";
}

