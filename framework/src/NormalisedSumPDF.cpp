/*!
 * @class NormalisedSumPDF
 *
 * An implementation of IPDF for adding the values of two other IPDFs, normalised relative to each other
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

//	RapidFit Headers
#include "ClassLookUp.h"
#include "NormalisedSumPDF.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <sstream>
#include <float.h>

PDF_CREATOR( NormalisedSumPDF );

using namespace::std;

NormalisedSumPDF::NormalisedSumPDF( const NormalisedSumPDF& input ) : BasePDF( (BasePDF) input ),
	prototypeDataPoint( input.prototypeDataPoint ), prototypeParameterSet( input.prototypeParameterSet ), doNotIntegrateList( input.doNotIntegrateList ),
	firstPDF( ClassLookUp::CopyPDF( input.firstPDF ) ), secondPDF( ClassLookUp::CopyPDF( input.secondPDF ) ),
	firstFraction( input.firstFraction ), firstIntegralCorrection( input.firstIntegralCorrection ), secondIntegralCorrection( input.secondIntegralCorrection ),
	fractionName( input.fractionName ), integrationBoundary(NULL)
{
	firstPDF->SetDebugMutex( this->DebugMutex(), false );
	secondPDF->SetDebugMutex( this->DebugMutex(), false );

	if( input.integrationBoundary != NULL ) integrationBoundary = new PhaseSpaceBoundary(*input.integrationBoundary);
}

//Constructor specifying fraction parameter name
//NormalisedSumPDF::NormalisedSumPDF( IPDF * FirstPDF, IPDF * SecondPDF, PhaseSpaceBoundary * InputBoundary, string FractionName ) : BasePDF(),
//	prototypeDataPoint(), prototypeParameterSet(), doNotIntegrateList(), firstPDF( ClassLookUp::CopyPDF( FirstPDF ) ), secondPDF( ClassLookUp::CopyPDF( SecondPDF ) ),
//	firstFraction(0.5), firstIntegralCorrection(), secondIntegralCorrection(), fractionName(FractionName), integrationBoundary(new PhaseSpaceBoundary(*InputBoundary) )

NormalisedSumPDF::NormalisedSumPDF( PDFConfigurator* config ) : BasePDF(), prototypeDataPoint(), prototypeParameterSet(), doNotIntegrateList(), firstPDF(NULL), secondPDF(NULL),
	firstFraction(0.5), firstIntegralCorrection(), secondIntegralCorrection(), fractionName(), integrationBoundary(NULL)
{

	if( config->GetDaughterPDFs().size() != 2 )
	{
		cerr << "NormalisedSumPDF requires ONLY 2 daughter PDFs" << endl;
		exit(-54263);
	}
	else
	{
		firstPDF = ClassLookUp::CopyPDF( config->GetDaughterPDFs()[0] );
		secondPDF = ClassLookUp::CopyPDF( config->GetDaughterPDFs()[1] );
		fractionName = config->getName( config->GetFractionName() );
		if( config->GetPhaseSpaceBoundary() != NULL ) integrationBoundary = new PhaseSpaceBoundary( *(config->GetPhaseSpaceBoundary()) );
	}

	this->SetName("NormalisedSumPDF");
	this->SetLabel( "NormalisedSumPDF_("+firstPDF->GetLabel()+")+("+secondPDF->GetLabel()+")" );

	cout << endl;
	cout << "Constructing NormalisedSum "<< this->GetLabel() << endl;
	cout << "FractionName:\t" << string(fractionName) << endl;
	cout << endl;

	firstPDF->SetDebugMutex( this->DebugMutex(), false );
	secondPDF->SetDebugMutex( this->DebugMutex(), false );

	MakePrototypes(integrationBoundary);

	//This by design will create a ParameterSet with the same structure as the prototypeParameterSet list
	allParameters.AddPhysicsParameters( firstPDF->GetPhysicsParameters(), false );
	allParameters.AddPhysicsParameters( secondPDF->GetPhysicsParameters(), false );
	PhysicsParameter* frac_param = new PhysicsParameter( fractionName );
	allParameters.AddPhysicsParameter( frac_param, false );

	//this->SetComponentStatus( firstPDF->GetComponentStatus() || secondPDF->GetComponentStatus() );
	//secondPDF->SetComponentStatus( firstPDF->GetComponentStatus() || secondPDF->GetComponentStatus() );
	//firstPDF->SetComponentStatus( firstPDF->GetComponentStatus() || secondPDF->GetComponentStatus() );
}

void NormalisedSumPDF::SetComponentStatus( const bool input )
{
	firstPDF->SetComponentStatus( input );
	secondPDF->SetComponentStatus( input );
	this->ReallySetComponentStatus( input );
}

bool NormalisedSumPDF::GetComponentStatus() const
{
	return this->ReallyGetComponentStatus();
}

vector<string> NormalisedSumPDF::PDFComponents()
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

	if( firstFraction > 0. && firstFraction < 1. )
	{
		for( unsigned int i=0; i< firstpdf_components.size(); ++i )
		{
			component_list.push_back( StringProcessing::AddNumberToLeft( firstpdf_components[i], 1 ) );
		}
		for( unsigned int j=0; j< secondpdf_components.size(); ++j )
		{
			component_list.push_back( StringProcessing::AddNumberToLeft( secondpdf_components[j], 2 ) );
		}
	}
	else if( firstFraction >= 1. )
	{
		component_list = vector<string>();
		for( unsigned int i=0; i< firstpdf_components.size(); ++i )
		{
			component_list.push_back( StringProcessing::AddNumberToLeft( firstpdf_components[i], 1 ) );
		}
	}
	else if( firstFraction <= 0. )
	{
		component_list = vector<string>();
		for( unsigned int j=0; j< secondpdf_components.size(); ++j )
		{
			component_list.push_back( StringProcessing::AddNumberToLeft( secondpdf_components[j], 2 ) );
		}
	}

	component_list = StringProcessing::MoveElementToStart( component_list, "0" );

	return component_list;
}

void NormalisedSumPDF::SetUpIntegrator( const RapidFitIntegratorConfig* input )
{
    firstPDF->SetUpIntegrator( input );
	secondPDF->SetUpIntegrator( input );
        this->GetPDFIntegrator()->SetUpIntegrator( input );
}

void NormalisedSumPDF::TurnCachingOff()
{
	firstPDF->TurnCachingOff();
	secondPDF->TurnCachingOff();
}

void NormalisedSumPDF::ChangePhaseSpace( PhaseSpaceBoundary * InputBoundary )
{
	if( integrationBoundary != NULL ) delete integrationBoundary;
	integrationBoundary = new PhaseSpaceBoundary( *InputBoundary );
	this->MakePrototypes( InputBoundary );
	firstPDF->ChangePhaseSpace( InputBoundary );
	secondPDF->ChangePhaseSpace( InputBoundary );
}

//Assemble the vectors of parameter/observable names needed
void NormalisedSumPDF::MakePrototypes( PhaseSpaceBoundary * InputBoundary )
{
	//Make the prototype parameter set
	prototypeParameterSet = StringProcessing::CombineUniques( firstPDF->GetPrototypeParameterSet(), secondPDF->GetPrototypeParameterSet() );

	//Add in the fraction Parameter as required
	string fractionStr = fractionName;
	if( StringProcessing::VectorContains( &prototypeParameterSet, &fractionStr ) == -1 ) prototypeParameterSet.push_back( fractionName );

	//Make the prototype data point
	vector<string> firstObservables = firstPDF->GetPrototypeDataPoint();
	vector<string> secondObservables = secondPDF->GetPrototypeDataPoint();
	prototypeDataPoint = StringProcessing::CombineUniques( firstObservables, secondObservables );

	//Make the do not integrate list
	doNotIntegrateList = StringProcessing::CombineUniques( firstPDF->GetDoNotIntegrateList(), secondPDF->GetDoNotIntegrateList() );

	//Make the corrections to the integrals for observables unused by only one PDF
	vector<string>::iterator observableIterator;
	IConstraint * inputConstraint=NULL;
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
NormalisedSumPDF::~NormalisedSumPDF()
{
	//cout << "Hello from Normalised destructor" << endl;
	if( firstPDF != NULL ) delete firstPDF;
	if( secondPDF != NULL ) delete secondPDF;
	if( integrationBoundary != NULL ) delete integrationBoundary;
}

//Set the function parameters
bool NormalisedSumPDF::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	PhysicsParameter * newFraction = NewParameterSet->GetPhysicsParameter(fractionName);
	if ( newFraction->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Parameter \"" << (string)fractionName << "\" expected but not found" << endl;
		return false;
	}
	else
	{
		double newFractionValue = newFraction->GetValue();

		//Stupidity check
		if ( newFractionValue > 1.0 || newFractionValue < 0.0 )
		{
			//cerr << "Requested impossible fraction: " << newFractionValue << endl;
			//return false;
		}
		else
		{
			firstFraction = newFractionValue;
			firstPDF->UpdatePhysicsParameters( NewParameterSet );
			secondPDF->UpdatePhysicsParameters( NewParameterSet );
			bool output = allParameters.SetPhysicsParameters( NewParameterSet );
			return output;
		}
	}
	return false;
}

//Return the integral of the function over the given boundary
double NormalisedSumPDF::Normalisation( DataPoint* NewDataPoint, PhaseSpaceBoundary * NewBoundary )
{
	(void)NewDataPoint;
	(void)NewBoundary;
	//The evaluate method already returns a normalised value
	return 1.0;
}

//Return the function value at the given point
double NormalisedSumPDF::Evaluate( DataPoint* NewDataPoint )
{
	if( firstFraction > 1.0 || firstFraction < 0.0 )
	{
		cerr << "Requested impossible fraction: " << firstFraction << endl;
		return DBL_MAX;
	}

	double termOne=0.;
	double termTwo=0.;

	double firstIntegral=0.;
	double secondIntegral=0.;
	if( firstFraction >= 1. )
	{
		firstIntegral = this->GetFirstIntegral( NewDataPoint );
		termOne = ( firstPDF->Evaluate( NewDataPoint ) ) / firstIntegral;
	}
	else if( firstFraction <= 0. )
	{
		secondIntegral = this->GetSecondIntegral( NewDataPoint );
		termTwo = ( secondPDF->Evaluate( NewDataPoint ) ) / secondIntegral;
	}
	else
	{
		//Calculate the integrals of the PDFs
		firstIntegral = this->GetFirstIntegral( NewDataPoint );
		secondIntegral = this->GetSecondIntegral( NewDataPoint );
		//Get the PDFs' values, normalised and weighted by firstFraction
		termOne = ( firstPDF->Evaluate( NewDataPoint ) * firstFraction ) / firstIntegral;
		termTwo = ( secondPDF->Evaluate( NewDataPoint ) * ( 1 - firstFraction ) ) / secondIntegral;
	}

	double sum=termOne + termTwo;

	if( std::isnan(sum) ||  sum <= 0. )
	{
		PDF_THREAD_LOCK

		cout << termOne*firstIntegral << "/" << firstIntegral << "\t+\t" << termTwo*secondIntegral << "/" << secondIntegral << endl;

		cout << firstPDF->GetLabel() << "\t\t\t\t\t" << secondPDF->GetLabel() << endl << endl;

		PDF_THREAD_UNLOCK
		throw(-653102);
	}

	//Return the sum
	return sum;
}

double NormalisedSumPDF::GetFirstIntegral( DataPoint* NewDataPoint )
{
	return firstPDF->Integral( NewDataPoint, integrationBoundary ) * firstIntegralCorrection;
}

double NormalisedSumPDF::GetSecondIntegral( DataPoint* NewDataPoint )
{
	return secondPDF->Integral( NewDataPoint, integrationBoundary ) * secondIntegralCorrection;
}

//Return the function value at the given point
double NormalisedSumPDF::EvaluateForNumericIntegral( DataPoint * NewDataPoint )
{
        double termOne=0.;
        double termTwo=0.;

        double firstIntegral=0.;
        double secondIntegral=0.;
        if( firstFraction >= 1. )
        {
                firstIntegral = this->GetFirstIntegral( NewDataPoint );
                termOne = ( firstPDF->EvaluateForNumericIntegral( NewDataPoint ) ) / firstIntegral;
        }
        else if( firstFraction <= 0. )
        {
                secondIntegral = this->GetSecondIntegral( NewDataPoint );
                termTwo = ( secondPDF->EvaluateForNumericIntegral( NewDataPoint ) ) / secondIntegral;
        }
        else
        {
                //Calculate the integrals of the PDFs
                firstIntegral = this->GetFirstIntegral( NewDataPoint );
                secondIntegral = this->GetSecondIntegral( NewDataPoint );
                //Get the PDFs' values, normalised and weighted by firstFraction
                termOne = ( firstPDF->EvaluateForNumericIntegral( NewDataPoint ) * firstFraction ) / firstIntegral;
                termTwo = ( secondPDF->EvaluateForNumericIntegral( NewDataPoint ) * ( 1 - firstFraction ) ) / secondIntegral;
        }
	//cout << "NSP numerical integrals " << firstIntegral << " , " << secondIntegral << endl;

	if( std::isnan(termOne) || std::isnan(termTwo) )
	{
		PDF_THREAD_LOCK
		cout << termOne*firstIntegral << "/" << firstIntegral << "\t+\t" << termTwo*secondIntegral << "/" << secondIntegral << endl;
		cout << firstPDF->GetLabel() << "\t\t\t\t\t" << secondPDF->GetLabel() << endl << endl;
		PDF_THREAD_UNLOCK
	}

	//Return the sum
	return termOne + termTwo;
}



//Return a prototype data point
vector<string> NormalisedSumPDF::GetPrototypeDataPoint()
{
	return prototypeDataPoint;
}

//Return a prototype set of physics parameters
vector<string> NormalisedSumPDF::GetPrototypeParameterSet()
{
	return prototypeParameterSet;
}

//Return a list of parameters not to be integrated
vector<string> NormalisedSumPDF::GetDoNotIntegrateList()
{
	return doNotIntegrateList;
}

bool NormalisedSumPDF::GetNumericalNormalisation() const
{
	return firstPDF->GetNumericalNormalisation() || secondPDF->GetNumericalNormalisation();
}

//Return the function value at the given point
double NormalisedSumPDF::EvaluateComponent( DataPoint* NewDataPoint, ComponentRef* componentIndexObj )
{
	if( componentIndexObj == NULL ) return this->Evaluate( NewDataPoint );
	string componentIndex = componentIndexObj->getComponentName();
	int component_num = componentIndexObj->getComponentNumber();

	//cout << endl << componentIndexObj->getComponentName() << endl;

	if( component_num == -1 || componentIndexObj->getSubComponent() == NULL )
	{
		component_num = StringProcessing::GetNumberOnLeft( componentIndex );
		componentIndexObj->setComponentNumber( component_num );
		componentIndexObj->addSubComponent( StringProcessing::RemoveFirstNumber( componentIndex ) );
	}
	//cout << "Eval Component: " << component_num << "\tAsking for " << StringProcessing::RemoveFirstNumber( componentIndex ) << endl << endl;

	double returnable_value = -1.;

	double termOne = -1.;
	double termTwo = -1.;

	//      All components are:
	//                              0       :       total of the whole NormalisedSumPDF
	//                              1xx     :       all of the components in the first PDF
	//                              2yy     :       all of the components in the second PDF

	//      !!!!!   NORMALISE THE VALUE     !!!!!
	switch( component_num )
	{
		case 1:
			//      We want a component from the first PDF, integrate over the second PDF
			termOne = firstPDF->EvaluateComponent( NewDataPoint, componentIndexObj->getSubComponent() ) / this->GetFirstIntegral( NewDataPoint );
			termTwo = 0.;
			break;

		case 2:
			//      We want a component from the second PDF, integrate over the first PDF
			termOne = 0.;
			termTwo = secondPDF->EvaluateComponent( NewDataPoint, componentIndexObj->getSubComponent() ) / this->GetSecondIntegral( NewDataPoint );
			break;

		case 0:
			termOne = firstPDF->EvaluateComponent( NewDataPoint, componentIndexObj->getSubComponent() ) / this->GetFirstIntegral( NewDataPoint );
			termTwo = secondPDF->EvaluateComponent( NewDataPoint, componentIndexObj->getSubComponent() ) / this->GetSecondIntegral( NewDataPoint );
			break;

		default:

			returnable_value = 0;//this->EvaluateForNumericIntegral( NewDataPoint );
			return returnable_value;
	}

	termOne = ( termOne * firstFraction );
	termTwo = ( termTwo * (1.-firstFraction) );


	//cout << "NormSum:\t" << setw(3) << componentIndex << " : " << component_num << " : " << StringProcessing::RemoveFirstNumber( componentIndex ) << "\t\t" << termOne+termTwo << " _:_ " << termOne << " , " << termTwo << "\t" << firstIntegral << " , " << secondIntegral << "\t" << this->GetLabel() << endl;

	returnable_value = termOne + termTwo;

	return returnable_value;
}

IPDF* NormalisedSumPDF::GetFirstPDF() const
{
	return firstPDF;
}

IPDF* NormalisedSumPDF::GetSecondPDF() const
{
	return secondPDF;
}

string NormalisedSumPDF::GetFractionName() const
{
	return fractionName;
}

void NormalisedSumPDF::SetCachingEnabled( bool input )
{
	firstPDF->SetCachingEnabled( input );
	secondPDF->SetCachingEnabled( input );
}

bool NormalisedSumPDF::GetCachingEnabled() const
{
	return firstPDF->GetCachingEnabled() && secondPDF->GetCachingEnabled();
}

string NormalisedSumPDF::XML() const
{
	stringstream xml;

	xml << "<NormalisedSumPDF>" << endl;
	xml << "<FractionName>" << (string)fractionName << "</FractionName>" << endl;
	xml << firstPDF->XML();
	xml << secondPDF->XML();
	xml << "</NormalisedSumPDF>" << endl;

	return xml.str();
}

void NormalisedSumPDF::SetDebugMutex( pthread_mutex_t* Input, bool can_remove )
{
	can_remove_mutex = can_remove;
	if( debug_mutex != NULL && can_remove_mutex ) delete debug_mutex;
	firstPDF->SetDebugMutex( Input, false );
	secondPDF->SetDebugMutex( Input, false );
	debug_mutex = Input;
}

void NormalisedSumPDF::SetDebug( DebugClass* input_debug )
{
	BasePDF::SetDebug( input_debug );
	firstPDF->SetDebug( input_debug );
	secondPDF->SetDebug( input_debug );
}

string NormalisedSumPDF::GetComponentName( ComponentRef* componentIndexObj )
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
				//	We wanr this PDF
				return this->GetLabel();

			default:
				return "Unknown-NormalisedSum";
		}
	}
	return "Invalid-NormalisedSum";
}

