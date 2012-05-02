/**
        @class InputParsing

        Collection of static functions for creating RapidFit data objects from string templates

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "StringProcessing.h"
#include "InputParsing.h"
#include "ObservableDiscreteConstraint.h"
#include "ObservableContinuousConstraint.h"
#include "ClassLookUp.h"
//	System Headers
#include <stdlib.h>
#include <iostream>

//Create a parameter set from an input string - assumes all parameters have minimum and maximum defined
ParameterSet * InputParsing::MakeParameterSet( string Input )
{
	vector<string> splitInput = StringProcessing::SplitString( Input, ',' );
	if ( splitInput.size() % 6 == 0 )
	{
		int numberParameters = (int)(splitInput.size() / 7);
		vector<string> names;
		vector<double> values;
		vector<double> minima;
		vector<double> maxima;
		vector<double> steps;
		vector<string> types;
		vector<string> units;

		//Retrieve all parameter components
		for ( int parameterIndex = 0; parameterIndex < numberParameters; ++parameterIndex )
		{
			names.push_back( splitInput[ unsigned(parameterIndex) ] );
			values.push_back( strtod( splitInput[ unsigned(parameterIndex + 1) ].c_str(), NULL ) );
			minima.push_back( strtod( splitInput[ unsigned(parameterIndex + 2) ].c_str(), NULL ) );
			maxima.push_back( strtod( splitInput[ unsigned(parameterIndex + 3) ].c_str(), NULL ) );
			steps.push_back( strtod( splitInput[ unsigned(parameterIndex + 4) ].c_str(), NULL ) );
			types.push_back( splitInput[ unsigned(parameterIndex + 5) ] );
			units.push_back( splitInput[ unsigned(parameterIndex + 6) ] );
		}

		ParameterSet * newParameters = new ParameterSet( names );

		//Populate the parameter set
		for ( int parameterIndex = 0; parameterIndex < numberParameters; ++parameterIndex )
		{
			newParameters->SetPhysicsParameter( names[unsigned(parameterIndex)], values[unsigned(parameterIndex)], minima[unsigned(parameterIndex)], maxima[unsigned(parameterIndex)], steps[unsigned(parameterIndex)], types[unsigned(parameterIndex)], units[unsigned(parameterIndex)] );
		}

		return newParameters;
	}
	else
	{
		cerr << "String describing parameters is incorrectly formed: number of components ( " << splitInput.size() << " ) is not a multiple of 6" << endl;
		return new ParameterSet( vector<string>() );
	}
}

//Create a parameter set from a vector of strings
ParameterSet * InputParsing::MakeParameterSet( vector<string> Input )
{
	vector<string> parameterNames, splitInput;
	vector<PhysicsParameter> physicsParameters;
	string name, type, unit;
	double value, minimum, maximum, steps;

	//Parse the input
	for (unsigned int parameterIndex = 0; parameterIndex < Input.size(); ++parameterIndex )
	{
		splitInput = StringProcessing::SplitString( Input[parameterIndex], ':' );
		if ( splitInput.size() == 7 )
		{
			//Parse the parameter components
			name = splitInput[0];
			value = strtod( splitInput[1].c_str(), NULL );
			minimum = strtod( splitInput[2].c_str(), NULL );
			maximum = strtod( splitInput[3].c_str(), NULL );
			steps = strtod( splitInput[4].c_str(), NULL );
			type = splitInput[5];
			unit = splitInput[6];

			//Make the parameter
			physicsParameters.push_back( PhysicsParameter( name, value, minimum, maximum, steps, type, unit ) );
			parameterNames.push_back(name);
		}
		else if ( splitInput.size() == 5 )
		{
			//Parse components of an unbounded parameter
			name = splitInput[0];
			value = strtod( splitInput[1].c_str(), NULL );
			steps =  strtod( splitInput[2].c_str(), NULL );
			type = splitInput[3];
			unit = splitInput[4];

			//Make the parameter
			physicsParameters.push_back( PhysicsParameter( name, value, steps, type, unit ) );
			parameterNames.push_back(name);
		}
		else
		{
			cerr << "Malformed parameter description: " << Input[parameterIndex] << " has " << Input.size() << " terms, not 4 or 6" << endl;
		}
	}

	//Create the set
	ParameterSet * newSet = new ParameterSet(parameterNames);

	//Populate the set
	for (unsigned int parameterIndex = 0; parameterIndex < physicsParameters.size(); ++parameterIndex )
	{
		newSet->SetPhysicsParameter( parameterNames[parameterIndex], &physicsParameters[parameterIndex] );
	}

	return newSet;
}

//Make a phase space boundary
PhaseSpaceBoundary * InputParsing::MakePhaseSpaceBoundary( string Input )
{
	vector<string> allNames, splitTemplate, splitMaxMin;
	vector<IConstraint*> newConstraints;

	//Get the templates for the individual boundaries
	vector<string> boundaryTemplates = StringProcessing::SplitString( Input, '{' );
	for (unsigned int templateIndex = 0; templateIndex < boundaryTemplates.size(); ++templateIndex )
	{
		StringProcessing::RemoveCharacter( boundaryTemplates[templateIndex], '}' );
		splitTemplate = StringProcessing::SplitString( boundaryTemplates[templateIndex], ':' );

		//Distinguish discrete and continuous boundaries
		if ( splitTemplate.size() == 3 )
		{
			//Either continuous, or discrete with one possible value
			splitMaxMin = StringProcessing::SplitString( splitTemplate[1], '<' );
			if ( splitMaxMin.size() == 1 )
			{
				//Discrete
				vector<double> allValues;
				allValues.push_back( strtod( splitTemplate[1].c_str(), NULL ) );
				ObservableDiscreteConstraint * newDiscrete = new ObservableDiscreteConstraint( splitTemplate[0], allValues,  splitTemplate[2] );
				newConstraints.push_back(newDiscrete);
				allNames.push_back( splitTemplate[0] );
			}
			else if (  splitMaxMin.size() == 2 )
			{
				//Continuous
				double min = strtod( splitMaxMin[0].c_str(), NULL );
				double max = strtod( splitMaxMin[1].c_str(), NULL );
				ObservableContinuousConstraint * newContinuous = new ObservableContinuousConstraint( splitTemplate[0], min, max, splitTemplate[2] );
				newConstraints.push_back(newContinuous);
				allNames.push_back( splitTemplate[0] );
			}
			else
			{
				cerr << "Incorrect input format: " << boundaryTemplates[templateIndex] << endl;
			}
		}
		else
		{
			//Discrete
			vector<double> allValues;
			for (unsigned int valueIndex = 1; valueIndex < splitTemplate.size() - 1; ++valueIndex )
			{
				allValues.push_back( strtod( splitTemplate[valueIndex].c_str(), NULL ) );
			}
			ObservableDiscreteConstraint * newDiscrete = new ObservableDiscreteConstraint( splitTemplate[0], allValues,  splitTemplate[2] );
                        newConstraints.push_back(newDiscrete);
			allNames.push_back( splitTemplate[0] );
		}
	}

	//Make the boundary
	PhaseSpaceBoundary * newBoundary = new PhaseSpaceBoundary(allNames);

	//Populate the boundary
	for (unsigned int constraintIndex = 0; constraintIndex < allNames.size(); ++constraintIndex )
	{
		newBoundary->SetConstraint( allNames[constraintIndex], newConstraints[constraintIndex] );
	}

	return newBoundary;
}

//Make a PDFWithData object
PDFWithData * InputParsing::MakePDFWithData( string PDFName, string DataSource, string PhaseSpacePrototype )
{
	PhaseSpaceBoundary * newBoundary = MakePhaseSpaceBoundary(PhaseSpacePrototype);

	//Assume for the moment we're just dealing with regular PDFs
	vector<string> pdfComponents = StringProcessing::SplitString( PDFName, '{' );
	IPDF * newPDF=NULL;
	if ( pdfComponents.size() == 1 )
	{
		//PELC 
		PDFConfigurator* dummy = new PDFConfigurator;
		cout << " In InputParsing::MakePDFWithData : WARNING : Calling ClassLookUp::LookUpPDFName with dummy configurator " << endl ;
		//~PELC
		//PELC newPDF = ClassLookUp::LookUpPDFName( pdfComponents[0], vector<string>(), vector<string>() );
		newPDF = ClassLookUp::LookUpPDFName( pdfComponents[0], dummy );
	}
	else if ( pdfComponents.size() == 3 )
	{
		string observables = pdfComponents[1];
		StringProcessing::RemoveCharacter( observables, '}' );
		vector<string> observableComponents = StringProcessing::SplitString( observables, ':' );
		string parameters = pdfComponents[2];
		StringProcessing::RemoveCharacter( parameters, '}' );
		vector<string> parameterComponents = StringProcessing::SplitString( parameters, ':' );

		//PELC 
		PDFConfigurator* dummy = new PDFConfigurator;
		cout << " In InputParsing::MakePDFWithData : WARNING : Calling ClassLookUp::LookUpPDFName with dummy configurator " << endl ;
		//~PELC
		
		//PELC newPDF = ClassLookUp::LookUpPDFName( pdfComponents[0], observableComponents, parameterComponents );
		newPDF = ClassLookUp::LookUpPDFName( pdfComponents[0], dummy );
	}
	else
	{
		cerr << "PDF prototype incorrectly formed" << endl;
	}

	vector<string> dataComponents = StringProcessing::SplitString( DataSource, ':' );
	if ( dataComponents.size() > 1 )
	{
		string dataSourceName = dataComponents[0];
		int dataAmount = atoi( dataComponents[1].c_str() );

		//Put the leftovers down as extra arguments
		vector<string> dataArguments, argumentNames;
		for (unsigned int argumentIndex = 2; argumentIndex < dataComponents.size(); ++argumentIndex )
		{
			dataArguments.push_back( dataComponents[argumentIndex] );

			if ( argumentIndex == 2 )
			{
				argumentNames.push_back("FileName");
			}
			else if ( argumentIndex == 3 )
			{
				argumentNames.push_back("NTuplePath");
			}
			else
			{
				argumentNames.push_back("Unknown");
			}
		}

		//Make the configuration objects
		string cutString = ""; // For now, just apply no cut
		vector< DataSetConfiguration* > singleDataSet;
		singleDataSet.push_back( new DataSetConfiguration( dataSourceName, dataAmount, cutString, dataArguments, argumentNames ) );
		return new PDFWithData( newPDF, newBoundary, singleDataSet );
	}
	else
	{
		cerr << "Incorrectly formed data source description: must be at least 2 terms" << endl;
		return NULL;
	}
}
