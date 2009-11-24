/**
        @class ToyStudyResult

        The result of a toy study.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "ToyStudyResult.h"
#include <iostream>

//Default constructor
ToyStudyResult::ToyStudyResult()
{
}

//Constructor with correct argument
ToyStudyResult::ToyStudyResult( vector<string> AllParameterNames ) : allNames(AllParameterNames)
{
	//Construct the result data structure
	for ( int nameIndex = 0; nameIndex < allNames.size(); nameIndex++ )
	{
		allValues.push_back( vector<double>() );
		allErrors.push_back( vector<double>() );
		allPulls.push_back( vector<double>() );
	}

	clock = new TStopwatch();
}

//Destructor
ToyStudyResult::~ToyStudyResult()
{
}

//Note the time the study starts
void ToyStudyResult::StartStopwatch()
{
	clock->Start();
}

//Add a new fit result
bool ToyStudyResult::AddFitResult( FitResult * NewResult )
{
	vector<double> newParameterValues, newParameterErrors, newParameterPulls;
	vector<string>::iterator nameIterator;
	ResultParameterSet * newSet = NewResult->GetResultParameterSet();

	//Check all expected parameters are found
	for ( nameIterator = allNames.begin(); nameIterator != allNames.end(); nameIterator++ )
	{
		ResultParameter * newResult = newSet->GetResultParameter( *nameIterator );
		if ( newResult->GetUnit() == "NameNotFoundError" )
		{
			//If any parameter is not found, fail
			cerr << "Expected fitted parameter \"" << *nameIterator << "\" not found" << endl;
			return false;
		}
		else
		{
			//Retrieve the parameter information
			newParameterValues.push_back( newResult->GetValue() );
			newParameterErrors.push_back( newResult->GetError() );
			newParameterPulls.push_back( newResult->GetPull() );
		}
	}

	//If you've got this far, all the parameters have been found, so add them to the record
	allResults.push_back(NewResult);
	for ( int nameIndex = 0; nameIndex < allNames.size(); nameIndex++ )
	{
		allValues[nameIndex].push_back( newParameterValues[nameIndex] );
		allErrors[nameIndex].push_back( newParameterErrors[nameIndex] );
		allPulls[nameIndex].push_back( newParameterPulls[nameIndex] );
	}

	//Store the duration
	clock->Stop();
	allRealTimes.push_back( clock->RealTime() );
	allCPUTimes.push_back( clock->CpuTime() );

	return true;
}

//Return vectors of values, errors and pulls for a particular parameter name
vector<double> ToyStudyResult::GetParameterValues( string ParameterName )
{
	for ( int nameIndex = 0; nameIndex < allNames.size(); nameIndex++ )
	{
		//If you find the parameter name, return the vector of values
		if ( ParameterName == allNames[nameIndex] )
		{
			return allValues[nameIndex];
		}
	}

	//If you get this far, the parameter name was not found
	cerr << "Result parameter name \"" << ParameterName << "\" not found" << endl;
	return vector<double>();
}
vector<double> ToyStudyResult::GetParameterErrors( string ParameterName )
{
	for ( int nameIndex = 0; nameIndex < allNames.size(); nameIndex++ )
	{
		//If you find the parameter name, return the vector of errors
		if ( ParameterName == allNames[nameIndex] )
		{
			return allErrors[nameIndex];
		}
	}

	//If you get this far, the parameter name was not found
	cerr << "Result parameter name \"" << ParameterName << "\" not found" << endl;
	return vector<double>();
}
vector<double> ToyStudyResult::GetParameterPulls( string ParameterName )
{
	for ( int nameIndex = 0; nameIndex < allNames.size(); nameIndex++ )
	{
		//If you find the parameter name, return the vector of pulls
		if ( ParameterName == allNames[nameIndex] )
		{
			return allPulls[nameIndex];
		}
	}

	//If you get this far, the parameter name was not found
	cerr << "Result parameter name \"" << ParameterName << "\" not found" << endl;
	return vector<double>();
}

//Allow access to the vector of results
int ToyStudyResult::NumberResults()
{
	return allResults.size();
}
FitResult * ToyStudyResult::GetFitResult( int Index )
{
	if ( Index < allResults.size() )
	{
		return allResults[Index];
	}
	else
	{
		cerr << "Index (" << Index << ") out of range" << endl;
		return new FitResult();
	}
}

//Return names of all variables
vector<string> ToyStudyResult::GetAllNames()
{
	return allNames;
}

//Return the time taken for each fit
vector<double> ToyStudyResult::GetAllRealTimes()
{
	return allRealTimes;
}
vector<double> ToyStudyResult::GetAllCPUTimes()
{
	return allCPUTimes;
}
