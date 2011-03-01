/**
  @class ToyStudyResult

  The result of a toy study.

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

#include "ToyStudyResult.h"
#include "LLscanResult.h"
#include <iostream>

//Default constructor
ToyStudyResult::ToyStudyResult()
{
}

//  Constructor to Return a single array from multiple arrays
ToyStudyResult::ToyStudyResult( vector<ToyStudyResult*> Result_Array )
{
	if( !Result_Array.empty() )
	{
		clock = new TStopwatch();
		allNames = Result_Array[0]->GetAllNames();
		//Construct the result data structure
		for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); nameIndex++ )
		{
			allValues.push_back( vector<double>() );
			allErrors.push_back( vector<double>() );
			allPulls.push_back( vector<double>() );
			allGenValues.push_back( vector<double>() );
		}
		for(unsigned int i=0; i < Result_Array.size(); i++ )
		{
			for( short int j=0; j < Result_Array[i]->NumberResults(); j++ )
			{
				AddFitResult(  Result_Array[i]->GetFitResult( j ), false );
			}

			vector<double> input_real_times = Result_Array[i]->GetAllRealTimes();
			vector<double> input_cpu_times = Result_Array[i]->GetAllCPUTimes();
			for(unsigned int j2=0; j2 < input_real_times.size(); j2++)
			{
				allRealTimes.push_back( input_real_times[j2] );
				allCPUTimes.push_back( input_cpu_times[j2] );
			}

		}
	}
}

//Constructor with correct argument
ToyStudyResult::ToyStudyResult( vector<string> AllParameterNames ) : allNames(AllParameterNames)
{
	//Construct the result data structure
	for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); nameIndex++ )
	{
		allValues.push_back( vector<double>() );
		allErrors.push_back( vector<double>() );
		allPulls.push_back( vector<double>() );
		allGenValues.push_back( vector<double>() );
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
bool ToyStudyResult::AddFitResult( FitResult * NewResult, bool with_clock )
{
	vector<double> newParameterValues, newParameterErrors, newParameterPulls, newParameterGenValues;
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
			newParameterGenValues.push_back( newResult->GetOriginalValue() );
		}
	}

	//If you've got this far, all the parameters have been found, so add them to the record
	allResults.push_back(NewResult);
	for (unsigned int nameIndex = 0; nameIndex < allNames.size(); nameIndex++ )
	{
		allValues[nameIndex].push_back( newParameterValues[nameIndex] );
		allErrors[nameIndex].push_back( newParameterErrors[nameIndex] );
		allPulls[nameIndex].push_back( newParameterPulls[nameIndex] );
		allGenValues[nameIndex].push_back( newParameterGenValues[nameIndex] );
	}

	if( with_clock )
	{
		//Store the duration
		clock->Stop();
		allRealTimes.push_back( clock->RealTime() );
		allCPUTimes.push_back( clock->CpuTime() );
	}

	return true;
}


//  Return an array of the MLL values of the fits
vector<double> ToyStudyResult::GetAllMLL()
{
	vector<double> output_MLL;
	for (unsigned short int i=0; i < allResults.size(); i++ )
	{
		if( allResults[i]->GetFitStatus() >0 ) 
		{
			output_MLL.push_back( allResults[i]->GetMinimumValue() );
		}
		else{
			output_MLL.push_back( LLSCAN_FIT_FAILURE_VALUE );
		}
	}
	return output_MLL;
}

//Return vectors of values, errors and pulls for a particular parameter name
vector<double> ToyStudyResult::GetParameterValues( string ParameterName )
{
	for (unsigned int nameIndex = 0; nameIndex < allNames.size(); nameIndex++ )
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
	for (unsigned int nameIndex = 0; nameIndex < allNames.size(); nameIndex++ )
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
	for (unsigned int nameIndex = 0; nameIndex < allNames.size(); nameIndex++ )
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
	if ( Index < int(allResults.size()) )
	{
		return allResults[Index];
	}
	else
	{
		cerr << "Index (" << Index << ") out of range" << endl;
		return new FitResult();
	}
}
double ToyStudyResult::GetRealTime( int Index )
{
	if( Index < int(allRealTimes.size()) ) return allRealTimes[ Index ];
	else return -1;
}
void ToyStudyResult::AddRealTimes( vector<double> input_times )
{
	for( unsigned int i=0; i< input_times.size(); i++)
	{
		allRealTimes.push_back( input_times[i] );
	}
}
void ToyStudyResult::SetRealTime( int Index, double input_time )
{
	if( int(allRealTimes.size()) < Index ) allRealTimes.resize( Index );
	allRealTimes[Index] = input_time;
}
double ToyStudyResult::GetCPUTime( int Index )
{
	if( Index < int(allCPUTimes.size()) ) return allCPUTimes[ Index ];
	else return -1;
}
void ToyStudyResult::AddCPUTimes( vector<double> input_times )
{
	for( unsigned int i=0; i< input_times.size(); i++)
	{
		allCPUTimes.push_back( input_times[i] );
	}
}
void ToyStudyResult::SetCPUTime( int Index, double input_time )
{
	if( Index < int(allCPUTimes.size()) ) allCPUTimes.resize( Index );
	allCPUTimes[Index] = input_time;  
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


vector<double> ToyStudyResult::GetFlatResult( int Index )
{
	vector<double> Flatresult;
	for(unsigned int i = 0; i<allNames.size(); i++)
	{
		Flatresult.push_back( (allValues[i][Index]) );
		Flatresult.push_back( (allErrors[i][Index]) );
		Flatresult.push_back( (allPulls[i][Index]) );
		Flatresult.push_back( (allGenValues[i][Index]) );

	}

	Flatresult.push_back(allRealTimes[Index]);
	Flatresult.push_back(allCPUTimes[Index]);
	Flatresult.push_back(allResults[Index]->GetFitStatus());
	Flatresult.push_back(allResults[Index]->GetMinimumValue());
	return Flatresult;
}
TString ToyStudyResult::GetFlatResultHeader()
{
	TString header = "";
	for(unsigned short int i = 0; i<allNames.size(); i++)
	{
		TString name = allNames[i];
		header += name + "_value:";
		header += name + "_error:";
		header += name + "_pull:";
		header += name + "_gen:";
	}
	header += "Fit_RealTime:Fit_CPUTime:Fit_Status:NLL";
	return header;
}
