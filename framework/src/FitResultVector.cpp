/**
  @class FitResultVector

  The result of a toy study.

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

//	RapidFit Headers
#include "FitResultVector.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <string>
#include <vector>

using namespace::std;

//  Constructor to Return a single array from multiple arrays
FitResultVector::FitResultVector( const vector<FitResultVector*> Result_Array ) :
	allResults(), allNames(), allValues(), allErrors(), allPulls(), allGenValues(),
	allRealTimes(), allCPUTimes(), clock(NULL)
	#ifdef RAPIDFIT_USETGLTIMER
	, gl_clock(NULL), allGLTimes()
	#endif
{
	if( !Result_Array.empty() )
	{
		clock = new TStopwatch();
		#ifdef RAPIDFIT_USETGLTIMER
		gl_clock = new TGLStopwatch();
		#endif
		allNames = Result_Array.empty()? vector<string>() : Result_Array[0]->GetAllNames();
		//Construct the result data structure
		for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex )
		{
			allValues.push_back( vector<double>() );
			allErrors.push_back( vector<double>() );
			allPulls.push_back( vector<double>() );
			allGenValues.push_back( vector<double>() );
		}
		for(unsigned int i=0; i < Result_Array.size(); ++i )
		{
			for( int j=0; j < Result_Array[i]->NumberResults(); ++j )
			{
				AddFitResult( Result_Array[i]->GetFitResult( j ), false );
			}

			vector<double> input_real_times = Result_Array[i]->GetAllRealTimes();
			vector<double> input_cpu_times = Result_Array[i]->GetAllCPUTimes();
			vector<double> input_gl_times = Result_Array[i]->GetAllGLTimes();
			for( unsigned int j2=0; j2 < input_real_times.size(); ++j2 )
			{
				allRealTimes.push_back( input_real_times[j2] );
				allCPUTimes.push_back( input_cpu_times[j2] );
				#ifdef RAPIDFIT_USETGLTIMER
				allGLTimes.push_back( input_gl_times[j2] );
				#endif
			}

			if( ((int)input_real_times.size() != Result_Array[i]->NumberResults()) || ((int)input_cpu_times.size() != Result_Array[i]->NumberResults()) )
			{
				cout << "ERROR: Time Data from Fits not consistent!!! Ignoring!" << endl;
			}
		}
	}
}

//Constructor with correct argument
FitResultVector::FitResultVector( const vector<string> AllParameterNames ) :
	allResults(), allNames(AllParameterNames), allValues(), allErrors(), allPulls(), allGenValues(),
	allRealTimes(), allCPUTimes(), clock(NULL)
	#ifdef RAPIDFIT_USETGLTIMER
	, gl_clock(NULL), allGLTimes()
	#endif
{
	vector<string> duplicates;
	allNames = StringProcessing::RemoveDuplicates( AllParameterNames, duplicates );
	if( allNames.size() != AllParameterNames.size() )
	{
		cerr << "WARNING: Cannot Generate a FitResultVector with 2 Occurances of the same name" << endl;
		for( vector<string>::iterator str_i = duplicates.begin(); str_i != duplicates.end(); ++str_i )
		{
			cout << *str_i << endl;
		}
	}


	//Construct the result data structure
	for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex )
	{
		allValues.push_back( vector<double>() );
		allErrors.push_back( vector<double>() );
		allPulls.push_back( vector<double>() );
		allGenValues.push_back( vector<double>() );
	}

	clock = new TStopwatch();
	#ifdef RAPIDFIT_USETGLTIMER
	gl_clock = new TGLStopwatch();
	#endif
}

//Destructor
FitResultVector::~FitResultVector()
{
	if( clock != NULL ) delete clock;
	#ifdef RAPIDFIT_USETGLTIMER
	if( gl_clock != NULL ) delete gl_clock;
	#endif
}

//Note the time the study starts
void FitResultVector::StartStopwatch()
{
	#ifdef RAPIDFIT_USETGLTIMER
	gl_clock->Start();
	#endif
	clock->Start();
}

//Add a new fit result
bool FitResultVector::AddFitResult( FitResult * NewResult, const bool with_clock )
{
	vector<double> newParameterValues, newParameterErrors, newParameterPulls, newParameterGenValues;
	vector<string>::iterator nameIterator;
	ResultParameterSet * newSet = NewResult->GetResultParameterSet();

	//Check all expected parameters are found
	for ( nameIterator = allNames.begin(); nameIterator != allNames.end(); ++nameIterator )
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
	for (unsigned int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex )
	{
		allValues[nameIndex].push_back( newParameterValues[nameIndex] );
		allErrors[nameIndex].push_back( newParameterErrors[nameIndex] );
		allPulls[nameIndex].push_back( newParameterPulls[nameIndex] );
		allGenValues[nameIndex].push_back( newParameterGenValues[nameIndex] );
	}

	if( with_clock && clock != NULL )
	{
		//Store the duration
		#ifdef RAPIDFIT_USETGLTIMER
		double thisTime = gl_clock->End();
		allGLTimes.push_back( thisTime );
		#endif
		clock->Stop();
		allRealTimes.push_back( clock->RealTime() );
		allCPUTimes.push_back( clock->CpuTime() );
	}

	return true;
}


//  Return an array of the MLL values of the fits
vector<double> FitResultVector::GetAllMLL() const
{
	vector<double> output_MLL;
	for (unsigned short int i=0; i < allResults.size(); ++i )
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
vector<double> FitResultVector::GetParameterValues( const string ParameterName ) const
{
	for (unsigned int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex )
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
vector<double> FitResultVector::GetParameterErrors( const string ParameterName ) const
{
	for (unsigned int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex )
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
vector<double> FitResultVector::GetParameterPulls( const string ParameterName ) const
{
	for (unsigned int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex )
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
int FitResultVector::NumberResults() const
{
	return (int)allResults.size();
}

unsigned int FitResultVector::size() const
{
	return (unsigned)allResults.size();
}

FitResult * FitResultVector::GetFitResult( const int Index ) const
{
	if ( Index < int(allResults.size()) )
	{
		return allResults[unsigned(Index)];
	}
	else
	{
		cerr << "Index (" << Index << ") out of range" << endl;
		return NULL;
	}
}

double FitResultVector::GetRealTime( const int Index ) const
{
	if( Index < int(allRealTimes.size()) ) return allRealTimes[ unsigned(Index) ];
	else return -1;
}

void FitResultVector::AddRealTimes( const vector<double> input_times )
{
	for( unsigned int i=0; i< input_times.size(); ++i)
	{
		allRealTimes.push_back( input_times[i] );
	}
}

void FitResultVector::AddRealTime( const double input_time )
{
	allRealTimes.push_back( input_time );
}

void FitResultVector::SetRealTime( const int Index, const double input_time )
{
	if( int(allRealTimes.size()) < Index ) allRealTimes.resize( unsigned(Index) );
	allRealTimes[unsigned(Index)] = input_time;
}

double FitResultVector::GetCPUTime( const int Index ) const
{
	if( Index < int(allCPUTimes.size()) ) return allCPUTimes[ unsigned(Index) ];
	else return -1;
}

void FitResultVector::AddCPUTimes( const vector<double> input_times )
{
	for( unsigned int i=0; i< input_times.size(); ++i)
	{
		allCPUTimes.push_back( input_times[i] );
	}
}

void FitResultVector::AddCPUTime( const double input_time )
{
	allCPUTimes.push_back( input_time );
}

void FitResultVector::AddGLTime( const double input_time )
{
	#ifdef RAPIDFIT_USETGLTIMER
        allGLTimes.push_back( input_time );
	#else
	(void) input_time;
	#endif
}

double FitResultVector::GetGLTime( const int Index ) const
{
	#ifdef RAPIDFIT_USETGLTIMER
        if( Index < int(allGLTimes.size()) ) return allGLTimes[ unsigned(Index) ];
        else return -1;
	#else
	(void) Index;
	return -1;
	#endif
}

void FitResultVector::SetCPUTime( const int Index, const double input_time )
{
	if( Index < int(allCPUTimes.size()) ) allCPUTimes.resize( unsigned(Index) );
	allCPUTimes[unsigned(Index)] = input_time;
}

//Return names of all variables
vector<string> FitResultVector::GetAllNames() const
{
	return allNames;
}

//Return the time taken for each fit
vector<double> FitResultVector::GetAllRealTimes() const
{
	return allRealTimes;
}

vector<double> FitResultVector::GetAllCPUTimes() const
{
	return allCPUTimes;
}

vector<double> FitResultVector::GetAllGLTimes() const
{
	#ifdef RAPIDFIT_USETGLTIMER
	return allGLTimes;
	#else
	return vector<double>( allCPUTimes.size(), -1. );
	#endif
}

vector<double> FitResultVector::GetFlatResult( const int Index ) const
{
	vector<double> Flatresult;
	for(unsigned int i = 0; i<allNames.size(); ++i)
	{
		Flatresult.push_back( (allValues[i][unsigned(Index)]) );
		if( (allResults[(unsigned)Index]->GetResultParameterSet()->GetResultParameter(allNames[i])->GetType() != "Fixed") || allResults[(unsigned)Index]->GetResultParameterSet()->GetResultParameter(allNames[i])->GetScanStatus() )
		{
			Flatresult.push_back( (allErrors[i][unsigned(Index)]) );
			Flatresult.push_back( (allPulls[i][unsigned(Index)]) );
			Flatresult.push_back( (allResults[(unsigned)Index]->GetResultParameterSet()->GetResultParameter(allNames[i])->GetMinimum() ) );
			Flatresult.push_back( (allResults[(unsigned)Index]->GetResultParameterSet()->GetResultParameter(allNames[i])->GetMaximum() ) );
			Flatresult.push_back( (allResults[(unsigned)Index]->GetResultParameterSet()->GetResultParameter(allNames[i])->GetStepSize() ) );
			Flatresult.push_back( (allResults[(unsigned)Index]->GetResultParameterSet()->GetResultParameter(allNames[i])->GetOriginalValue() ) );
		}
		Flatresult.push_back( ( allResults[(unsigned)Index]->GetResultParameterSet()->GetResultParameter(allNames[i])->GetScanStatus() ? 1.0 : 0.0 ) );
	}

	Flatresult.push_back(allRealTimes[unsigned(Index)]);
	Flatresult.push_back(allCPUTimes[unsigned(Index)]);
	Flatresult.push_back(allResults[unsigned(Index)]->GetFitStatus());
	Flatresult.push_back(allResults[unsigned(Index)]->GetMinimumValue());
	return Flatresult;
}

TString FitResultVector::GetFlatResultHeader() const
{
	TString header = "";
	for(unsigned short int i = 0; i<allNames.size(); ++i)
	{
		TString name = allNames[i];
		header += name + "_value:";
		if( (allResults[ 0 ]->GetResultParameterSet()->GetResultParameter(allNames[i])->GetType() != "Fixed" ) || allResults[ 0 ]->GetResultParameterSet()->GetResultParameter(allNames[i])->GetScanStatus() )
		{
			header += name + "_error:";
			header += name + "_pull:";
			header += name + "_min:";
			header += name + "_max:";
			header += name + "_step:";
			header += name + "_gen:";
		}
		header += name + "_scan:";
	}
	header += "Fit_RealTime:Fit_CPUTime:Fit_Status:NLL";
	return header;
}

void FitResultVector::Print() const
{
	cout << "FitResultVector:" << endl;
	for( vector<FitResult*>::const_iterator result_i = allResults.begin(); result_i != allResults.end(); ++result_i )
	{
		(*result_i)->Print();
	}
}

