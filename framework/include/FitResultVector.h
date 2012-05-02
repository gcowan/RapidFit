/*!
 * @class FitResultVector
 *
 * The result of a toy study.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
*/

#pragma once
#ifndef TOY_STUDY_RESULT_H
#define TOY_STUDY_RESULT_H

#define LLSCAN_FIT_FAILURE_VALUE -9999

//	ROOT Headers
#include "TStopwatch.h"
//	RapidFit Headers
#include "FitResult.h"
//	System Headers
#include <vector>
#include <string>

using namespace::std;

class FitResultVector
{
	public:
		/*!
		 * Create a FitResultVector from a vector of FitResults
		 */
		FitResultVector( vector<FitResultVector*> );

		/*!
		 * Create a FitResultVector based on the Input ResultParameter Names
		 */
		FitResultVector( vector<string> );

		/*!
		 * Destructor
		 */
		~FitResultVector();

		/*!
		 * Use the operator[] in the same way as GetFitResult[i], for lazy programmers ;)
		 */
		const FitResult* operator []( unsigned int index )const
		{
			return allResults[index];
		}

		/*!
		 * Add a new FitResult to the FitResultVector
		 * By default this looks to calculate the internal clock difference due to the length of time to perform a fit
		 * To Stop this behaviour, pass (FitResult*, false)
		 */
		bool AddFitResult( FitResult*, bool=true );

		/*!
		 * The Number of FitResults in this FitResultVector
		 */
		unsigned int size();

		/*!
		 * The Number of FitResults in this FitResultVector
		 */
		int NumberResults();

		/*!
		 * Get a requested FitResult from the Vector
		 */
		FitResult * GetFitResult(int);

		/*!
		 * Set the Stat of the 'stopwatch' to time how long a fit takes
		 */
		void StartStopwatch();

		/*!
		 * Get the Names of all of the Parameters in the FitResultVector
		 * This includes all of the Parameters in the FitResult as well as some information on elapsed real and cpu time taken to converge
		 */
		vector<string> GetAllNames();

		/*!
		 * Contains the length of 'real' time taken for all of the fits to converge
		 */
		vector<double> GetAllRealTimes();

		/*!
		 * Contains the length of 'CPU' time taken for all of the fits to converge
		 */
		vector<double> GetAllCPUTimes();

		/*!
		 * Contains the final value 'NLL' of the Function after Minimisation
		 */
		vector<double> GetAllMLL();

		/*!
		 * Get a vector containing the value from all of the FitResults for the Requested ResultParameter
		 */
		vector<double> GetParameterValues(string);

		/*!
		 * Get a vector containing the error from all of the FitResults for the Requested ResultParameter
		 */
		vector<double> GetParameterErrors(string);

		/*!
		 * Get a vector containing the pull from all of the FitResults for the Requested ResultParameter
		 */
		vector<double> GetParameterPulls(string);

		/*!
		 * Get the 'Flat' result of all values/errors/pulls from the requested FitResult
		 * This is useful when populating a TTree with all of the fit data
		 */
		vector<double> GetFlatResult(int);

		/*!
		 * Names of the Branches that a TTree will have in order to store all of this data, seperated by ':'
		 */
		TString GetFlatResultHeader();

		/*!
		 * Get the CPU time for the requested FitResult
		 */
		double GetCPUTime( int );

		/*!
		 * Get the Real CPU time for the requested FitResult
		 */
		double GetRealTime( int );

		/*!
		 * Fill the CPU times for all results with the given vector
		 */
		void AddCPUTimes( vector<double> );

		/*!
		 * Fill the Real times for all results with the given vector
		 */
		void AddRealTimes( vector<double> );

		/*!
		 * Add a CPU time to the last FitResult
		 */
		void AddCPUTime( double );

		/*!
		 * Add a Real time to the last FitResult
		 */
		void AddRealTime( double );

		/*!
		 * Add a CPU time to the requested FitResult
		 */
		void SetCPUTime( int, double );

		/*!
		 * Add a Real time to the requested FitResult
		 */
		void SetRealTime( int, double );

		/*!
		 * Output some debugging info
		 */
		void Print() const;

	private:
		/*!
		 * Don't Copy the class this way!
		 */
		FitResultVector ( const FitResultVector& );

		/*!
		 * Don't Copy the class this way!
		 */
		FitResultVector& operator = ( const FitResultVector& );

		/*!
		 * A vector of pointers to all of the FitResults associated with this FitResultVector
		 */
		vector< FitResult* > allResults;

		/*!
		 * A vector containing the names of all of the ResultParameters present in all FitResults
		 */
		vector<string> allNames;

		/*!
		 * Internal copy of all of the Result Parameter values, errors, pulls and true values
		 *
		 * Was this wise storing a COPY of ALL of this data in the class?
		 *
		 * Some of these have been replaced with functions to assemble the objects on request by polling the internal FitResults, but this hasn't been completed for all objects yet
		 * This allows you to work with just vectors not nested vectors which can confuse things
		 */
		vector<vector<double> > allValues, allErrors, allPulls, allGenValues;

		vector<double> allRealTimes;		/*!	All of the Real times for all FitResults			*/
		vector<double> allCPUTimes;		/*!	All of the CPU times for all FitResults				*/
		TStopwatch * clock;			/*!	Stop Watch to time the length of time to generate a FitResult	*/
};

#endif

