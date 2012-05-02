//	MCStudy Class to sequentially step through a data file and perform fits on subsets

#pragma once
#ifndef MCStudy_H
#define MCStudy_H

//	RapidFit Headers
#include "IStudy.h"
#include "XMLConfigReader.h"
#include "FitResultVector.h"
//	System Headers
#include <vector>
#include <string>

using namespace::std;

class MCStudy	:	public IStudy
{
	public:
		//	Public Constructors
		MCStudy( XMLConfigReader* );			//	Read defaults from XML
		MCStudy( XMLConfigReader* , vector<int>, int, vector<int> );	//	Read defaults, apart from events_to_step_over, num_repeats, starting_entry
		//	Destructor
		~MCStudy();

		void DoWholeStudy( int = -999 );				//	Perform the Study
		FitResultVector* GetStudyResult();		//	Get the output of the study

		void SetNumRepeats( int );			//	Set number of Repeats
		void SetNumEvents( vector<int> );		//	Set number of events_to_step_over
		void SetStartingEntry( vector<int> );		//	Set starting entry of study
		void SetPhysParams( vector<string> );		//	Set runtime Physics Parameters
		void SetCommandLineParams( vector<string> );	//	Set Command Line Physics Parameters


	private:
		//	Don't allow copying as we have pointers to objects and these would likely not copy well!
		MCStudy ( const MCStudy& );
		MCStudy& operator = ( const MCStudy& );

		//	Internal Parameters
		XMLConfigReader* input_xml;
		vector<int> events_to_step_over;
		vector<string> CommandLineParams;
		vector<int> StartingEntries;

		vector<string> PhysParams;
};

#endif

