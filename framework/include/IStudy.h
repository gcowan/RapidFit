
//	This class is intended as a simple interface for studies that people will use RapidFit for
//	It provides a template of 'essential' functions that every study should implement without being too strict

//	This is more of a template of how to do something systematically within RapidFit and interface it to the rest of the framwork

#ifndef IStudy_H
#define IStudy_H

//	RapidFit Headers
#include "XMLConfigReader.h"
#include "ToyStudyResult.h"
//	System Headers
#include <vector>
#include <string>

using namespace::std;

class IStudy
{
	public:

		//	Can't have virtual public Constructors, and it doesn't make sense to either
		virtual ~IStudy() {};

		virtual void DoWholeStudy() = 0;				//	Perform the Study
		virtual ToyStudyResult* GetStudyResult() = 0;			//	Get the output of the study

		virtual void SetNumRepeats( int ) = 0;				//	Set number of Repeats
		virtual void SetCommandLineParams( vector<string> ) = 0;	//	Set Command Line Physics Parameters

};

#endif
