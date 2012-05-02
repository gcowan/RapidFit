#ifndef RAPIDFIT_TOY_H
#define RAPIDFIT_TOY_H

#include "TString.h"
#include "TRandom3.h"
#include "TTree.h"

#include <vector>

using namespace::std;

class ToyStudyAnalysis
{
	public:
		static int Toy_Study( TTree* input_trees, TRandom3* rand_gen, vector<string> OtherOptions );

};

#endif

