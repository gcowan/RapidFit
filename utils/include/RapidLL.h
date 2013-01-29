
#ifndef RAPID_LL_H
#define RAPID_LL_H

#include "TString.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TGraph.h"

#include <vector>
#include <string>

using namespace::std;

class RapidLL
{
	public:
		static unsigned int GetFunctionLineWidth();

		static unsigned int GetAxisWidth();

		static void Help();

		static TGraph* PlotRapidLL( TString controlled_parameter, TTree* input_trees, TRandom3* rand_gen, vector<string> other_params );

		//	Do the work to extract a pair of vectors
		//	X data is stored in output.first
		//	Y data is stored in output.second
		static pair<vector<double>,vector<double> > LL_Plot_Histo( TTree* input_TTree, TString Cut_String, double Global_Best_NLL, TString NLL, TString param );

		//	Do everyting required to return a 1DLL plot from the given 'param' branch within 'input_TTree'
		static TGraph* LL_Plot( TTree* input_TTree, TString Cut_String, double Global_Best_NLL, TString NLL, TString param, TRandom3* rand );

		static void OverlayMutliplePlots( TMultiGraph* graph );
};

#endif

