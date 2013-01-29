
#ifndef RAPID_2DLL_H
#define RAPID_2DLL_H

#include "TString.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH2D.h"
#include "TGraph2D.h"

#include <vector>
#include <string>

using namespace::std;

class Rapid2DLL
{
	public:

		static void Help();

		//	This class has been completely re-written again from the ground up as the old codebase was simply impossible to folllow/maintain

		static unsigned int GetFunctionLineWidth();

		static unsigned int GetAxisWidth();

		static double GetLegTextSize();

		static unsigned int GetColors( unsigned int input );

		static unsigned int GetStyle( unsigned int input );

		static vector<pair<double,TString> > GetContour( TString input  );

		static int PlotRapidFit2DLL( TString controlled_parameter1, TString controlled_parameter2, TTree* input_tree, TRandom3* rand_gen, vector<string> other_params );

		static void Plot_Free_Parameters( TTree* input_tree, TString controlled_parameter1, TString controlled_parameter2, TRandom3* rand );

		static void Plot_Contours( TString controlled_parameter1, TString controlled_parameter2, TH1* nll_hist,
					vector<pair<TMultiGraph*,TString> > nll_contours, TString filename, TRandom* rand, vector<string> other_options );

	private:

		static void get_Plotting_Data( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom* rand, vector<vector<Double_t> >& nll_data,
					vector<vector<Double_t> >& nll_data_rotated, vector<vector<Double_t> >& coords, vector<vector<Double_t> >& );

};

#endif

