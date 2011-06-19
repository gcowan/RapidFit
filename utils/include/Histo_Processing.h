#ifndef _HISTO_PROCESS
#define _HISTO_PROCESS

//	ROOT Headers
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TPolyMarker3D.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraph2D.h"
//	System Headers

#include <vector>

using namespace::std;

//  This has been adapted from the original code in RapidFits Statistics code
//  It is intended to take a histogram and automatically rebin according to this function
//  As such it uses as much inbuilt functionality in root as possible
//Return the ideal number of bins for a histogram of a vector of doubles
////Uses D. Scott's method, published 1979
int OptimumBinNumber( TH1* input_hist, int axis=1 );

//	Return the optimal number of bins for a given axis
//	1) X	2) Y	3) Z
unsigned int GetOptimalBins( TH1* input_hist, int axis=1 );

//	Function which returns the Unique corrdinates contained within a TPolyMarker3D* set of 3d points
vector<vector<Float_t> > Unique_Coords( TPolyMarker3D *pm );

//	Main plotting algorithm for plotting a TGraph2D from the input_tree based on the other parameters passed to it
//	This removes degenerate datapoints as ROOT throws it's toys out of the pram if you don't...
TGraph2D* Plotter( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom3* random = new TRandom3() );

//      Produce Plotting Histograms from an input TTree
TH2D* Plot_From_Cut( TTree* wanted_tree, TString Draw_String, TString Cut_String, TRandom3* random, TString param1="", TString param2="" );


//      Produce plots for all Physics Parameters stored within the given TTree
//      In theory could also extract the parameters from the input TTree, but as the user likely knows what they want, ask-em
void Physics_Plots( vector<TString> all_parameter_values, Float_t* best_fit_values, TTree* input_tree, TRandom3* rand_gen, TString Param1_Param2, bool CV_Drift, TH2** Physics_Param_Plots, TString Cut_String);

//	Output the Nuisence parameters to GRaphs
void Finalize_Physics_Plots( TH2* All_Physics_Plots[], vector<TString> all_parameter_values, TString param1string, TString param2string, TString outputdir, bool CV_Drift );

//	Plot a the unique coordinates on a grid to see where there is data to be analysed
void LL2D_Grid( TTree* input_tree, TString Cut_String, TString param1_val, TString param2_val, TRandom3* random, TString Suffix, TString outputdir );

//	Plot a TH2 object with various contours and styles
void Plot_Styled_Contour( TH2* input_hist, int cont_num, double* input_conts, double* confs, TString outputdir, TString Name );

//	DANGEROUS
void Plot_Styled_Contour2( TGraph2D* input_graph, int cont_num, double* input_conts, double* confs, TString outputdir, TString Name );

//	Plot 2 TH2 object contours on the same TCanvas
void Plot_Both( TH2* pllhist, TH2* FC_Plot, int nconts, double* fcconts, double *llconts, double* confs, TString outputdir, TString Legend_Name_1="NLL", TString Legend_Name_2="FC" );

//	Perform the full FC analysis on a RapidFit dataset
TH2D* FC_TOYS( TTree* input_tree, TString Fit_Cut_String, TString param1, TString param2, TString NLL, TString Fit_Cut, double NLL_Global_Best, TTree* FC_Output, TString Double_Tolerance, TRandom3* random );

//	Analyse the output TTree from the FC analysis in FC_TOYS
void FC_Stats( TTree* FC_Output, TString param1, TString param2, TRandom3* rand, TString outputdir );

//	return the *UNIQUE* corrdinates contained in the input_tree from the Draw_String after applying the Cut_String
vector<vector<Float_t> > Plotter_Data( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom3* random );

//	Create a ttree from a vector of vectors
TTree* vecvec2TTree( vector<vector<Float_t> > input_vec );

//	
TGraph2D* Plot_From_Cut_lo( TTree* wanted_tree, TString Draw_String, TString Cut_String, TRandom3* random, TString param1, TString param2 );

#endif
