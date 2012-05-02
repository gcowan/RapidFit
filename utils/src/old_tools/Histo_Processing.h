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
#include "TCanvas.h"
//	Utils Headers
#include "Template_Functions.h"
//	System Headers
#include <vector>


using namespace::std;

bool double_equals_test(Double_t first, Double_t second);

//  This has been adapted from the original code in RapidFits Statistics code
//  It is intended to take a histogram and automatically rebin according to this function
//  As such it uses as much inbuilt functionality in root as possible
//Return the ideal number of bins for a histogram of a vector of doubles
////Uses D. Scott's method, published 1979
int OptimumBinNumber( TH1* input_hist, int axis=1 );

//	Return the optimal number of bins for a given axis
//	1) X	2) Y	3) Z
unsigned int GetOptimalBins( TH1* input_hist, int axis=1 );

//	Optimally rebins the histogram you provide
//	By default this operates on the X axis assuming that you have a 1D histo but can rebin each axis independently upto 3D plots
void OptimallyRebin( TH1* input_hist, int axis=1 );

//      This will find all of the unique coordinates stored in a 1D vector of doubles
vector<vector<Double_t> > Unique_Coords( vector<Double_t> input );
//      This will find all of the unique coordinates stored in a 2D vector of pairs of doubles
vector<vector<Double_t> > Unique_Coords( vector<pair<Double_t,Double_t> > input );
//      This will find all of the unique coordinates stored in a n-D vector of doubles		UNTESTED UNTESTED UNTESTED UNTESTED
vector<vector<Double_t> > Unique_Coords( vector<pair<Double_t,Double_t> > input );
//	Function which returns the Unique corrdinates contained within a TPolyMarker3D* set of 3d points
vector<vector<Double_t> > Unique_Coords( TPolyMarker3D *pm );

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

//	DANGEROUS	untested and unfinished code due to other work
//void Plot_Styled_Contour2( TGraph2D* input_graph, int cont_num, double* input_conts, double* confs, TString outputdir, TString Name );

//	Plot 2 TH2 object contours on the same TCanvas
void Plot_Both( TH2* pllhist, TH2* FC_Plot, int nconts, double* fcconts, double *llconts, double* confs, TString outputdir, TString Legend_Name_1="NLL", TString Legend_Name_2="FC" );

//	Perform the full FC analysis on a RapidFit dataset
TH2D* FC_TOYS( TTree* input_tree, TString Fit_Cut_String, TString param1, TString param2, TString NLL, TString Fit_Cut, double NLL_Global_Best, TTree* FC_Output, TString Double_Tolerance, TRandom3* random );

//	Analyse the output TTree from the FC analysis in FC_TOYS
void FC_Stats( TTree* FC_Output, TString param1, TString param2, TRandom3* rand, TString outputdir );

//	return the *UNIQUE* corrdinates contained in the input_tree from the Draw_String after applying the Cut_String
//
//	ATTENTION: THIS WILL MANGLE DATA COMPARED TO HOW IT APPEARS IN THE INPUT TREE, BUT IT GIVEC UNIQUE SORTED DATA SO I DONT CARE!
vector<vector<Double_t> > Plotter_Data( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom3* random );

//	Create a ttree from a vector of vectors
TTree* vecvec2TTree( vector<vector<Float_t> > input_vec );

//	
TGraph2D* Plot_From_Cut_lo( TTree* wanted_tree, TString Draw_String, TString Cut_String, TRandom3* random, TString param1, TString param2 );

//	Algorithm for comparing a 2D point and returning a boolean comparison
//	used for stl::unique
bool Unique_2D_Double( pair<double,double> one_pair, pair<double,double> two_pair );

//	Algorithm for comparing a 2D point with another and returning a decison on the biggest one
//	used for stl::sort
bool Sort_first_Double( pair<double,double> one_pair, pair<double,double> two_pair );
bool Sort_second_Double( pair<double,double> one_pair, pair<double,double> two_pair );
//	Will replace these with templates eventually

//	Do the work to extract a pair of vectors
//	X data is stored in output.first
//	Y data is stored in output.second
pair<vector<double>,vector<double> > LL_Plot_Histo( TTree* input_TTree, TString Cut_String, double Global_Best_NLL, TString NLL, TString param );

//	Do everyting required to return a 1DLL plot from the given 'param' branch within 'input_TTree'
TGraph* LL_Plot( TTree* input_TTree, TString Cut_String, double Global_Best_NLL, TString NLL, TString param, TRandom3* rand );


//  More correct to fit to the gamma distribution than the landau function for those that require it
TF1* gamma_func( int OutputLevel=-1 );
//  More correct to fit to the landau dist for some other dists from toys
TF1* landau_func();

//	Return the string which corresponds to the best fit function for the dataset by the best chi2
//	Evaluates and compares the gaus, gamma & landau functions from ROOT
TString Best_Fit_Function( TH1* input, int OutputLevel=-1 );

//	Fit the TH1 whilst catching a lot of unwanted output
void Silent_Fit( TH1* input_histo, TString fit_type, int OutputLevel=-1 );

//	Draw something whilst catching all of the root output to the standard streams
void Silent_Draw( TCanvas* c1, TH1* input_histo, TString options="", int OutputLevel=-1 );

//	Print whilst catching all of the root output to the standard streams
void Silent_Print( TCanvas* c1, TString Print_String, int OutputLevel=-1 );

//	Get a histogram from your draw string, after weighting from the input_tree
TH1* Get_Histo( TTree* input_tree, TString draw_str, TString weight_str, TRandom3* rand );

//	Get a graph from your input draw string, after weighting from the input_tree
TGraph* Get_Graph( TTree* input_tree, TString draw_str, TString weight_str, TRandom3* rand );

//	Get a TH1 from a vector of Doubles (could even write this to be type agnostric but I have the recast function in the template header
TH1* Get_TH1( vector<Double_t> input, TRandom3* rand, int bins=100 );

//	Get a TH2 from a vector of vectors
TH2* Get_TH2( vector<vector<Double_t> > input, TRandom3* rand, int bins1=100, int bins2=100 );

//	Get a TH3 from a vector of vectors
TH3* Get_TH3( vector<vector<Double_t> > input, TRandom3* rand, int bins1=100, int bins2=100, int bins3=100 );

#endif

