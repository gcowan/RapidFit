#ifndef DOFC_ANALYSIS_H
#define DOFC_ANALYSIS_H

#include "TString.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "OutputPlots.h"

#include <vector>
#include <string>

using namespace::std;

class FeldmanCousinsAnalysis
{
	public:

		/*!
		 * @brief Perform the full analysis of a FC output file from RapidFit
		 */
		static void DoFCAnalysis( TTree* input_tree, vector<string>& controlled_parameter_name, TRandom3* rand, vector<string>& OtherOptions );

		/*!
		 * @brief Return a pair of vectors first=CV of each control parameter second=CV_err of each control parameter
		 */
		static pair<vector<double>*, vector<double>* > Get_Best_CV( TTree* input_tree, vector<string>& controlled_param_name );

		/*!
		 * @brief Print some info about the toys at this coordinate
		 */
		static void Print_Info_At_Grid_ij( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, vector<vector<Double_t> > this_Grid_Coordinate );

		/*!
		 * @brief Print some useful output about the number of toys at this coordinate
		 */
		static void Print_Toy_Info( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, int Fixed, int Free );

		/*!
		 * @brief Construct a string which describes this unique coordinate with no spaces
		 */
		static TString Flat_Coord_String( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate );

		/*!
		 * @brief Construct multiple plots for the nuicence parameter distributions at this point
		 */
		static void Plot_Nuisance_Parameters( TTree* input_tree, TString& free_toys, TString& fixed_toys, vector<TString>& Free_Params, vector<string> controlled_parameter_name, vector<Double_t> coordinate, TRandom3* rand, double DLL_Data );

		/*!
		 * @brief Print some useful information to the screen&file about this grid point
		 */
		static void Output_GridPoint( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, double LOCAL_DATA_DLL, vector<Double_t>& Toy_DLL_Dist );

		/*!
		 * @brief Construct and write out the nll distributions of the toys to root file and as plots
		 */
		static void NLL_dists( vector<Double_t>& param_data_free, vector<Double_t>& param_data_fixed, vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, TRandom3* rand );

		/*!
		 * @brief Plot the data in a 1D output format
		 */
		static vector<OutputPlots*> Plot_1D( vector<pair<double, double> > ALL_CL_FROM_FC, vector<double> DATA_DLL, double Global_CV, double Global_CV_err, TRandom3* rand );

		/*!
		 * @brief Plot the data in a 2D output format
		 */
		static vector<OutputPlots*> Plot_2D( vector<pair<vector<double>,double> > ALL_CL_FROM_FC, vector<pair<vector<double>,double> > DATA_DLL, vector<double> Global_CV, vector<double> Global_CV_err, TRandom3* rand );

		/*!
		 * @brief Construct an OutputPlot based on the vectors for FC data
		 */
		static OutputPlots* Plot_1D_Vectors( pair<vector<double>,vector<double> > theory, pair<vector<double>,vector<double> > data, pair<vector<double>,vector<double> > FC, TRandom3* rand );

};

#endif

