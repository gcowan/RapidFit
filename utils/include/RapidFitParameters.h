#ifndef RAPIDFITPARAMETERS_H
#define RAPIDFITPARAMETERS_H

#include "TString.h"
#include "TTree.h"
#include "TRandom3.h"

#include <vector>
#include <string>

using namespace::std;

//	Collection of functions for handling the fit parameters

const TString Fit_Status = "Fit_Status";
const TString NLL = "NLL";
const TString notgen = "-9999.";
const TString error_suffix = "_error";
const TString value_suffix = "_value";
const TString pull_suffix = "_pull";
const TString gen_suffix = "_gen";
const TString ScanStatus_suffix = "_scan";
const TString double_tolerance = "0.0000001";

//	Does the input tree contain toys in the fit or only results from fits to data
bool HasToys( TTree* input_tree, vector<string> controlled_parameter_name, TRandom3* rand );

//      Construct a draw string which will only draw the control parameter(s)
TString Construct_Draw_String( vector<string>& controlled_parameter_name );

//      Construct a Cut String which will only accept fixed scan parameters and will ignore the CV
TString Construct_Cut_String( TTree* input_tree, vector<string>& controlled_parameter_name );

//      Return the fixed/free conditions which uses all control parameters
TString Construct_Fixed_condition( vector<string>& controlled_parameter_name, bool Fixed );



//      Construct a string which will only give the fits to a data file at this given coordinate in a scan
TString Data_At_Grid_ij( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate );

//      Construct a string which applys a cut to only get all good fits at a given coordinate in n fixed parameters
TString Fits_At_Grid_ij( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, bool generate_toys );

//	Get all Free Parameters in a RapidFit output file that DO NOT include any scanned parameters These do not contain the '_value' extention
vector<TString> get_free_non_scanned_parameters( TTree* input_tree, vector<string>& controlled_parameter_name );

//	Get ALL Free Parameters in a RapidFit output file without the extention '_value'
vector<TString> get_free_parameters( TTree* local_tree );

//	Get All Free Parameters with the extention '_value' on the end of the parameter name
vector<TString> get_free_parameter_values( TTree* local_tree );


#endif
