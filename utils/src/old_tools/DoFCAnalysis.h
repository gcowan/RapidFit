#ifndef DOFC_ANALYSIS_H
#define DOFC_ANALYSIS_H

#include <vector>

using namespace::std;

//      Construct a draw string which will only draw the control parameter(s)
TString Construct_Draw_String( vector<string>& controlled_parameter_name );

//      Construct a Cut String which will only accept fixed scan parameters and will ignore the CV
TString Construct_Cut_String( TTree* input_tree, vector<string>& controlled_parameter_name );

void DoFCAnalysis( TTree* input_tree, vector<string>& controlled_parameter_name, TRandom3* rand );

#endif
