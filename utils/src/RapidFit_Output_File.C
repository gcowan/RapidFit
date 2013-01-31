#include "TTree.h"
#include "TString.h"

#include "EdStyle.h"
#include "RapidFit_Output_File.h"
#include "Histo_Processing.h"
#include "TTree_Processing.h"
#include "StringOperations.h"
#include "Template_Functions.h"

#include <vector>
#include <string>
#include <cstdlib>

using namespace::std;

double RapidFit_Output_File::GetCV_val( TTree* input_tree, TString param_name, bool append )
{
	input_tree->SetEstimate(1);
	TString DrawStr( param_name ); if( append ) DrawStr.Append("_value");
	input_tree->Draw( DrawStr, "", "goff", 1, 0 );
	double returnable = input_tree->GetV1()[0];
	input_tree->SetEstimate(input_tree->GetEntries() );
	return returnable;
}
	

pair<double,double> RapidFit_Output_File::GetCV_dat( TTree* input_tree, TString param_name )
{
	input_tree->SetEstimate(1);
	input_tree->Draw( param_name+"_value:"+param_name+"_error", "", "goff", 1, 0 );

	pair<double,double> returnable = make_pair( input_tree->GetV1()[0], input_tree->GetV2()[0] );
	input_tree->SetEstimate( input_tree->GetEntries() );
	return returnable;
}

bool RapidFit_Output_File::HasToys( TTree* input_tree, vector<string> controlled_parameter_name, TRandom3* rand )
{
	TString Grid_Draw_String = Construct_Draw_String( controlled_parameter_name );

	TString Grid_Cut_String = Construct_Cut_String( input_tree, controlled_parameter_name, false ); 

	//	Get 10 unique points that correspond to data in this file, don't need to read the whole file as 10 toys is evidence of toys and quicker
	vector<vector<double> > toyz = TTree_Processing::Plotter_Data( input_tree, Grid_Draw_String, Grid_Cut_String, rand, 10 );

	bool returnable=true;

	if( toyz.empty() ) returnable = false;
	if( !toyz.empty() && toyz[0].empty() ) returnable = false;

	return returnable;
}

//      Construct a draw string which will only draw the control parameter(s)
TString RapidFit_Output_File::Construct_Draw_String( vector<string>& controlled_parameter_name )
{
	vector<string>::iterator name_i = controlled_parameter_name.begin();
	TString Grid_Draw_String = controlled_parameter_name[0].c_str();
	Grid_Draw_String.Append("_value");
	name_i = controlled_parameter_name.begin();
	++name_i;
	for( ; name_i != controlled_parameter_name.end(); ++name_i )
	{
		Grid_Draw_String.Append(":");
		Grid_Draw_String.Append( name_i->c_str() );
		Grid_Draw_String.Append("_value");
	}
	//Grid_Draw_String.Append(":NLL");
	return Grid_Draw_String;
}

//      Construct a Cut String which will only accept fixed scan parameters and will ignore the CV
TString RapidFit_Output_File::Construct_Cut_String( TTree* input_tree, vector<string>& controlled_parameter_name, bool option )
{
	vector<TString> global_param_cuts;

	vector<string>::iterator name_i = controlled_parameter_name.begin();

	//print( controlled_parameter_name );
	//exit(0);

	for( unsigned int counter=0; name_i != controlled_parameter_name.end(); ++name_i, ++counter )
	{
		TString plot_val = *name_i;
		input_tree->Draw( plot_val+"_gen", "", "goff", 1, 0 );
		double global_gen = input_tree->GetV1()[0];
		TString CV_gen; CV_gen +=global_gen;
		TString param_cut = "( abs( " + TString(name_i->c_str()) + "_gen - "+CV_gen+" ) >" + double_tolerance + ")";
		if( option )
		{
			param_cut.Append( "&&(abs("+TString(name_i->c_str())+"_gen-"+TString(name_i->c_str())+"_value)<1E-5)");
		}
		else
		{
			param_cut.Append( "&&(abs("+TString(name_i->c_str())+"_gen-"+TString(name_i->c_str())+"_value)>1E-5)");
		}
		global_param_cuts.push_back( param_cut );
	}

	TString Grid_Cut_String = global_param_cuts[0];
	vector<TString>::iterator cut_i = global_param_cuts.begin();
	++cut_i;
	for( ; cut_i != global_param_cuts.end(); ++cut_i )
	{
		Grid_Cut_String.Append( "&&" );
		Grid_Cut_String.Append( *cut_i );
	}
	//Grid_Cut_String.Append("&&(Fit_Status==3)");

	return Grid_Cut_String;
}

//	Construct a string which will only give the fits to a data file at this given coordinate in a scan
TString RapidFit_Output_File::Data_At_Grid_ij( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate )
{
	vector<string>::iterator param_i = controlled_parameter_name.begin();
	vector<TString> data_cuts;
	vector<Double_t>::iterator grid_ij = coordinate.begin();
	for( ; param_i != controlled_parameter_name.end(); ++param_i, ++grid_ij )
	{
		TString grid_str; grid_str+=*grid_ij;
		TString temp_data_cut = "( abs( " + TString(*param_i) + "_value - " + grid_str + " ) < " + double_tolerance + ")&&(" + TString(*param_i) + "_scan == 1)";
		data_cuts.push_back( temp_data_cut );
	}
	TString data_cut = data_cuts[0];
	vector<TString>::iterator data_cut_i = data_cuts.begin();
	++data_cut_i;
	for( ; data_cut_i != data_cuts.end(); ++data_cut_i )
	{
		data_cut.Append("&&");
		data_cut.Append( *data_cut_i );
	}
	return data_cut;
}

//	Construct a string which applys a cut to only get all good fits at a given coordinate in n fixed parameters
TString RapidFit_Output_File::Fits_At_Grid_ij( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, bool generate_toys )
{
	TString conditional;
	if( generate_toys == true ) conditional = "<";
	else conditional = ">";
	vector<string>::iterator param_i = controlled_parameter_name.begin();
	vector<TString> gen_cuts;

	vector<Double_t>::iterator grid_ij = coordinate.begin();

	for( ; param_i != controlled_parameter_name.end(); ++param_i, ++grid_ij )
	{
		TString grid_str; grid_str+=*grid_ij;
		TString temp_gen_cut = "( abs( " + TString(param_i->c_str()) + "_gen - " + grid_str + ") "+conditional+" " + double_tolerance + ")";
		gen_cuts.push_back( temp_gen_cut );
	}
	TString gen_cut = gen_cuts[0];
	vector<TString>::iterator cut_i = gen_cuts.begin();
	++cut_i;
	for( ; cut_i != gen_cuts.end(); ++cut_i )
	{
		gen_cut.Append("&&");
		gen_cut.Append( *cut_i );
	}
	gen_cut.Append("&&(Fit_Status==3)");
	return gen_cut;
}

//	Return the fixed/free conditions which uses all control parameters
TString RapidFit_Output_File::Construct_Fixed_condition( vector<string>& controlled_parameter_name, bool Fixed )
{
	//TString conditional, conditional2;
	//if( Fixed ) { conditional = "=="; conditional2 = "<"; }
	//else { conditional = "!="; conditional2 = ">"; }
	vector<string>::iterator param_i = controlled_parameter_name.begin();
	vector<TString> fix_conditions;

	for( ; param_i != controlled_parameter_name.end(); ++param_i )
	{
		//TString temp_fixed_str = "( abs( " + TString( *param_i ) + "_value - " + TString( *param_i ) + "_gen ) " + conditional + " " + double_tolerance + " )";
		TString temp_fixed_str;
		//TString conditional;
		if( Fixed )
		{
			temp_fixed_str.Append( "(" + /* "abs(" + */ TString( *param_i ) + "_value == " + TString( *param_i ) + "_gen ) " ); //< " + double_tolerance + " )" );
			temp_fixed_str.Append( "&&(" + TString( *param_i ) + "_scan != 1 )" );
			//conditional.Append( "<=" );
		}
		else
		{
			//temp_fixed_str.Append( "( abs( " + TString( *param_i ) + "_value - " + TString( *param_i ) + "_gen ) >= " + double_tolerance + ")" );
			temp_fixed_str.Append( "( " + TString( *param_i ) + "_value != " + TString( *param_i ) + "_gen )" );
			//conditional.Append( ">" );
		}
		//temp_fixed_str.Append( "&&( abs( " + TString( *param_i ) + "_error ) " + conditional + " " + double_tolerance + ")" );
		fix_conditions.push_back( temp_fixed_str );
	}
	vector<TString>::iterator fix_i = fix_conditions.begin();
	TString fix_condition = fix_conditions[0];
	++fix_i;
	for( ; fix_i != fix_conditions.end(); ++fix_i )
	{
		fix_condition.Append("&&");
		fix_condition.Append( *fix_i );
	}
	fix_condition.Append("&&(Fit_Status==3)");
	return fix_condition;
}

//	Return a vector of the Free non scanned Parameters in the input_tree
vector<TString> RapidFit_Output_File::get_free_non_scanned_parameters( TTree* input_tree, vector<string>& controlled_parameter_name )
{
	vector<TString> Free_Params = get_free_parameters( input_tree );
	vector<vector<TString>::iterator> remove_list;
	for( vector<TString>::iterator param_i = Free_Params.begin(); param_i != Free_Params.end(); ++param_i )
	{
		for( vector<string>::iterator c_param_i = controlled_parameter_name.begin(); c_param_i != controlled_parameter_name.end(); ++c_param_i )
		{
			if( c_param_i->compare( param_i->Data() ) == 0 ) remove_list.push_back( param_i );
		}
	}

	for( vector<vector<TString>::iterator>::iterator r_i = remove_list.begin(); r_i != remove_list.end(); ++r_i )
	{
		Free_Params.erase( *r_i );
	}
	return Free_Params;
}

//	Return a vector of only the name of the Parameters in the TTree
vector<TString> RapidFit_Output_File::get_free_parameters( TTree* local_tree )
{
	vector<TString> all_free_values = get_free_parameter_values( local_tree );
	for( vector<TString>::iterator i=all_free_values.begin(); i != all_free_values.end(); ++i )
	{
		(*i) = EdStyle::Remove_Suffix( *i );
	}
	return all_free_values;
}


//	Return a vector of parameter_value strings for only the free parameters in a TTree
vector<TString> RapidFit_Output_File::get_free_parameter_values( TTree* local_tree )
{
	//      Get a list of all branches in allresults
	vector<TString> all_parameters = TTree_Processing::get_branch_names( local_tree );
	//      Get a list of all branches in allresults with '_value' in their name
	vector<TString> all_parameter_values = StringOperations::filter_names( all_parameters, string(value_suffix.Data()) );
	vector<TString> all_parameter_errors = StringOperations::filter_names( all_parameters, string(error_suffix.Data()) );

	vector<TString> to_be_removed;
	for( vector<TString>::iterator j=all_parameter_values.begin(); j!= all_parameter_values.end(); ++j )
	{
		TString temp = EdStyle::Remove_Suffix( *j );
		temp+=error_suffix;
		//      If I have a value, but NOT an error then it was a fixed parameter
		if( StringOperations::VectorContains( all_parameter_errors, temp ) == -1 )
		{
			to_be_removed.push_back( *j );
		}
		else
		{
			//	If the maximum error for this parameter in the whole study is <= 0
			//	then the parameter was either fixed or was ignored in fitting so we should ignore it here
			double max_error = local_tree->GetMaximum( temp );
			if( max_error <= 0. )
			{
				to_be_removed.push_back( *j );
			}
		}
	}

	for( vector<TString>::iterator j= to_be_removed.begin(); j!= to_be_removed.end(); ++j )
	{
		int position = StringOperations::VectorContains( all_parameter_values, *j );
		if( position != -1 ) all_parameter_values.erase( all_parameter_values.begin() + position );
	}
	return all_parameter_values;
}

//      Get All Control Parameters from this TTree
vector<string> RapidFit_Output_File::get_control_parameters( TTree* input_tree )
{
	vector<string> controlled_parameters;

	//      Get a list of all branches in allresults
	vector<TString> all_parameters = TTree_Processing::get_branch_names( input_tree );
	//      Get a list of all branches containing information on the controlled parameters in the scan
	vector<TString> all_parameter_status = StringOperations::filter_names( all_parameters, string(ScanStatus_suffix.Data()) );

	vector<TString>::iterator status_i = all_parameter_status.begin();

	for( ; status_i != all_parameter_status.end(); ++status_i )
	{
		//cout << *status_i << endl;
		input_tree->SetEstimate( input_tree->GetEntries() );
		input_tree->Draw( *status_i, "", "goff" );
		double* temp_p = input_tree->GetV1();
		double running_total=0.;
		for( unsigned int j=0; j< input_tree->GetEntries(); ++j )
		{
			if( fabs( temp_p[j] - 1. ) < DOUBLE_TOLERANCE ) running_total+=temp_p[j];
		}
		if( fabs( running_total ) > 0.00001 )
		{
			controlled_parameters.push_back( status_i->Data() );
		}
	}
	return controlled_parameters;
}

//      Check that the Global Minima as defined at entry 0 within the file is actually the best minima or not
void RapidFit_Output_File::Check_Minima( TTree* input_tree, TString Cut_String, Float_t* Global_Best_NLL, TString NLL, TString param1_val, TString param2_val )
{
	TTree* wanted = input_tree->CopyTree( Cut_String, "fast", input_tree->GetEntries() );

	Float_t Global_Min_NLL = Float_t( wanted->GetMinimum( NLL ) );

	//      Only if we have a better minima Minuit did NOT find
	if( (Global_Min_NLL-*Global_Best_NLL) < 0 )
	{
		TString Global_Min_NLL_Str;
		Global_Min_NLL_Str+=Global_Min_NLL;
		TString Catch( "abs(" + NLL + "-" + Global_Min_NLL_Str + ")<" + double_tolerance);
		TTree* local_best = wanted->CopyTree( Catch, "fast", wanted->GetEntries() );
		//      GetMinimum == GetMaximum == Get 0th event
		double true_X = local_best->GetMinimum(param1_val);
		double true_Y = local_best->GetMinimum(param2_val);
		cout << "\n\t\tWARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING\n"<<endl;
		cout << "\t\tTRUE MINIMUM NLL = " << Global_Min_NLL << "\t\tAt:\t\tX:" << true_X << "\tY:\t" << true_Y <<endl<<endl;
		//              cout << "\tNEW MINIMA FOUND AT :\tX:\t" << X_true_min << "\tY:\t" << Y_true_min << endl<<endl;
		cout << "\t\tWARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING\n"<<endl;
		*Global_Best_NLL=Global_Min_NLL;
	}
}

