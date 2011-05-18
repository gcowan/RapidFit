//	This is the 3rd complete re-write of the plotting code used in the analysis of 2DLL and FC plots from the Edinburgh RapidFit Fitter
//	The reasoning behind this re-write is complex but the complexity is driven by shortcomings in the root framework,
//	Whilst the speed and actual plots are a credit to the things that root does well

//	ROOT Headers
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TString.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph2D.h"
#include "TPolyMarker3D.h"
#include "TList.h"
#include "TRandom3.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TString.h"
//	RapidFit Headers
#include "EdStyle.h"
//	System Headers
#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace::std;

inline TString prettyPrint(Double_t value){
	char pretty[20];
	TString prettyString;
	sprintf (pretty, "%1.3g",value);
	prettyString = pretty;
	return prettyString;
}

inline TString prettyPrint(Double_t val, Double_t err){
	TString outstr = "$";
	char errstr[20];
	char valstr[20];
	Int_t n = sprintf (errstr, "%1.1g",err);
	sprintf(valstr,"%1.*f",n-2,val);
	outstr += valstr;
	outstr += "\\pm";
	outstr += errstr;
	outstr+= "$";
	return outstr;
}

TPaveText* addLHCbLabel(TString footer){
	//				
	TPaveText * label = new TPaveText(0.18, 0.73, 0.18, 0.88,"BRNDC");
	label->SetFillStyle(0);		//Transparent i.e. Opacity of 0 :D
	label->SetBorderSize(0);
	label->SetTextAlign(11);
	label->SetTextSize(Float_t(0.04));
	TText * labeltext = 0;
	labeltext = label->AddText("LHC#font[12]{b} 2010 Data");
	labeltext = label->AddText("#sqrt{s} = 7TeV");
	labeltext = label->AddText(footer);
	return label;
}


//	Pass this the pointer to the ntuple you're handing and it will give you a list
//	of branches that it contains in an easy to handle manner
vector<TString> get_branch_names( TTree* local_tree )
{
	//	To be populated and returned to the user
	vector<TString> temp_branch_names;

	//	Get the list of branches within the TTree
	TObjArray* branch_obj_array = local_tree->GetListOfBranches();

	//	Loop over all found branch objects and request their names
	for( unsigned short int i=0; i < branch_obj_array->GetEntries() ; ++i )
	{
		TObject* branch_object = (*branch_obj_array)[i];
		temp_branch_names.push_back((const char*) branch_object->GetName());
	}

	//	Return the vector of names I have found
	return temp_branch_names;
}

//	Pass this the name of the ntuple you're handing and it will give you a list
//	of branches that it contains in an easy to handle manner
vector<TString> get_branch_names( TString absolute_path )
{
	//	Get the tree that has been defined by the user
	TTree* local_tree = (TTree*) gDirectory->Get((const char*) absolute_path);
	//	Return a vector of Branch names as defined by the user
	return get_branch_names( local_tree );
}

//	Pass this an array of strings and it will find all strings matching a substring
vector<TString> filter_names( vector<TString> all_names, string substring )
{
	//	To be populated and returned to the user
	vector<TString> returnable_names;

	//	Address of the end of the string object
	size_t found=string::npos;

	//	Loop over all strings in the vector
	for( unsigned short int i=0; i<all_names.size(); ++i )
	{	//	Again using STL functions :)
		string temp_str = (all_names[i].Data());
		//	Attempt to find the coordinate of the substring within the string
		found = temp_str.find( substring );
		if( found!=string::npos )	//	If the substring is found
		{
			returnable_names.push_back(all_names[i]);
		}
	}
	//	Return all strings found containing the substring
	return returnable_names;
}

//	Pass this an array of strings and it will find all strings matching a substring
vector<TString> filter_names( vector<TString> all_names, TString substring )
{
	return filter_names( all_names, string(substring.Data()) );
}

//	This will find all of the unique corrdinates stored in a TPolyMarker3D object
//	NB this was written due to the fact that ROOT thows away the contents of {3/4}D TTree->Draw() objects... (God only knows why)
vector<vector<Float_t> > Unique_Coords( TPolyMarker3D *pm )
{
	//	Tolerance of the DOUBLE
	double DT = 1E-5;

	//	Get the number of coordinates in the TPolyMarker3D
	int number_of_points = pm->GetN();
	//	Get the data contained in the TPolyMarker3D
	Float_t* Coord_Data_pointer = pm->GetP();

	//	Will populate and return to the user
	vector<vector<Float_t> > Returnable_Coord_Data;

	//	Used within a for loop in logic, externally created/destroyed
	bool add_point = true;
	bool temp_decision_1 = true;
	bool temp_decision_2 = true;
	bool temp_decision_3 = true;

	//	Run over all of the points contained within the TPolyMarker3D object
	for( int i=0; i< number_of_points; ++i )
	{
		//	Assume we haven't seen this point yet
		add_point = true;
		//	Check the anzats
		for( unsigned int j=0; j< Returnable_Coord_Data.size(); ++j )
		{
			temp_decision_1 = true;
			temp_decision_2 = true;
			temp_decision_3 = true;
			if( fabs( Coord_Data_pointer[i*3] - Returnable_Coord_Data[j][0] ) < DT ) temp_decision_1 = false;
			if( fabs( Coord_Data_pointer[i*3+1] - Returnable_Coord_Data[j][1] ) < DT ) temp_decision_2 = false;
			if( fabs( Coord_Data_pointer[i*3+2] - Returnable_Coord_Data[j][2] ) < DT ) temp_decision_3 = false;

			//	If all of the information is the same do NOT add the point to the new vector of data
			if( ( ( temp_decision_1 == temp_decision_2 ) && ( temp_decision_2 == temp_decision_3 ) ) && ( temp_decision_1 == false ) )
				add_point = false;
		}
		//	If we haven't seen this point yet add it to the array of points
		if( add_point )
		{
			//cout << Coord_Data_pointer[i*3] << "\t" << Coord_Data_pointer[i*3+1] << "\t" << Coord_Data_pointer[i*3+2] << endl;
			vector<Float_t> temp_vector;
			temp_vector.push_back( Coord_Data_pointer[i*3] );
			temp_vector.push_back( Coord_Data_pointer[i*3+1] );
			temp_vector.push_back( Coord_Data_pointer[i*3+2] );
			Returnable_Coord_Data.push_back( temp_vector );
		}
	}

	//	Return a 2D vector of unique 3D corrdinates of size npoints*3
	return Returnable_Coord_Data;
}

TGraph2D* Plotter( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom3* random )
{
	//	Plot the graph using TTree->Draw()
	//	The resulting graph (TH3) contains empty axis and a TPolyMarker3D object
	input_tree->SetEstimate(input_tree->GetEntries());  // Fix the size of the array of doubles to be created (There will never be more than this
	input_tree->Draw( Draw_String, Cut_String );

	//	Get the Points that have been plotted (TPolyMarker3D object named "TPolyMarker3D", see ROOTtalk)
	TPolyMarker3D *pm = (TPolyMarker3D*)gPad->FindObject("TPolyMarker3D");
	double temp = random->Rndm();
	TString Name = "TPoly3_";
	Name+=temp;
	pm->SetName(Name);

	//	Get a list of ONLY unique coordinates due to the short comings of the interpolation within TGraph2D
	vector<vector<Float_t> > Coord_Data = Unique_Coords( pm );

	//	Make a new EMPTY TGraph2D with a unique name
	TGraph2D* new_plot = new TGraph2D();
	TString Plot="Plot_";
	temp = random->Rndm();
	Plot+=temp;
	new_plot->SetName(Plot);
	new_plot->SetTitle(Plot);
	//	Set Binning of Plot
	new_plot->SetNpx( int( ceil( sqrt( double(Coord_Data.size()) ) ) ) );
	new_plot->SetNpy( int( ceil( sqrt( double(Coord_Data.size()) ) ) ) );

	//	Add the data to the TGraph2D object as points
	for( unsigned int i=0; i< Coord_Data.size(); ++i )
	{
		//cout << i << "\t" << Coord_Data[i][0] << "\t" << Coord_Data[i][1] << "\t" << Coord_Data[i][2] << endl;
		new_plot->SetPoint( int(i), Coord_Data[i][0], Coord_Data[i][1], Coord_Data[i][2] );
	}

	//	Return the TGraph2D object
	return new_plot;
}

//	Perform a cut on a TTree object and return a new TTree with a unique name
TTree* Cut( TTree* input_tree, TString Cut_String, TRandom3* random )
{
	//	Tree to be returned
	TTree* wanted_tree = NULL;
	//	Perform The Cut
	wanted_tree = input_tree->CopyTree( Cut_String, "fast", input_tree->GetEntries() );

	//	Generate a Unique Name
	TString Name="Tree_";
	double rand = random->Rndm();
	Name+=rand;

	//	Set Unique Name
	wanted_tree->SetName(Name);
	wanted_tree->SetTitle(Name);

	//	Return the TTree
	return wanted_tree;
}

vector<Double_t> Get_Data( TTree* input_tree, TString Cut_String, TString Data_String )
{
	input_tree->Draw( Data_String, Cut_String );
	int Data_Num = int( input_tree->GetSelectedRows() );
	Double_t* temp_pointer = input_tree->GetV1();
	vector<Double_t> Returnable_Data;
	for( int i=0; i < Data_Num; ++i )
		Returnable_Data.push_back( double(temp_pointer[i]) );
	return Returnable_Data;
}

vector<vector<Double_t> > Get_All_Data( TTree* input_tree, TString Cut_String, TString Data_String )
{
	int dimention_num = 1;
	string input_str = Data_String.Data();
	for( size_t found=0; found!=string::npos; )
	{
		found = input_str.find( ":", found+1, string::npos );
		if( found!=string::npos ) ++dimention_num;
	}
	cout << dimention_num << endl;
	input_tree->Draw( Data_String, Cut_String );
	int Data_Num = int( input_tree->GetSelectedRows() );
	Double_t* temp_pointer = input_tree->GetV1();
	vector<vector<Double_t> > Returnable_Data;
	for( int i=0; i < Data_Num; ++i )
	{
		vector<Double_t> temp_vec;
		temp_vec.push_back( double(temp_pointer[i]) );
		Returnable_Data.push_back( temp_vec );
	}
	return Returnable_Data;
}

//	Produce Plotting Histograms from an input TTree
TH2D* Plot_From_Cut( TTree* wanted_tree, TString Draw_String, TString Cut_String, TRandom3* random, TString param1="", TString param2="" )
{
	//	Create a canvas
	TString Canvas_Name("Canvas_");
	double rand = random->Rndm();
	Canvas_Name+=rand;
	TCanvas* temp_canvas = new TCanvas( Canvas_Name, Canvas_Name );
	//	Plot the Graph
	TGraph2D* new_graph = Plotter( wanted_tree, Draw_String, Cut_String, random );
	//	Update the Canvas to initialize anything that requires this step... this is very much a ROOT thing
	temp_canvas->Update();
	//	Return the Histogram from within this graph for plotting
	TH2D* Returnable_Hist = new_graph->GetHistogram();
	Returnable_Hist->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param2 ) );
	Returnable_Hist->GetYaxis()->SetTitle( EdStyle::GetParamRootName( param1 ) );
	return Returnable_Hist;
}

//	Check that the Global Minima as defined at entry 0 within the file is actually the best minima or not
void Check_Minima( TTree* input_tree, TString Cut_String, Float_t* Global_Best_NLL, TString NLL, TString Double_Tolerance, TString param1_val, TString param2_val )
{
	TTree* wanted = input_tree->CopyTree( Cut_String, "fast", input_tree->GetEntries() );

	Float_t Global_Min_NLL = Float_t( wanted->GetMinimum( NLL ) );

	//	Only if we have a better minima Minuit did NOT find
	if( (Global_Min_NLL-*Global_Best_NLL) < 0 )
	{
		TString Global_Min_NLL_Str;
		Global_Min_NLL_Str+=Global_Min_NLL;
		TString Catch( "abs(" + NLL + "-" + Global_Min_NLL_Str + ")<" + Double_Tolerance);
		TTree* local_best = wanted->CopyTree( Catch );
		//	GetMinimum == GetMaximum == Get 0th event
		double true_X = local_best->GetMinimum(param1_val);
		double true_Y = local_best->GetMinimum(param2_val);
		cout << "\n\t\tWARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING\n"<<endl;
		cout << "\t\tTRUE MINIMUM NLL = " << Global_Min_NLL << "\t\tAt:\t\tX:" << true_X << "\tY:\t" << true_Y <<endl<<endl;
//		cout << "\tNEW MINIMA FOUND AT :\tX:\t" << X_true_min << "\tY:\t" << Y_true_min << endl<<endl;
		cout << "\t\tWARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING\n"<<endl;
		*Global_Best_NLL=Global_Min_NLL;
	}

}

//	Produce plots for all Physics Parameters stored within the given TTree
//	In theory could also extract the parameters from the input TTree, but as the user likely knows what they want, ask-em
void Physics_Plots( vector<TString> all_parameter_values, Float_t* best_fit_values, TTree* input_tree, TRandom3* rand_gen, TString Param1_Param2, bool CV_Drift, TH2** Physics_Param_Plots, TString Cut_String)
{
	//	Construct a plot string for the physics parameters and plot them
	cout << endl << "STARTING PLOTS SHOWING THE VARIATION OF PHYSICS PARAMETERS THROUGHOUT THE SCANS" << endl;

	//	Store all of the plotting strings
	TString* Physics_Param_DrawString = new TString[unsigned(all_parameter_values.size())];
	//	Loop over all of the wanted parameters that have been passed
	for( unsigned int i=0; i<all_parameter_values.size(); ++i )
	{
		cout << endl << "PLOTTING: " << all_parameter_values[i] << "\t" << i+1 << " of: " << all_parameter_values.size() << endl;

		//	Construct Plotting String based on input
		Physics_Param_DrawString[i] = "(" + all_parameter_values[i];
		if( CV_Drift ) {
			TString Best_Value;
			Best_Value+=best_fit_values[i];
			Physics_Param_DrawString[i] += "-" + Best_Value;
		}
		Physics_Param_DrawString[i] += ")";
		Physics_Param_DrawString[i] += Param1_Param2;

		//	Actually perform the plots using the same tool as before
		cout << Physics_Param_DrawString[i] << endl;
		cout << Cut_String << endl;
		Physics_Param_Plots[i] = Plot_From_Cut( input_tree, Physics_Param_DrawString[i], Cut_String, rand_gen );

		//	Plot the graph
		TString Canvas_Name("Plot_Me_");
		double rand_canv = rand_gen->Rndm();
		Canvas_Name+=rand_canv;
		TCanvas* plot_me = new TCanvas( Canvas_Name, Canvas_Name, 1680, 1050 );
		plot_me->SetTitle( "" );
		plot_me->SetName( all_parameter_values[i] );
		Physics_Param_Plots[i]->Draw("colz");
		plot_me->Update();
	}
}

void Finalize_Physics_Plots( TH2* All_Physics_Plots[], vector<TString> all_parameter_values, TString param1string, TString param2string, TString outputdir, bool CV_Drift )
{
	for( unsigned short int i=0; i < all_parameter_values.size(); ++i )
	{
		TString Name("Physics_Param_");
		Name+=i;

		TCanvas* final_physics_canvas = new TCanvas( Name, Name, 1680, 1050 );

		All_Physics_Plots[i]->SetTitle("");
		All_Physics_Plots[i]->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param2string ) );
		All_Physics_Plots[i]->GetYaxis()->SetTitle( EdStyle::GetParamRootName( param1string ) );

		if( CV_Drift ) 	all_parameter_values[i].Append("_var");

		All_Physics_Plots[i]->Draw("colz");
		addLHCbLabel( all_parameter_values[i] )->Draw();

		final_physics_canvas->Update();

		//	Output Filename for the Plot
		TString Output_File( outputdir );
		Output_File+= "/" + all_parameter_values[i] + ".png";

		final_physics_canvas->Print( Output_File );
	}
	return;
}

bool Toys_Check( TTree* input_tree, TString Cut_String, TTree* Toy_Tree )
{
	Toy_Tree = input_tree->CopyTree( Cut_String, "fast", input_tree->GetEntries() );
	if( Toy_Tree->GetEntries() > 0 ) return true;
	return false;
}

TH2D* FC_TOYS( TTree* input_tree, TString Fit_Cut_String, TString param1, TString param2, TString NLL, TString Fit_Cut, double NLL_Global_Best, TTree* FC_Output, TString Double_Tolerance, TRandom3* random )
{

	TString param1_gen = param1 + "_gen";
	TString param2_gen = param2 + "_gen";
	TString param1_val = param1 + "_value";
	TString param2_val = param2 + "_value";
	TString param1_err = param1 + "_error";
	TString param2_err = param2 + "_error";

	TTree* Grid = Cut( input_tree, Fit_Cut_String, random );
	//	The TTree object from the PLL fit contains information on the corrdinates of the fit
	Float_t Param_1_Coord=0, Param_2_Coord=0, NLL_Local_Best=0;
	Grid->SetBranchAddress( param1_val, &Param_1_Coord );
	Grid->SetBranchAddress( param2_val, &Param_2_Coord );
	Grid->SetBranchAddress( NLL, &NLL_Local_Best );

	//	Easiest to store the output in a new TTree
	//	Want to store the CL calculated at every point and the efficiency of Toys at this coordinate
	Float_t CL=0., Toy_Num=0., Successful_Toys=0., Processed_Toys=0., Fit_Status=3.;
	TString CL_Branch = "CL";
	TString Success_Branch = "Success";
	TString Total_Branch = "Total";
	TString Processed_Toys_Branch = "Processed_Toys";
	FC_Output->Branch( CL_Branch, &CL );
	FC_Output->Branch( Success_Branch, &Successful_Toys );
	FC_Output->Branch( Total_Branch, &Toy_Num );
	FC_Output->Branch( param1_val, &Param_1_Coord );
	FC_Output->Branch( param2_val, &Param_2_Coord );
	FC_Output->Branch( Processed_Toys_Branch, &Processed_Toys );
	FC_Output->Branch( "Fit_Status", &Fit_Status );

	vector<vector<Float_t> > Used_Coordinate;

	
	bool Add_Point=true;
	bool decision_1=true;
	bool decision_2=true;

	//double rand=0;

	//Find the toys generated at the previously discovered gridpoints
	for( int i = 0; i < Grid->GetEntries(); ++i )
	{
		vector<Double_t> Floated_Toys_NLL;
		vector<Double_t> Fixed_Toys_NLL;

		//	Read in information about this Grid_Coordinate, coordinate and the NLL from the local best fit
		Grid->GetEntry( i );

		Add_Point = true;
		for( unsigned int j=0; j < Used_Coordinate.size(); ++j )
		{
			decision_1 = true;
			decision_2 = true;
			if( fabs( Param_1_Coord  - Used_Coordinate[j][0] ) < 1E-5 ) { decision_1 = false; }
			if( fabs( Param_2_Coord  - Used_Coordinate[j][1] ) < 1E-5 ) { decision_2 = false; }
			if( ( decision_1 == false ) && ( decision_2 == false ) ) { Add_Point = false; }
		}
		if( Add_Point ){
			vector<Float_t> temp_Coord;
			temp_Coord.push_back( Param_1_Coord );
			temp_Coord.push_back( Param_2_Coord );
			Used_Coordinate.push_back( temp_Coord );
		} else { continue; }

		//	Cut on only toys at this Grid coord
		TString Param_1_Coord_Str;
		Param_1_Coord_Str+=Param_1_Coord;
		TString Param_2_Coord_Str;
		Param_2_Coord_Str+=Param_2_Coord;
		TString Param_1_Grid = "(abs(" + param1_gen + "-" + Param_1_Coord_Str + ")<" + Double_Tolerance + ")";
		TString Param_2_Grid = "(abs(" + param2_gen + "-" + Param_2_Coord_Str + ")<" + Double_Tolerance + ")";

		TString Toys_At_Grid_Point = Param_1_Grid + "&&" + Param_2_Grid;

		Floated_Toys_NLL = Get_Data( input_tree, Toys_At_Grid_Point, NLL );
		Toy_Num = Float_t( Floated_Toys_NLL.size() );

		Toys_At_Grid_Point += "&&" + Fit_Cut;

		Floated_Toys_NLL = Get_Data( input_tree, Toys_At_Grid_Point, NLL );
		Successful_Toys = Float_t( Floated_Toys_NLL.size() );

		TString float_param_1 = "(abs(" + param1_err + ")>" + Double_Tolerance + ")";
		TString float_param_2 = "(abs(" + param2_err + ")>" + Double_Tolerance + ")";

		TString Floated_Toys_At_Grid_Point = Toys_At_Grid_Point + "&&" + float_param_1 + "&&" + float_param_2;


		TString fixed_param_1 = "(abs(" + param1_err + ")<" + Double_Tolerance + ")";
		TString fixed_param_2 = "(abs(" + param2_err + ")<" + Double_Tolerance + ")";

		TString Fixed_Toys_At_Grid_Point = Toys_At_Grid_Point + "&&" + fixed_param_1 + "&&" + fixed_param_2;

		//	Get all toys at this grid point that are Floated

		Floated_Toys_NLL = Get_Data( input_tree, Floated_Toys_At_Grid_Point, NLL );
		//TTree* Floated_Toys = input_tree->CopyTree( Floated_Toys_At_Grid_Point, "fast", input_tree->GetEntries() );
		//TString Floated_Name("Floated_");
		//rand = random->Rndm();
		//Floated_Name+=rand;
		//Floated_Toys->SetName(Floated_Name);

		//	Get all toys at this grid point that are Fixed
		Fixed_Toys_NLL = Get_Data( input_tree, Fixed_Toys_At_Grid_Point, NLL );
		//TTree* Fixed_Toys = input_tree->CopyTree( Fixed_Toys_At_Grid_Point, "fast", input_tree->GetEntries() );
		//TString Fixed_Name("Fixed_");
		//rand = random->Rndm();
		//Fixed_Name+=rand;
		//Fixed_Toys->SetName(Fixed_Name);

		int Floated_Toy_Num = int( Floated_Toys_NLL.size() );
		int Fixed_Toy_Num = int( Fixed_Toys_NLL.size() );

		if( Floated_Toy_Num != Fixed_Toy_Num )
		{
			cerr << endl << "NUMBER OF FIXED AND FLOATED TOYS AT COORDINATE:\t" << Param_1_Coord << ":" << Param_2_Coord <<endl;
			cerr << "ARE DIFFERENT, USING THE SAME SAMBLE SIZE OF EACH WHICH REDUCES THE ACCURACY" << endl << endl;
			if( Floated_Toy_Num < Fixed_Toy_Num ) Fixed_Toy_Num = Floated_Toy_Num;
			else Floated_Toy_Num = Fixed_Toy_Num;
		}

		//	By definition here!
		Processed_Toys = Float_t(Floated_Toy_Num);

		//	Remember coordinates stored in param1gridpoints and param2gridpoints
		//	With NLL values stored in NLL_Local_Best
		//	This again would be nicer&easier if associated data was stored in the one object


		//	Using GetEntry to look at each toy at each point so tell it where to store the value we get
		//Float_t Floated_NLL=0, Fixed_NLL=0;
		//Floated_Toys->SetBranchAddress( NLL, &Floated_NLL );
		//Fixed_Toys->SetBranchAddress( NLL, &Fixed_NLL );

		if( Fixed_Toy_Num != 0 ){

			//Loop over the toys, pulling out the NLL ratio
			UInt_t smaller_Toys = 0;

			Double_t Ratio = NLL_Local_Best - NLL_Global_Best;
			for(unsigned short int j = 0; j < Fixed_Toy_Num; ++j){

			//	VERY MEMORY INTENSIVE!!!
			//	Floated_Toys->GetEntry(j);
			//	Fixed_Toys->GetEntry(j);

				//THE LINE BELOW IS THE FELDMAN-COUSINS ORDERING METHOD USED BY CDF/HEIDELBERG: 
				//if the toyratio is smaller than the data ratio at this point, increment:
				if ( (Fixed_Toys_NLL[j]-Floated_Toys_NLL[j]) < Ratio )
				{
					++smaller_Toys;
				}
			}

			//The C.L. is the percentage of toys that were smaller
			CL = Float_t(smaller_Toys)/Float_t(Fixed_Toy_Num);

		} else {
			//	THIS IS SERIOUS, tell the user about it!
			cerr << "\t" << Param_1_Coord << ":" << Param_2_Coord << "\tWARNING: NO TOYS FOUND HERE! " << endl;
			CL = +9999.;	//	This plots a spike here which shows on contours
		}
	cout << Processed_Toys << "\tTOYS PROCESSED AT:\t" << setprecision(4) << Param_1_Coord << "\t:\t" << Param_2_Coord << endl;
	//	Store the relevent information for plotting in the FC_Output TTree
	FC_Output->Fill();
	}

	TString FCName="FC_Plot_";
	double temp = random->Rndm();
	FCName+=temp;

	TCanvas* FC_Plot = new TCanvas( FCName, FCName, 1680, 1050 );
	TString Param1_Param2 = ":" + param1_val + ":" + param2_val;

	FC_Output->Draw( "CL" + Param1_Param2 );
	FC_Plot->Update();
	TPolyMarker3D *pm = (TPolyMarker3D*)gPad->FindObject("TPolyMarker3D");
	temp = random->Rndm();
	TString Name = "TPoly3_";
	Name+=temp;
	pm->SetName(Name);

	vector<vector<Float_t> > Coord_Data = Unique_Coords( pm );

	TGraph2D* new_plot = new TGraph2D();
	TString Plot="FC_Plot_";
	temp = random->Rndm();
	Plot+=temp;
	new_plot->SetName(Plot);
	new_plot->SetTitle(Plot);
	//	Set Binning of Plot
	new_plot->SetNpx( int( ceil( sqrt( double(Coord_Data.size()) ) ) ) );
	new_plot->SetNpy( int( ceil( sqrt( double(Coord_Data.size()) ) ) ) );

	for( unsigned int i=0; i< Coord_Data.size(); ++i )
	{
		//cout << i << "\t" << Coord_Data[i][0] << "\t" << Coord_Data[i][1] << "\t" << Coord_Data[i][2] << endl;
		new_plot->SetPoint( int(i), Coord_Data[i][0], Coord_Data[i][1], Coord_Data[i][2] );
	}

	TH2D* Returnable_Hist = new_plot->GetHistogram();

	Returnable_Hist->GetXaxis()->SetTitle( EdStyle::GetParamRootName( param2 ) );
        Returnable_Hist->GetYaxis()->SetTitle( EdStyle::GetParamRootName( param1 ) );

	return Returnable_Hist;
}

TH1* LL2D_Grid( TTree* input_tree, TString Cut_String, TString param1_val, TString param2_val, TRandom3* random, TString Suffix )
{
	TString Name("Canvas");
	double rand = random->Rndm();
	Name+=rand;
	TCanvas* GRID = new TCanvas(Name,Name,1680,1050);
	input_tree->SetEstimate(input_tree->GetEntries());  // Fix the size of the array of doubles to be created (There will never be more than this
	TString Draw_Str = param1_val + ":" + param2_val;
	input_tree->Draw( Draw_Str, Cut_String );
	GRID->Print("Coordinate_Grid"+Suffix+".png");
	return input_tree->GetHistogram();
}

void Plot_Both( TH2* pllhist, TH2* FC_Plot, int nconts, double* fcconts, double *llconts, double* confs, TString outputdir )
{
	TCanvas* Temp_1 = new TCanvas("Temp","Temp",1680,1050);  

	TList* contLevel = NULL;
	TGraph* Line     = NULL;
	vector<TGraph*> Contour_Lines;

	//	Construct the Legend
	TLegend *leg = new TLegend(0.75,0.89,0.95,0.7);
	leg->SetHeader("Conf. Levels");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);

	//	Construct the Contours in the LL plot
	pllhist->SetContour( nconts, llconts );
	pllhist->Draw("cont LIST");
	Temp_1->Update();

	//	Get the Contours in the LL Plot
	TObjArray *LL_Contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

	//	Loop over all contours
	for( int i = 0; i < LL_Contours->GetSize(); ++i )
	{
		//	Name that contour...
		TString confname = "";
		confname += confs[i];
		confname += "% C.L. (PLL)";

		//	Get the List of lines making up this contour
		contLevel = (TList*) LL_Contours->At(i);

		//	Loop over all lines constructing this contour (max 4)
		for(int j =0; j < contLevel->GetSize(); ++j)
		{
			//	Current line
			Line = (TGraph*) contLevel->At(j);

			//	Change the contour line Style
			TGraph *gc = (TGraph*) Line->Clone();
			gc->SetLineStyle( Style_t(i+1) );

			//	Store this line
			Contour_Lines.push_back(gc);

			//	Add an entry in the Legend for it
			if( j==0 )  leg->AddEntry( gc, confname, "L");
		}
	}

	//	Construct the Legend
	FC_Plot->SetContour(nconts,fcconts);
	FC_Plot->Draw("cont LIST");
	Temp_1->Update();

	//	Get the Contours in the FC Plot
	TObjArray *FC_Contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	Temp_1->Update();

	//	Loop over all contours
	for(int i = 0; i < FC_Contours->GetSize(); ++i)
	{
		//	Name that contour...
		TString confname = "";
		confname += confs[i];
		confname += "% C.L. (FC)";

		//	Get the List of lines making up this contour
		contLevel = (TList*) FC_Contours->At(i);

		//	Loop over all lines constructing this contour (max 4)
		for(int j =0; j<contLevel->GetSize(); ++j)
		{
			//	Current line
			Line = (TGraph*) contLevel->At(j);
		
			//	Set the line Color
			TGraph *gc = (TGraph*) Line->Clone();
			gc->SetLineColor( Color_t(i+2) );

			//	Store this line
			Contour_Lines.push_back(gc);

			//	Add an entry in the Legend for it
			if(j==0)  leg->AddEntry( gc, confname, "L");
		}
	}

	TCanvas* Overlay_Output = new TCanvas( "Output_Overlay", "Output_Overlay", 1680, 1050);

	//	First construct the Axis
	pllhist->Draw("AXIS");

	//	Now plot all of the lines from all of the contours
	for( unsigned int i = 0; i < Contour_Lines.size(); ++i )
		Contour_Lines[i]->Draw("L SAME");

	//	Add the rest of the details to the plot
	addLHCbLabel( "LL & FC Overlay" )->Draw();
	leg->Draw();
	Overlay_Output->Update();

	TString Overlay_FileName( outputdir + "/LL_FC_Overlay");

	Overlay_Output->Print( Overlay_FileName + ".png" );
	Overlay_Output->Print( Overlay_FileName + ".pdf" );

	//      Return the plots back to their input state
	pllhist->SetContour(20);
	FC_Plot->SetContour(20);
}

void Plot_Styled_Contour( TH2* input_hist, int cont_num, double* input_conts, double* confs, TString outputdir, TString Name )
{
	vector<TString> Plot_Type;
	Plot_Type.push_back( "Temp" );
	Plot_Type.push_back( "Cont" );
	Plot_Type.push_back( "Conf" );
	Plot_Type.push_back( "Pub" );

	TString Draw_String;
	TString Canvas_Name("Styled_Canvas_");
	
	TString Base_Name = outputdir + "/" + Name + "_";

	for( unsigned int i = 0; i < Plot_Type.size(); ++i )
	{
	  
		if( i == 0 )	Draw_String = "colz";
		if( i == 1 )	Draw_String = "cont1z";

		if( i == 2 || i == 3 )
		{
			Draw_String = "cont LIST";
			input_hist->SetContour( cont_num, input_conts );
		}
		
		Canvas_Name+=i;
		TCanvas* Styled_Output_Canvas = new TCanvas( Canvas_Name, Canvas_Name, 1680, 1050 );

		input_hist->Draw( Draw_String );
		Styled_Output_Canvas->Update();

		if( i == 0 )
		{
			input_hist->Draw( Draw_String );
			Styled_Output_Canvas->Update();
		}

		if( i == 2 || i == 3 )
		{
			TObjArray *contObjArr = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

			TList* contLevel = NULL;
			TGraph* curv = NULL;
			TGraph* gc = NULL;
			double cl=0;
			TString confname;

			int TotalConts = contObjArr->GetSize();

			TLegend *leg = new TLegend(0.80,0.89,0.95,0.7);
			leg->SetHeader("Conf. Levels");
			leg->SetBorderSize(0);
			leg->SetFillStyle(0);

			input_hist->Draw("AXIS");

			for(int j = 0; j < TotalConts; ++j )
			{
				confname = "";
				cl = confs[j];
				confname +=cl;
				confname += "% C.L.";
				contLevel = (TList*)contObjArr->At(j);
				for(int k =0; k < contLevel->GetSize(); ++k)
				{
					curv = (TGraph*)contLevel->At(k);
					gc = (TGraph*)curv->Clone();
					if( Plot_Type[i] == "Pub"  )	gc->SetLineStyle( Style_t(j+1) );
					if( Plot_Type[i] == "Conf" )	gc->SetLineColor( Color_t(j+2) );
					gc->Draw("L");
				}
				leg->AddEntry( gc, confname, "L");
			}

			leg->Draw();
		}

		TString Output_Name = Base_Name + Plot_Type[i];

		addLHCbLabel( Name )->Draw();
		input_hist->SetTitle("");
		Styled_Output_Canvas->Update();

		Styled_Output_Canvas->Print( Output_Name + ".png" );
		Styled_Output_Canvas->Print( Output_Name + ".pdf" );

	}

	//	Lets conserve the efforts of this plotting tool :D
	cout << endl << "Writing TH2D object for comparisons with other scans" << endl << endl;

	TString Hist_FileName = outputdir+"/"+Name+".root";

        TFile * output = new TFile( Hist_FileName, "RECREATE" );

        input_hist->Write();
        output->Close();

	// Return to default
	input_hist->SetContour(20);
	return;
}

int main( int argc, char* argv[] )
{
	if( (argc != 5) && (argc != 6) ) exit(-1);
	//	Use UUID based seed from ROOT, just used for unique identification of ROOT objects
	TRandom3* rand_gen = new TRandom3(0);

	//	Setup the Canvas and such
	EdStyle* RapidFit_Style = new EdStyle();
	RapidFit_Style->SetStyle();

	//	By Design
	TString outputdir = argv[4];
	TString param1string = argv[2];
	TString param2string = argv[3];

	//	Make OutputDir
	gSystem->mkdir( outputdir );

	//	Open the File
	TFile * input = TFile::Open( argv[1] );
	TNtuple* allresults;
	input->GetObject("RapidFitResult", allresults);
	if(!allresults){
		input->GetObject("RapidFitResult/RapidFitResult",allresults);
		if(!allresults){
			cout << "Couldn't find ntuple RapidFitResult in TFile" << endl;
			exit(1);
		}
	}

	//	Switches for different plotting
	bool CV_Drift = false;
	bool Want_Physics_Params = false;
	bool want_FC = true;

	if( argc == 6 )
	{
		Want_Physics_Params = ( TString(argv[5]) == "CVPlots" );
	}

	//	Strings that are universally defined
	TString Double_Tolerance="0.000001";
	TString Fit_Status = "Fit_Status";
	TString NLL = "NLL";
	TString notgen = "-9999.";
	TString error_suffix = "_error";
	TString value_suffix = "_value";
	TString gen_suffix = "_gen";
	TString Copy_Option = "fast";


	//	Strings that are relevent to this plot
	//	This is currently given by the user, but in time we can write an algorithm to determine this
	TString param1_gen = param1string + gen_suffix;
	TString param1_val = param1string + value_suffix;
	TString param1_err = param1string + error_suffix;
	TString param2_gen = param2string + gen_suffix;
	TString param2_val = param2string + value_suffix;
	TString param2_err = param2string + error_suffix;


	//	Get a list of ALL parameters within the file, this is given for free in the algorithms I wrote above

	//	Get a list of all branches in allresults
	vector<TString> all_parameters = get_branch_names( allresults );
	//	Get a list of all branches in allresults with '_value' in their name
	vector<TString> all_parameter_values = filter_names( all_parameters, value_suffix );



	//	Now to read in the Global Minima

	//	Store the Global fit corrdinate and value
	Float_t nll=0, Global_Best_NLL=0, x_point=0, y_point=0;

	//	Store the value of the 'nuisence' parameters as well
	Float_t* best_fit_values = new Float_t[all_parameter_values.size()];
	//	Construct an array to hold the best values for the rest of the parameters
	Float_t* best_fit_temp_values = new Float_t[all_parameter_values.size()];
	allresults->SetBranchAddress(NLL,&nll);

	//	Tell ROOT where to store the values
	for( unsigned short int i=0; i < all_parameter_values.size(); ++i )
	{
		allresults->SetBranchAddress(all_parameter_values[i],&best_fit_temp_values[i]);
	}


	//	Actually store the values in resident memory of the program
	cout << endl << "BEST FIT RESULTS AS DETERMINED FROM THE STANDARD FIT:" << endl;
	allresults->GetEntry(0);
	Global_Best_NLL = nll;
	for( unsigned short int i=0; i < all_parameter_values.size(); ++i )
	{
		cout << all_parameter_values[i] << ":\t" << best_fit_temp_values[i] << endl;
		best_fit_values[i] = best_fit_temp_values[i];
		if( string(all_parameter_values[i].Data()).compare(param1_val.Data()) == 0 ) x_point = best_fit_temp_values[i];
		if( string(all_parameter_values[i].Data()).compare(param2_val.Data()) == 0 ) y_point = best_fit_temp_values[i];
	}
	cout << endl;



	//	General Cuts to be applied for various plots

	//	Fit_Status == 3
	TString Fit_Cut = "(abs(" + Fit_Status + "-3.0)<"+Double_Tolerance+")";


	//	param1_gen == notgen	&&	param1_err == 0
	TString Param_1_Cut = "(abs(" + param1_gen + "-" + notgen +")<" + Double_Tolerance + ")&&(abs("+ param1_err + "-0.0)<"+ Double_Tolerance + ")";
	//	param2_gen == notgen	&&	param2_err == 0
	TString Param_2_Cut = "(abs(" + param2_gen + "-" + notgen + ")<"+ Double_Tolerance + ")&&(abs(" +param2_err +"-0.0)<" + Double_Tolerance + ")";

	//	Combine the individual Cuts
	TString Fit_Cut_String = Param_1_Cut + "&&" + Param_2_Cut + "&&" + Fit_Cut;



	//	Toys have a defined generation Value, Fits to Data DO NOT

	//	param1_gen != notgen
	TString p1isatoy = "(abs(" + param1_gen + "-" + notgen + ")>" + Double_Tolerance + ")";
	//	param2_gen != notgen
	TString p2isatoy = "(abs(" + param2_gen + "-" + notgen + ")>" + Double_Tolerance + ")";

	//	Combine the individual Cuts
	TString Toy_Cut_String = p1isatoy + "&&" + p2isatoy;




	//	Check for Toys in the file and wether I should run the FC code

	allresults->Draw( NLL, Toy_Cut_String, "goff" );
	bool Has_Toys = allresults->GetSelectedRows() > 0;

	cout << endl << "NUMBER OF TOYS IN FILE:\t" << allresults->GetSelectedRows() << endl;




	//	Fit values for the global fit are now stored in:
	//
	//	Global_Best_NLL,	best_fit_values,	x_point,	y_point

	//	Tell the user
	cout << "GLOBAL DATA BEST FIT NLL:\t" << setprecision(10) << Global_Best_NLL << "\tAT:\tX:" << setprecision(10) << x_point << "\tY:\t" <<setprecision(10)<< y_point << endl;

	//	Check wether the minima as defined from the Global fit is the true minima within the phase-space

	Check_Minima( allresults, Fit_Cut_String, &Global_Best_NLL, NLL, Double_Tolerance, param1_val, param2_val );


	TString NLL_Min;		//	Of course ROOT doesn't have USEFUL constructors!
	NLL_Min+=Global_Best_NLL;

	//	Plot the distribution of successfully fitted grid points for the PLL scan
	//	NB:	For FC this will likely saturate due to multiple layers of fits
	cout << endl << "FOUND UNIQUE GRID POINTS, PLOTTING" << endl;
	LL2D_Grid( allresults, Fit_Cut, param1_val, param2_val, rand_gen, "LL" );

	//	Construct a plot string for the NLL plot and plot it

	//	PLL part of Draw_String
	TString PLL = "(" + NLL + "-" + NLL_Min + ")";
	//	Gridding part of Draw_String
	TString Param1_Param2 = ":" + param1_val + ":" + param2_val;


	//	Draw String for PLL
	TString NLL_DrawString = PLL + Param1_Param2;

	cout << endl << "PLOTTING NLL VARIATION" << endl;
	//	PLL plot
	TH2* pllhist = Plot_From_Cut( allresults, NLL_DrawString, Fit_Cut_String, rand_gen, param1string, param2string );


	//	Array storing the addresses of all of the Physics Plots still in memory
	TH2** All_Physics_Plots = new TH2*[all_parameter_values.size()];
	//	Making use of switch
	if( Want_Physics_Params )
	{
		//	Physics Plots
		Physics_Plots( all_parameter_values, best_fit_values, allresults, rand_gen, Param1_Param2, CV_Drift, All_Physics_Plots, Fit_Cut_String);
		cout << endl << "Finalising CV Plots" << endl << endl;
		Finalize_Physics_Plots( All_Physics_Plots, all_parameter_values, param1string, param2string, outputdir, CV_Drift );
	}



	TTree* FC_Output = new TTree( "FC_Output", "FC_Output" );;
	TH2* FC_Plot = NULL;
	//	Making use of switch
	if( Has_Toys && want_FC )
	{
		//	FC Plot
		cout << "FOUND TOYS IN FILE, PLOTTING FC" <<endl;
		FC_Plot = FC_TOYS( allresults, Fit_Cut_String, param1string, param2string, NLL, Fit_Cut, Global_Best_NLL, FC_Output, Double_Tolerance, rand_gen );

		//LL2D_Grid( FC_Output, Fit_Cut, param1_val, param2_val, rand_gen, "FC" );
	}

	//	Contours to be used in plotting
	int cont_num = 3;

	double* pllconts = new double[unsigned(cont_num)];
	pllconts[0] = 1.15; pllconts[1] = 2.36;
	pllconts[2] = 3.0;//  pllconts[3] = 4.61;

	double* fcconts = new double[unsigned(cont_num)];
	fcconts[0] = 0.68; fcconts[1] = 0.9;
	fcconts[2] = 0.95;// fcconts[3] = 0.99;

	double* confs = new double[unsigned(cont_num)];
	confs[0] = 68.0; confs[1] = 90.0;
	confs[2] = 95.0;// confs[3] = 99.0;

	cout <<endl<< "SAVING GRAPHS" << endl;

	Plot_Styled_Contour( pllhist, cont_num, pllconts, confs, outputdir, "Likelihood Profile" );

	if( Has_Toys && want_FC )	//	If FC was generated
	{
		Plot_Styled_Contour( FC_Plot, cont_num, fcconts, confs, outputdir, "FeldmanCousins Profile" );

		TString FC_Tuple_File( outputdir+ "FC_Tuple.root" );
		TFile* new_FC_Output = new TFile( FC_Tuple_File, "RECREATE" );
		FC_Output->Write();
		new_FC_Output->Close();

		Plot_Both( pllhist, FC_Plot, cont_num, fcconts, pllconts, confs, outputdir );
	}

	return 0;
}
