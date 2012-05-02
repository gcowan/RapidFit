/* stackergen: Part of the simpletools package
 * (c) Conor Fitzpatrick, 2008
 *
 * If you find this program useful in whole or in part 
 * please cite this paper: 
 *
 * Feel free to send bugreports, feature requests, patches etc to:
 * conor.fitzpatrick@cern.ch
 *
 */

#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCint.h"
#include "TMath.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "Riostream.h"
#include "TAxis.h"
#include "TRandom3.h"
#include "TArrow.h"
#include "TEntryList.h"
#include "TLegend.h"
#include "THStack.h"
#include "TColor.h"
#include "Rtypes.h"
#include "TObjArray.h"
#include "TKey.h"
#include "TVirtualHistPainter.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TPaveText.h"
#include "TMarker.h"
#include <list>
#include <cctype>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <stdio.h>

using std::cout;
using std::endl;

typedef vector<Float_t> array1f;
typedef vector<TGraph*> array1g;
typedef vector<TString> array1s;

bool unique_xy_coord( pair< pair<double,double>, vector<double> > first, pair< pair<double,double>, vector<double> > second )
{
	bool x_unique = (fabs(first.first.first - second.first.first) > 1E-9 );
	bool y_unique = (fabs(first.first.second - second.first.second) > 1E-9 );
	return !(x_unique || y_unique);
}

inline TString prettyPrint(Double_t value){
	char pretty[20];
	TString prettyString;
	sprintf (pretty, "%1.3g",value);
	prettyString = pretty;
	return prettyString;
}

//	Pass this the pointer to the ntuple you're handing and it will give you a list
//	of branches that it contains in an easy to handle manner
vector<TString> get_branch_names( TTree* local_tree )
{
	vector<TString> temp_branch_names;
	TObjArray* branch_obj_array = local_tree->GetListOfBranches();
	for(unsigned short int i=0; i<branch_obj_array->GetEntries();++i)
	{
		TObject* branch_object = (*branch_obj_array)[i];
		temp_branch_names.push_back((const char*) branch_object->GetName());
	}
	return temp_branch_names;
}

//	Pass this the name of the ntuple you're handing and it will give you a list
//	of branches that it contains in an easy to handle manner
vector<TString> get_branch_names( TString absolute_path )
{
	TTree* local_tree = (TTree*) gDirectory->Get((const char*) absolute_path);
	return get_branch_names( local_tree );
}

//	Pass this an array of strings and it will find all strings matching a substring
vector<TString> filter_names( vector<TString> all_names, string substring )
{
	vector<TString> returnable_names;
	size_t found=string::npos;
	for( unsigned short int i=0; i<all_names.size(); ++i )
	{	//	Again using STL functions :)
		string temp_str = (all_names[i].Data());
		found = temp_str.find( substring );
		if( found!=string::npos )	//	If the substring is found
		{
			returnable_names.push_back(all_names[i]);
		}
	}
	return returnable_names;
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
	TPaveText * label = new TPaveText(0.14, 0.73, 0.14, 0.88,"BRNDC");
	//TPaveText * label = new TPaveText(0.12, 0.58, 0.12, 0.43,"BRNDC");
	label->SetFillStyle(0);		//Transparent i.e. Opacity of 0 :D
	//	label->SetFillColor(0);
	label->SetBorderSize(0);
	label->SetTextAlign(11);
	label->SetTextSize(Float_t(0.04));
	TText * labeltext = 0;
	labeltext = label->AddText("LHC#font[12]{b} 2011 Data");
	labeltext = label->AddText("#sqrt{s} = 7TeV");
	labeltext = label->AddText(footer);
	return label;
}

TCanvas *makeTempCanvas(TGraph2D* graph, TString labelname){
	TCanvas *tempCanvas = new TCanvas(labelname+" temperature",labelname+" temperature",1024,768);
	tempCanvas->SetRightMargin(0.15);
	graph->Draw("colz");
	addLHCbLabel(labelname)->Draw();
	return tempCanvas;
}

TCanvas *makeTempCanvas(TH2* hist, TString labelname){
	TCanvas *tempCanvas = new TCanvas(labelname+" temperature",labelname+" temperature",1024,768);

	tempCanvas->SetRightMargin(0.15);
	hist->Draw("colz");
	addLHCbLabel(labelname)->Draw();
	return tempCanvas;
}
TCanvas *makePointCanvas(TGraph2D* graph, TString labelname){
	TCanvas *pointCanvas = new TCanvas(labelname+" points",labelname+" points",1024,768);
	pointCanvas->SetRightMargin(0.15);
	pointCanvas->SetPhi(0.0);
	pointCanvas->SetTheta(90.0);
	graph->Draw("PCOL");
	addLHCbLabel(labelname)->Draw();
	return pointCanvas;
}
TCanvas *makeContCanvas(TGraph2D* graph, TString labelname){
	TCanvas *contourCanvas = new TCanvas(labelname+" contour",labelname+" contour",1024,768);
	contourCanvas->SetRightMargin(0.15);
	graph->Draw("cont1z");
	addLHCbLabel(labelname)->Draw();
	return contourCanvas;
}

TCanvas *makeContCanvas(TH2* hist, TString labelname){
	TCanvas *contourCanvas = new TCanvas(labelname+" contour",labelname+" contour",1024,768);
	contourCanvas->SetRightMargin(0.15);
	hist->Draw("cont1z");
	addLHCbLabel(labelname)->Draw();
	return contourCanvas;
}

void scaleContours( TGraph * gc, double x_point, double y_point, double x_scale, double y_scale ) {
	double x = 0., y = 0.;
	double radius = 0.;
	double sinTheta = 0.;
	double cosTheta = 0.;
	double scale = 1.;
	// y_point is actually for the x-axis...
	double x_min = y_point;
	double y_min = x_point;
	// check which quadrant the global minimum is in	
	bool sm = false;
	bool flipped = false;
	if ( x_min > -1. && y_min > 0. ) sm = true;
	for ( int p = 0; p < gc->GetN(); p++ ) {
		gc->GetPoint( p, x, y );
		// need to flip the minimum depending upon which contour
		if ( !flipped && ( (sm && (x < -1. && y < 0.)) || (!sm && (x > -1 && y > 0.) )) ) {
			flipped = true;
			x_min = TMath::Pi() - x_min - 2*TMath::Pi(); // also move by 2pi
			y_min = -y_min;
		}
		radius = sqrt((x_min - x)*(x_min - x) + (y_min - y)*(y_min - y));	
		sinTheta = (y_min - y)/radius;
		cosTheta = (x_min - x)/radius;
		scale = sqrt( (x_scale * cosTheta)*(x_scale * cosTheta) + (y_scale * sinTheta)*(y_scale * sinTheta) );
		// set the new coordinate after scaling the radius and doing a bit of trig
		gc->SetPoint( p, -scale*radius*cosTheta + x_min, -scale*radius*sinTheta + y_min);
	}
}

TCanvas *makeConfCanvas(TH2* hist, TString labelname,UInt_t nconts, double* conts, double* confs, bool pub, double x_point, double y_point){
	//TH2* hist = (TH2*)_hist->Clone("confhist");
	TString pname = "color";
	if(pub){pname = "pub";}
	TCanvas *confCanvas = new TCanvas(pname + " " + labelname+" CL contours",pname + " " + labelname+" CL contours",1024,768);
	confCanvas->SetGrid();
	//confCanvas->SetRightMargin(0.15);
	hist->SetContour(nconts,conts);
	hist->Draw("cont LIST");
	confCanvas->Update();
	TObjArray *contObjArr = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TList* contLevel = NULL;
	TGraph* curv     = NULL;
	TGraph* gc    = NULL;
	double phi[1] = {-0.036};
	double phiE[1] = {0.002};
	double dg[1] = {0.087};
	double dgE[1] = {0.021}; 
	TGraphErrors* sm = new TGraphErrors( 1, phi, dg, phiE, dgE);
	sm->SetMarkerStyle(21);
	int TotalConts = contObjArr->GetSize();
	TLegend *leg = new TLegend(0.68,0.39,0.83,0.2);
	leg->SetHeader("Conf. Levels");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	hist->GetYaxis()->SetNdivisions(505);
	hist->GetXaxis()->SetTitle("#phi_{s} [rad]");
	hist->GetYaxis()->SetTitle("#Delta#Gamma_{s} [ps^{-1}]");
	hist->GetYaxis()->SetTitleOffset(1.2);
	hist->Draw("AXIS");
	for(int i = 0; i < TotalConts; ++i){
		TString confname = "";
		double cl = confs[i];
		confname +=cl;
		confname += "\% C.L.";
		contLevel = (TList*)contObjArr->At(i);
		for(int j = 0; j < contLevel->GetSize(); ++j){
			curv = (TGraph*)contLevel->At(j);
			gc = (TGraph*)curv->Clone();
			scaleContours( gc, x_point, y_point, 1.08, 1.011 ); // If you want to scale the contours to include systematics 
			if(pub){	
				gc->SetLineStyle(Style_t(i+1));
			}else{
				int c = 0;
				int s = 0;
				if ( i == 0 ) { c = 2; s = 1;} 
				if ( i == 1 ) { c = kViolet; s = 9;}
				if ( i == 2 ) { c = 4; s = 2;}
				gc->SetLineColor(Color_t(c));
				gc->SetLineStyle(s);
				gc->SetLineWidth(2.5);
			}
			gc->Draw("L");
		}
		leg->AddEntry(gc,confname, "L");
	}
	addLHCbLabel(labelname)->Draw();
	leg->Draw();
	sm->Draw("PLSAME");
	confCanvas->Update();
	//hist->Delete();
	hist->SetContour(20);
	return confCanvas;


}

TCanvas *makeBothCanvas(TH2* histfc, TH2* histll, TString labelname, UInt_t nconts, double* fcconts, double *llconts, double* confs){
	//TH2* hist = (TH2*)_hist->Clone("confhist");
	TString pname = "both";
	TCanvas *confCanvas = new TCanvas(pname + " " + labelname+" CL contours",pname + " " + labelname+" CL contours",1024,768);
	confCanvas->SetRightMargin(0.15);
	histll->SetContour(nconts,llconts);
	histll->Draw("cont LIST");
	confCanvas->Update();
	TObjArray *llcontObjArr = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TList* contLevel = NULL;
	TGraph* curv     = NULL;
	array1g gcarr;
	int llTotalConts = llcontObjArr->GetSize();


	TLegend *leg = new TLegend(0.7,0.89,0.85,0.7);
	leg->SetHeader("Conf. Levels");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	for(int i = 0; i < llTotalConts; ++i){
		TString confname = "";
		double cl = confs[i];
		confname +=cl;
		confname += "\% C.L. (PLL)";
		contLevel = (TList*)llcontObjArr->At(i);
		for(int j =0; j<contLevel->GetSize(); ++j){
			curv = (TGraph*)contLevel->At(j);
			TGraph *gc = (TGraph*)curv->Clone();
			gc->SetLineStyle(Style_t(i+1));
			gcarr.push_back(gc);

			if(j==0){leg->AddEntry(gc,confname, "L");}
		}

	}

	histfc->SetContour(nconts,fcconts);
	histfc->Draw("cont LIST");
	confCanvas->Update();
	TObjArray *fccontObjArr = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	confCanvas->Update();
	int fcTotalConts = fccontObjArr->GetSize();
	for(int i = 0; i < fcTotalConts; ++i){
		TString confname = "";
		double cl = confs[i];
		confname +=cl;
		confname += "\% C.L. (FC)";
		contLevel = (TList*)fccontObjArr->At(i);
		for(int j =0; j<contLevel->GetSize(); ++j){
			curv = (TGraph*)contLevel->At(j);
			TGraph *gc = (TGraph*)curv->Clone();
			gc->SetLineColor(Color_t(i+2));
			gcarr.push_back(gc);
			if(j==0){leg->AddEntry(gc,confname, "L");}
		}
	}
	histll->Draw("AXIS");
	for(UInt_t i =0; i<gcarr.size(); ++i){
		gcarr[i]->Draw("L SAME");
	}
	addLHCbLabel(labelname)->Draw();
	leg->Draw();
	confCanvas->Update();
	//hist->Delete();
	histll->SetContour(20);
	histfc->SetContour(20);
	return confCanvas;


}


int main(int argc, char *argv[]){

	if( (argc!=5) && (argc!=6) ){
		cout << "Plots 2D FCScans from rapidfit flat ntuples. Usage:" << endl;
		cout <<  argv[0] << " <flatntuple.root> <param1> <param2> <outputdir>"	<< endl;
		exit(1);
	}
	gStyle->SetPalette(1);
	gROOT->SetStyle("Plain");
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetStatColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetFillColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetLineWidth(2);
	bool plot_point=false;
	bool plot_CV=false;
	bool AbsCV=false;
	if( (argc>5) && ( string(argv[5])=="Points" ) ) { plot_point=true; }
	if( (argc>5) && ( string(argv[5])=="CV" ) ) { plot_CV=true; }
	if( (argc>5) && ( string(argv[5])=="AbsCV" ) ) { plot_CV=true; AbsCV=true; }
	if( (argc>5) && ( string(argv[5])=="CVPoints" ) ) { plot_point=true; plot_CV=true; }



	TString outputdir = argv[4];
	TString param1string = argv[2];
	TString param2string = argv[3];
	gSystem->mkdir( outputdir );
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
	vector<TString> all_parameters = get_branch_names( allresults );

	vector<TString> all_parameter_values = filter_names( all_parameters, "_value" );
	vector<TString> all_parameter_errors = filter_names( all_parameters, "_error" );
	all_parameter_values.insert(all_parameter_values.end(),all_parameter_errors.begin(),all_parameter_errors.end());

	TFile * output = new TFile(outputdir+"/llscanresults.root","RECREATE");

	Long64_t entries = allresults->GetEntries();
	cout << "INPUT NTUPLE CONTAINS: " << entries << " ENTRIES" << endl;

	TString param1valstr = param1string + "_value";
	TString param1genstr = param1string + "_gen";
	TString param1errstr = param1string + "_error";

	TString param2valstr = param2string + "_value";
	TString param2genstr = param2string + "_gen";
	TString param2errstr = param2string + "_error";

	TString FRstr = "Fit_Status";
	TString NLLstr = "NLL";


	//The first entry in every file will be the data result floated globally
	//Gen values will be -9999, 
	//Par values will be floated
	//Par errors !=0
	//The second entry in every file is the data result fixed at the current gridpoint. 
	//Gen values will be -9999, 
	//Par values will be set to the current gridpoint
	//Par errors will be 0
	//Entries 3-> Ntoyspersj+3 will be the fixed toy values
	//Gen values will be set to the current gridpoint
	//Par values will be set to the current gridpoint
	//Par errors will be 0
	//Entries Ntoyspersj+3->2*Ntoyspersj+3 will be the floated toy values
	//Gen values will be set to the current gridpoint
	//Par values will be floated
	//Par errors !=0

	Float_t shift = Float_t(0.000001);
	TString shiftstr = "";
	shiftstr += shift;

	TString notgen = "-9999.";
	Float_t nll, nlldatabest, x_point, y_point, global_x,global_y;
	array1f param1gridpoints, param2gridpoints, dataRatiogridpoints, clgridpoints, toysgridpoints;
	TTree* first_entry = allresults->CopyTree("","",1,0);
	first_entry->SetBranchAddress(NLLstr,&nll);
	first_entry->SetBranchAddress(param2valstr,&global_x);
	first_entry->SetBranchAddress(param1valstr,&global_y);
	Float_t best_fit_temp_values[all_parameter_values.size()];
	for( unsigned short int i=0; i < all_parameter_values.size(); ++i )
	{
		first_entry->SetBranchAddress(all_parameter_values[i],&best_fit_temp_values[i]);
	}
	// Find the global minimum from the first entry:

	//	THIS SHOULD BE THE ONLY USE OF GetEntry it's slow and doesn't emply something 'intelligent'
	//	nlldatabest = allresults->CopyTree("","",0)->GetMinimum("NLL");
	//	param1databest = allresults->CopyTree("","",0)->GetMinimum(param1valstr);
	//	param2databest = allresults->CopyTree("","",0)->GetMinimum(param2valstr);
	first_entry->GetEntry(0);
	nlldatabest = nll;
	Float_t best_fit_values[all_parameter_values.size()];
	for( unsigned short int i=0; i < all_parameter_values.size(); ++i )
	{
		best_fit_values[i] = best_fit_temp_values[i];
	}
	TString Best_Catch("NLL==");
	Best_Catch+=nlldatabest;
	TTree* local_best = allresults->CopyTree( Best_Catch );
	x_point = (Float_t)local_best->GetMinimum(param1valstr);
	y_point = (Float_t)local_best->GetMinimum(param2valstr);
	//	x_point = global_x;
	//	y_point = global_y;
	cout << "GLOBAL DATA MINIMUM NLL:\t" << setprecision(10)<< nlldatabest << "\tAT:\tX:" << x_point << "\tY:\t" << y_point << endl;

	// Find the positions of the gridpoints we scanned, and find the profileLL at each point for data:
	TString datafixedstr = param1genstr + "==" + notgen + "&&" + param2genstr + "==" + notgen + "&&"+ param1errstr + "==0.0&&" +param2errstr +"==0.0&&" + FRstr + "==3.0";
	//	TTree* datafixed = allresults->CopyTree(datafixedstr);


	//	40 x 40 x 200toys translates into 320000 seperate read statements...
	//	Use root's internal macros which are faster at processing root files esp for reading in bulk objects
	TString Plot_Str(param1valstr+":"+param2valstr+":"+NLLstr);
	TString PlotNLL_Str(NLLstr);
	allresults->SetEstimate(allresults->GetEntries());	//	Just incase you ever have more than 1E6 entries **
	int numberOfEventsAfterCut = int(allresults->Draw(Plot_Str,datafixedstr,"goff"));//	Draw without plotting
	//	That was 3 large read statements which is A LOT faster.
	Double_t* param1val_pointer = allresults->GetV1();	//	array of size	datafixed->GetEntries() **
	Double_t* param2val_pointer = allresults->GetV2();	//	array of size	datafixed->GetEntries()
	Double_t* nllval_pointer    = allresults->GetV3();	//	array of size	datafixed->GetEntries()
	//	param1gridpoints.reserve( numberOfEventsAfterCut );
	//	param2gridpoints.reserve( numberOfEventsAfterCut );
	//	dataRatiogridpoints.reserve( numberOfEventsAfterCut );
	//	Temp store all data in vectors

	list<pair<pair<double,double>,double> > searchable_list_data;

	for( short int i=0; i<numberOfEventsAfterCut; ++i )
	{
		param1gridpoints.push_back( float(param1val_pointer[i]) );
		param2gridpoints.push_back( float(param2val_pointer[i]) );
		dataRatiogridpoints.push_back( float(nllval_pointer[i]) );
		pair<pair<double,double>,double> new_coordinate_value;
		new_coordinate_value.first = make_pair( param1val_pointer[i], param2val_pointer[i] );
		new_coordinate_value.second = nllval_pointer[i];
		searchable_list_data.push_back( new_coordinate_value );
	}

	//------------------------------------------------------------------------------------------------------
	//	This section deals with fit minima and the duality in deltaGamma Phi_s

	double true_Z = allresults->CopyTree(datafixedstr)->GetMinimum("NLL");
	//	Move values into vectors
	
	if( (true_Z-nlldatabest) < 0 )
	{
		//	I know this is comparing doubles exactly, however this is also doing maths in Loki... you be the judge
		TString Catch(datafixedstr);
		Catch.Append("&&NLL==");
		Catch+=true_Z;
		TTree* local_best = allresults->CopyTree( Catch );
		double true_X = local_best->GetMinimum(param1valstr);
		double true_Y = local_best->GetMinimum(param2valstr);
		cout << "\n\t\tWARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING\n"<<endl;
		cout << "\t\tTRUE MINIMUM NLL = " << true_Z << "\t\tAt:\t\tX:" << true_X << "\tY:\t" << true_Y <<endl<<endl;
		//		cout << "\tNEW MINIMA FOUND AT :\tX:\t" << X_true_min << "\tY:\t" << Y_true_min << endl<<endl;
		cout << "\t\tWARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING\n"<<endl;
		nlldatabest=(Float_t)true_Z;
	}
	
	TString Condition(datafixedstr); Condition.Append("&&NLL==");	Condition+=nlldatabest;
	TTree* new_best = allresults->CopyTree( Condition );
	if( ( ( allresults->GetEntries() /2 ) - new_best->GetEntries()) == 0 )
	{
		TString Condition_Orig(datafixedstr);
		Condition_Orig.Append("&&(NLL!=");
		Condition_Orig+=true_Z;	Condition_Orig.Append(")");
		new_best = allresults->CopyTree( Condition_Orig );
		double new_Z = new_best->GetMinimum( "NLL" );
		TString Condition2(datafixedstr);
		Condition2.Append("NLL==");
		Condition2+=new_Z;
		new_best = allresults->CopyTree(Condition2);
		double new_X = new_best->GetMinimum(param1valstr);
		double new_Y = new_best->GetMinimum(param2valstr);
		cout << "SECOND MINIMA:\t\t" << setprecision(10) << new_Z << "\tFOUND AT:\tX:"<< new_X << "\tY:\t" <<new_Y<<endl;
		//	fabs doubles based on difference between them and best_X (deltaGamma)
		//	sort based on difference
	} else {
		cout << new_best->GetEntries() << "\tPOINTS FOUND WHICH CORRESPOND TO MINIMA"<<endl;
	}



	//------------------------------------------------------------------------------------------------------

	int numberOfFloatedPhysicsParams=0;
	vector<TString> Floated_Parameters;
	vector<vector<double> > Floated_values;
	//	Read in the CV of all floated parameters
	vector<vector<double> > all_CV_values;
	//	all_CV_values.reserve( all_parameter_values.size() );
	if( plot_CV )
	{
		//  Hold the Data in a temp object
		TString PlotString;

		for ( unsigned int obsIndex = 0; obsIndex < all_parameter_values.size(); obsIndex+=3 )
		{
			//	Un-comment the next 3 lines for Absolute Central Values in the plots
			if( AbsCV )
			{
				best_fit_values[obsIndex] = 0.0;
				if((obsIndex+1)<=all_parameter_values.size())best_fit_values[obsIndex+1] = 0.0;
				if((obsIndex+2)<=all_parameter_values.size())best_fit_values[obsIndex+2] = 0.0;
			}

			//			data_array.reserve(3);
			vector<Double_t *> data_array;
			//  Construct a Plot String to use the TTree->Draw Method
			PlotString="(";
			PlotString.Append(all_parameter_values[obsIndex]);
			PlotString.Append("-");
			PlotString+=best_fit_values[obsIndex];
			PlotString.Append(")");
			for( unsigned short int i=1; ((obsIndex+i)<all_parameter_values.size())&&((i<3)); ++i )
			{
				PlotString.Append(":(");
				PlotString.Append(all_parameter_values[obsIndex+i]);
				PlotString.Append("-");
				PlotString+=best_fit_values[obsIndex+i];
				PlotString.Append(")");
			}
			//  Draw 3 observables at a time in some large plot
			//  use the 'goff' option to turn graphical output (and annoying text output from con/de-structors) off
			//  (it doesn't matter what this looks like and we can throw it away from here :)
			int num2 = (int)allresults->Draw( PlotString , datafixedstr, "goff" );

			//  Store pointers to the objects for ease of access
			TH1* temp_hist = (TH1*)allresults->GetHistogram()->Clone();
			TString Unique("");
			Unique+=obsIndex;
			temp_hist->SetName(Unique);
			//			temp_hist->SetTitle(Unique);

			bool flag1=true, flag2=true, flag3=true;


			//			double X_min_range = temp_hist->GetXaxis()->GetXmin();
			//			double X_max_range = temp_hist->GetXaxis()->GetXmax();
			//			double Y_min_range = temp_hist->GetYaxis()->GetXmin();
			//			double Y_max_range = temp_hist->GetYaxis()->GetXmax();
			//			double Z_min_range = temp_hist->GetZaxis()->GetXmin();
			//			double Z_max_range = temp_hist->GetZaxis()->GetXmax();
			//			double bin_width_X = temp_hist->GetXaxis()->GetBinWidth(0);
			//			double bin_width_Y = temp_hist->GetYaxis()->GetBinWidth(0);
			//			double bin_width_Z = temp_hist->GetZaxis()->GetBinWidth(0);
			//			if( fabs(fabs((X_min_range + X_max_range) - X_max_range*2.)/X_max_range-2) > 1E-3 ) flag1=true;
			//			if( fabs(fabs((Y_min_range + Y_max_range) - Y_max_range*2.)/Y_max_range-2) > 1E-3 ) flag2=true;
			//			if( fabs(fabs((Z_min_range + Z_max_range) - Z_max_range*2.)/Z_max_range-2) > 1E-3 ) flag3=true;
			//			flag1=flag2=flag3=true;
			//			if( (num1==num1) && num1>1E-9 )	flag1=true;
			//			if( (num2==num2) && num2>1E-9 )	flag2=true;
			//			if( (num3==num3) && num3>1E-9 )	flag3=true;


			if( flag1 ){	data_array.push_back( allresults->GetV1() ); ++numberOfFloatedPhysicsParams; Floated_Parameters.push_back(all_parameter_values[obsIndex]);}
			//  GetV1 returns a pointer to an array of doubles for the first ('X') axis in the plot we've just drawn
			if( ((obsIndex+1) < all_parameter_values.size()) && (flag2) ) {  ++numberOfFloatedPhysicsParams; data_array.push_back( allresults->GetV2() ); Floated_Parameters.push_back(all_parameter_values[obsIndex+1]);}
			if( ((obsIndex+2) < all_parameter_values.size()) && (flag3) ) {  ++numberOfFloatedPhysicsParams; data_array.push_back( allresults->GetV3() ); Floated_Parameters.push_back(all_parameter_values[obsIndex+2]);}

			//  As You can only plot 3 at a time using this mechanism and Root will overide the results between drawing the plots
			//  Save the actual data to somewhere protected in our memory
			for(unsigned short int i=0; i < data_array.size(); ++i )
			{
				vector<double> temp_vector;
				temp_vector.reserve( num2 );
				for(int j=0; j < num2; ++j )
				{
					temp_vector.push_back( Double_t(data_array[i][j]) );	//	ith array and jth event in
				}
				all_CV_values.push_back( temp_vector );
			}
		}
		//		Floated_values.resize( numberOfFloatedPhysicsParams );
		cout << "Number of Additional Floated Physics Parameters: " << numberOfFloatedPhysicsParams <<endl;
	}

	list<pair<pair<double, double>, vector<double> > >  Coordinates_NLL_CV;
	//Coordinates_NLL_CV.reserve(numberOfEventsAfterCut);		//	Doesn't apply to lists

	cout << param1gridpoints.size() <<endl;
	for(unsigned short int i = 0; i < param1gridpoints.size(); ++i)
	{
		//	Store the coordinate we are at
		pair<double,double> Coords;
		Coords.first=param1gridpoints[i];
		Coords.second=param2gridpoints[i];
		//	Store all interesting information from the fit
		vector<double> NLL_CV_values;//( all_CV_values.size()+1 );		//	Store all CV values+NLL
		NLL_CV_values.push_back( dataRatiogridpoints[i] );		//	Always store NLL as first value
		for( unsigned short int j=0; j < all_CV_values.size() ; ++j )
		{
			NLL_CV_values.push_back( all_CV_values[j][i] );	//	Store CV at point i for param j
		}

		//	Construct an object which holds all of this data at this point
		pair<pair<double,double>,vector<double> > Coordinate_Element;
		Coordinate_Element.first = Coords;
		Coordinate_Element.second = NLL_CV_values;

		//	Store all of the information at this point
		Coordinates_NLL_CV.push_back(Coordinate_Element);
	}

	cout << "FINDING UNIQUE COORDINATES ONLY" <<endl;
	//	Using the STL syntax for finding unique objects in a list seems nicer
	Coordinates_NLL_CV.unique( unique_xy_coord );
	int numberOfUniqueEventsAfterCut = (int)Coordinates_NLL_CV.size();
	cout << "FOUND " << numberOfUniqueEventsAfterCut << " UNIQUE GRIDPOINTS" << endl;

	//	Some cleanup so we can continue
	while( !param1gridpoints.empty() )
	{
		param1gridpoints.pop_back();
		param2gridpoints.pop_back();
		dataRatiogridpoints.pop_back();
	}

	//	Store all values corresponding to unique data points in X and Y
	list<pair<pair<double,double>,vector<double> > >::iterator coord_iter=Coordinates_NLL_CV.begin();
	Floated_values.resize( all_CV_values.size() );
	for( unsigned short int i=0; i<Floated_values.size(); ++i ){
		Floated_values[i].resize( numberOfUniqueEventsAfterCut ); }

	for( short int i=0; coord_iter != Coordinates_NLL_CV.end(); ++coord_iter, ++i )
	{
		param1gridpoints.push_back( float(coord_iter->first.first) );
		param2gridpoints.push_back( float(coord_iter->first.second) );
		dataRatiogridpoints.push_back( float(coord_iter->second[0] - nlldatabest) );
		for( unsigned short int j=0; j<Floated_values.size(); ++j )
		{
			Floated_values[j][i] = coord_iter->second[j+1];
		}
	}

	bool hastoys = false;

	TString p1isatoy = "(abs("+param1genstr+"-"+notgen+")>"+shiftstr+")";
	TString p2isatoy = "(abs("+param2genstr+"-"+notgen+")>"+shiftstr+")";
	TString isatoy = p1isatoy + "&&" + p2isatoy;

	TString all_toys_cut = datafixedstr+"&&"+isatoy;

	//	TTree* toys = datafixed->CopyTree(isatoy);
	int numberOfToyEventsAfterCut = int(allresults->Draw(Plot_Str,all_toys_cut,"goff"));//	Draw without plotting
	cout << "FOUND " << numberOfToyEventsAfterCut << " TOYS" << endl;

	//Is this just an LLscan, or were toys for FC generated? 
	if(numberOfToyEventsAfterCut > 0){hastoys = true;}

	if(hastoys){
		//Find the toys generated at the previously discovered gridpoints
		for(int i = 0; i<numberOfEventsAfterCut; ++i){
			//cout << "step " << i << " of " <<  param1gridpoints.size() << endl;

			// Build up cutstrings for the toys
			TString toyscutstr = FRstr + "==3.0&&(abs(" +param1genstr + "-";
			toyscutstr += param1gridpoints[i];
			toyscutstr += ")<"+shiftstr+")&&(abs(" + param2genstr + "-";
			toyscutstr += param2gridpoints[i];
			toyscutstr += ")<"+shiftstr+")&&";

			//Reduce the ntuple to only toys generated at this gridpoint, and only toys that had floated params
			TString floatedtoyscutstr = toyscutstr;
			floatedtoyscutstr += "(abs("+param1errstr +")>"+shiftstr + ")&&(abs(" + +param1errstr +")>"+shiftstr + ")";


			//	Extract the Floated Toys from the complete dataset
			TString Floated_Toys=all_toys_cut+"&&"+floatedtoyscutstr;
			int numberOfFloatedToys = int(allresults->Draw(PlotNLL_Str,Floated_Toys,"goff"));
			Double_t* NLLfloatedtoy_pointer = allresults->GetV1();
			vector<double> NLLtoyfloat(numberOfFloatedToys);
			for( unsigned short int j=0; j < numberOfFloatedToys; ++j )
			{
				NLLtoyfloat.push_back( NLLfloatedtoy_pointer[j] );
			}


			//Reduce the ntuple to only toys generated at this gridpoint, and only toys that had fixed params
			TString fixedtoyscutstr = toyscutstr;
			fixedtoyscutstr += "(abs("+param1errstr +")<"+shiftstr + ")&&(abs(" + +param1errstr +")<"+shiftstr + ")";

			//Do the same for toys that had fixed params
			TString Fixed_Toys=all_toys_cut+"&&"+fixedtoyscutstr;
			int numberOfFixedToys = int(allresults->Draw(PlotNLL_Str,Fixed_Toys,"goff"));
			Double_t* NLLfixedtoy_pointer = allresults->GetV1();
			vector<double> NLLtoyfixed(numberOfFixedToys);
			for( unsigned short int j=0; j < numberOfFixedToys; ++j )
			{
				NLLtoyfixed.push_back( NLLfixedtoy_pointer[j] );
			}



			if( numberOfFixedToys != numberOfFloatedToys ) {
				//We've got serious problems! 
				cout << "Different number of toy fits found for gridpoint " << param1gridpoints[i] << "," << param2gridpoints[i] << " (" << numberOfFixedToys << "," << numberOfFloatedToys << ") Can't continue!" << endl;
				exit(1);
			}

			if( numberOfFixedToys != 0 ){

				//Loop over the toys, pulling out the NLL ratio
				UInt_t toyNLLsmaller = 0;
				toysgridpoints.push_back((Float_t)numberOfFixedToys);
				for(unsigned short int j = 0; j<numberOfFixedToys; ++j){
					//THE LINE BELOW IS THE FELDMAN-COUSINS ORDERING METHOD USED BY CDF/HEIDELBERG: 
					//if the toyratio is smaller than the data ratio at this point, increment:
					if((NLLtoyfixed[j]-NLLtoyfloat[j])<dataRatiogridpoints[i]){++toyNLLsmaller;}
				}
				//The C.L. is the percentage of toys that were smaller
				double cl = ((Double_t)toyNLLsmaller/(Double_t)numberOfFixedToys);
				clgridpoints.push_back(float(cl));
			}else{
				cout << param1gridpoints[i] << " " << param2gridpoints[i] << " WARNING: Found no toys here. " << endl;
				clgridpoints.push_back(-9999.0);
			}

		}
	}
	//We now have 4 vectors: The gridpoints in x,y and the conf. limits in z or the profile LL in z. We need to make TGraphs from these bad boys. 

	//Copy vectors to arrays: 
	int npoints = int(param1gridpoints.size());
	Double_t* p2points = new Double_t [npoints];
	copy( param2gridpoints.begin(), param2gridpoints.end(),p2points);
	Double_t* p1points = new Double_t [npoints];
	copy( param1gridpoints.begin(), param1gridpoints.end(),p1points);
	Double_t* pllpoints = new Double_t [npoints];
	copy( dataRatiogridpoints.begin(), dataRatiogridpoints.end(),pllpoints);

	int np = (int)sqrt((double)npoints);
	//We get the data profile likelihood free, so let's plot it: 
	TGraph2D *pllgraph = new TGraph2D(npoints, p2points, p1points, pllpoints);
	pllgraph->SetName("pllgraph");
	pllgraph->SetNpx(np);
	pllgraph->SetNpy(np);

	TH2D* pllhist = pllgraph->GetHistogram();
	pllhist->SetName("pllhist");
	pllhist->SetTitle("");
	pllhist->GetXaxis()->SetTitle(param2string);
	pllhist->GetYaxis()->SetTitle(param1string);
	pllhist->Write();

	if( plot_CV )
	{
		cout << "CONSTRUCTING CV PLOTS, NB For a large number of events, this is CPU intensive!\n"<<endl;
		for( unsigned short int i=0; i<Floated_Parameters.size(); ++i )
		{
			cout << "Constructing CV Plot " << i+1 << " of " << Floated_Parameters.size() <<endl;
			Double_t* CV_values = new Double_t [Floated_values[i].size()];
			copy( Floated_values[i].begin(), Floated_values[i].end(), CV_values );
			cerr << Floated_Parameters.size() << " " << Floated_values[i].size() << " " << npoints << endl;
			TGraph2D* Local_Graph = new TGraph2D((int)Floated_values[i].size(), p2points, p1points, CV_values);
			Local_Graph->SetName(Floated_Parameters[i]);
			Local_Graph->SetNpx(np);
			Local_Graph->SetNpy(np);
			TH2D* Local_Hist = Local_Graph->GetHistogram();
			Local_Hist->SetTitle("");
			Local_Hist->GetXaxis()->SetTitle(param2string);
			Local_Hist->GetYaxis()->SetTitle(param1string);
			TCanvas *Local_Canvas = makeTempCanvas(Local_Hist,Floated_Parameters[i]);
			Local_Canvas->Print(outputdir+"/"+param1string+"_"+param2string+"_"+Floated_Parameters[i]+".png");
			Local_Canvas->Print(outputdir+"/"+param1string+"_"+param2string+"_"+Floated_Parameters[i]+".pdf");
			cout << "Finished Constructing This CV Plot" <<endl;
		}
	}


	double pllconts[3] = {1.15,2.31,3.0};
	double fcconts[3] = {0.68,0.9,0.95};
	double confs[3] = {68.0,90.0,95.0};

	TMarker* global_minima = new TMarker( x_point, y_point, 20 );
	global_minima->SetMarkerSize(1);
	global_minima->SetMarkerColor(1);


	TCanvas *pllPoints = makePointCanvas(pllgraph, "PLL Gridpoints");

	pllPoints->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_points.pdf");
	pllPoints->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_points.png");

	TCanvas *plltemp = makeTempCanvas(pllhist,"Profile Likelihood");
	if( plot_point )  {
		global_minima->Draw("SAME");
		plltemp->Update();
	}
	plltemp->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_temp.pdf");
	plltemp->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_temp.png");

	TCanvas *pllcont = makeContCanvas(pllhist,"Profile Likelihood");
	if( plot_point )  {
		global_minima->Draw("SAME");
		pllcont->Update();
	}
	pllcont->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_cont.pdf");
	pllcont->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_cont.png");
	TCanvas *pllconf = makeConfCanvas(pllhist,"Profile Likelihood",3,pllconts,confs,false, x_point, y_point);
	if( plot_point )  {
		global_minima->Draw("SAME");
		pllconf->Update();
	}
	pllconf->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_conf.pdf");
	pllconf->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_conf.png");
	TCanvas *pllpub = makeConfCanvas(pllhist,"Profile Likelihood",3,pllconts,confs,true, x_point, y_point);
	if( plot_point )  {
		global_minima->Draw("SAME");
		pllcont->Update();
	}
	pllpub->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_pub.pdf");
	pllpub->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_pub.png");

	// Now let's plot the feldman cousins in the same way: 
	if(hastoys){
		Double_t* clpoints = new Double_t [npoints];
		copy( clgridpoints.begin(), clgridpoints.end(),clpoints);
		Double_t* toyspoints = new Double_t [npoints];
		copy( toysgridpoints.begin(), toysgridpoints.end(),toyspoints);

		TGraph2D *fcgraph = new TGraph2D(npoints, p2points, p1points, clpoints);
		fcgraph->SetName("fcgraph");
		fcgraph->SetNpx(np);
		fcgraph->SetNpy(np);
		TH2D* fchist = fcgraph->GetHistogram();
		fchist->SetName("fchist");
		fchist->SetTitle("");
		fchist->GetXaxis()->SetTitle(param2string);
		fchist->GetYaxis()->SetTitle(param1string);
		fchist->Write();

		TCanvas *fctemp = makeTempCanvas(fchist,"Feldman Cousins");
		fctemp->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_temp.pdf");
		fctemp->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_temp.png");
		TCanvas *fccont = makeContCanvas(fchist,"Feldman Cousins");
		fccont->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_cont.pdf");
		fccont->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_cont.png");
		TCanvas *fcconf = makeConfCanvas(fchist,"Feldman Cousins",3,fcconts,confs,false, x_point, y_point);
		fcconf->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_conf.pdf");
		fcconf->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_conf.png");
		TCanvas *fcpub = makeConfCanvas(fchist,"Feldman Cousins",3,fcconts,confs,true, x_point, y_point);
		fcpub->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_pub.pdf");
		fcpub->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_pub.png");
		TCanvas *both = makeBothCanvas(fchist,pllhist,"FC and PLL",3,fcconts,pllconts,confs);
		both->Print(outputdir+"/"+param1string+"_"+param2string+"_both.pdf");
		both->Print(outputdir+"/"+param1string+"_"+param2string+"_both.png");

		TGraph2D *toygraph = new TGraph2D(npoints, p2points, p1points, toyspoints);
		toygraph->SetName("ntoysgraph");
		toygraph->SetNpx(np);
		toygraph->SetNpy(np);
		TH2D* toyhist = toygraph->GetHistogram();
		toyhist->SetName("ntoyshist");
		toyhist->SetTitle("");
		toyhist->GetXaxis()->SetTitle(param2string);
		toyhist->GetYaxis()->SetTitle(param1string);
		toyhist->Write();

		TCanvas *toyPoints = makePointCanvas(toygraph, "FC Gridpoints");
		toyPoints->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_points.pdf");
		toyPoints->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_points.png");

		TCanvas *toytemp = makeTempCanvas(toyhist,"Toys per gridpoint");
		toytemp->Print(outputdir+"/"+param1string+"_"+param2string+"_ntoys_temp.pdf");
		toytemp->Print(outputdir+"/"+param1string+"_"+param2string+"_ntoys_temp.png");


		output->Write();
		output->Close();
	}

}

