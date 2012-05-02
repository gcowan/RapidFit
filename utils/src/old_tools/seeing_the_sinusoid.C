
#include "TString.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

#include "Histo_Processing.h"
#include "NTuple_Processing.h"

#include "EdStyle.h"

#include <vector>
#include <string>
#include <stdio.h>
#include <ctype.h>

using namespace::std;

const TString format_str = "PE9";

struct FitStruct{
	explicit FitStruct() : fitFunction(NULL), min(0.), max(0.), fitOption("WL"), fitDrawOption("") {}

	TF1* fitFunction;
	double min;
	double max;
	TString fitOption;
	TString fitDrawOption;

	private:
	FitStruct(const FitStruct&);
	FitStruct& operator=(const FitStruct&);
};

void DrawWrite( TH1* histo_p, TString title, TString name )
{
	histo_p->SetTitle( title );
	histo_p->SetName( name );
	histo_p->Draw( format_str );
	histo_p->Write();
	return;
}

TH1* Generate_Histo( TTree* decay_tree, TString DrawStr, TString WeightStr, TRandom3* rand, TString title, TString name, FitStruct* fit_struct, bool rebin=false, int rebinVal=-1 )
{
	TH1* histo_p = Get_Histo( decay_tree, DrawStr, WeightStr, rand );
	if( rebin==true )
	{
		if( rebinVal = -1 ) OptimallyRebin( histo_p );
		else histo_p->Rebin( rebinVal );
	}
	if( fit_struct != NULL ) histo_p->Fit( fit_struct->fitFunction, fit_struct->fitOption, fit_struct->fitDrawOption, fit_struct->min, fit_struct->max );
	DrawWrite( histo_p, title, name );
	return histo_p;
}

TH1* Generate_Histo( TTree* decay_tree, TString DrawStr, TString WeightStr, TRandom3* rand, TString title, TString name, bool rebin, int rebinVal )
{
	return Generate_Histo( decay_tree, DrawStr, WeightStr, rand, title, name, NULL, rebin, rebinVal );
}

TH1* Generate_Histo( TTree* decay_tree, TString DrawStr, TString WeightStr, TRandom3* rand, TString title, TString name )
{
	return Generate_Histo( decay_tree, DrawStr, WeightStr, rand, title, name, NULL );
}

TString Sanitize_Text( TString input )
{
	TString output;
	string internal( input.Data() );
	for( unsigned int i=0; i< internal.size(); ++i )
	{
		if( isalnum( internal.at( i ) ) )
		{
			output.Append( internal.at(i) );
		}
		else
		{
			output.Append("_");
		}
	}
	return output;
}

vector<TH1*> Plot_Tagged_Combinations( TTree* decayTree, TString X_var, TString selection_cuts, TString tag_str, TString weight_0, TString weight_1, TString weight_2, TRandom3* rand,
		TString Weight0_Name, TString Weight1_Name, TString Weight2_Name, bool rebin=false, int rebinVal=-1, FitStruct* func_form1=NULL, FitStruct* func_form2=NULL, bool smooth=false )
{
	TCanvas* c2 = new TCanvas( "new", "new", 1680, 1050);
	c2->Divide(3,3);
	TString rand_str="";rand_str+=rand->Rndm();
	c2->SetName(rand_str);

	int framenum=0;

	TString Weight0_Name_ROOT = Sanitize_Text( Weight0_Name );
	TString Weight1_Name_ROOT = Sanitize_Text( Weight1_Name );
	TString Weight2_Name_ROOT = Sanitize_Text( Weight2_Name );

	vector<TH1*> histos;
	TH1* temp = NULL;

	c2->cd(++framenum);
	temp = Generate_Histo( decayTree, X_var, selection_cuts+"*("+tag_str+"==1)*"+weight_0, rand, Weight0_Name + " tag==1", Weight0_Name_ROOT+"_tag1", func_form1, rebin, rebinVal );
	if( smooth ) temp->Smooth();
	histos.push_back( temp );

	c2->cd(++framenum);
	temp = Generate_Histo( decayTree, X_var, selection_cuts+"*("+tag_str+"==1)*"+weight_1, rand, Weight1_Name + " tag==1", Weight1_Name_ROOT+"_tag1", func_form1, rebin, rebinVal );
	if( smooth ) temp->Smooth();
	histos.push_back( temp );

	c2->cd(++framenum);
	temp = Generate_Histo( decayTree, X_var, selection_cuts+"*("+tag_str+"==1)*"+weight_2, rand, Weight2_Name + " tag==1", Weight2_Name_ROOT+"_tag1", func_form1, rebin, rebinVal );
	if( smooth ) temp->Smooth();
	histos.push_back( temp );

	c2->cd(++framenum);
	temp = Generate_Histo( decayTree, X_var, selection_cuts+"*(abs("+tag_str+"+1)<1E-5)*"+weight_0, rand, Weight0_Name + " tag==-1", Weight0_Name_ROOT+"_tag-1", func_form1, rebin, rebinVal );
	if( smooth ) temp->Smooth();
	histos.push_back( temp );

	c2->cd(++framenum);
	temp = Generate_Histo( decayTree, X_var, selection_cuts+"*(abs("+tag_str+"+1)<1E-5)*"+weight_1, rand, Weight1_Name + " tag==-1", Weight1_Name_ROOT+"_tag-1", func_form1, rebin, rebinVal );
	if( smooth ) temp->Smooth();
	histos.push_back( temp );

	c2->cd(++framenum);
	temp = Generate_Histo( decayTree, X_var, selection_cuts+"*(abs("+tag_str+"+1)<1E-5)*"+weight_2, rand, Weight2_Name + " tag==-1", Weight2_Name_ROOT+"_tag-1", func_form1, rebin, rebinVal );
	if( smooth ) temp->Smooth();
	histos.push_back( temp );

	c2->cd(++framenum);
	temp = Generate_Histo( decayTree, X_var, selection_cuts+"*(("+tag_str+"==1)-("+tag_str+"==-1))*"+weight_0, rand, Weight0_Name + " tag==1 - tag==-1", Weight0_Name_ROOT+"_diff", func_form2, rebin, rebinVal );
	if( smooth ) temp->Smooth();
	histos.push_back( temp );

	c2->cd(++framenum);
	temp = Generate_Histo( decayTree, X_var, selection_cuts+"*(("+tag_str+"==1)-("+tag_str+"==-1))*"+weight_1, rand, Weight1_Name + " tag==1 - tag==-1", Weight1_Name_ROOT+"_diff", func_form2, rebin, rebinVal );
	if( smooth ) temp->Smooth();
	histos.push_back( temp );

	c2->cd(++framenum);
	temp = Generate_Histo( decayTree, X_var, selection_cuts+"*(("+tag_str+"==1)-("+tag_str+"==-1))*"+weight_2, rand, Weight2_Name + " tag==1 - tag==-1", Weight2_Name_ROOT+"_diff", func_form2, rebin, rebinVal );
	if( smooth ) temp->Smooth();
	histos.push_back( temp );

	c2->Print("output_"+rand_str+".png");

	return histos;
}


int main( int argc, char* argv[] )
{
	(void) argc; (void) argv;
	EdStyle* myStyle = new EdStyle();
	myStyle->SetStyle();
	gStyle->SetOptStat(0);
	TRandom3* rand = new TRandom3();
	TString rand_str;

	double dms = 17.63;//17.716;
	double period = 2.0*3.14159/dms;
	double lower_time = 0.5;//1.0;
	double upper_time = 8.0;

	TString deltaM_str; deltaM_str+=dms;
	TString lower_time_str; lower_time_str+=lower_time;
	TString upper_time_str; upper_time_str+=upper_time;
	TString period_str; period_str+=period;

	bool realData = false;

	string fileName = argv[1];

	//      Open a file ane get the requested TTree from the file
	TTree* decayTree = GetFirstTree( fileName );

	if( decayTree->GetName() == "DecayTree" )
	{
		realData = true;
	}


	TString time_str = "time";
	TString tag_str = "tag";
	TString mistag_str = "mistag";

	if( realData )
	{
		tag_str+="_os";
	}


	TString time_cuts = "(time>=" +lower_time_str + ")*(time<=" + upper_time_str + ")";

	TString analysis_cut = "(sel==1)*(triggerDecision==1)*(12>abs(mdau2-1020))";

	TString tagged_cut = "("+tag_str+"*"+tag_str+"==1)";

	TString mistag_cut = "("+mistag_str+"<=0.5)";

	TString mass_cut = "(mass>5200)*(mass<5550)*((mdau2-1020)*(mdau2-1020)<400)";

	TString tac_cut = "(timeAcceptanceCategory<2)*(timeAcceptanceCategory>-1)";

	TString selection_cuts = time_cuts + "*" + tagged_cut;

	if( realData )
	{
		selection_cuts += "*" + mistag_cut + "*" + analysis_cut + "*" + tac_cut + "*" + mass_cut;
	}

	//	From 1-angle analysis
	TString weight_1_1angle = "(5.*cosTheta*cosTheta - 1.0)";		// CP even, +sin terms
	TString weight_2_1angle = "(2.- 5.*cosTheta*cosTheta )";		// CP odd, -sin terms

	TString even_weight_1angle = "(("+tag_str+"==+1)*" + weight_1_1angle + "-" + "("+tag_str+"==-1)*" + weight_1_1angle + ")";
	TString odd_weight_1angle =  "(("+tag_str+"==+1)*" + weight_2_1angle + "-" + "("+tag_str+"==-1)*" + weight_2_1angle + ")";


	//	From full-angular analysis
	TString A0_weight = "((5*cosTheta*cosTheta-sin(acos(cosTheta))*sin(acos(cosTheta))*cos(2.*acos(cosPsi)))-1)/2.";
	TString Apara_weight = "((5*cosTheta*cosTheta+sin(acos(cosTheta))*sin(acos(cosTheta))*cos(2.*acos(cosPsi)))-1)/2.";
	TString Aperp_weight = "(2-5*cosTheta*cosTheta)";

	TString Re_int = "(5./sqrt(2.))*(sin(2*acos(cosPsi))*sin(2*acos(cosPsi)))";
	TString Im_A0Ap = "(-5./2.)*(sin(2*acos(cosTheta))*sin(phi))";
	TString Im_ApAp = "(25./4.)/sqrt(2.)*(sin(2*acos(cosPsi))*sin(2*acos(cosTheta))*cos(phi))";

	//TString weight_1 = "((" + A0_weight + ")+(" + Apara_weight + ")+(" + Re_int +"))";		//	+sin terms
	//TString weight_2 = "((" + Aperp_weight + ")+(" + Im_A0Ap + ")+(" + Im_ApAp + "))";		//	-sin terms

	//TString weight_1 = "(" + Re_int +")";
	//TString weight_2 = "((" + Im_A0Ap + ")+(" + Im_ApAp + "))";

	TString weight_1 = "((" + A0_weight + ")+(" + Apara_weight + "))";
	TString weight_2 = "(" + Aperp_weight + ")";


	TString even_weight = "(("+tag_str+"==+1)*" + weight_1 + "-" + "("+tag_str+"==-1)*" + weight_1 + ")";
	TString odd_weight =  "(("+tag_str+"==+1)*" + weight_2 + "-" + "("+tag_str+"==-1)*" + weight_2 + ")";

	TString full_weight;
	TH1* histo_p=NULL;

	if( realData )
	{
		TString background_weight = "(1. - 0.73*(timeAcceptanceCategory==0)*exp(-0.145*time) - 0.62*(timeAcceptanceCategory==1)*exp(-0.366*time))";

		even_weight += "*" + background_weight;
		odd_weight += "*" + background_weight;

		even_weight_1angle += "*" + background_weight;
		odd_weight_1angle += "*" + background_weight;
	}

	TF1 *_osc = new TF1( "decay_osc", "[0]*sin([1]*x+[2])", lower_time, upper_time );
	_osc->SetParameter(0,50.0 );
	_osc->SetParName(0,"Amplitude");
	_osc->SetParameter(1,2*3.14159);
	_osc->SetParName(1,"Period");
	_osc->SetParameter(2,0);
	_osc->SetParName(2,"Offset");

        TF1 *_osc2 = new TF1( "flat_osc", "[0]*sin([1]*x+[2])+[3]+[4]*x", lower_time, upper_time );
        _osc2->SetParameter(0,50.0 );
        _osc2->SetParName(0,"Amplitude");
        _osc2->SetParameter(1,2*3.14159);
        _osc2->SetParName(1,"Period");
        _osc2->SetParameter(2,0);
        _osc2->SetParName(2,"Offset");
	_osc2->SetParameter(3,100);
	_osc2->SetParName(3,"Const");
	_osc2->SetParameter(4,-0.1);
	_osc2->SetParName(4,"Grad");

	TF1 *decay_osc = new TF1( "decay_osc2", "[0]*exp([3]*x)*sin(x*[1]+[2])", lower_time, upper_time );
	decay_osc->SetParameter(0,100.0 );
	decay_osc->SetParName(0,"Amplitude_");
	decay_osc->SetParameter(1,dms);
	decay_osc->SetParName(1,"Period_");
	decay_osc->SetParameter(2,0);
	decay_osc->SetParName(2,"Offset_");
	decay_osc->SetParameter(3,-0.66);
	decay_osc->SetParName(3,"Gamma_");

	TF1 *exp_osc = new TF1( "exp_osc", "[0]*exp([4]*x)+[1]*sin(x*[2]+[3])", lower_time, upper_time );
        exp_osc->SetParameter(0,40.0 );
        exp_osc->SetParName(0,"Amplitude");
        exp_osc->SetParameter(1,dms);
        exp_osc->SetParName(1,"Period");
        exp_osc->SetParameter(2,0);
        exp_osc->SetParName(2,"Offset");
	exp_osc->SetParameter(3,0.);
	exp_osc->SetParName(3,"Amp2");
        exp_osc->SetParameter(4,-0.657);
        exp_osc->SetParName(4,"Gamma");

	TFile* output = new TFile( "sinusoid_output.root", "UPDATE");

	output->cd();
	gDirectory->mkdir("1angle");
	gDirectory->cd("1angle");
	TCanvas* c1 = new TCanvas("new","new",1680,1050);
	c1->Divide(3,3);
	rand_str="";rand_str+=rand->Rndm();
	c1->SetName(rand_str);

	int framenum=0;

	struct FitStruct* fit_struct_0 = new FitStruct();
	fit_struct_0->fitFunction = _osc;
	fit_struct_0->min = 0.;
	fit_struct_0->max = 1.;

	struct FitStruct* fit_struct_1 = new FitStruct();
	fit_struct_1->fitFunction = decay_osc;
	fit_struct_1->min = lower_time;
	fit_struct_1->max = upper_time;

	struct FitStruct* fit_struct_2 = new FitStruct();
	fit_struct_2->fitFunction = _osc;
	fit_struct_2->min = 0.;			//by construction
	fit_struct_2->max = 1.;

	struct FitStruct* fit_struct_3 = new FitStruct();
	fit_struct_3->fitFunction = _osc2;
	fit_struct_3->min = 0.;
        fit_struct_3->max = 1.;

        struct FitStruct* exp_osc_struct = new FitStruct();
        exp_osc_struct->fitFunction = exp_osc;
        exp_osc_struct->min = lower_time;
        exp_osc_struct->max = upper_time;

	vector<TH1*> first;

	first = Plot_Tagged_Combinations( decayTree, "time", selection_cuts, tag_str, "1", weight_1_1angle, weight_2_1angle, rand,
			"All Data", "+sin(#Deltam) 'even' Terms 1-angle", "-sin(#Deltam) 'odd' Terms 1-angle",
			true, 10, exp_osc_struct, fit_struct_1, false);

	TString folded_time = "((time-" + lower_time_str + ")/" + "(2.0*3.1415926535897/"+deltaM_str+")" + ")";

	folded_time = "(" + folded_time + "-" + "int(" + folded_time + ")" + ")";

        TCanvas* c2 = new TCanvas("new2","new2",1680,1050);
        c2->Divide(3,3);
        rand_str="";rand_str+=rand->Rndm();
        c2->SetName(rand_str);

        vector<TH1*> plots;

	plots = Plot_Tagged_Combinations( decayTree, folded_time, selection_cuts, tag_str, "1", weight_1_1angle, weight_2_1angle, rand,
                        "All Data folded", "+sin(#Deltam) 'even' Terms 1-angle folded", "-sin(#Deltam) 'odd' Terms 1-angle folded",
                        true, 10, fit_struct_3, fit_struct_0, false );

	TCanvas* c3 = new TCanvas("total_proj","total_proj",1680, 1050);

	Double_t* temp_arr = NULL;

	TH1* new_hist = new TH1D( "total_proj_h", "total_proj_h", plots[7]->GetNbinsX(), temp_arr );

	new_hist->Add( plots[7], plots[8], 1., -1. );

	//new_hist->Smooth();
	//OptimallyRebin( new_hist );
	//new_hist->Smooth();

	new_hist->Draw( "EP" );
	c3->Update();

	c3->Print("output_proj.png");

	//Plot_Tagged_Combinations( decayTree, "time", selection_cuts, tag_str, A0_weight, Apara_weight, Aperp_weight, rand,
	//		"A_0", "A_#parallel", "A_#perp", true, exp_osc_struct, fit_struct_1 );

	//Plot_Tagged_Combinations( decayTree, "time", selection_cuts, tag_str, Re_int, Im_A0Ap, Im_ApAp, rand,
	//		"Re_int", "Im(A_{0}A_#perp)", "Im(A_{#para}A_#perp)", true, exp_osc_struct, fit_struct_1 );

	return 0;
}

