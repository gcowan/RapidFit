// $Id: EdStyle.cpp,v 1.1 2009/11/10 10:35:46 gcowan Exp $
/**
  @class EdStyle

  Default style for RapidFit plots

  @author Greig A Cowan greig.alan.cowan@cern.ch
  @date 2009-10-02
 */

// Include files 



// local
#include "EdStyle.h"
#include "TLatex.h"
#include "TText.h"
#include "TPaveText.h"
#include "TROOT.h"
#include <iostream>
using namespace std;

//-----------------------------------------------------------------------------
// Implementation file for class : EdStyle
//
// 2007-12-07 : Greig Cowan
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================

EdStyle::EdStyle( ) {

	// use helvetica-bold-r-normal, precision 2 (rotatable)
	Int_t lhcbFont = 62;
	// line thickness
	Double_t lhcbWidth = 1.75;

	// use plain black on white colors
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetStatColor(0);
	gStyle->SetPalette(1);
	//gStyle->SetTitleColor(0);
	//gStyle->SetFillColor(0);

	// set the paper & margin sizes
	gStyle->SetPaperSize(20,26);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.05); // increase for colz plots!!
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetPadLeftMargin(0.14);

	// use large fonts
	gStyle->SetTextFont(lhcbFont);
	gStyle->SetTextSize(0.08);
	gStyle->SetLabelFont(lhcbFont,"x");
	gStyle->SetLabelFont(lhcbFont,"y");
	gStyle->SetLabelFont(lhcbFont,"z");
	gStyle->SetLabelSize(0.05,"x");
	gStyle->SetLabelSize(0.05,"y");
	gStyle->SetLabelSize(0.05,"z");
	gStyle->SetTitleFont(lhcbFont);
	gStyle->SetTitleSize(0.06,"x");
	gStyle->SetTitleSize(0.06,"y");
	gStyle->SetTitleSize(0.06,"z");

	// use bold lines and markers
	gStyle->SetLineWidth(lhcbWidth);
	gStyle->SetFrameLineWidth(3.);
	gStyle->SetHistLineWidth(lhcbWidth);
	gStyle->SetFuncWidth(lhcbWidth);
	gStyle->SetGridWidth(lhcbWidth);
	gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
	//gStyle->SetMarkerStyle(15);
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(1.5);

	// label offsets
	gStyle->SetLabelOffset(0.015);

	// by default, do not display histogram decorations:
	gStyle->SetOptStat(0);  
	gStyle->SetOptStat(1110);  // show only nent, mean, rms
	//gStyle->SetOptTitle(0);
	gStyle->SetOptFit(0);
	//gStyle->SetOptFit(1011); // show probability, parameters and errors

	// look of the statistics box:
	gStyle->SetStatBorderSize(1);
	gStyle->SetStatFont(lhcbFont);
	gStyle->SetStatFontSize(0.05);
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	gStyle->SetStatW(0.25);
	gStyle->SetStatH(0.15);

	// put tick marks on top and RHS of plots
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

	// histogram divisions: only 5 in x to avoid label overlaps
	gStyle->SetNdivisions(505,"x");
	gStyle->SetNdivisions(510,"y");

	TPaveText *lhcbName = new TPaveText(0.65,0.8,0.9,0.9,"BRNDC");
	lhcbName->SetFillColor(0);
	lhcbName->SetTextAlign(12);
	lhcbName->SetBorderSize(0);
	lhcbName->AddText("LHCb");

	TPaveText *lhcbPrelimR = new TPaveText(0.70 - gStyle->GetPadRightMargin(),
			0.80 - gStyle->GetPadTopMargin(),
			0.95 - gStyle->GetPadRightMargin(),
			0.85 - gStyle->GetPadTopMargin(),
			"BRNDC");
	lhcbPrelimR->SetFillColor(0);
	lhcbPrelimR->SetTextAlign(12);
	lhcbPrelimR->SetBorderSize(0);
	lhcbPrelimR->AddText("#splitline{LHCb}{#scale[1.0]{Preliminary}}");

	TPaveText *lhcbPrelimL = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
			0.87 - gStyle->GetPadTopMargin(),
			gStyle->GetPadLeftMargin() + 0.30,
			0.95 - gStyle->GetPadTopMargin(),
			"BRNDC");
	lhcbPrelimL->SetFillColor(0);
	lhcbPrelimL->SetTextAlign(12);
	lhcbPrelimL->SetBorderSize(0);
	lhcbPrelimL->AddText("#splitline{LHCb}{#scale[1.0]{Preliminary}}");

	/*
	TPaveText *lhcb7TeVPrelimR = new TPaveText(0.70 - gStyle->GetPadRightMargin(),
			0.75 - gStyle->SetPadTopMargin(0.05),
			0.95 - gStyle->GetPadRightMargin(),
			0.85 - gStyle->SetPadTopMargin(0.05),
			"BRNDC");
	lhcb7TeVPrelimR->SetFillColor(0);
	lhcb7TeVPrelimR->SetTextAlign(12);
	lhcb7TeVPrelimR->SetBorderSize(0);
	lhcb7TeVPrelimR->AddText("#splitline{#splitline{LHCb}{Preliminary}}{#scale[0.7]{#sqrt{s} = 7 TeV Data}}");

	TPaveText *lhcb7TeVPrelimL = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
			0.78 - gStyle->SetPadTopMargin(0.05),
			gStyle->GetPadLeftMargin() + 0.30,
			0.88 - gStyle->SetPadTopMargin(0.05),
			"BRNDC");
	lhcb7TeVPrelimL->SetFillColor(0);
	lhcb7TeVPrelimL->SetTextAlign(12);
	lhcb7TeVPrelimL->SetBorderSize(0);
	lhcb7TeVPrelimL->SetTextSize(0.06);
	lhcb7TeVPrelimL->AddText("#splitline{#splitline{LHCb}{Preliminary}}{#scale[0.7]{#sqrt{s} = 7 TeV Data}}");

	TPaveText *lhcb0_9TeVPrelimR = new TPaveText(0.70 - gStyle->GetPadRightMargin(),
			0.75 - gStyle->SetPadTopMargin(0.05),
			0.95 - gStyle->GetPadRightMargin(),
			0.85 - gStyle->SetPadTopMargin(0.05),
			"BRNDC");
	lhcb7TeVPrelimR->SetFillColor(0);
	lhcb7TeVPrelimR->SetTextAlign(12);
	lhcb7TeVPrelimR->SetBorderSize(0);
	lhcb7TeVPrelimR->AddText("#splitline{#splitline{LHCb}{Preliminary}}{#scale[0.7]{#sqrt{s} = 900 eV Data}}");

	TPaveText *lhcb0_9TeVPrelimL = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
			0.78 - gStyle->SetPadTopMargin(0.05),
			gStyle->GetPadLeftMargin() + 0.30,
			0.88 - gStyle->SetPadTopMargin(0.05),
			"BRNDC");
	lhcb0_9TeVPrelimL->SetFillColor(0);
	lhcb0_9TeVPrelimL->SetTextAlign(12);
	lhcb0_9TeVPrelimL->SetBorderSize(0);
	lhcb0_9TeVPrelimL->SetTextSize(0.06);
	lhcb0_9TeVPrelimL->AddText("#splitline{#splitline{LHCb}{Preliminary}}{#scale[0.7]{#sqrt{s} = 900 GeV Data}}");
	*/

	TText *lhcbLabel = new TText();
	lhcbLabel->SetTextFont(lhcbFont);
	lhcbLabel->SetTextColor(1);
	lhcbLabel->SetTextSize(0.04);
	lhcbLabel->SetTextAlign(12);

	TLatex *lhcbLatex = new TLatex();
	lhcbLatex->SetTextFont(lhcbFont);
	lhcbLatex->SetTextColor(1);
	lhcbLatex->SetTextSize(0.04);
	lhcbLatex->SetTextAlign(12);

	//gROOT->SetStyle("gStyle");
	//gROOT->ForceStyle();

}


/*
   EdStyle::EdStyle(  ) {
//
// based on a style file from BaBar
//

//..BABAR style from RooLogon.C in workdir

// use plain black on white colors
Int_t icol=0;
gStyle->SetFrameBorderMode(icol);
gStyle->SetCanvasBorderMode(icol);
gStyle->SetPadBorderMode(icol);
gStyle->SetPadColor(icol);
gStyle->SetCanvasColor(icol);
gStyle->SetStatColor(icol);
//gStyle->SetFillColor(icol);

// set the paper & margin sizes
//   gStyle->SetPaperSize(20,26);
//   gStyle->SetPadTopMargin(0.05);
//   gStyle->SetPadRightMargin(0.05);
//   gStyle->SetPadBottomMargin(0.16);
//   gStyle->SetPadLeftMargin(0.12);
/*   
// use large fonts
//Int_t font=72;
Int_t font=42;
Double_t tsize=0.045;
gStyle->SetTextFont(font);
gStyle->SetTextSize(tsize);
gStyle->SetLabelFont(font,"x");
gStyle->SetTitleFont(font,"x");
gStyle->SetLabelFont(font,"y");
gStyle->SetTitleFont(font,"y");
gStyle->SetLabelFont(font,"z");
gStyle->SetTitleFont(font,"z");

gStyle->SetLabelSize(tsize,"x");
gStyle->SetTitleSize(tsize,"x");
gStyle->SetLabelSize(tsize,"y");
gStyle->SetTitleSize(tsize,"y");
gStyle->SetLabelSize(tsize,"z");
gStyle->SetTitleSize(tsize,"z");
gStyle->SetTitleSize(0.05);


//use bold lines and markers
gStyle->SetMarkerStyle(20);
gStyle->SetMarkerSize(1.2);
gStyle->SetHistLineWidth(2.0);
gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

//get rid of X error bars and y error bar caps
//gStyle->SetErrorX(0.001);

//do not display any of the standard histogram decorations
gStyle->SetOptTitle(0);

gStyle->SetOptStat(0);

gStyle->SetOptFit(0);

// put tick marks on top and RHS of plots
gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);

//   gROOT->SetStyle("GREIG");

//gStyle->SetPadTickX(1);
//gStyle->SetPadTickY(1);
 */
//gStyle->SetOptStat(1111);
/*
   gStyle->SetOptFit(1111);
   gStyle->SetMarkerStyle(20);
   gStyle->SetMarkerSize(1.);
   }
 */

//=============================================================================
// Destructor
//=============================================================================
EdStyle::~EdStyle() {} 

//=============================================================================
