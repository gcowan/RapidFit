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
 gStyle->SetOptStat(1111);
 gStyle->SetOptFit(1111);
 gStyle->SetMarkerStyle(20);
 gStyle->SetMarkerSize(1.);
}
//=============================================================================
// Destructor
//=============================================================================
EdStyle::~EdStyle() {} 

//=============================================================================
