#ifndef _OUTPUT_PLOTS_H
#define _OUTPUT_PLOTS_H

#include "TString.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace::std;

class OutputPlots
{
	public:
		OutputPlots( TString Plot_Type, TCanvas* inputCanv, TObject* inputGraph, TLegend* inputLegend, TPaveText* inputLabel, TString inputString="APL" );
		void SetDrawString( TString DrawOptions );
		void Draw( TString DrawOptions );
		void Write( TString fileName="Output.pdf" );
		void SetResolution( int x, int y );
		vector<TObject*> GetInternalGraphs();
		TString GetPlotType();

	private:
		OutputPlots( const OutputPlots& );
		OutputPlots& operator= ( const OutputPlots& );

		TCanvas* canvas;
		TObject* internal_graphs;
		TLegend* legend;
		TPaveText* label;
		TString DrawString;
		TString plotType;
};

#endif
