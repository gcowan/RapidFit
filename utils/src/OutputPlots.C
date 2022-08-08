#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "OutputPlots.h"

using namespace::std;

OutputPlots::OutputPlots( TString Plot_Type, TCanvas* inputCanv, TObject* inputGraph, TLegend* inputLegend, TPaveText* inputLabel, TString inputString )
	: plotType( Plot_Type ), canvas( inputCanv ), internal_graphs( inputGraph ), legend( inputLegend ), label( inputLabel ), DrawString( inputString )
{
}

void OutputPlots::SetDrawString( TString DrawOptions )
{
(void) DrawOptions;	
}

void OutputPlots::Draw( TString DrawOptions )
{
(void) DrawOptions;
}

void OutputPlots::Write( TString fileName )
{
(void) fileName;
}

void OutputPlots::SetResolution( int x, int y )
{
(void) x; (void) y;
}

vector<TObject*> OutputPlots::GetInternalGraphs()
{
return vector<TObject*>();
}

TString OutputPlots::GetPlotType()
{
	return plotType;
}
