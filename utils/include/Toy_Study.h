#ifndef RAPIDFIT_TOY_H
#define RAPIDFIT_TOY_H

#include "TString.h"
#include "TRandom3.h"
#include "TTree.h"

#include <vector>

using namespace::std;

class ToyStudyAnalysis
{
	public:
		static int Toy_Study( TTree* input_trees, TRandom3* rand_gen, vector<string> OtherOptions );

		static void Help();

	private:

		static void LatexDocHeader( stringstream& latex );

		static void LatexDocFooter( stringstream& latex );

		static void TableHeader( stringstream& latex );

		static void TableFooter( stringstream& latex );

		static void ImageSplit( stringstream& latex );

		static void ProcessImageContent( stringstream& latex, vector<TString> all_parameter_plots );

		static void ProcessTableContent( stringstream& latex, vector<TString> all_parameter_plots, vector<pair<pair<double,double>,pair<double,double> > > table_content );
};

#endif

