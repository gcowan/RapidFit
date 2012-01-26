#ifndef AccCorr_H
#define AccCorr_H

#include <iostream>
#include "PDFConfigurator.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH3D.h"
#include "math.h"
#include "TMath.h"
#include "RooMath.h"
#include "TAxis.h"
#include "TH1.h"
#include "Mathematics.h"
#include "XMLConfigReader.h"
#include <vector>
#include "Bs2PhiPhi_AccCorr.h"

class AccCorr
{

	public:
		AccCorr();
		AccCorr( const AccCorr& );
		virtual void AccCorrInit(PDFConfigurator*);
		virtual ~AccCorr();
		virtual double AccEval(vector<double>);
		void MakeWeights(PDFConfigurator*);
		vector<double> ReturnWeights();
	
	private:
		AccCorr& operator = ( const AccCorr& );		
		bool useAccFit;
		TH3D* histo;
                int nxbins, nybins, nzbins;
                double xmin, xmax, ymin, ymax, zmin, zmax, deltax, deltay, deltaz, total_num_entries;
		double angNorm;
		double w1,w2,w3,w4,w5,w6;

};

#endif
