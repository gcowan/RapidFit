// $Id: DPHistoBackground.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class DPHistoBackground DPHistoBackground_withAngDist.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution + non-trivial angular distribution realised by a histogram
 *
 *  @author Greig Cowan 
 *  @date 2012-04-03
 */

#ifndef DPHistoBackground_H
#define DPHistoBackground_H

#include "TROOT.h"
#include "TFile.h"
#include "THnSparse.h"

#include "BasePDF.h"

class DPHistoBackground : public BasePDF
{
	public:
		//DPHistoBackground();
		DPHistoBackground( PDFConfigurator* );
		DPHistoBackground( const DPHistoBackground& );
		~DPHistoBackground();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		DPHistoBackground& operator=( const DPHistoBackground& );
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);
		double buildPDFnumerator();
		double buildPDFdenominator();
		double angleMassFactor();

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF.
		ObservableRef massName;		
		ObservableRef cosTheta1Name;	
		ObservableRef phiName;			
		ObservableRef cosTheta2Name;	

		double mass;
		double cos1;
		double cos2;
		double phi;

		//Additions to deal with 3-D angular distribution via a histogram
                string fullFileName;
		TFile * histogramFile;
		THnSparse * histo;
		TAxis *xaxis, *yaxis, *zaxis, *maxis;
		int nxbins, nybins, nzbins, nmbins;
		double xmin, xmax, ymin, ymax, zmin, zmax, mmin, mmax, deltax, deltay, deltaz, deltam;
		double total_num_entries;
		bool useFlatAngularDistribution;
};

#endif
