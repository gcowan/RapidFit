// $Id: PerEventMistagHistogram.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class PerEventMistagHistogram PerEventMistagHistogram_withAngDist.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution + non-trivial angular distribution realised by a histogram
 *
 *  @author Greig Cowan 
 *  @date 2012-04-03
 */

#ifndef PerEventMistagHistogram_H
#define PerEventMistagHistogram_H

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"

#include "BasePDF.h"

class PerEventMistagHistogram : public BasePDF
{
	public:
		//PerEventMistagHistogram();
		PerEventMistagHistogram( PDFConfigurator* );
		//PerEventMistagHistogram( const PerEventMistagHistogram& );
		~PerEventMistagHistogram();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		//PerEventMistagHistogram& operator=( const PerEventMistagHistogram& );
		void MakePrototypes();
		vector<string> GetDoNotIntegrateList();
		bool SetPhysicsParameters(ParameterSet*);
		double buildPDFnumerator();
		double buildPDFdenominator();
		double timeMassFactor();

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF.
		ObservableRef mistagName;	

		double mistag;

		//Additions to deal with 3-D angular distribution via a histogram
		TH1D *histo;
		TAxis *xaxis;
		int nxbins;
		double xmin, xmax, deltax;
		double total_num_entries;
};

#endif
