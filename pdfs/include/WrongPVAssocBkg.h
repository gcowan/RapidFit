// $Id: WrongPVAssocBkg.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class WrongPVAssocBkg WrongPVAssocBkg_withAngDist.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution + non-trivial angular distribution realised by a histogram
 *
 *  @author Greig Cowan 
 *  @date 2012-04-03
 */

#ifndef WrongPVAssocBkg_H
#define WrongPVAssocBkg_H

#include "TROOT.h"
#include "TFile.h"
#include "TH2D.h"

#include "BasePDF.h"

class WrongPVAssocBkg : public BasePDF
{
	public:
		//WrongPVAssocBkg();
		WrongPVAssocBkg( PDFConfigurator* );
		WrongPVAssocBkg( const WrongPVAssocBkg& );
		~WrongPVAssocBkg();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		WrongPVAssocBkg& operator=( const WrongPVAssocBkg& );
		void MakePrototypes();
		vector<string> GetDoNotIntegrateList();
		bool SetPhysicsParameters(ParameterSet*);
		double buildPDFnumerator();
		double buildPDFdenominator();
		double timeMassFactor();

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF.
		ObservableRef massName;		
		ObservableRef timeName;	
		ObservableRef eventResolutionName;	

		double mass;
		double time;

		//Additions to deal with 3-D angular distribution via a histogram
		TH2D *histo;
		TAxis *xaxis, *yaxis;
		int nxbins, nybins;
		double xmin, xmax, ymin, ymax, deltax, deltay;
		double total_num_entries;
		double total_num_entries_phase_space;
		bool normalisationDone ;
		bool _makeFlat ;
};

#endif
