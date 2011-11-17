// $Id: Bs2JpsiPhiMassBkgLL.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhiMassBkgLL Bs2JpsiPhiMassBkgLL.h
 *
 *  RapidFit PDF for Bs2JpsiPhi long lived background
 *
 *  @author Pete Clarke
 *  @date 2011 10 23
 */

#ifndef Bs2DsPi_lowmassbkg_H
#define Bs2DsPi_lowmassbkg_H

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"


#include "BasePDF.h"


class Bs2DsPi_lowmassbkg : public BasePDF
{
	public:
		Bs2DsPi_lowmassbkg(PDFConfigurator*);
		~Bs2DsPi_lowmassbkg();
		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);
		//virtual bool SetPhysicsParameters(ParameterSet*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);



	private:	

		void MakePrototypes();

		// Physics parameters

		// Observables
		ObservableRef massName;
		ObservableRef constraint_massName;

		double mass;

		//Additions to deal with 1-D mass distribution via a histogram
		TH1D *histo;
		TAxis *xaxis;
		int nxbins;
		double xmin, xmax, deltax;
		double total_num_entries;



};

#endif
