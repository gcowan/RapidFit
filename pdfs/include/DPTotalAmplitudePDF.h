// $Id: DPTotalAmplitudePDF.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class DPTotalAmplitudePDF DPTotalAmplitudePDF.h
 *
 *  PDF for Bd2JpsiKpi dalitz
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef DPTotalAmplitudePDF_H
#define DPTotalAmplitudePDF_H

#include "BasePDF.h"
#include "../dalitz/include/DPComponent.hh"
#include "../dalitz/include/DPWignerFunctionJ1.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TLorentzVector.h"

class DPTotalAmplitudePDF : public BasePDF
{
	public:
		DPTotalAmplitudePDF( PDFConfigurator* );
		DPTotalAmplitudePDF( const DPTotalAmplitudePDF & );
		~DPTotalAmplitudePDF();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);
		double EvaluateComponent(DataPoint * measurement, ComponentRef* Component);
		vector<string> PDFComponents();

	protected:
                virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);

		// Experimental observables
		ObservableRef m23Name;
		ObservableRef cosTheta1Name;
		ObservableRef cosTheta2Name;
		ObservableRef phiName;
        ObservableRef pionIDName;

		int componentIndex;

		double m23;
		double cosTheta1; // This is the member variable used in the "builder" functions
		double cosTheta2; // These are the physics parameters varied in the fit and passed from the XML;
		double phi;
        int pionID;

		// These contain the ObservableRefs that correspond
		// to the observable names that are used in the
		// PDF.
		ObservableRef magA0ZplusName;
		ObservableRef magApZplusName;
		ObservableRef magAmZplusName;
		ObservableRef phaseA0ZplusName;
		ObservableRef phaseApZplusName;
		ObservableRef phaseAmZplusName;

		ObservableRef magA0Kst892Name;
		ObservableRef magApKst892Name;
		ObservableRef magAmKst892Name;
		ObservableRef phaseA0Kst892Name;
		ObservableRef phaseApKst892Name;
		ObservableRef phaseAmKst892Name;

		ObservableRef magA0Kst1410Name;
		ObservableRef magApKst1410Name;
		ObservableRef magAmKst1410Name;
		ObservableRef phaseA0Kst1410Name;
		ObservableRef phaseApKst1410Name;
		ObservableRef phaseAmKst1410Name;

		ObservableRef magA0Kst1680Name;
		ObservableRef magApKst1680Name;
		ObservableRef magAmKst1680Name;
		ObservableRef phaseA0Kst1680Name;
		ObservableRef phaseApKst1680Name;
		ObservableRef phaseAmKst1680Name;

		ObservableRef magA0K01430Name;
		ObservableRef phaseA0K01430Name;

		ObservableRef magA0K21430Name;
		ObservableRef magApK21430Name;
		ObservableRef magAmK21430Name;
		ObservableRef phaseA0K21430Name;
		ObservableRef phaseApK21430Name;
		ObservableRef phaseAmK21430Name;

		ObservableRef magA0K31780Name;
		ObservableRef magApK31780Name;
		ObservableRef magAmK31780Name;
		ObservableRef phaseA0K31780Name;
		ObservableRef phaseApK31780Name;
		ObservableRef phaseAmK31780Name;

		ObservableRef magA0K800Name;
		ObservableRef phaseA0K800Name;

		ObservableRef magA0NRName;
		ObservableRef phaseA0NRName;

		ObservableRef massZplusName;
		ObservableRef widthZplusName;
		ObservableRef massKst892Name;
		ObservableRef widthKst892Name;
		ObservableRef massKst1410Name;
		ObservableRef widthKst1410Name;
		ObservableRef massKst1680Name;
		ObservableRef widthKst1680Name;
		ObservableRef massK01430Name;
		ObservableRef widthK01430Name;
		ObservableRef massK21430Name;
		ObservableRef widthK21430Name;
		ObservableRef massK31780Name;
		ObservableRef widthK31780Name;
		ObservableRef massK800Name;
		ObservableRef widthK800Name;

		ObservableRef mag_LASSName;
		ObservableRef phase_LASSName;
		ObservableRef a_LASSName;
		ObservableRef r_LASSName;

		double magA0Zplus;
		double magApZplus;
		double magAmZplus;
		double phaseA0Zplus;
		double phaseApZplus;
		double phaseAmZplus;

		double magA0Kst892;
		double magApKst892;
		double magAmKst892;
		double phaseA0Kst892;
		double phaseApKst892;
		double phaseAmKst892;

		double magA0Kst1410;
		double magApKst1410;
		double magAmKst1410;
		double phaseA0Kst1410;
		double phaseApKst1410;
		double phaseAmKst1410;

		double magA0Kst1680;
		double magApKst1680;
		double magAmKst1680;
		double phaseA0Kst1680;
		double phaseApKst1680;
		double phaseAmKst1680;

		double magA0K01430;
		double phaseA0K01430;

		double magA0K21430;
		double magApK21430;
		double magAmK21430;
		double phaseA0K21430;
		double phaseApK21430;
		double phaseAmK21430;

		double magA0K31780;
		double magApK31780;
		double magAmK31780;
		double phaseA0K31780;
		double phaseApK31780;
		double phaseAmK31780;

		double magA0K800;
		double phaseA0K800;

		double magA0NR;
		double phaseA0NR;

		double massZplus;
		double widthZplus;
		double massKst892;
		double widthKst892;
		double massKst1410;
		double widthKst1410;
		double massKst1680;
		double widthKst1680;
		double massK01430;
		double widthK01430;
		double massK21430;
		double widthK21430;
		double massK31780;
		double widthK31780;
		double massK800;
		double widthK800;

		//LASS
		double mag_LASS;
		double phase_LASS;
		double a_LASS;
		double r_LASS;

		// These values are cached since they are calculated for each event
		TLorentzVector pMuPlus;
		TLorentzVector pMuMinus;
		TLorentzVector pPi;
		TLorentzVector pK;
	 	TLorentzVector pB;
                // Cos of the angle between psi reference axis
                double cosARefs;

		std::vector<DPComponent*> KpiComponents;
		std::vector<DPComponent*> ZComponents;

		DPWignerFunctionJ1 wigner;

		double massPsi;

		bool useAngularAcceptance;
		bool useFourDHistogram;
		string fullFileName;
		TFile * histogramFile;
		TH1D * angularAccHistCosTheta1;
		TH1D * angularAccHistPhi;
		TH2D * angularAccHistMassCosTheta2;
                //THnSparse * histo;
                TH3D * histo;
                TAxis *xaxis, *yaxis, *zaxis, *maxis;
                int nxbins, nybins, nzbins, nmbins;
};

#endif
