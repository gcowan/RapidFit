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

		int componentIndex;

		double m23;
		double cosTheta1; // This is the member variable used in the "builder" functions 
		double cosTheta2; // These are the physics parameters varied in the fit and passed from the XML;
		double phi;

		// These contain the ObservableRefs that correspond
		// to the observable names that are used in the
		// PDF. 
		ObservableRef fracA0sqZplusName;
		ObservableRef fracApsqZplusName;
		ObservableRef fracZplusName;
		ObservableRef phaseA0ZplusName;
		ObservableRef phaseApZplusName;
		ObservableRef phaseAmZplusName;
		
		ObservableRef fracA0sqKst892Name;
		ObservableRef fracApsqKst892Name;
		ObservableRef phaseA0Kst892Name;
		ObservableRef phaseApKst892Name;
		ObservableRef phaseAmKst892Name;
		
		ObservableRef fracA0sqKst1410Name;
		ObservableRef fracApsqKst1410Name;
		ObservableRef fracKst1410Name;
		ObservableRef phaseA0Kst1410Name;
		ObservableRef phaseApKst1410Name;
		ObservableRef phaseAmKst1410Name;
		
		ObservableRef fracA0sqKst1680Name;
		ObservableRef fracApsqKst1680Name;
		ObservableRef fracKst1680Name;
		ObservableRef phaseA0Kst1680Name;
		ObservableRef phaseApKst1680Name;
		ObservableRef phaseAmKst1680Name;
		
		ObservableRef fracK01430Name;
		ObservableRef phaseA0K01430Name;
		ObservableRef fracK21430Name;
		ObservableRef phaseA0K21430Name;
		
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

		ObservableRef frac_LASSName;
		ObservableRef a_LASSName;
		ObservableRef r_LASSName;

		double fracA0sqZplus;
		double fracApsqZplus;
		double fracZplus;
		double phaseA0Zplus;
		double phaseApZplus;
		double phaseAmZplus;
		
		double fracA0sqKst892;
		double fracApsqKst892;
		double phaseA0Kst892;
		double phaseApKst892;
		double phaseAmKst892;
		
		double fracA0sqKst1410;
		double fracApsqKst1410;
		double fracKst1410;
		double phaseA0Kst1410;
		double phaseApKst1410;
		double phaseAmKst1410;
		
		double fracA0sqKst1680;
		double fracApsqKst1680;
		double fracKst1680;
		double phaseA0Kst1680;
		double phaseApKst1680;
		double phaseAmKst1680;
		
		double fracK01430;
		double phaseA0K01430;
		double fracK21430;
		double phaseA0K21430;
		
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


		//LASS
		double frac_LASS;
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

		double mJpsi;

		bool useAngularAcceptance;
		string fullFileName;
		TFile * histogramFile;
		//TH1D * angularAccHistCosTheta1;
		//TH1D * angularAccHistPhi;
		//TH2D * angularAccHistMassCosTheta2;
                THnSparse * histo;
                TAxis *xaxis, *yaxis, *zaxis, *maxis;
                int nxbins, nybins, nzbins, nmbins;
};

#endif
