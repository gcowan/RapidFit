// $Id: DPTotalAmplitudePDF_withAcc_withBkg.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class DPTotalAmplitudePDF_withAcc_withBkg DPTotalAmplitudePDF_withAcc_withBkg.h
 *
 *  PDF for Bd2JpsiKpi dalitz
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef DPTotalAmplitudePDF_withAcc_withBkg_H
#define DPTotalAmplitudePDF_withAcc_withBkg_H

#include "BasePDF.h"
#include "../dalitz/include/DPComponent.hh"
#include "../dalitz/include/DPWignerFunctionJ1over2.hh"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TLorentzVector.h"

class DPTotalAmplitudePDF_withAcc_withBkg : public BasePDF
{
	public:
		DPTotalAmplitudePDF_withAcc_withBkg( PDFConfigurator* );
		DPTotalAmplitudePDF_withAcc_withBkg( const DPTotalAmplitudePDF_withAcc_withBkg & );
		~DPTotalAmplitudePDF_withAcc_withBkg();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);
		double EvaluateComponent(DataPoint * measurement, ComponentRef* Component);
		vector<string> PDFComponents();

	protected:
                virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);
        vector<string> GetDoNotIntegrateList();

		// Experimental observables
		ObservableRef m23Name;
		ObservableRef cosTheta1Name;
		ObservableRef cosTheta2Name;
		ObservableRef phiName;
        ObservableRef pionIDName;
		ObservableRef m13Name;
		ObservableRef cosZName;
		ObservableRef cosPsiZName;
		ObservableRef phiZName;
        ObservableRef alphaName;

		int componentIndex;

		double m23;
		double cosTheta1;
		double cosTheta2;
		double phi;
        int pionID;
		double m13;
		double cosZ;
		double cosPsiZ;
		double phiZ;
        double alpha;

		// These contain the ObservableRefs that correspond
		// to the observable names that are used in the
		// PDF.
		ObservableRef fractionName;
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

		ObservableRef magA0K42045Name;
		ObservableRef magApK42045Name;
		ObservableRef magAmK42045Name;
		ObservableRef phaseA0K42045Name;
		ObservableRef phaseApK42045Name;
		ObservableRef phaseAmK42045Name;

		ObservableRef magA0K52380Name;
		ObservableRef magApK52380Name;
		ObservableRef magAmK52380Name;
		ObservableRef phaseA0K52380Name;
		ObservableRef phaseApK52380Name;
		ObservableRef phaseAmK52380Name;

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
		ObservableRef massK42045Name;
		ObservableRef widthK42045Name;
		ObservableRef massK52380Name;
		ObservableRef widthK52380Name;
		ObservableRef massK800Name;
		ObservableRef widthK800Name;

		ObservableRef mag_LASSName;
		ObservableRef phase_LASSName;
		ObservableRef a_LASSName;
		ObservableRef r_LASSName;

		double fraction;
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

		double magA0K42045;
		double magApK42045;
		double magAmK42045;
		double phaseA0K42045;
		double phaseApK42045;
		double phaseAmK42045;

		double magA0K52380;
		double magApK52380;
		double magAmK52380;
		double phaseA0K52380;
		double phaseApK52380;
		double phaseAmK52380;

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
		double massK42045;
		double widthK42045;
		double massK52380;
		double widthK52380;
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

		DPWignerFunctionJ1over2 wigner;

		double massPsi;
		double massB;

        // Acceptance parameterisation
        bool useAngularAcceptance;
        static const int l_max = 6;
        static const int i_max = 4;
        static const int k_max = 2;
        static const int j_max = 2;
        double c[l_max+1][i_max+1][k_max+1][j_max+1];
        static const int l_max_b = 6;
        static const int i_max_b = 2;
        static const int k_max_b = 1;
        static const int j_max_b = 2;
        double b[l_max_b+1][i_max_b+1][k_max_b+1][j_max_b+1];
};

#endif
