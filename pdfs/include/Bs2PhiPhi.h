/** @class Bs2PhiPhi Bs2PhiPhi.h
 *
 *  RapidFit PDF for Bs2PhiPhi
 *
 *  @author Young Min Kim
 *  @date 12 Nov 2009
 */

#ifndef Bs2PhiPhi_H
#define Bs2PhiPhi_H

#include "BasePDF.h"

class Bs2PhiPhi : public BasePDF
{
    public:
        Bs2PhiPhi( PDFConfigurator* );
        ~Bs2PhiPhi();

        //Calculate the PDF value
        virtual double Evaluate(DataPoint*);
        virtual bool SetPhysicsParameters(ParameterSet*);
        //Return a list of parameters not to be integrated
        virtual vector<string> GetDoNotIntegrateList();

    protected:
        //Calculate the PDF normalisation
        virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

    private:
        void MakePrototypes();

        //Cached values
        double cach_v1, cach_v2;
        double cach_Azero, cach_Apara, cach_Aperp;
        double cach_sinPhis, cach_cosPhis;
        double cach_sinDelta_[2];
        double cach_cosDelta_[2];
        double cach_cosDeltaSubtr;
        bool normalisationCacheValid, evaluationCacheValid;

        // These contain the ObservableRefs that correspond
        // to the physics parameter names that will be
        // used in the minimiser.
        ObservableRef gamma_sName;       // gamma_s
        ObservableRef gamma_lName;       // gamma_L
        ObservableRef gamma_hName;       // gamma_H
        ObservableRef deltaMName;        // delta mass
        ObservableRef Phi_sName;         // what we want to measure!
        ObservableRef Azero_sqName;      // amplitude
        ObservableRef Apara_sqName;      // amplitude
        ObservableRef Aperp_sqName;      // amplitude
        ObservableRef delta_1Name;       // strong phase
        ObservableRef delta_2Name;       // strong phase

        // These contain the ObservableRefs that correspond
        // to the observable names that are used in the
        // PDF. 
        ObservableRef timeName;    // proper time
        ObservableRef theta_1Name;  // angle between Phi1 p vector and its K+ daughter's p vector
        ObservableRef theta_2Name;  // angle between Phi2 p vector and its K+ daughter's p vector
        ObservableRef phiName;     // angle between the decay planes of the 2 Phi resonances
        ObservableRef tagName;     // B tag
        ObservableRef mistagName;  // B mistag

	ObservableRef timeconstraintName;	
    
        void getPhysicsParameters( double&, double&, double&, double&, double&
                                 , double&, double&, double&, double&, double&);
        // 10 parameters: gamma_{s,l,h}, deltaM, Phi_s, A{zero,para,perp}_sq, delta_{1,2}
        
        void getAngularFunctions( double&, double&, double&, double&, double&, double&, DataPoint*);
        // 6 output parameters: f_{1~6}
        
        void getTimeDependentAmplitudes( double&, double&, double&, double&, double&, double&, DataPoint*, int);
        // double AzeroAzero, AparaApara, AperpAperp, ImAparaAperp, ReAzeroApara, ImAzeroAperp
        // int Btype
        
        void getTimeAmplitudeIntegrals(double&, double&, double&, PhaseSpaceBoundary*, int);
        // double AzeroAzeroInt, AparaAparaInt, AperpAperpInt
        // int Btype

        //!! Pending re-write
        inline double getMainIntAnswer(double, double, double, double, double, double, double, int, int);
        // double tmin, tmax, gamma_s, gamma_l, gamma_h, Dms, phis
        // int Btype, Swap (for K3 which has a few sign changes compared to K1,2)
    
            // Work out what interference terms these correspond to.        
        //inline double A4def(const double, const double, double, double, double, double, double, double, double, double, int);
        //inline double A5def(const double, const double, double, double, double, double, double, double, double, double, int);
        // double tmin, tmax, k0, tauL, tauH, tauBar, Dms, phis, tphase1, tphase2
        // int Btype
};

#endif
