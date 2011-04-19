// $Id: Bs2JpsiPhi_SignalAlt_MO_v1.h,v 1.1 2009/12/06  Pete Clarke Exp $
/** @class Bs2JpsiPhi_SignalAlt_MO_v1 
 *
 *  Bs2JpsiPhi_SignalAlt series with mistag as observable
 *
 *  @author Pete Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-13
 */

#ifndef Bs2JpsiPhi_SignalAlt_MO_v1_H
#define Bs2JpsiPhi_SignalAlt_MO_v1_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif

#include "Bs2JpsiPhi_SignalAlt_BaseClass.h"
#include "RooComplex.h"

#include <exception>

class Bs2JpsiPhi_SignalAlt_MO_v1 : public BasePDF,  public Bs2JpsiPhi_SignalAlt_BaseClass
{
	public:
		Bs2JpsiPhi_SignalAlt_MO_v1();
		~Bs2JpsiPhi_SignalAlt_MO_v1();

		//Mandatory method to evaluate the PDF value:
		virtual double Evaluate(DataPoint*);

		//Other operating methods
		virtual bool SetPhysicsParameters(ParameterSet*);
		virtual vector<string> GetDoNotIntegrateList();

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:

		void MakePrototypes();

		bool normalisationCacheValid ;
		double normalisationCacheValueRes1[3] ;
		double normalisationCacheValueRes2[3] ;
		
		pair<string,int> mistagScaleName ;
		pair<string,int> mistagOffsetName ;
		double mistagScale ;
		double mistagOffset ;
};

#endif
