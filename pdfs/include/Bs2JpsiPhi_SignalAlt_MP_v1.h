// $Id: Bs2JpsiPhi_SignalAlt_MP_v1.h,v 1.1 2009/12/06  Pete Clarke Exp $
/** @class Bs2JpsiPhi_SignalAlt_MP_v1.h
 *
 *  PDF for Bs2JpsiPhi Signal - new release for feb 2011
 *
 *  @author Pete Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-13
 *
 */

#ifndef Bs2JpsiPhi_SignalAlt_MP_v1_H
#define Bs2JpsiPhi_SignalAlt_MP_v1_H

#include "BasePDF.h"

#include "Bs2JpsiPhi_SignalAlt_BaseClass.h"
#include "RooComplex.h"

#include <exception>
#include <iostream>



class Bs2JpsiPhi_SignalAlt_MP_v1 : public BasePDF,  public Bs2JpsiPhi_SignalAlt_BaseClass
{
	public:
		Bs2JpsiPhi_SignalAlt_MP_v1( PDFConfigurator* );
		~Bs2JpsiPhi_SignalAlt_MP_v1();

		//Mandatory method to evaluate the PDF value:
		virtual double Evaluate(DataPoint*) ;

		//Other operating methods
		virtual bool SetPhysicsParameters(ParameterSet*);
		virtual vector<string> GetDoNotIntegrateList();

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*)  ;

	private:
	
		void MakePrototypes();

		bool normalisationCacheValid ;
		double normalisationCacheValueRes1[3] ;
		double normalisationCacheValueRes2[3] ;

};

#endif
