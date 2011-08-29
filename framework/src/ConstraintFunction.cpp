/**
  @class ConstraintFunction

  Where external, experimental constraints on PhysicsParameters are calculated

  @author Benjamin M Wynne bwynne@cern.ch
  @date 21-01-10
  */

//	RapidFit Headers
#include "ConstraintFunction.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>

using namespace std;

//Default constructor
ConstraintFunction::ConstraintFunction() : Found_Position(), allConstraints()
{
}

//Constructor with correct arguments
ConstraintFunction::ConstraintFunction( vector< ExternalConstraint* > NewConstraints ) : Found_Position(), allConstraints(NewConstraints)
{
}

//Destructor
ConstraintFunction::~ConstraintFunction()
{
}

//Perform the constraint calculation
double ConstraintFunction::Evaluate( ParameterSet * NewParameters )
{
	vector<string> parameterNames = NewParameters->GetAllNames();
	double constraintValue = 0.0;

	if( Found_Position.size() != allConstraints.size() )
	{
		for( unsigned int constraintIndex = 0; constraintIndex < allConstraints.size(); ++constraintIndex )
		{
			string name = allConstraints[constraintIndex]->GetName();
			Found_Position.push_back( StringProcessing::VectorContains( &parameterNames, &name ) );
		}
	}

	//Loop over all ExternalConstraints
	for (unsigned int constraintIndex = 0; constraintIndex < allConstraints.size(); ++constraintIndex )
	{
		string name = allConstraints[constraintIndex]->GetName();
		if ( name == "GammaL" )
		{
			//GammaL = Gamma + ( deltaGamma / 2 )
			// Get gamma and delta gamma
			double gamma = NewParameters->GetPhysicsParameter("gamma")->GetValue();
			double dgam =  NewParameters->GetPhysicsParameter("deltaGamma")->GetValue();
			double gaml_fit = gamma + (dgam/2.0);
			double gaml_con = allConstraints[constraintIndex]->GetValue();
			double gaussSqrt = ( gaml_fit -  gaml_con ) / allConstraints[constraintIndex]->GetError();
			constraintValue += gaussSqrt * gaussSqrt;
		}
		else if ( name == "GammaObs" )
		{
			//GammaObs^-1 = Gamma^-1 * ( 1 + (deltaGamma/2*Gamma)^2 ) / ( 1 - (deltaGamma/2*Gamma)^2 )
			// Get gamma and delta gamma
			double gamma = NewParameters->GetPhysicsParameter("gamma")->GetValue();
			double dgam =  NewParameters->GetPhysicsParameter("deltaGamma")->GetValue();
			double gamobs_fit = gamma * (1-(dgam/2/gamma)*(dgam/2/gamma)) / (1+(dgam/2/gamma)*(dgam/2/gamma));
			double gamobs_con = allConstraints[constraintIndex]->GetValue();
			double gaussSqrt = ( gamobs_fit -  gamobs_con ) / allConstraints[constraintIndex]->GetError();
			constraintValue += gaussSqrt * gaussSqrt;
		}
		else if ( name == "SpecialForGreig" )
		{
			//GammaObs^-1 = Gamma^-1 * ( 1 + (deltaGamma/2*Gamma)^2 ) / ( 1 - (deltaGamma/2*Gamma)^2 )
			// Get gamma and delta gamma
			double gamma = NewParameters->GetPhysicsParameter("gamma")->GetValue();
			double dgam =  NewParameters->GetPhysicsParameter("deltaGamma")->GetValue();
			double gamma_con = 0.7 ;
			double dgam_con = 0.07 ;
			double v[2] ;
			v[0] = gamma-gamma_con ;
			v[1] = dgam-dgam_con ;
			double E[2][2] ;
			E[0][0] = 1 ;
			E[1][1] = 2 ;
			E[1][0] = 3 ;
			E[0][1] = E[1][0] ;
			double chisq = 0;
			for( int ix=0; ix<2; ++ix ) {
				for( int iy=0; iy<2; ++iy ) {
					chisq += v[ix]*E[ix][iy]*v[iy] ;
				}
			}
			constraintValue += chisq ;
		}
		else if ( name == "ATOTAL" )
		{
			// This is a special one to contrain the sum of the amplitudes to 1
			double Aperpsq = NewParameters->GetPhysicsParameter("Aperp_sq")->GetValue();
			double Azerosq =  NewParameters->GetPhysicsParameter("Azero_sq")->GetValue();
			double Assq =  NewParameters->GetPhysicsParameter("As_sq")->GetValue();
			double Atot_val = allConstraints[constraintIndex]->GetValue();
			double Atot_constraint = allConstraints[constraintIndex]->GetError();
			double excess = (Atot_val-Aperpsq-Azerosq-Assq) ;
			double penalty ;
			if( excess >= 0 ) penalty = 0 ;
			else penalty = (excess*excess) / (Atot_constraint*Atot_constraint) ;
			constraintValue += penalty;
		}
		else if ( name == "GLandGH" )
		{
			// This is a special one to stop GL and GH going negative
			double gamma = NewParameters->GetPhysicsParameter("gamma")->GetValue();
			double dgam =  NewParameters->GetPhysicsParameter("deltaGamma")->GetValue();
			double G1 = gamma-dgam/2.0 ;
			double G2 = gamma+dgam/2.0 ;
			//double val = allConstraints[constraintIndex]->GetValue();
			double constraint = allConstraints[constraintIndex]->GetError();
			double penalty = 0 ;
			if( G1 < 0. ) penalty += (G1*G1)/(constraint*constraint) ;
			if( G2 < 0. ) penalty += (G2*G2)/(constraint*constraint) ;
			if( penalty > 0. ) cout << " GLandGH constraint being applied " << endl ;
			constraintValue += penalty;
		}
		else if ( name == "CosSqPlusSinSq" )
		{
			// This is a special one to contrain cosphis+sinphis to 1
			double cosphis = NewParameters->GetPhysicsParameter("cosphis")->GetValue();
			double sinphis =  NewParameters->GetPhysicsParameter("sinphis")->GetValue();
			double val = allConstraints[constraintIndex]->GetValue();
			double constraint = allConstraints[constraintIndex]->GetError();
			double excess = (val - cosphis*cosphis - sinphis*sinphis ) ;
			double penalty = (excess*excess) / (constraint*constraint) ;
			constraintValue += penalty;
		}
		else if ( Found_Position[constraintIndex] >= 0 )
		{
			//Do standard gaussian constraint calculation
			double parameterValue = NewParameters->GetPhysicsParameter(name)->GetValue();
			double gaussSqrt = ( parameterValue - allConstraints[constraintIndex]->GetValue() ) / allConstraints[constraintIndex]->GetError();
			constraintValue += gaussSqrt * gaussSqrt;
		}
	}

	return 0.5 * constraintValue;
}

