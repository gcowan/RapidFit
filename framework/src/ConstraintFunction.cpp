/*!
 * @class ConstraintFunction
 *
 * Where external, experimental constraints on PhysicsParameters are calculated
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern,ch
 */

//	RapidFit Headers
#include "ConstraintFunction.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <vector>
#include <sstream>

using namespace::std;

//Constructor with correct arguments
ConstraintFunction::ConstraintFunction( const vector< ExternalConstraint* > NewConstraints ) : Found_Position(), allConstraints()
{
	for( vector<ExternalConstraint*>::const_iterator const_i = NewConstraints.begin(); const_i != NewConstraints.end(); ++const_i )
	{
		allConstraints.push_back( new ExternalConstraint( *(*const_i) ) );
	}
}

ConstraintFunction::ConstraintFunction( const ConstraintFunction& input ) : Found_Position( input.Found_Position ), allConstraints()
{
	for( vector<ExternalConstraint*>::const_iterator const_i = input.allConstraints.begin(); const_i != input.allConstraints.end(); ++const_i )
	{
		allConstraints.push_back( new ExternalConstraint( *(*const_i) ) );
	}
}


//Destructor
ConstraintFunction::~ConstraintFunction()
{
	while( !allConstraints.empty() )
	{
		if( allConstraints.back() != NULL ) delete allConstraints.back();
		allConstraints.pop_back();
	}
}

//Perform the constraint calculation
double ConstraintFunction::Evaluate( const ParameterSet * NewParameters )
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

	unsigned int num_i=0;
	vector< vector<ExternalConstraint*>::iterator > to_be_removed;
	//Loop over all ExternalConstraints
	for( vector<ExternalConstraint*>::iterator constraintIndex = allConstraints.begin(); constraintIndex != allConstraints.end(); ++constraintIndex, ++num_i )
	{
		string name = (*constraintIndex)->GetName();
		if ( name == "GammaL" )
		{
			//GammaL = Gamma + ( deltaGamma / 2 )
			// Get gamma and delta gamma
			string gammaName("gamma");
			string deltaGammaName("deltaGamma");
			int gamma_i = StringProcessing::VectorContains(&parameterNames,&gammaName);
			int dG_i = StringProcessing::VectorContains(&parameterNames,&deltaGammaName);
			if( (gamma_i == -1) || (dG_i == -1) )
			{
				cerr << "Cannot Apply " << name << " constraint" << endl;
				to_be_removed.push_back(constraintIndex);
				continue;
			}
			double gamma = NewParameters->GetPhysicsParameter(gammaName)->GetValue();
			double dgam =  NewParameters->GetPhysicsParameter(deltaGammaName)->GetValue();
			double gaml_fit = gamma + (dgam/2.0);
			double gaml_con = (*constraintIndex)->GetValue();
			double gaussSqrt = ( gaml_fit -  gaml_con ) / (*constraintIndex)->GetError();
			constraintValue += gaussSqrt * gaussSqrt;
		}
		else if ( name == "GammaH" )
		{
			//GammaL = Gamma + ( deltaGamma / 2 )
			// Get gamma and delta gamma
			string gammaName("gamma");
			string deltaGammaName("deltaGamma");
			int gamma_i = StringProcessing::VectorContains(&parameterNames,&gammaName);
			int dG_i = StringProcessing::VectorContains(&parameterNames,&deltaGammaName);
			if( (gamma_i == -1) || (dG_i == -1) )
			{
				cerr << "Cannot Apply " << name << " constraint" << endl;
				to_be_removed.push_back(constraintIndex);
				continue;
			}
			double gamma = NewParameters->GetPhysicsParameter(gammaName)->GetValue();
			double dgam =  NewParameters->GetPhysicsParameter(deltaGammaName)->GetValue();
			double gaml_fit = gamma - (dgam/2.0);
			double gaml_con = (*constraintIndex)->GetValue();
			double gaussSqrt = ( gaml_fit -  gaml_con ) / (*constraintIndex)->GetError();
			constraintValue += gaussSqrt * gaussSqrt;
		}
		else if ( name == "GammaObs" )
		{
			//GammaObs^-1 = Gamma^-1 * ( 1 + (deltaGamma/2*Gamma)^2 ) / ( 1 - (deltaGamma/2*Gamma)^2 )
			// Get gamma and delta gamma
			string gammaName("gamma");
			string deltaGammaName("deltaGamma");
			int gamma_i = StringProcessing::VectorContains(&parameterNames,&gammaName);
			int dG_i = StringProcessing::VectorContains(&parameterNames,&deltaGammaName);
			if( (gamma_i == -1) || (dG_i == -1) )
			{
				cerr << "Cannot Apply " << name << " constraint" << endl;
				to_be_removed.push_back(constraintIndex);
				continue;
			}
			double gamma = NewParameters->GetPhysicsParameter(gammaName)->GetValue();
			double dgam =  NewParameters->GetPhysicsParameter(deltaGammaName)->GetValue();
			double gamobs_fit = gamma * (1-(dgam/2/gamma)*(dgam/2/gamma)) / (1+(dgam/2/gamma)*(dgam/2/gamma));
			double gamobs_con = (*constraintIndex)->GetValue();
			double gaussSqrt = ( gamobs_fit -  gamobs_con ) / (*constraintIndex)->GetError();
			constraintValue += gaussSqrt * gaussSqrt;
		}
		else if ( name == "SpecialForGreig" )
		{
			//GammaObs^-1 = Gamma^-1 * ( 1 + (deltaGamma/2*Gamma)^2 ) / ( 1 - (deltaGamma/2*Gamma)^2 )
			// Get gamma and delta gamma
			string gammaName("gamma");
			string deltaGammaName("deltaGamma");
			int gamma_i = StringProcessing::VectorContains(&parameterNames,&gammaName);
			int dG_i = StringProcessing::VectorContains(&parameterNames,&deltaGammaName);
			if( (gamma_i == -1) || (dG_i == -1) )
			{
				cerr << "Cannot Apply " << name << " constraint" << endl;
				to_be_removed.push_back(constraintIndex);
				continue;
			}
			double gamma = NewParameters->GetPhysicsParameter(gammaName)->GetValue();
			double dgam =  NewParameters->GetPhysicsParameter(deltaGammaName)->GetValue();
			double gamma_con = 0.657 ;
			double dgam_con = 0.123 ;
			double v[2] ;
			v[0] = gamma-gamma_con ;
			v[1] = dgam-dgam_con ;
			double E[2][2] ;
			E[0][0] = 7213.6;//7213.6;  //13566.7
			E[1][1] = 1087.3;//1087.3;  //1221.
			E[1][0] = 587.1;//587.1;   //1221.
			E[0][1] = E[1][0] ;
			double chisq = 0;
			for( int ix=0; ix<2; ++ix ) {
				for( int iy=0; iy<2; ++iy ) {
					chisq += v[ix]*E[ix][iy]*v[iy] ;
				}
			}
			constraintValue += chisq ;
		}
		else if ( name == "Gammas-DeltaGamma-1fbPaper" )
		{
			//GammaObs^-1 = Gamma^-1 * ( 1 + (deltaGamma/2*Gamma)^2 ) / ( 1 - (deltaGamma/2*Gamma)^2 )
			// Get gamma and delta gamma
			string gammaName("gamma");
			string deltaGammaName("deltaGamma");
			int gamma_i = StringProcessing::VectorContains(&parameterNames,&gammaName);
			int dG_i = StringProcessing::VectorContains(&parameterNames,&deltaGammaName);
			if( (gamma_i == -1) || (dG_i == -1) )
			{
				cerr << "Cannot Apply " << name << " constraint" << endl;
				to_be_removed.push_back(constraintIndex);
				continue;
			}
			double gamma = NewParameters->GetPhysicsParameter(gammaName)->GetValue();
			double dgam =  NewParameters->GetPhysicsParameter(deltaGammaName)->GetValue();
			double gamma_con = 0.6631 ;
			double dgam_con = 0.100 ;
			double v[2] ;
			v[0] = gamma-gamma_con ;
			v[1] = dgam-dgam_con ;
			double E[2][2] ;
			E[0][0] = 16850.7 ;
			E[1][1] = 3988.85  ;
			E[1][0] = 1904.58  ;
			E[0][1] = E[1][0] ;
			double chisq = 0;
			for( int ix=0; ix<2; ++ix ) {
				for( int iy=0; iy<2; ++iy ) {
					chisq += v[ix]*E[ix][iy]*v[iy] ;
				}
			}
			constraintValue += chisq ;
		}
		else if ( name == "Gammas-DeltaGamma-1fbPaper-combined" )
		{
			//GammaObs^-1 = Gamma^-1 * ( 1 + (deltaGamma/2*Gamma)^2 ) / ( 1 - (deltaGamma/2*Gamma)^2 )
			// Get gamma and delta gamma
			string gammaName("gamma");
			string deltaGammaName("deltaGamma");
			int gamma_i = StringProcessing::VectorContains(&parameterNames,&gammaName);
			int dG_i = StringProcessing::VectorContains(&parameterNames,&deltaGammaName);
			if( (gamma_i == -1) || (dG_i == -1) )
			{
				cerr << "Cannot Apply " << name << " constraint" << endl;
				to_be_removed.push_back(constraintIndex);
				continue;
			}
			double gamma = NewParameters->GetPhysicsParameter(gammaName)->GetValue();
			double dgam =  NewParameters->GetPhysicsParameter(deltaGammaName)->GetValue();
			double gamma_con = 0.661 ;
			double dgam_con = 0.106 ;
			double v[2] ;
			v[0] = gamma-gamma_con ;
			v[1] = dgam-dgam_con ;
			double E[2][2] ;
			E[0][0] = 19273. ;
			E[1][1] = 5895.26  ;
			E[1][0] = -498.83  ;
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
			string AperpName("Aperp_sq");
			string AzeroName("Azero_sq");
			//string AsName("As_sq");
			int Aperp_i = StringProcessing::VectorContains(&parameterNames,&AperpName);
			int Azero_i = StringProcessing::VectorContains(&parameterNames,&AzeroName);
			//int As_i = StringProcessing::VectorContains(&parameterNames,&AsName);
			if( ((Aperp_i == -1) || (Azero_i == -1)) )//|| (As_i == -1) )
			{
				cerr << "Cannot Apply " << name << " constraint" << endl;
				to_be_removed.push_back(constraintIndex);
				continue;
			}
			double Aperpsq = NewParameters->GetPhysicsParameter(AperpName)->GetValue();
			double Azerosq =  NewParameters->GetPhysicsParameter(AzeroName)->GetValue();
			//double Assq =  NewParameters->GetPhysicsParameter(AsName)->GetValue();
			double Atot_val = (*constraintIndex)->GetValue();
			double Atot_constraint = (*constraintIndex)->GetError();
			double excess = (Atot_val-Aperpsq-Azerosq);//-Assq) ;
			double penalty ;
			if( excess >= 0 ) penalty = 0;
			else penalty = (excess*excess) / (Atot_constraint*Atot_constraint) ;
			constraintValue += penalty;
		}
		else if ( name == "GLandGH" )
		{
			// This is a special one to stop GL and GH going negative
			string gammaName("gamma");
			string deltaGammaName("deltaGamma");
			int gamma_i = StringProcessing::VectorContains(&parameterNames,&gammaName);
			int dG_i = StringProcessing::VectorContains(&parameterNames,&deltaGammaName);
			if( (gamma_i == -1) || (dG_i == -1) )
			{
				cerr << "Cannot Apply " << name << " constraint" << endl;
				to_be_removed.push_back(constraintIndex);
				continue;
			}
			double gamma = NewParameters->GetPhysicsParameter("gamma")->GetValue();
			double dgam =  NewParameters->GetPhysicsParameter("deltaGamma")->GetValue();
			double G1 = gamma-dgam/2.0 ;
			double G2 = gamma+dgam/2.0 ;
			//double val = allConstraints[constraintIndex]->GetValue();
			double constraint = (*constraintIndex)->GetError();
			double penalty = 0 ;
			if( G1 < 0. ) penalty += (G1*G1)/(constraint*constraint) ;
			if( G2 < 0. ) penalty += (G2*G2)/(constraint*constraint) ;
			if( penalty > 0. ) cout << " GLandGH constraint being applied " << endl ;
			constraintValue += penalty;
		}
		else if ( name == "CosSqPlusSinSq" )
		{
			// This is a special one to contrain cosphis+sinphis to 1
			string cosphiName("cosphis");
			string sinphiName("sinphis");
			int cos_i = StringProcessing::VectorContains(&parameterNames,&cosphiName);
			int sin_i = StringProcessing::VectorContains(&parameterNames,&sinphiName);
			if( (cos_i == -1) || (sin_i == -1) )
			{
				cerr << "Cannot Apply " << name << " constraint" << endl;
				to_be_removed.push_back(constraintIndex);
				continue;
			}
			double cosphis = NewParameters->GetPhysicsParameter(cosphiName)->GetValue();
			double sinphis =  NewParameters->GetPhysicsParameter(sinphiName)->GetValue();
			double val = (*constraintIndex)->GetValue();
			double constraint = (*constraintIndex)->GetError();
			double excess = (val - cosphis*cosphis - sinphis*sinphis ) ;
			double penalty = (excess*excess) / (constraint*constraint) ;
			constraintValue += penalty;
		}
		else if ( Found_Position[num_i] >= 0 )
		{
			int dummy_i = StringProcessing::VectorContains(&parameterNames,&name);
			if( dummy_i == -1 )
			{
				cerr << "Cannot Apply " << name << " constraint" << endl;
				to_be_removed.push_back(constraintIndex);
				continue;
			}
			//Do standard gaussian constraint calculation
			double parameterValue = NewParameters->GetPhysicsParameter(name)->GetValue();
			double gaussSqrt = ( parameterValue - (*constraintIndex)->GetValue() ) / (*constraintIndex)->GetError();
			constraintValue += gaussSqrt * gaussSqrt;
		}
	}

	for( unsigned int i=0; i< to_be_removed.size(); ++i )
	{
		allConstraints.erase( to_be_removed[i] );
	}

	return 0.5 * constraintValue;
}

void ConstraintFunction::Print() const
{
	cout << "ConstraintFunction:" << endl;
	cout << "This is to be coded up when needed" << endl;
}

string ConstraintFunction::XML() const
{
	stringstream xml;

	for( vector<ExternalConstraint*>::const_iterator const_i = allConstraints.begin(); const_i != allConstraints.end(); ++const_i )
	{
		xml << (*const_i)->XML() << endl;
		xml << endl;
	}

	return xml.str();
}

vector<string> ConstraintFunction::ConstrainedParameter() const
{
	vector<string> constnames;
	for( unsigned int i=0; i< allConstraints.size(); ++i )
	{
		constnames = StringProcessing::CombineUniques( constnames, vector<string>(1,allConstraints[i]->GetName()) );
	}
	return constnames;
}

