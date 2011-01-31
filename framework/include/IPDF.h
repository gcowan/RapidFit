/**
        @interface IPDF

        Common interface for all PDFs

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef IPDF_H
#define IPDF_H

#include "DataPoint.h"
#include "PhaseSpaceBoundary.h"
#include "ParameterSet.h"
#include <vector>
#include <string>

using namespace std;

class IPDF
{
	public:
		//Indicate whether the function has been set up correctly
		virtual bool IsValid() = 0;

		//Set the function parameters
		virtual bool SetPhysicsParameters( ParameterSet* ) = 0;

		//Return the integral of the function over the given boundary
		virtual double Integral( DataPoint*, PhaseSpaceBoundary* ) = 0;

		//Return the function value at the given point
		virtual double Evaluate( DataPoint* ) = 0;

		//Return the components of the function value at the given point
		virtual vector<double> EvaluateComponents( DataPoint* ) = 0;
	
		//Return a prototype data point
		virtual vector<string> GetPrototypeDataPoint() = 0;

		//Return a prototype set of physics parameters
		virtual vector<string> GetPrototypeParameterSet() = 0;

		//Return a list of parameters not to be integrated
		virtual vector<string> GetDoNotIntegrateList() = 0;

		//Update the integral cache
		virtual void UpdateIntegralCache() = 0;

		//  Get a pointer to the seed function
		//  Using a pointer so we have one seed per normal study
		TRandom3 * GetRandomFunction()
		{
			if( seed_function.empty() )
			{
//				cout << "Seed not found in PDF, generating TRandom3(0) random function." << endl;
				seed_function.push_back( new TRandom3(0) );
			}
			return seed_function.back();
		}

		//  Set the Random Generator to be some externally defined instance
		void SetRandomFunction( TRandom3 * new_function )
		{
			while( !seed_function.empty() ){  seed_function.pop_back();  }  //Get rid of old seed(s)

			seed_function.push_back(  new_function  );       //Insert reference to new seed
		}

			//  Seed the Random Number Generator correctly
		void SetRandomFunction( int new_seed )
		{
			while( !seed_function.empty() ){  seed_function.pop_back();  }    //Get rid of old seed(s)

			seed_function.push_back(  new TRandom3( new_seed ) ); //Insert reference to new seed

			while( !seed_num.empty() ){  seed_num.pop_back();  }  //  Get rid of old seed(s)

			seed_num.push_back( new_seed );             //  Store the seed for internal reference
		}

		//  Return the numerical seed
		int GetSeedNum()
		{
			if( !seed_num.empty() )return seed_num.back();
			else return 0;
		}

		private:

		vector<TRandom3 *> seed_function;
		vector<int> seed_num;

};

#endif
