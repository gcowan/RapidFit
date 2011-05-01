
//	ROOT Headers
#include "TString.h"
#include "TRandom3.h"
//	RapidFit Headers
#include "IPDF.h"
//	System Headers
#include <string.h>

IPDF::IPDF() : cached_files(), stored_ID(), hasCachedMCGenerator(false), seed_function(), seed_num()
{
}

//  Get a pointer to the seed function
//  Using a pointer so we have one seed per normal study
TRandom3 * IPDF::GetRandomFunction()
{
	if( seed_function.empty() )
	{
		//				cout << "Seed not found in PDF, generating TRandom3(0) random function." << endl;
		seed_function.push_back( new TRandom3(0) );
	}
	return seed_function.back();
}

//  Set the Random Generator to be some externally defined instance
void IPDF::SetRandomFunction( TRandom3 * new_function )
{
	while( !seed_function.empty() ){  seed_function.pop_back();  }  //Get rid of old seed(s)
	seed_function.push_back(  new_function  );       //Insert reference to new seed
}

//  Seed the Random Number Generator correctly
void IPDF::SetRandomFunction( int new_seed )
{
	while( !seed_function.empty() ){  seed_function.pop_back();  }    //Get rid of old seed(s)

	seed_function.push_back(  new TRandom3( unsigned(new_seed) ) ); //Insert reference to new seed

	while( !seed_num.empty() ){  seed_num.pop_back();  }  //  Get rid of old seed(s)

	seed_num.push_back( new_seed );             //  Store the seed for internal reference
}

//  Return the numerical seed
int IPDF::GetSeedNum()
{
	if( !seed_num.empty() )return seed_num.back();
	else return 0;
}

//	Set the ID of the PDF (Used internally as a way of identifiying this PDF
void IPDF::SET_ID( TString id_string )
{
	stored_ID = id_string.Data();
}
void IPDF::SET_ID( string id_string )
{
	stored_ID = id_string;
}

//	Get the ID of this particular PDF
string IPDF::GET_ID()
{
	return stored_ID;
}

//	Set the Status of a cache for the MC generator associated with this PDF
void IPDF::SetMCCacheStatus( bool newStatus)	
{
	if( (newStatus == false) && hasCachedMCGenerator )
	{
		cout << GET_ID() << ":\tRemoving Cache" << endl;
		Remove_Cache();
	}
	hasCachedMCGenerator = newStatus;
}

void IPDF::Remove_Cache()
{
	while( !cached_files.empty() )
	{
		remove ( cached_files.back().c_str() );
		cached_files.pop_back();
	}
}

//	Get the Status of the MC generator for this PDF
bool IPDF::GetMCCacheStatus()
{
	return hasCachedMCGenerator;
}

//	Add an on-disk object which exists for the lifetime of this PDF
void IPDF::AddCacheObject( TString obj_name )
{
	cached_files.push_back( obj_name.Data() );
}
void IPDF::AddCacheObject( string obj_name )
{
	cached_files.push_back( obj_name );
}
