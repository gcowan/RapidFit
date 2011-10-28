/**
        @class BasePDF

        Class that provides a general implementation of IPDF.
        Can inherit from this to make a PDF without worrying about the details.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "BasePDF.h"
#include "ObservableRef.h"
//	System Headers
#include <iostream>

//Constructor
BasePDF::BasePDF() : cachedIntegral(-1.0), cacheValid(false), allParameters(), allObservables(), valid(false), observables(),
	cached_files(), stored_ID(), hasCachedMCGenerator(false), seed_function(), seed_num(), PDFName("Base")
{
}

//Destructor
BasePDF::~BasePDF()
{
	Remove_Cache();
}

//Indicate whether the function has been set up correctly
bool BasePDF::IsValid()
{
	return valid;
}

//Set the function parameters
bool BasePDF::SetPhysicsParameters(ParameterSet * NewParameterSet)
{
	bool success = allParameters.SetPhysicsParameters(NewParameterSet);

	//Invalidate the cache
	cacheValid = false;
	cachedIntegral = -1.0;

	return success;
}

//Return the integral of the function over the given boundary
double BasePDF::Integral(DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary)
{
	if( observables.size() != this->GetPrototypeDataPoint().size() )
	{
		observables=ObservableRef( this->GetPrototypeDataPoint() );
	}
	//Check the boundary is within the correct phase space
	for( unsigned int nameIndex=0; nameIndex < allObservables.size(); ++nameIndex )
	{
		IConstraint * testConstraint = NewBoundary->GetConstraint( observables[nameIndex] );
		if (testConstraint->GetUnit() == "NameNotFoundError")
		{
			cerr << "PDF cannot integrate over phase space: observable \"" << observables[nameIndex].Name().c_str() << "\" not found" << endl;
			return -1.0;
		}
	}

	//Use the cached integral, so long as it makes sense
	if (cacheValid)
	{
		return cachedIntegral;
	}
	else
	{
		//If the PDF normalisation does not depend on the DataPoint, cache the value to save time
		cachedIntegral = this->Normalisation(NewBoundary);

		//Check the cache is valid
		if ( cachedIntegral > 0.0 )
		{
			cacheValid = true;
			return cachedIntegral;
		}
		else
		{
			//Integral must depend on the DataPoint
			return this->Normalisation(NewDataPoint, NewBoundary);
		}
	}
	return -1.;
}

//Do the integration
double BasePDF::Normalisation(PhaseSpaceBoundary * NewBoundary)
{
	//	Stupid gcc
	(void)NewBoundary;
	//Just a default value
	return -1.0;
}

//Do the integration
double BasePDF::Normalisation(DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary)
{
	//	Stupid gcc
	(void)NewDataPoint;
	(void)NewBoundary;
	//Just a default value
	return -1.0;
}

/*//Return the function value at the given point
double BasePDF::Evaluate(DataPoint * NewDataPoint)
{
	//Check the datapoint is within the PDF boundary
	if ( xRange.Evaluate(NewDataPoint) > 0.0 )
	{
		return  Value(NewDataPoint);
	}
	else
	{
		cerr << "Datapoint is outside PDF boundary" << endl;
		return 0.0;
	}
}*/

//Calculate the function value
double BasePDF::Evaluate(DataPoint * NewDataPoint)
{
	//	Stupid gcc
	(void)NewDataPoint;
	//Just a default value
	return  -1.0;
}

//Calculate the function value for numerical integration
double BasePDF::EvaluateForNumericIntegral(DataPoint * NewDataPoint)
{
	return this->Evaluate( NewDataPoint ) ;
}

//Calculate the function value
vector<double> BasePDF::EvaluateComponents(DataPoint * NewDataPoint)
{
	//Just assume a single component equal to the ordinary evaluate method
	vector<double> components ;
	components.push_back( this->Evaluate( NewDataPoint ) ) ;
	return  components ;
}

//Return a prototype data point
vector<string> BasePDF::GetPrototypeDataPoint()
{
	return allObservables;
}

//Return a prototype set of physics parameters
vector<string> BasePDF::GetPrototypeParameterSet()
{
	return allParameters.GetAllNames();
}

ParameterSet* BasePDF::GetActualParameterSet()
{
	return &allParameters;
}

//Return a list of parameters not to be integrated
vector<string> BasePDF::GetDoNotIntegrateList()
{
	vector<string> emptyVector;
	return emptyVector;
}

// Update integral cache
void BasePDF::UpdateIntegralCache()
{
}


//  Get a pointer to the seed function
//  Using a pointer so we have one seed per normal study
TRandom3 * BasePDF::GetRandomFunction()
{
	if( seed_function.empty() )
	{
		//				cout << "Seed not found in PDF, generating TRandom3(0) random function." << endl;
		seed_function.push_back( new TRandom3(0) );
	}
	return seed_function.back();
}

//  Set the Random Generator to be some externally defined instance
void BasePDF::SetRandomFunction( TRandom3 * new_function )
{
	while( !seed_function.empty() ){  seed_function.pop_back();  }  //Get rid of old seed(s)
	seed_function.push_back(  new_function  );       //Insert reference to new seed
}

//  Seed the Random Number Generator correctly
void BasePDF::SetRandomFunction( int new_seed )
{
	while( !seed_function.empty() ){  seed_function.pop_back();  }    //Get rid of old seed(s)

	seed_function.push_back(  new TRandom3( unsigned(new_seed) ) ); //Insert reference to new seed

	while( !seed_num.empty() ){  seed_num.pop_back();  }  //  Get rid of old seed(s)

	seed_num.push_back( new_seed );             //  Store the seed for internal reference
}

//  Return the numerical seed
int BasePDF::GetSeedNum()
{
	if( !seed_num.empty() )return seed_num.back();
	else return 0;
}

//	Set the ID of the PDF (Used internally as a way of identifiying this PDF
void BasePDF::SET_ID( TString id_string )
{
	stored_ID = id_string.Data();
}
void BasePDF::SET_ID( string id_string )
{
	stored_ID = id_string;
}

//	Get the ID of this particular PDF
string BasePDF::GET_ID()
{
	return stored_ID;
}

//	Set the Status of a cache for the MC generator associated with this PDF
void BasePDF::SetMCCacheStatus( bool newStatus)	
{
	if( (newStatus == false) && hasCachedMCGenerator )
	{
		//cout << GET_Name() << ":\tRemoving Cache" << endl;
		Remove_Cache();
	}
	hasCachedMCGenerator = newStatus;
}

void BasePDF::Remove_Cache( bool leave_on_disk )
{
	while( !cached_files.empty() )
	{
		if( !leave_on_disk ) remove ( cached_files.back().c_str() );
		cached_files.pop_back();
	}
}

//	Get the Status of the MC generator for this PDF
bool BasePDF::GetMCCacheStatus()
{
	return hasCachedMCGenerator;
}

//	Add an on-disk object which exists for the lifetime of this PDF
void BasePDF::AddCacheObject( TString obj_name )
{
	cached_files.push_back( obj_name.Data() );
}
void BasePDF::AddCacheObject( string obj_name )
{
	cached_files.push_back( obj_name );
}

string BasePDF::GetName()
{
	return PDFName;
}

void BasePDF::SetName( string input )
{
	PDFName = input;
	TString Random;	Random+=this->GetRandomFunction()->Rndm()*1000;
	string new_name = input+"_"+Random.Data();
	this->SET_ID( new_name );
}

