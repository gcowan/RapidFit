/**
        @class RootFileDataSet

        A collection of data points stored in a Root NTuple

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	ROOT Headers
#include "TFile.h"
#include "TKey.h"
#include "TLeaf.h"
#include "Rtypes.h"
//	RapidFit Headers
#include "RootFileDataSet.h"
//	System Headers
#include <iostream>

//Default constructor
RootFileDataSet::RootFileDataSet()
{
}
//Contructor with file name as argument
RootFileDataSet::RootFileDataSet( string FileName, PhaseSpaceBoundary * Boundary ) :  dataBoundary( Boundary ), branches( Boundary->GetAllNames().size() ), observableValues( Boundary->GetAllNames().size() )
{
	//Initialise the output data point
	vector<string> observableNames = Boundary->GetAllNames();
        outputDataPoint = new DataPoint(observableNames);

	//Find the first NTuple in the file
	inputFile = new TFile( FileName.c_str(), "UPDATE" );
	TKey * testKey;
	TIter fileKeys( inputFile->GetListOfKeys() );
	while ( (testKey = (TKey*)fileKeys()) )
	{
		string className = testKey->GetClassName();
		if ( className == "TNtuple" )
		{
			cout << "Found TNtuple: " << testKey->GetName() << endl;
			rootNTuple = (TNtuple*)inputFile->Get( testKey->GetName() );

			if ( CheckTNtupleWithBoundary( rootNTuple, dataBoundary ) )
			{
				SetBranches();
				return;
			}
		}
	}

	//If you get this far, no compatible ntuple found
	cerr << "No compatible ntuple found, did you give the correct path?" << endl;
}

//Constructor that specifies an NTuple
RootFileDataSet::RootFileDataSet( string FileName, string TuplePath, PhaseSpaceBoundary * Boundary ) : dataBoundary( Boundary )
{
	//Initialise the output data point
	vector<string> observableNames = Boundary->GetAllNames();
        outputDataPoint = new DataPoint(observableNames);

	//Create the objects to access the file
	inputFile = new TFile( FileName.c_str(), "UPDATE" );
	rootNTuple = (TNtuple*)( inputFile->Get( TuplePath.c_str() ) );

	//Check the NTuple found is in the expected format
	if ( CheckTNtupleWithBoundary( rootNTuple, dataBoundary ) )
	{
		SetBranches();
	}
	else
	{
		cerr << "NTuple is incompatible with boundary, do you have all the columns in this ntuple?" << endl;
	}
}

//Destructor
RootFileDataSet::~RootFileDataSet()
{
	inputFile->Close();
	delete inputFile;
	delete rootNTuple;
}

// set up the branches which will be used to read in the data
bool RootFileDataSet::SetBranches()
{
	cout << "Setting branches" << endl;
	vector<string> observableNames = dataBoundary->GetAllNames();
        for ( unsigned short int obsIndex = 0; obsIndex < observableNames.size(); ++obsIndex )
        {
		rootNTuple->SetBranchAddress(observableNames[obsIndex].c_str(), &(observableValues[obsIndex]), &(branches[obsIndex]));
	}

	return true;
}

//Check if an NTuple is compatible with a boundary
bool RootFileDataSet::CheckTNtupleWithBoundary( TNtuple * TestTuple, PhaseSpaceBoundary * TestBoundary )
{
	bool compatible = true;

	vector<string> allNames = TestBoundary->GetAllNames();
	string lastName;
	for ( unsigned short int observableIndex = 0; observableIndex < allNames.size(); ++observableIndex )
	{
		bool found = false;
		//Find the requested observable name in the ntuple
		TIter observableIterator( TestTuple->GetListOfLeaves() );
		TLeaf * observableLeaf;
		while ( ( observableLeaf = (TLeaf*)observableIterator() ) )
		{
			string leafName = observableLeaf->GetName();
			lastName = leafName;
			if ( leafName == allNames[observableIndex] )
			{
				found = true;
				break;
			}
		}

		//If the name isn't present, the ntuple is incompatible
		if ( !found )
		{
			compatible = false;
			cerr << "\n" << allNames[observableIndex] << "\tNOT FOUND in Ntuple!!\n" << endl;
			break;
		}
	}

	return output;
}

//Return a data point from the NTuple
DataPoint * RootFileDataSet::GetDataPoint( int DataIndex )
{
	if ( DataIndex < rootNTuple->GetEntries() )
	{
		rootNTuple->GetEntry( DataIndex );	

		//Delete the old reference data point and make a new one
		vector<string> observableNames = dataBoundary->GetAllNames();
		delete outputDataPoint;
		outputDataPoint = new DataPoint(observableNames);

		//Populate the output data point
		for ( unsigned short int observableIndex = 0; observableIndex < observableNames.size(); ++observableIndex )
		{
			string name = observableNames[observableIndex];
		        string unit = dataBoundary->GetConstraint(name)->GetUnit();
			outputDataPoint->SetObservable( name, observableValues[observableIndex], 0.0, unit );
		}
	}
	else
	{
		cerr << "Requested data point index (" << DataIndex << ") is out of range (" << rootNTuple->GetEntries() << "). DataPoint reference not updated" << endl;
	}

	return outputDataPoint;
}

//Add a data point to the NTuple
bool RootFileDataSet::AddDataPoint( DataPoint * NewDataPoint )
{
	if ( dataBoundary->IsPointInBoundary(NewDataPoint) )
	{
		//Make an array of the observable values
		vector<string> observableNames = dataBoundary->GetAllNames();
		Float_t* dataPointArray = new Float_t[ observableNames.size() ];
		for ( unsigned short int observableIndex = 0; observableIndex < observableNames.size(); ++observableIndex )
		{
			string name = observableNames[observableIndex];
			Observable * inputObservable = NewDataPoint->GetObservable(name);
			dataPointArray[observableIndex] = Float_t(inputObservable->GetValue());
		}

		//Put the array into the NTuple
		rootNTuple->Fill(dataPointArray);
		return true;
	}
	else
	{
		cerr << "DataPoint is incompatible with NTuple" << endl;
		return false;
	}
}

//Return the number of data points
int RootFileDataSet::GetDataNumber()
{
	return int(rootNTuple->GetEntries());
}

//Return the phase space boundary
PhaseSpaceBoundary * RootFileDataSet::GetBoundary()
{
	return dataBoundary;
}
