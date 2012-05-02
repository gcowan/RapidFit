
//	RapidFit Headers
#include "FeldmanCousins_Study.h"
#include "FitResultVector.h"
#include "OutputConfiguration.h"
#include "MinimiserConfiguration.h"
#include "FitFunctionConfiguration.h"
#include "XMLConfigReader.h"
#include "PDFWithData.h"
#include "StringProcessing.h"
#include "FitAssembler.h"
#include "ResultFormatter.h"
//	System Headers
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace::std;

//	This is the Main loop of Code which performs the Numerical work needed for the Feldman-Cousins Fit in RapidFit
//	The output of this code is to be analysed by Conor's plotting tool
//	I make no claims as to how FC plots are extracted from this I will just tell you what this code Does
//
//	Step1:		Find Global Minima
//	Step2:		Find Minima at point i,j
//	Step3:		Setup Framework to, and Generate Data for toy study with control paramater for i,j Fixed
//			You can either set nuisence parameters to their central values when fit at Step 1 or 2
//			At the time of writing the nuisence parameters are still under investigation
//	Step4:		Repeat Step 3 as many times as requested, or do 100 toys as default
//	Step5:		Format output into standard block format. This gives us 2 things.
FitResultVector* FeldmanCousins_Study::FeldmanCousins( FitResultVector* GlobalResult, FitResultVector* _2DResultForFC, vector<unsigned int> numberRepeats, unsigned int NuisenceModel, bool FC_Debug_Flag, OutputConfiguration* makeOutput, MinimiserConfiguration* theMinimiser, FitFunctionConfiguration* theFunction, XMLConfigReader* xmlFile, vector< PDFWithData* > pdfsAndData, int OutputLevel )
{
	double Random_Seed = pdfsAndData[0]->GetPDF()->GetRandomFunction()->Rndm();
	FitResult* GlobalFitResult = GlobalResult->GetFitResult( 0 );
	vector<FitResultVector*> AllResults;

	//		Want to loop over all points on a 2D Plot that have been allocated to this instance
	for( int iFC=0; iFC < _2DResultForFC->NumberResults(); ++iFC )
	{
		TRandom3* new_rand = new TRandom3(UInt_t(Random_Seed));
		//		GET INPUT Data from fit Results
		//
		//  Get a ParameterSet that contains the input parameters from the output of a fit
		vector<pair<string, string> > _2DLLscanList = makeOutput->Get2DScanList();
		string name1 = _2DLLscanList[0].first;
		string name2 = _2DLLscanList[0].second;

		ParameterSet* InputParamSet = NULL;
		ParameterSet* InputFreeSet = NULL;
		ParameterSet* ControlParamSet = NULL;

		//	True by definition:
		double lim1 = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetResultParameter( name1 )->GetValue();
		double lim2 = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetResultParameter( name2 )->GetValue();

		if( NuisenceModel == 1 )
		{
			//		Use the inputs from Step 1
			InputParamSet = GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet();
			InputFreeSet = GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet();
			ControlParamSet = GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet();
		}
		else if( NuisenceModel == 2 )
		{
			//		Use the inputs from Step 2
			InputParamSet = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet();
			InputFreeSet = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet();
			ControlParamSet = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet();
		}


		//		Just for clarity
		//		also when using values from Step1 these are not set correctly for the study
		ControlParamSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
		ControlParamSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
		InputParamSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
		InputParamSet->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
		InputParamSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
		InputParamSet->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
		InputParamSet->GetPhysicsParameter( name1 )->SetType( "Fixed" );
		InputParamSet->GetPhysicsParameter( name2 )->SetType( "Fixed" );
		InputFreeSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
		InputFreeSet->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
		InputFreeSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
		InputFreeSet->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
		InputFreeSet->GetPhysicsParameter( name1 )->SetType( "Free" );
		InputFreeSet->GetPhysicsParameter( name2 )->SetType( "Free" );

		cout << "Here" << endl;

		//	Collect all of the relevent Data from the XML
		//	Note: most of these had to be written for FCscans
		MinimiserConfiguration * ToyStudyMinimiser = theMinimiser;
		theMinimiser->SetOutputLevel(-999);
		FitFunctionConfiguration* ToyStudyFunction = theFunction;
		vector<PhaseSpaceBoundary*> PhaseSpaceForToys = xmlFile->GetPhaseSpaceBoundaries();
		vector<ConstraintFunction*> ConstraintsForToys = xmlFile->GetConstraints();

		//	Read the number of events from the existing DataSets
		//	I think there needs be an equivalent function drafted to use sWeights
		vector<unsigned int> EventsPerPDF;
		vector<double> sweight_error;

		for( unsigned short int pdf_num=0; pdf_num < pdfsAndData.size(); ++pdf_num )
		{
			string sWeightObs="Nsig_sw";
			bool sWeighted=false;
			vector<string> phase_names = PhaseSpaceForToys[pdf_num]->GetAllNames();
			int sweight_num = StringProcessing::VectorContains( &phase_names, &sWeightObs );
			if( sweight_num != -1 ) sWeighted=true;
			if( sWeighted )
			{
				for(unsigned short int i=0; i< PhaseSpaceForToys.size(); ++i )
				{
					vector<double> new_constraint(1,1.0);
					PhaseSpaceForToys[i]->SetConstraint( "Nsig_sw", new_constraint, " " );
				}
				IDataSet* file_input = pdfsAndData[pdf_num]->GetDataSet();
				int point_number = file_input->GetDataNumber();
				double data_num=0; double data_err=0;
				for( int i=0; i < point_number; ++i)
				{
					double sweight_val = file_input->GetDataPoint( i )->GetObservable(sWeightObs)->GetValue();
					data_num+=sweight_val;
					data_err+=sweight_val*sweight_val;
				}
				data_err=sqrt(data_err);
				sweight_error.push_back( unsigned(int(data_err)) );
				EventsPerPDF.push_back( unsigned(int(data_num)) );
			}
			else	{
				EventsPerPDF.push_back( unsigned(pdfsAndData[pdf_num]->GetDataSet()->GetDataNumber()) );
				sweight_error.push_back( 0. );
			}
		}



		//		Number of datasets wanted i.e. Number of Toys
		unsigned int wanted_number_of_toys = 100;
		if( !numberRepeats.empty() ) wanted_number_of_toys = unsigned(numberRepeats[0]);



		//		GENERATE DATA
		//	We want to now generate a new DataSet for a Toy Study
		//	I choose a Foam DataSet with no Cuts or extra arguments as these are not required
		//
		//	Data is now only generated and held in memory for as long as it's needed as this
		//	has a lighter memory footprint
		//
		//	Generate Once, Fit twice
		//	This was a bit of a pain to code up.
		//	Although I like the idea of coding it up into a GenerateToyData object in future to avoid this

		vector<PDFWithData*> PDFsWithDataForToys;
		vector<vector<IDataSet*> > Memory_Data;
		Memory_Data.resize( pdfsAndData.size() );

		//	NB:
		//	Generating data the first time is more complex as we need to construct the correct PDFWithData
		//	The object setup here will be used later in the code as a template with cached data being replaced
		//	This reduces the amount of code and makes this procedure easier to understand as a whole

		cout << "\n\n\t\tGenerating Data For First Toy\n" << endl;
		//	We may have multiple PDFs in the XML
		for( unsigned short int pdf_num=0; pdf_num < pdfsAndData.size(); ++pdf_num )
		{
			//	Generate the Data for this Study using the given the XML PhaseSpace and PDF
			IPDF* PDF_from_XML = pdfsAndData[pdf_num]->GetPDF();

			//	Create a DataSetConfiguration to make the datasets
			vector<string> empty_args;
			vector<string> empty_arg_names;
			unsigned int wanted_events = EventsPerPDF[pdf_num];
			if( sweight_error[pdf_num] > 0. )  wanted_events *= unsigned(int(new_rand->Gaus()*sweight_error[pdf_num]));
			DataSetConfiguration* Toy_Foam_DataSet= new DataSetConfiguration( "Foam", wanted_events, "", empty_args, empty_arg_names, PDF_from_XML );

			//	Set the Input Physics Parameters to be as defined above ControlSet
			//	This is to seperate changing bottles for Physics from that used at Generation
			//	as we want the data to have a wider scope than 1 fit with the data
			Toy_Foam_DataSet->SetPhysicsParameters( ControlParamSet );

			//	Let the user know what's going on
			cout << "generating data for pdf: " << pdf_num << "\n";

			//	Make the data
			IDataSet* new_dataset = Toy_Foam_DataSet->MakeDataSet( PhaseSpaceForToys[pdf_num], PDF_from_XML );
			//	Store the Data Object in Memory
			Memory_Data[pdf_num].push_back( new_dataset );

			//	Store the information for each PDFWithData
			vector<DataSetConfiguration*> DataSetConfigForToys;
			DataSetConfigForToys.push_back( Toy_Foam_DataSet );
			PDFWithData* ToyPDFWithData = new PDFWithData( PDF_from_XML, PhaseSpaceForToys[pdf_num], DataSetConfigForToys );

			//	Store this PDFWithData
			ToyPDFWithData->SetPhysicsParameters( InputParamSet );
			PDFsWithDataForToys.push_back( ToyPDFWithData );
		}

		//		Result Vectors for the Fit for clarity and Output Formatting
		FitResultVector* study1Results = new FitResultVector( GlobalFitResult->GetResultParameterSet()->GetAllNames() );
		FitResultVector* study2Results = new FitResultVector( GlobalFitResult->GetResultParameterSet()->GetAllNames() );

		//	This Forms the main sequence of the fit

		cout << "\n\n\t\tPerforming Fits to Toys in FC\n"<<endl;
		//		Perform Fits
		//	This will record ONLY Data from working fits i.e. FitStatus==3 unless the Debug flag is on
		for( unsigned short int dataset_num=0; dataset_num < wanted_number_of_toys; ++dataset_num )
		{

			bool toy_failed = false;

			FitResult* fit1Result = NULL;
			FitResult* fit2Result = NULL;

			ParameterSet* LocalInputFreeSet = NULL;
			ParameterSet* LocalInputFixedSet = NULL;

			if( NuisenceModel == 1 )
			{
				//	Assuming Nuisence Parameters set to Global Minima
				//	We need to set this for EVERY FIT in order to have correct generation/pull values
				LocalInputFreeSet = GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet();
				LocalInputFixedSet = GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet();
			}
			else if( NuisenceModel == 2 )
			{
				//	Assuming Nuisence Parameters set to Local Minima
				//
				LocalInputFreeSet = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet();
				LocalInputFixedSet = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet();
			}
			else if( NuisenceModel == 3 )
			{
				//	Generate a ParameterSet with random numbers based on
				//	Fit Results and assign it to the ControlParamSet which is used for generation
				//	Also need to Setup Local ParameterSets as with standard case
			}

			ControlParamSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
			ControlParamSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
			LocalInputFreeSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
			LocalInputFreeSet->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
			LocalInputFreeSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
			LocalInputFreeSet->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
			LocalInputFreeSet->GetPhysicsParameter( name1 )->SetType( "Free" );
			LocalInputFreeSet->GetPhysicsParameter( name2 )->SetType( "Free" );

			LocalInputFixedSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
			LocalInputFixedSet->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
			LocalInputFixedSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
			LocalInputFixedSet->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
			LocalInputFixedSet->GetPhysicsParameter( name1 )->SetType( "Fixed" );
			LocalInputFixedSet->GetPhysicsParameter( name2 )->SetType( "Fixed" );

			//	We need to set some factors before we perform the fit
			for( unsigned short int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); ++pdf_num )
			{
				vector<IDataSet*> wanted_set;
				wanted_set.push_back( Memory_Data[pdf_num][dataset_num] );
				PDFsWithDataForToys[pdf_num]->AddCachedData( wanted_set );
				PDFsWithDataForToys[pdf_num]->SetPhysicsParameters( LocalInputFreeSet );
			}

			cout << "\n\n\t\tPerforming Fit To Toy: "<< (dataset_num+1) <<" of " << wanted_number_of_toys << endl<<endl;

			streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;
			cout << "\n\tStarting Fit:" <<endl;
			//	If the user wanted silence we point the Std Output Streams to /dev/null
			if( OutputLevel <= -1 )
			{
				cout_bak = cout.rdbuf();
				cerr_bak = cerr.rdbuf();
				clog_bak = clog.rdbuf();
				cout.rdbuf(0);
				cerr.rdbuf(0);
				clog.rdbuf(0);
			}
			//	Fit once with control parameters Free
			fit1Result = FitAssembler::DoSafeFit( ToyStudyMinimiser, ToyStudyFunction, LocalInputFreeSet, PDFsWithDataForToys, ConstraintsForToys, OutputLevel );
			//	Reset Std Output Streams
			if( OutputLevel <= -1 )
			{
				cout.rdbuf(cout_bak);
				cerr.rdbuf(cerr_bak);
				clog.rdbuf(clog_bak);
			}
			cout << "\tFit Finished!\n" <<endl;

			//	Only Fit again to this dataset if it fits well with +2 dof
			//	This has the obvious savings in CPU resources
			//	We may want to keep this information so leaving it configurable.
			if( ( fit1Result->GetFitStatus() == 3 ) || FC_Debug_Flag )
			{
				if( FC_Debug_Flag )	cout << "\n\nYou are aware your requesting all output?\n"<<endl;

				ResultFormatter::ReviewOutput( fit1Result );

				//  Fit secondGlobalResult with control parameters Fixed
				for( unsigned short int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); ++pdf_num )
				{
					vector<IDataSet*> wanted_set;
					wanted_set.push_back( Memory_Data[pdf_num][dataset_num] );
					PDFsWithDataForToys[pdf_num]->AddCachedData( wanted_set );
					PDFsWithDataForToys[pdf_num]->SetPhysicsParameters( LocalInputFixedSet );
				}

				cout << "\n\t\tFirst Fit Successful, Performing the Second Fit " << (dataset_num+1) << " of " << wanted_number_of_toys <<endl;

				cout << "\n\tStarting Fit:" <<endl;
				//	If the user wanted silence we point the Std Output Streams to /dev/null
				if( OutputLevel <= -1 )
				{
					cout_bak = cout.rdbuf();
					cerr_bak = cerr.rdbuf();
					clog_bak = clog.rdbuf();
					cout.rdbuf(0);
					cerr.rdbuf(0);
					clog.rdbuf(0);
				}
				//	Use the SafeFit as this always returns something when a PDF has been written to throw not exit
				fit2Result = FitAssembler::DoSafeFit( ToyStudyMinimiser, ToyStudyFunction, LocalInputFixedSet, PDFsWithDataForToys, ConstraintsForToys, OutputLevel );
				//	Reset Std Output Streams
				if( OutputLevel <= -1 )
				{
					cout.rdbuf(cout_bak);
					cerr.rdbuf(cerr_bak);
					clog.rdbuf(clog_bak);
				}
				cout << "\tFit Finished!\n" <<endl;

				//	If either Fit Failed we want to 'dump the results' and run an extra Fit.
				if( (fit1Result->GetFitStatus() != 3) || (fit2Result->GetFitStatus() != 3) ) toy_failed = true;
			}
			else	toy_failed = false;

			if( !toy_failed || FC_Debug_Flag )
				ResultFormatter::ReviewOutput( fit2Result );

			//	Do we want to store the Data OR run another toy to get a better Fit
			if( toy_failed )
			{
				cerr << "\n\t\tA Single Toy Study Failed... Requesting More Data for another pass.\n" << endl;
				//  Increment counter so we can guarantee we get 'wanted_number_of_toys' which have fitted correctly
				++wanted_number_of_toys;
			}

			//	In my opinion the fact this is so easy and doesn't cause root to throw up is a testiment to the other authors in RapidFit :D
			cout << "\n\n\t\tDeleting Used Data:\n"<<endl;
			for( unsigned int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); ++pdf_num )
			{
				cout << "deleting data for pdf: " << pdf_num << "\n";
				//delete Memory_Data[pdf_num].back();
				Memory_Data[pdf_num].back() = NULL;
			}

			// Only generate data if I'm going to fit to it
			if( dataset_num < (wanted_number_of_toys-1) )
			{
				cout << "\n\n\t\tGenerating Data For Next Toy:\n" <<endl;

				for( unsigned short int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); ++pdf_num )
				{
					// See above for more detailed description
					cout << "generating data for pdf: " << pdf_num << "\n";
					PDFsWithDataForToys[pdf_num]->SetPhysicsParameters( ControlParamSet );
					unsigned int wanted_events = EventsPerPDF[pdf_num];
					if( sweight_error[pdf_num] > 0. )  wanted_events *= unsigned(int(new_rand->Gaus()*sweight_error[pdf_num]));
					IDataSet* new_dataset = PDFsWithDataForToys[pdf_num]->GetDataSetConfig()->MakeDataSet( PhaseSpaceForToys[pdf_num], PDFsWithDataForToys[pdf_num]->GetPDF(), int(wanted_events) );
					Memory_Data[pdf_num].push_back( new_dataset );
				}
			}

			if( FC_Debug_Flag || !toy_failed ){
				if( toy_failed )
				{
					fit1Result->ForceFitStatus(-2);
					fit2Result->ForceFitStatus(-2);
				}
				study1Results->AddFitResult( fit1Result );
				study2Results->AddFitResult( fit2Result );
			}
		}


		//			STEP 6
		//
		//		Standard Output Format which makes the results the same running either
		//		the whole scan on one machine or running the whole set on a batch system
		FitResultVector* ThisStudy = new FitResultVector( GlobalFitResult->GetResultParameterSet()->GetAllNames() );
		ThisStudy->AddFitResult( _2DResultForFC->GetFitResult( iFC ), false );
		ThisStudy->AddCPUTimes( _2DResultForFC->GetAllCPUTimes() );
		ThisStudy->AddRealTimes( _2DResultForFC->GetAllRealTimes() );
		//	The Generated Value for the Global and Local fit are best defined as -9999 as a sensible default
		for( unsigned short int num=0; num < GlobalFitResult->GetResultParameterSet()->GetAllNames().size(); ++num )
		{
			string name = GlobalFitResult->GetResultParameterSet()->GetAllNames()[num];
			GlobalFitResult->GetResultParameterSet()->GetResultParameter( name )->ForcePullValue( -9999 );
			GlobalFitResult->GetResultParameterSet()->GetResultParameter( name )->ForceOriginalValue( -9999 );
			ThisStudy->GetFitResult(0)->GetResultParameterSet()->GetResultParameter( name )->ForcePullValue( -9999 );
			ThisStudy->GetFitResult(0)->GetResultParameterSet()->GetResultParameter( name )->ForceOriginalValue( -9999 );
		}
		AllResults.push_back( GlobalResult );
		AllResults.push_back( ThisStudy );
		AllResults.push_back( study1Results );
		AllResults.push_back( study2Results );

		delete new_rand;
	}

	cout << "\n\nNumerical part of FeldMan-Cousins Scan Complete,\n\n Now process the data through the plotting tool :D\n\n";
	return new FitResultVector( AllResults );
}

