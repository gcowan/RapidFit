
#include "VectoredFeldmanCousins.h"
#include "ClassLookUp.h"
#include "FitAssembler.h"
#include "ParameterSet.h"
#include <iomanip>

using namespace::std;

VectoredFeldmanCousins::VectoredFeldmanCousins( FitResultVector* input_GlobalResult, FitResultVector* ResultsForFC, unsigned int inputNuisenceModel, OutputConfiguration* new_makeOutput, MinimiserConfiguration* newMinimiser, FitFunctionConfiguration* newFunction, XMLConfigReader* new_xmlFile, vector< PDFWithData* > new_pdfsAndData ) : 
	GlobalFitResult(), GlobalFitPhysicsParameters(), FitAtGridPoints(), cout_bak(NULL), cerr_bak(NULL), clog_bak(NULL), allPhaseSpaces(),
	input_pdfsAndData(), controlled_parameters(), nuisenceModel(0), stored_pdfs(), stored_dataconfigs(), sWeighted_study(), sweight_error(), generate_n_events(0),
	ParameterSetWithFreeParameters(), ParameterSetWithFixedParameters()
{
	cout_bak = cout.rdbuf();
	cerr_bak = cerr.rdbuf();
	clog_bak = clog.rdbuf();

	//	External knowledge from the previous fits in RapidFit
	GlobalFitResult = input_GlobalResult;
	GlobalFitPhysicsParameters = GlobalFitResult->GetFitResult(0)->GetResultParameterSet()->GetDummyParameterSet();
	FitAtGridPoints = ResultsForFC;

	//	Setup the parameters we want to use throughout
	input_pdfsAndData = new_pdfsAndData;
	xmlConfig = new_xmlFile;
	numberStudies = xmlConfig->GetNumberRepeats();
	theMinimiser = newMinimiser;
	theFunction = newFunction;
	theFunction->SetIntegratorTest( false );
	allConstraints = xmlConfig->GetConstraints();
	delete_objects = false;

	nuisenceModel = inputNuisenceModel;

	sWeighted_study = newFunction->GetWeightsWereUsed();

	controlled_parameters = new_makeOutput->GetScanList();
	vector<pair<string, string> > _2d_param_names = new_makeOutput->Get2DScanList();
	for( unsigned int i=0; i< _2d_param_names.size(); ++i )
	{
		controlled_parameters.push_back( _2d_param_names[i].first );
		controlled_parameters.push_back( _2d_param_names[i].second );
	}

	this->InitalizeNewPDFWithData();
}

VectoredFeldmanCousins::~VectoredFeldmanCousins()
{
	if( GlobalFitPhysicsParameters != NULL ) delete GlobalFitPhysicsParameters;
	while( !allPhaseSpaces.empty() ) { if( allPhaseSpaces.back() != NULL ) { delete allPhaseSpaces.back(); } allPhaseSpaces.pop_back(); }
	while( !stored_pdfs.empty() ) { if( stored_pdfs.back() != NULL ) { delete stored_pdfs.back(); } stored_pdfs.pop_back(); }
	while( !stored_dataconfigs.empty() ) { if( stored_dataconfigs.back() != NULL ) { delete stored_dataconfigs.back(); } stored_dataconfigs.pop_back(); }
}

void VectoredFeldmanCousins::InitalizeNewPDFWithData()
{
	for( vector<PDFWithData*>::iterator pdfdat_i = input_pdfsAndData.begin(); pdfdat_i != input_pdfsAndData.end(); ++pdfdat_i )
	{
		IDataSet* file_input = (*pdfdat_i)->GetDataSet();
		//	Read in the real data and culculate the sWeight error and dataset size
		allPhaseSpaces.push_back( new PhaseSpaceBoundary( *(file_input->GetBoundary()) ) );

		//	sWeighted FC studies are performed to compare the sensitivity between FC and the NLL from data
		//	In order to do this we need to account for the alpha parameter
		//	The NLL function is scaled by alpha in the fit to determine the correct sensitivity to a parameter
		//	Because of this if toys are generated with a scaled yield to have a comprable sensitivity

		if( sWeighted_study )
		{
			string sWeightName = theFunction->GetWeightName();
			ObservableRef sWeightObsRef ( sWeightName );

			string alphaName;
			ObservableRef alphaNameRef;

			//	When using a custom alpha with the result stored on an event by event basis then we need to retrieve it this way

			if( theFunction->GetUseCustomAlpha() )
			{
				alphaName = theFunction->GetAlphaName();
				alphaNameRef = ObservableRef( alphaName );
			}

			//Only generate toys with an sWeight of 1
			vector<double> new_constraint(1,1.0);
			allPhaseSpaces.back()->SetConstraint( sWeightName, new_constraint, " " );
			// 	If we have a per-event alpha we need to modify the constraint to remove it from the toys
			if( theFunction->GetUseCustomAlpha() )
			{
				allPhaseSpaces.back()->SetConstraint( alphaName, new_constraint, " " );
			}

			int datasetsize = file_input->GetDataNumber();
			double data_num=0.; double data_err=0.;
			for( int i=0; i < datasetsize; ++i)
			{
				double sweight_val = file_input->GetDataPoint( i )->GetObservable( sWeightObsRef )->GetValue();
				if( theFunction->GetUseCustomAlpha() )
				{
					double this_alpha = file_input->GetDataPoint( i )->GetObservable( alphaNameRef )->GetValue();
					data_num+=sweight_val*this_alpha;
					data_err+=sweight_val*sweight_val*this_alpha*this_alpha;
				}
				else
				{
					data_num+=sweight_val;
					data_err+=sweight_val*sweight_val;
				}
			}
			data_err=sqrt(data_err);

			//	When alpha has been calculated at runTime it isn't stored in a new Observable
			//	Because of this we need to call the GetAlpha to determine what was calculated for the dataset
			//	This is not done on a per-event basis and alpha is initialized to 1. by default so this should always be safe

			if( !theFunction->GetUseCustomAlpha() )
			{
				if( theFunction->GetNormaliseWeights() || theFunction->GetSingleNormaliseWeights() )
				{
					double this_alpha = file_input->GetAlpha();
					data_err*=this_alpha;
					data_num*=this_alpha;
				}
			}

			sweight_error.push_back( data_err );
			generate_n_events.push_back( data_num );
		}
		else
		{
			generate_n_events.push_back( file_input->GetDataNumber() );
		}

		vector<string> empty_args, empty_arg_names;	//Not used for data generators

		IPDF* PDFfromFit = ClassLookUp::CopyPDF( (*pdfdat_i)->GetPDF() );
		stored_pdfs.push_back( PDFfromFit );

		vector<DataSetConfiguration*> temp_datasetconfig;
		DataSetConfiguration* configForToys = new DataSetConfiguration( "Foam", -1, "", empty_args, empty_arg_names, stored_pdfs.back() );
		stored_dataconfigs.push_back( configForToys );
		temp_datasetconfig.push_back( configForToys );

		PDFWithData* ToyPDFWithData = new PDFWithData( stored_pdfs.back(), allPhaseSpaces.back(), temp_datasetconfig );
		pdfsAndData.push_back( ToyPDFWithData );
	}
}

//	Initialize ParameterSetWithFreeParameters and ParameterSetWithFixedParameters as copies of the input dataset
void VectoredFeldmanCousins::InitializePhysicsParameters( ParameterSet* inputParameters )
{
	ParameterSet* temp_freeParam = new ParameterSet( *inputParameters );
	ParameterSet* temp_fixedParam = new ParameterSet( *inputParameters );

	vector<string>::iterator fixed_param_i = controlled_parameters.begin();
	for( ; fixed_param_i != controlled_parameters.end(); ++fixed_param_i )
	{
		temp_fixedParam->GetPhysicsParameter( *fixed_param_i )->SetType( "Fixed" );
		temp_freeParam->GetPhysicsParameter( *fixed_param_i )->SetType( "Free" );//Just to be explicit, even though this is expected to not be needed
	}

	if( ParameterSetWithFreeParameters != NULL ) delete ParameterSetWithFreeParameters;
	ParameterSetWithFreeParameters = temp_freeParam;
	if( ParameterSetWithFixedParameters != NULL ) delete ParameterSetWithFixedParameters;
	ParameterSetWithFixedParameters = temp_fixedParam;
}

vector<IDataSet*> VectoredFeldmanCousins::GetNewDataSets( ParameterSet* input_params )
{
	vector<IDataSet*> output_datasets;

	vector<double>::iterator sWeight_errors = sweight_error.begin();
	vector<double>::iterator wanted_events = generate_n_events.begin();
	vector<PDFWithData*>::iterator pdfdat_i = pdfsAndData.begin();
	vector<PhaseSpaceBoundary*>::iterator phaseSpace_i = allPhaseSpaces.begin();

	for( ; pdfdat_i != pdfsAndData.end(); ++pdfdat_i, ++wanted_events, ++phaseSpace_i )
	{
		//vector<ParameterSet*> temp_vec( 1, input_params );
		(*pdfdat_i)->SetPhysicsParameters( input_params );

		TRandom3* new_rand = (*pdfdat_i)->GetPDF()->GetRandomFunction();

		double num_wanted_events=-1;
		if( sWeighted_study )
		{
			num_wanted_events = new_rand->Gaus( (*wanted_events), (*sWeight_errors) );
			++sWeight_errors;
		}
		else
		{
			num_wanted_events = (*wanted_events);
		}

		IDataSet* new_dataset=NULL;

		//if( (*pdfdat_i)->GetCacheList().empty() )
		//{
		new_dataset = (*pdfdat_i)->GetDataSetConfig()->MakeDataSet( (*phaseSpace_i), (*pdfdat_i)->GetPDF(), (int)num_wanted_events );
		//	(*pdfdat_i)->AddCachedData( new_dataset );
		//}
		//else
		//{
		//	new_dataset = (*pdfdat_i)->GetDataSet();
		//}

		output_datasets.push_back( new_dataset );
	}
	return output_datasets;
}

void VectoredFeldmanCousins::DoWholeStudy( int OutputLevel )
{
	//cout << "Doing Study: " << OutputLevel << endl;
	theMinimiser->SetOutputLevel( OutputLevel );
	cout << "\n\tStarting Fit:" <<endl;

	for( unsigned int i=0; i< pdfsAndData.size(); ++i )
	{
		pdfsAndData[i]->ClearCache();
	}

	vector<FitResultVector*> temp_complete_vec;

	for( unsigned int result_i=0; result_i < (unsigned)FitAtGridPoints->NumberResults(); ++result_i )
	{
		cout << "Running at Point: " << result_i+1 << " of: " << FitAtGridPoints->NumberResults() << endl;
		FitResult* InputResult = FitAtGridPoints->GetFitResult( (int)result_i );

		cout << endl << "Setting Up Physics Parameters" << endl << endl;
		ParameterSet* gridPointSet = new ParameterSet( *(InputResult->GetResultParameterSet()->GetDummyParameterSet()) );
		this->InitializePhysicsParameters( gridPointSet );

		for( vector<PDFWithData*>::iterator pdfdat_i=pdfsAndData.begin(); pdfdat_i != pdfsAndData.end(); ++pdfdat_i )
		{
			(*pdfdat_i)->SetPhysicsParameters( ParameterSetWithFreeParameters );
		}

		vector<FitResultVector*> grid_pointResultVector;

		grid_pointResultVector.push_back( GlobalFitResult );
		FitResultVector* temp_gridpoint = new FitResultVector( GlobalFitResult->GetAllNames() );
		temp_gridpoint->AddFitResult( InputResult, false );
		temp_gridpoint->AddCPUTime( FitAtGridPoints->GetCPUTime( (int)result_i ) );
		temp_gridpoint->AddRealTime( FitAtGridPoints->GetRealTime( (int)result_i ) );
		temp_gridpoint->AddGLTime( FitAtGridPoints->GetGLTime( (int)result_i ) );
		grid_pointResultVector.push_back( temp_gridpoint );

		unsigned int this_study = (unsigned)numberStudies;

		for( unsigned short int dataset_num=0; dataset_num < this_study; ++dataset_num )
		{
			cout << "Generating and Fitting to Toy DataSet: " << dataset_num+1 << " of: " << this_study << endl;

			cout << endl << "Generating Toy Dataset at this Coordinate" << endl;
			cout << "Original DataSet was ";
			if( !sWeighted_study ) cout << "NOT ";
			cout << "an sWeighted Study" << endl;

			ParameterSet* FittingParameterSetWithFreeParameters =  this->getParameterSet( ParameterSetWithFreeParameters, InputResult->GetResultParameterSet() );
			ParameterSet* FittingParameterSetWithFixedParameters = this->getParameterSet( ParameterSetWithFixedParameters, InputResult->GetResultParameterSet() );

			//	Generate and cache a set of data
			this->SetOutput( -1 );//OutputLevel );
			cout << "Generating Data With:" << endl;
			FittingParameterSetWithFreeParameters->Print(); cout << endl;
			vector<IDataSet*> dataset_p = this->GetNewDataSets( FittingParameterSetWithFreeParameters );
			this->ResetOutput();
			cout << "Storing The DataSet at this Coordinate" << endl;
			this->SetOutput( -1 );//OutputLevel );
			for( unsigned int i=0; i< pdfsAndData.size(); ++i )
			{
				pdfsAndData[i]->AddCachedData( dataset_p[i] );
				pdfsAndData[i]->SetUseCache( true );
			}
			this->ResetOutput();

			cout << endl << "Control Parameter(s) are:" << endl;
			for( vector<string>::iterator param_i = controlled_parameters.begin(); param_i != controlled_parameters.end(); ++param_i )
			{
				cout << *param_i << ", " ;
			}
			cout << endl;

			grid_pointResultVector.push_back( GlobalFitResult );
			FitResultVector* temp_vec = new FitResultVector( GlobalFitResult->GetAllNames() );
			temp_vec->StartStopwatch();
			cout << endl << "Fitting to the Dataset with Control Parameter(s) Fixed" << endl;
			this->SetOutput( -1 );//OutputLevel );
			cout << "Fitting With:" << endl;
			FittingParameterSetWithFixedParameters->Print(); cout << endl;
			FitResult* fit1Result = FitAssembler::DoSafeFit( theMinimiser, theFunction, FittingParameterSetWithFixedParameters, pdfsAndData, allConstraints, true, OutputLevel );
			for( vector<string>::iterator param_i = controlled_parameters.begin(); param_i != controlled_parameters.end(); ++param_i )
			{
				fit1Result->GetResultParameterSet()->GetResultParameter( *param_i )->ForceType( "Fixed" );
			}

			this->ResetOutput();
			if( fit1Result->GetFitStatus() != 3 )
			{
				cout << "Fit FAILED!!!!" << endl;
				cout << endl << "Requesting additional toy dataset and not re-fitting to this one!" << endl;
				for( unsigned int i=0; i< pdfsAndData.size(); ++i )
				{
					pdfsAndData[i]->ClearCache();
				}
				delete fit1Result; delete temp_vec;
				//Remove the global result from the vector containing all results
				grid_pointResultVector.pop_back();
				++this_study;
				if( FittingParameterSetWithFreeParameters != NULL ) delete FittingParameterSetWithFreeParameters;
				if( FittingParameterSetWithFixedParameters != NULL ) delete FittingParameterSetWithFixedParameters;
				continue;
			}
			temp_vec->AddFitResult( fit1Result );
			grid_pointResultVector.push_back( temp_vec );
			cout << "Fit Finished" << endl;

			grid_pointResultVector.push_back( GlobalFitResult );
			FitResultVector* temp_vec2 = new FitResultVector( GlobalFitResult->GetAllNames() );
			temp_vec2->StartStopwatch();
			cout << endl << "Fitting to the Dataset with Control Parameter(s) Free" << endl;
			this->SetOutput( -1 );//OutputLevel );
			FitResult* fit2Result = FitAssembler::DoSafeFit( theMinimiser, theFunction, FittingParameterSetWithFreeParameters, pdfsAndData, allConstraints, true, OutputLevel );
			for( vector<string>::iterator param_i = controlled_parameters.begin(); param_i != controlled_parameters.end(); ++param_i )
			{
				fit2Result->GetResultParameterSet()->GetResultParameter( *param_i )->ForceType( "Free" );
			}

			if( fit2Result->GetFitStatus() != 3 )
			{
				this->ResetOutput();
				cout << "Fit FAILED!!!!!" << endl;
				cout << endl << "Requesting additional toy dataset" << endl;
				for( unsigned int i=0; i< pdfsAndData.size(); ++i )
				{
					pdfsAndData[i]->ClearCache();
				}
				delete fit1Result; delete temp_vec;
				delete fit2Result; delete temp_vec2;
				//Remove the global result(s) and the fixed point result from this vector of all of the results
				grid_pointResultVector.pop_back();
				grid_pointResultVector.pop_back();
				grid_pointResultVector.pop_back();
				++this_study;
				if( FittingParameterSetWithFreeParameters != NULL ) delete FittingParameterSetWithFreeParameters;
				if( FittingParameterSetWithFixedParameters != NULL ) delete FittingParameterSetWithFixedParameters;
				continue;
			}
			temp_vec2->AddFitResult( fit2Result );
			grid_pointResultVector.push_back( temp_vec2 );
			this->ResetOutput();
			cout << "Fit Finished" << endl;

			cout << endl << "Removing cached DataSet" << endl << endl;
			for( unsigned int i=0; i< pdfsAndData.size(); ++i )
			{
				pdfsAndData[i]->ClearCache();
			}

			if( FittingParameterSetWithFreeParameters != NULL ) delete FittingParameterSetWithFreeParameters;
			if( FittingParameterSetWithFixedParameters != NULL ) delete FittingParameterSetWithFixedParameters;
		}

		cout << "Storing the Result for All toys at this grid point" << endl;
		FitResultVector* allGridPointResult = new FitResultVector( grid_pointResultVector );
		temp_complete_vec.push_back( allGridPointResult );
	}
	cout << endl << "Finalizing all FC Results" << endl;
	allResults = new FitResultVector( temp_complete_vec );
}

void VectoredFeldmanCousins::SetOutput( int OutputLevel )
{
	//	If the user wanted silence we point the Std Output Streams to /dev/null
	if( OutputLevel <= 0 )
	{
		cout.rdbuf(NULL);
		cerr.rdbuf(NULL);
		clog.rdbuf(NULL);
	}
}

void VectoredFeldmanCousins::ResetOutput()
{
	cout.rdbuf( cout_bak );
	cerr.rdbuf( cerr_bak );
	clog.rdbuf( clog_bak );
}

FitResultVector* VectoredFeldmanCousins::GetStudyResult()
{
	return allResults;
}

void VectoredFeldmanCousins::SetNumRepeats( int input )
{
	numberStudies = input;
}

void VectoredFeldmanCousins::SetCommandLineParams( vector<string> input )
{
	(void) input;
}

ParameterSet* VectoredFeldmanCousins::getParameterSet( ParameterSet* inputSet, ResultParameterSet* inputResult )
{
	// Model 1 nuisence parameters are not changed
	if( nuisenceModel == 0 )
	{
		return new ParameterSet( *inputSet );
	}
	// Model 2 nuisence parameters are varied within 1 sigma of their true value
	else if( nuisenceModel == 1 )
	{
		ParameterSet* tempSet = new ParameterSet( *inputSet );
		vector<string> all_params = tempSet->GetAllNames();
		for( vector<string>::iterator param_i = all_params.begin(); param_i != all_params.end(); ++param_i )
		{
			PhysicsParameter* thisParameter = tempSet->GetPhysicsParameter( *param_i );
			ResultParameter* thisResult = inputResult->GetResultParameter( *param_i );
			if( thisParameter->GetType() != "Fixed" && !thisResult->GetScanStatus() )
			{
				TRandom3* rand_gen = input_pdfsAndData[0]->GetPDF()->GetRandomFunction();
				double new_value = rand_gen->Gaus( thisResult->GetValue(), thisResult->GetError() );
				thisParameter->SetBlindedValue( new_value );
				thisParameter->SetStepSize( thisResult->GetError() );
			}
		}
		return tempSet;
	}
	else
	{
		cout << endl << "\t\tUNKNOWN NUICENCE MODEL, NOT DOING ANTYTHING!!!!" << endl << endl;
		nuisenceModel = 0;
		return this->getParameterSet( inputSet, inputResult );
	}
}

