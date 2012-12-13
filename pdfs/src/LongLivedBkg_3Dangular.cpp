// $Id: LongLivedBkg_3Dangular.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class LongLivedBkg_3Dangular LongLivedBkg_3Dangular.cpp
 *
 *  PDF for  long lived background with 3D histogram angular description. A 3D root Histogram required provided by user and input in the xml file. The transversity angles in the order (x,y,z) --> (cosPsi, cosTheta, phi)
 *
 *  @author Ailsa Sparkes
 *  @date 2011-05-30
 */

#include "LongLivedBkg_3Dangular.h"
#include "Mathematics.h"
#include <iostream>
#include <fstream>
#include "math.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH3D.h"
#include "TAxis.h"
#include "TH1.h"

#include "StringProcessing.h"

using namespace::std;

PDF_CREATOR( LongLivedBkg_3Dangular );

//.....................................................
//Constructor
LongLivedBkg_3Dangular::LongLivedBkg_3Dangular( PDFConfigurator* config ) :

	// Physics parameters
	f_LL1Name		( config->getName("f_LL1")  )
	, tauLL1Name		( config->getName("tau_LL1") )
	, tauLL2Name		( config->getName("tau_LL2") )
	//Detector parameters
	, timeResLL1FracName( config->getName("timeResLL1Frac") )
	, sigmaLL1Name		( config->getName("sigma_LL1") )
	, sigmaLL2Name		( config->getName("sigma_LL2") )
	// Observables
	, timeName			( config->getName("time") )
	, cosThetaName		( config->getName("cosTheta") )
	, phiName			( config->getName("phi") )
	, cosPsiName		( config->getName("cosPsi") )
	, cthetakName 		( config->getName("helcosthetaK") )
	, cthetalName		( config->getName("helcosthetaL") )
	, phihName			( config->getName("helphi") )
	, timeconstName		( config->getName("time") )
	, resScaleName                  ( config->getName("timeResolutionScale") )
	, eventResolutionName	( config->getName("eventResolution") )
	, timeOffsetName                ( config->getName("timeOffset") )
	, tauLL1(), tauLL2(), f_LL1(), sigmaLL(), sigmaLL1(), sigmaLL2(), timeResLL1Frac(), tlow(), thigh(), time(), cos1(),
	cos2(), phi(), histo(), xaxis(), yaxis(), zaxis(), nxbins(), nybins(), nzbins(), xmin(), xmax(), ymin(),
	ymax(), zmin(), zmax(), deltax(), deltay(), deltaz(), total_num_entries(), useFlatAngularDistribution(true),
	_useHelicityBasis(false), _usePerEvent(false), _usePunziSigmat(false)
	, _useTimeAcceptance(false), timeAcc(NULL)
{

	cout << "Constructing PDF: LongLivedBkg_3Dangular  " << endl;

	//Configure
	_useHelicityBasis = config->isTrue( "UseHelicityBasis" );
	_usePerEvent = config->isTrue( "UseEventResolution" );
	_usePunziSigmat = config->isTrue( "UsePunziSigmat" );
	_useTimeAcceptance = config->isTrue( "UseTimeAcceptance" );

	//Make prototypes
	this->MakePrototypes();

	//...........................................
	// Configure to use time acceptance machinery
	if( _useTimeAcceptance )
	{
		if( config->hasConfigurationValue( "TimeAcceptanceType", "Upper" ) )
		{
			timeAcc = new SlicedAcceptance( 0., 14.0, /*0.0157*/ 0.0112);
			cout << "LongLivedBkg_3Dangular:: Constructing timeAcc: Upper time acceptance beta=0.0112 [0 < t < 14] " << endl;
		}
		else if( config->getConfigurationValue( "TimeAcceptanceFile" ) != "" )
		{
			timeAcc = new SlicedAcceptance( "File" , config->getConfigurationValue( "TimeAcceptanceFile" ) );
			cout << "LongLivedBkg_3Dangular:: Constructing timeAcc: using file: " << config->getConfigurationValue( "TimeAcceptanceFile" ) << endl;
		}
	}
	else
	{
		timeAcc = new SlicedAcceptance( 0., 14. );
		cout << "LongLivedBkg_3Dangular:: Constructing timeAcc: DEFAULT FLAT [0 < t < 14]  " << endl;
	}

	this->InitializeCaching();	//	Depends on Acceptance File

	//Find name of histogram needed to define 3-D angular distribution
	string fileName = config->getConfigurationValue( "AngularDistributionHistogram" );

	//Initialise depending upon whether configuration parameter was found
	if( fileName == "" )
	{
		cout << " No AngularDistributionHistogram found: using flat background " << endl;
		useFlatAngularDistribution = true ;
	}
	else
	{
		cout << " AngularDistributionHistogram requested: " << fileName << endl ;
		useFlatAngularDistribution = false ;

		//File location
		string fullFileName = StringProcessing::FindFileName( fileName );

		//Read in histo
		TFile* f =  TFile::Open(fullFileName.c_str());
		if( _useHelicityBasis )
		{
			histo = new TH3D( *(TH3D*) f->Get("histoHel") ); //(fileName.c_str())));
		}
		else
		{
			histo = new TH3D( *(TH3D*) f->Get("histo") ); //(fileName.c_str())));
		}

		TString Name = histo->GetName();
		TString Title = histo->GetTitle();
		TString ext = "_";
		ext+=this->GetRandomFunction()->Rndm();
		Name.Append(ext);
		Title.Append(ext);

		xaxis = histo->GetXaxis();
		xmin = xaxis->GetXmin();
		xmax = xaxis->GetXmax();
		nxbins = histo->GetNbinsX();
		deltax = (xmax-xmin)/nxbins;
		cout << " X axis Name: " << xaxis->GetName() << "\tTitle: " << xaxis->GetTitle() << "\t\t" << "X axis Min: " << xmin << "\tMax: " << xmax << "\tBins: " << nxbins << endl;

		yaxis = histo->GetYaxis();
		ymin = yaxis->GetXmin();
		ymax = yaxis->GetXmax();
		nybins = histo->GetNbinsY();
		deltay = (ymax-ymin)/nybins;
		cout << " Y axis Name: " << yaxis->GetName() << "\tTitle: " << yaxis->GetTitle() << "\t\t" << "Y axis Min: " << ymin << "\tMax: " << ymax << "\tBins: " << nybins << endl;

		zaxis = histo->GetZaxis();
		zmin = zaxis->GetXmin();
		zmax = zaxis->GetXmax();
		nzbins = histo->GetNbinsZ();
		deltaz = (zmax-zmin)/nzbins;
		cout << " Z axis Name: " << zaxis->GetName() << "\tTitle: " << zaxis->GetTitle() << "\t\t" << "Z axis Min: " << zmin << "\tMax: " << zmax << "\tBins: " << nzbins << "\t\t\t";

		//method for Checking whether histogram is sensible

		int total_num_bins = nxbins * nybins * nzbins;
		int sum = 0;

		vector<int> zero_bins;

		//loop over each bin in histogram and print out how many zero bins there are
		for (int i=1; i < nxbins+1; ++i)
		{
			for (int j=1; j < nybins+1; ++j)
			{
				for (int k=1; k < nzbins+1; ++k)
				{
					double bin_content = histo->GetBinContent(i,j,k);
					//cout << "Bin content: " << bin_content << endl;
					if(bin_content<=0)
					{
						zero_bins.push_back(0);
					}
					//cout << " Zero bins " << zero_bins.size() << endl;}
					else if (bin_content>0)
					{
						sum += (int) bin_content;
					}
			}
		}
	}

	total_num_entries = histo->GetEntries();
	int average_bin_content = sum / total_num_bins;

	cout << "Total Bins: " << total_num_bins << "\tEmpty Bins: " << zero_bins.size() << "\tAvg Bin Content: " << average_bin_content << endl;

	//cout << "\n\t\t" << "****" << "For total number of bins " << total_num_bins << " there are " << zero_bins.size() << " bins containing zero events " << "****" << endl;
	//cout <<  "\t\t\t" << "***" << "Average number of entries of non-zero bins: " << average_bin_content << "***" << endl;
	//cout << endl;

	// Check.  This order works for both bases since phi is always the third one.
	if ((xmax-xmin) < 2. || (ymax-ymin) < 2. || (zmax-zmin) < 2.*Mathematics::Pi() )
	{
		cout << endl << "In LongLivedBkg_3Dangular::LongLivedBkg_3Dangular: The full angular range is not used in this histogram - the PDF does not support this case" << endl;
		exit(1);
	}

	//cout << "Finishing processing histo" << endl;
}

}

//................................................................
//Destructor
LongLivedBkg_3Dangular::~LongLivedBkg_3Dangular()
{
}

void LongLivedBkg_3Dangular::InitializeCaching()
{
	if( _useTimeAcceptance )
	{
		for( unsigned int i=0; i< timeAcc->numberOfSlices(); ++i )
		{
			TString LL1a_int_vec_name="LLBKG_LL1a_time_"; LL1a_int_vec_name+=i;
			LL1a_int_vec.push_back( PseudoObservable( LL1a_int_vec_name.Data() ) );
			LL1a_int_vec.back().AddFunction( Mathematics::ExpInt_Wrapper );

			TString LL2a_int_vec_name="LLBKG_LL2a_time_"; LL2a_int_vec_name+=i;
			LL2a_int_vec.push_back( PseudoObservable( LL2a_int_vec_name.Data() ) );
			LL2a_int_vec.back().AddFunction( Mathematics::ExpInt_Wrapper );

			TString LL1b_int_vec_name="LLBKG_LL1b_time_"; LL1b_int_vec_name+=i;
			LL1b_int_vec.push_back( PseudoObservable( LL1b_int_vec_name.Data() ) );
			LL1b_int_vec.back().AddFunction( Mathematics::ExpInt_Wrapper );

			TString LL2b_int_vec_name="LLBKG_LL2b_time_"; LL2b_int_vec_name+=i;
			LL2b_int_vec.push_back( PseudoObservable( LL2b_int_vec_name.Data() ) );
			LL2b_int_vec.back().AddFunction( Mathematics::ExpInt_Wrapper );
		}
	}
	else
	{
		TString LL1a_Cache_name="LLBKG_LL1a_time";
		LL1a_int = PseudoObservable( LL1a_Cache_name.Data() );
		LL1a_int.AddFunction( Mathematics::ExpInt_Wrapper );

		TString LL2a_Cache_name="LLBKG_LL2a_time";
		LL2a_int = PseudoObservable( LL2a_Cache_name.Data() );
		LL2a_int.AddFunction( Mathematics::ExpInt_Wrapper );		

		TString LL1b_Cache_name="LLBKG_LL1b_time";
		LL1b_int = PseudoObservable( LL1b_Cache_name.Data() );
		LL1b_int.AddFunction( Mathematics::ExpInt_Wrapper );

		TString LL2b_Cache_name="LLBKG_LL2b_time";
		LL2b_int = PseudoObservable( LL2b_Cache_name.Data() );
		LL2b_int.AddFunction( Mathematics::ExpInt_Wrapper );
	}
}

vector<string> LongLivedBkg_3Dangular::GetDoNotIntegrateList()
{
	vector<string> doNotList;
	if( _usePerEvent && !_usePunziSigmat ) doNotList.push_back( eventResolutionName );
	return doNotList;
}

//..................................................................
//Make the data point and parameter set
void LongLivedBkg_3Dangular::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	if( _useHelicityBasis ) {
		allObservables.push_back( cthetakName );
		allObservables.push_back( cthetalName );
		allObservables.push_back( phihName );
	}
	else {
		allObservables.push_back( cosThetaName );
		allObservables.push_back( phiName );
		allObservables.push_back( cosPsiName );
	}

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( f_LL1Name );
	parameterNames.push_back( tauLL1Name );
	parameterNames.push_back( tauLL2Name );

	if( !_usePerEvent )
	{
		parameterNames.push_back( timeResLL1FracName );
		parameterNames.push_back( sigmaLL1Name );
		parameterNames.push_back( sigmaLL2Name );
	}
	else
	{
		parameterNames.push_back( resScaleName );
		allObservables.push_back( eventResolutionName );
	}

	if( _useTimeAcceptance )
	{
		parameterNames.push_back( timeOffsetName );
	}

	allParameters = ParameterSet(parameterNames);
}

//.................................................................
bool LongLivedBkg_3Dangular::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	f_LL1       = allParameters.GetPhysicsParameter( f_LL1Name )->GetValue();
	tauLL1      = allParameters.GetPhysicsParameter( tauLL1Name )->GetValue();
	tauLL2      = allParameters.GetPhysicsParameter( tauLL2Name )->GetValue();

	if( !_usePerEvent )
	{
		timeResLL1Frac = allParameters.GetPhysicsParameter( timeResLL1FracName )->GetValue();
		sigmaLL1    = allParameters.GetPhysicsParameter( sigmaLL1Name )->GetValue();
		sigmaLL2    = allParameters.GetPhysicsParameter( sigmaLL2Name )->GetValue();
	}
	else
	{
		resolutionScale = allParameters.GetPhysicsParameter( resScaleName )->GetValue();
	}

	if( _useTimeAcceptance )
	{
		timeOffset = allParameters.GetPhysicsParameter( timeOffsetName )->GetValue();
	}

	return isOK;
}

//..............................................................
//Main method to build the PDF return value
double LongLivedBkg_3Dangular::Evaluate(DataPoint * measurement)
{
	_datapoint=measurement;		//	Needed for AngAcc caching

	// Observable
	Observable* timeObs = measurement->GetObservable( timeName );
	time = timeObs->GetValue();

	IConstraint* timeConst = measurement->GetPhaseSpaceBoundary()->GetConstraint( timeName );

	tlow = timeConst->GetMinimum();
	thigh = timeConst->GetMaximum();

	if( _useHelicityBasis )
	{
		cos1   = measurement->GetObservable( cthetakName )->GetValue();
		cos2   = measurement->GetObservable( cthetalName )->GetValue();
		phi    = measurement->GetObservable( phihName )->GetValue();  // Pi offset is difference between angle calculator and "Our Paper"
	}
	else
	{
		cos1   = measurement->GetObservable( cosPsiName )->GetValue();
		cos2   = measurement->GetObservable( cosThetaName )->GetValue();
		phi    = measurement->GetObservable( phiName )->GetValue();
	}

	double returnValue = 0;
	double val1=-1., val2=-1.;

	//Deal with propertime resolution
	if( _usePerEvent )
	{
		sigmaLL = measurement->GetObservable( eventResolutionName )->GetValue() * resolutionScale;

		returnValue =  this->buildPDFnumerator( 1 );
	}
	else
	{
		if( timeResLL1Frac >= 0.9999 )
		{
			// Set the member variable for time resolution to the first value and calculate
			sigmaLL = sigmaLL1;
			returnValue =  this->buildPDFnumerator( 1 );
		}
		else
		{
			// Set the member variable for time resolution to the first value and calculate
			sigmaLL = sigmaLL1;
			val1 = this->buildPDFnumerator( 1 );

			sigmaLL = sigmaLL2;
			val2 = this->buildPDFnumerator( 2 );

			returnValue = (timeResLL1Frac*val1 + (1. - timeResLL1Frac)*val2);
		}
	}


	if (returnValue <= 0)
	{
		PDF_THREAD_LOCK
		//cout << "3D Background PDF returns zero!" << endl;
		returnValue=1E-99;
		//measurement->Print();
		PDF_THREAD_UNLOCK
	}

	if( _useTimeAcceptance )
	{
		returnValue *= timeAcc->getValue( timeObs, timeOffset );
	}

	return returnValue;
}


//.............................................................
// Core calculation of PDF value
double LongLivedBkg_3Dangular::buildPDFnumerator( unsigned int sigmaNum )
{
	// Sum of two exponentials, using the time resolution functions

	double returnValue = 0.;

	double val1=-1., val2=-1., norm1=-1., norm2=-1.;
	if( _usePerEvent )
	{
		double inv_tau=1./tauLL1;

		val1 = Mathematics::Exp(time, inv_tau, sigmaLL);

		norm1= this->buildPDFdenominator( tauLL1, sigmaLL, tlow, thigh, 1, 1 );

		returnValue = val1 / norm1;
	}
	else
	{
		if( f_LL1 >= 0.9999 )
		{
			if( tauLL1 <= 0 )
			{
				PDF_THREAD_LOCK
				cout << " In LongLivedBkg_3Dangular() you gave a negative or zero lifetime for tauLL1 " << endl ;
				PDF_THREAD_UNLOCK
				throw(-7632);
			}

			double inv_tau=1./tauLL1;

			val1 = Mathematics::Exp(time, inv_tau, sigmaLL);

			norm1= this->buildPDFdenominator( tauLL1, sigmaLL, tlow, thigh, 1, sigmaNum );

			returnValue = val1 /norm1;
		}
		else
		{
			if( (tauLL1 <= 0) ||  (tauLL2 <= 0) )
			{
				PDF_THREAD_LOCK
				cout << " In LongLivedBkg_3Dangular() you gave a negative or zero lifetime for tauLL1/2 " << endl ;
				PDF_THREAD_UNLOCK
				throw(7632);
			}

			double inv_tau_LL1 = 1./ tauLL1;
			double inv_tau_LL2 = 1./ tauLL2;


			val1 = Mathematics::Exp(time, inv_tau_LL1, sigmaLL );
			val2 = Mathematics::Exp(time, inv_tau_LL2, sigmaLL );


			norm1= this->buildPDFdenominator( tauLL1, sigmaLL, tlow, thigh, 1, sigmaNum );
			norm2= this->buildPDFdenominator( tauLL2, sigmaLL, tlow, thigh, 2, sigmaNum );

			returnValue = f_LL1 * val1/norm1 + (1. - f_LL1) * val2/norm2;
		}
	}

	returnValue *= angularFactor();

	return returnValue;
}


//..............................................................
// Normlisation
double LongLivedBkg_3Dangular::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void) boundary;
	return 1.;
}

//.............................................................
//
double LongLivedBkg_3Dangular::buildPDFdenominator( double input_tauLL, double input_sigmaLL, double t_min, double t_max, unsigned int LLnum, unsigned int sigmaNum )
{
	double inv_tau=1./input_tauLL;

	double tlo_boundary = t_min;
	double thi_boundary = t_max;

	double returnValue=0.;
	double tlow_temp=0., thigh_temp=0.;

	vector<double> ExpInt_Input( 4, 0. );

	if( _useTimeAcceptance )
	{
		ExpInt_Input[0] = tlow_temp; ExpInt_Input[1] = thigh_temp;  ExpInt_Input[2] = inv_tau; ExpInt_Input[3] = input_sigmaLL;
		//This loops over each time slice, does the normalisation between the limits, and accumulates
		for( unsigned int islice = 0; islice < timeAcc->numberOfSlices(); ++islice )
		{
			tlow_temp = tlo_boundary > timeAcc->getSlice(islice)->tlow() ? tlo_boundary : timeAcc->getSlice(islice)->tlow();
			thigh_temp = thi_boundary < timeAcc->getSlice(islice)->thigh() ? thi_boundary : timeAcc->getSlice(islice)->thigh();

			if( thigh_temp > tlow_temp )
			{
				ExpInt_Input[0] = tlow_temp; ExpInt_Input[1] = thigh_temp;

				if( LLnum == 1 )
				{
					if( sigmaNum == 1)
					{
						returnValue += _datapoint->GetPseudoObservable( LL1a_int_vec[islice], ExpInt_Input ) * timeAcc->getSlice(islice)->height();
					}
					else
					{
						returnValue += _datapoint->GetPseudoObservable( LL1b_int_vec[islice], ExpInt_Input ) * timeAcc->getSlice(islice)->height();
					}
				}
				else if ( LLnum == 2 )
				{
					if( sigmaNum == 1 )
					{
						returnValue += _datapoint->GetPseudoObservable( LL2a_int_vec[islice], ExpInt_Input ) * timeAcc->getSlice(islice)->height();
					}
					else
					{
						returnValue += _datapoint->GetPseudoObservable( LL2b_int_vec[islice], ExpInt_Input ) * timeAcc->getSlice(islice)->height();
					}
				}
				else
				{
					returnValue += Mathematics::ExpInt( tlow_temp, thigh_temp, inv_tau, input_sigmaLL ) * timeAcc->getSlice(islice)->height();
				}
			}

		}
	}
	else
	{
		ExpInt_Input[0] = t_min; ExpInt_Input[1] = t_max;  ExpInt_Input[2] = inv_tau; ExpInt_Input[3] = input_sigmaLL;
		if( LLnum == 1 )
		{
			if( sigmaNum == 1 )
			{
				returnValue = _datapoint->GetPseudoObservable( LL1a_int, ExpInt_Input );
			}
			else
			{
				returnValue = _datapoint->GetPseudoObservable( LL1b_int, ExpInt_Input );
			}
		}
		else if( LLnum == 2 )
		{
			if( sigmaNum == 1 )
			{
				returnValue = _datapoint->GetPseudoObservable( LL2a_int, ExpInt_Input );
			}
			else
			{
				returnValue = _datapoint->GetPseudoObservable( LL2b_int, ExpInt_Input );
			}
		}
		else
		{
			returnValue = Mathematics::ExpInt( t_min, t_min, inv_tau, input_sigmaLL );
		}
	}

	return returnValue;
}


//................................................................
//Angular distribution function
double LongLivedBkg_3Dangular::angularFactor( )
{
	if( useFlatAngularDistribution ) return 0.125 * Mathematics::_Over_PI();

	double returnValue=0.;

	int globalbin=-1;
	int xbin=-1, ybin=-1, zbin=-1;
	double num_entries_bin=-1.;

	Observable* cos1Obs = NULL;

	if( !_useHelicityBasis ) cos1Obs = _datapoint->GetObservable( cosPsiName );
	else cos1Obs = _datapoint->GetObservable( cthetakName );

	if( cos1Obs->GetBkgBinNumber() != -1 )
	{
		return cos1Obs->GetBkgAcceptance();
	}
	else
	{
		Observable* cos2Obs = NULL;
		Observable* phiObs = NULL;
		if( _useHelicityBasis ) {
			cos2Obs   = _datapoint->GetObservable( cthetalName );
			phiObs    = _datapoint->GetObservable( phihName );  // Pi offset is difference between angle calculator and "Our Paper"
		}
		else
		{
			cos2Obs   = _datapoint->GetObservable( cosThetaName );
			phiObs    = _datapoint->GetObservable( phiName );
		}

		//Find global bin number for values of angles, find number of entries per bin, divide by volume per bin and normalise with total number of entries in the histogram
		xbin = xaxis->FindFixBin( cos1 ); if( xbin > nxbins ) xbin = nxbins; if( xbin == 0 ) xbin = 1;
		ybin = yaxis->FindFixBin( cos2 ); if( ybin > nybins ) ybin = nybins; if( ybin == 0 ) ybin = 1;
		zbin = zaxis->FindFixBin( phi ); if( zbin > nzbins ) zbin = nzbins; if( zbin == 0 ) zbin = 1;

		globalbin = histo->GetBin( xbin, ybin, zbin );
		num_entries_bin = (double)histo->GetBinContent(globalbin);

		//Angular factor normalized with phase space of histogram and total number of entries in the histogram
		returnValue = (double)num_entries_bin / (deltax * deltay * deltaz) / (double)total_num_entries;
		//returnValue = (double)num_entries_bin / histo->Integral();

		cos1Obs->SetBkgBinNumber( xbin ); cos1Obs->SetBkgAcceptance( returnValue );
		cos2Obs->SetBkgBinNumber( ybin ); cos2Obs->SetBkgAcceptance( returnValue );
		phiObs->SetBkgBinNumber( zbin ); phiObs->SetBkgAcceptance( returnValue );
	}

	return returnValue;
}

