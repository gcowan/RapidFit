/** @class Bd2JpsiKstar_sWave_Fs Bd2JpsiKstar_sWave_Fs.cpp
 *
 *  RapidFit PDF for Bd2JpsiPhi with average angular acceptance input as a paramter
 *
 *  @author Ailsa Sparkes asparkes@cern.ch
 *  @date 2011-01-26
 */

#include "TMath.h"
#include <cmath>

#include "Bd2JpsiKstar_sWave_Fs.h"
#include "Mathematics.h"
#include "SlicedAcceptance.h"
#include <iostream>
#include <fstream>
//#include "math.h"
#include "TFile.h"
#include "TH3D.h"
#include "TROOT.h"

PDF_CREATOR( Bd2JpsiKstar_sWave_Fs );

Bd2JpsiKstar_sWave_Fs::Bd2JpsiKstar_sWave_Fs( const Bd2JpsiKstar_sWave_Fs& input ) : BasePDF( (BasePDF&) input ),
	cachedAzeroAzeroIntB(input.cachedAzeroAzeroIntB), cachedAparaAparaIntB(input.cachedAparaAparaIntB), cachedAperpAperpIntB(input.cachedAperpAperpIntB),
	cachedAparaAperpIntB(input.cachedAparaAperpIntB), cachedAzeroAparaIntB(input.cachedAzeroAparaIntB), cachedAzeroAperpIntB(input.cachedAzeroAperpIntB),
	cachedAsAsIntB(input.cachedAsAsIntB), cachedAparaAsIntB(input.cachedAparaAsIntB), cachedAperpAsIntB(input.cachedAperpAsIntB), cachedAzeroAsIntB(input.cachedAzeroAsIntB),
	AzeroAzeroB(input.AzeroAzeroB), AparaAparaB(input.AparaAparaB), AperpAperpB(input.AperpAperpB), AsAsB(input.AsAsB), ImAparaAperpB(input.ImAparaAperpB),
	ReAzeroAparaB(input.ReAzeroAparaB), ImAzeroAperpB(input.ImAzeroAperpB), ReAparaAsB(input.ReAparaAsB), ImAperpAsB(input.ImAperpAsB), ReAzeroAsB(input.ReAzeroAsB),
	cachedSinDeltaPerpPara(input.cachedSinDeltaPerpPara), cachedCosDeltaPara(input.cachedCosDeltaPara), cachedSinDeltaPerp(input.cachedSinDeltaPerp),
	cachedCosDeltaParaS(input.cachedCosDeltaParaS), cachedSinDeltaPerpS(input.cachedSinDeltaPerpS), cachedCosDeltaS(input.cachedCosDeltaS), cachedAzero(input.cachedAzero),
	cachedApara(input.cachedApara), cachedAperp(input.cachedAperp), cachedAs(input.cachedAs), timeAcc(NULL), gammaName(input.gammaName), Rzero_sqName(input.Rzero_sqName),
	Rpara_sqName(input.Rpara_sqName), Rperp_sqName(input.Rperp_sqName), As_sqName(input.As_sqName), delta_zeroName(input.delta_zeroName), delta_paraName(input.delta_paraName),
	delta_perpName(input.delta_perpName), CspName(input.CspName), angAccI1Name(input.angAccI1Name), angAccI2Name(input.angAccI2Name), angAccI3Name(input.angAccI3Name),
	angAccI4Name(input.angAccI4Name), angAccI5Name(input.angAccI5Name), angAccI6Name(input.angAccI6Name), angAccI7Name(input.angAccI7Name), angAccI8Name(input.angAccI8Name),
	angAccI9Name(input.angAccI9Name), angAccI10Name(input.angAccI10Name), timeRes1Name(input.timeRes1Name), timeRes2Name(input.timeRes2Name), timeRes1FractionName(input.timeRes1FractionName),
	_useTimeAcceptance(input._useTimeAcceptance), normalisationCacheValid(input.normalisationCacheValid), evaluationCacheValid(input.evaluationCacheValid), timeName(input.timeName),
	cosThetaName(input.cosThetaName), phiName(input.phiName), cosPsiName(input.cosPsiName), KstarFlavourName(input.KstarFlavourName), timeconstraintName(input.timeconstraintName),
	gamma(input.gamma), Rzero_sq(input.Rzero_sq), Rpara_sq(input.Rpara_sq), Rperp_sq(input.Rperp_sq), As_sq(input.As_sq), AzeroApara(input.AzeroApara), AzeroAperp(input.AzeroAperp),
	AparaAperp(input.AparaAperp), AparaAs(input.AparaAs), AperpAs(input.AperpAs), AzeroAs(input.AzeroAs), delta_zero(input.delta_zero), delta_para(input.delta_para), delta_perp(input.delta_perp),
	delta_s(input.delta_s), omega(input.omega), timeRes(input.timeRes), timeRes1(input.timeRes1), timeRes2(input.timeRes2), timeRes1Frac(input.timeRes1Frac), angAccI1(input.angAccI1), angAccI2(input.angAccI2),
	angAccI3(input.angAccI3), angAccI4(input.angAccI4), angAccI5(input.angAccI5), angAccI6(input.angAccI6), angAccI7(input.angAccI7), angAccI8(input.angAccI8), angAccI9(input.angAccI9),
	angAccI10(input.angAccI10), Ap_sq(input.Ap_sq), Ap(input.Ap), time(input.time), cosTheta(input.cosTheta), phi(input.phi), cosPsi(input.cosPsi), KstarFlavour(input.KstarFlavour), tlo(input.tlo),
	thi(input.thi), useFlatAngularDistribution(input.useFlatAngularDistribution), componentIndex(input.componentIndex), delta_sName(input.delta_sName), Azero_sq(input.Azero_sq),
	Apara_sq(input.Apara_sq), Aperp_sq(input.Aperp_sq), Csp(input.Csp), CspAs(input.CspAs), _plotAllComponents(input._plotAllComponents),
    histo(input.histo), xaxis(NULL), yaxis(NULL), zaxis(NULL),
    nxbins(input.nxbins), nybins(input.nybins), nzbins(input.nzbins), xmin(input.xmin), ymin(input.ymin), zmin(input.zmin),
	xmax(input.xmax), ymax(input.ymax), zmax(input.zmax), deltax(input.deltax), deltay(input.deltay), deltaz(input.deltaz),
    total_num_entries(input.total_num_entries), _datapoint(NULL),
	f(NULL), fileName(input.fileName)
{
	if( input.timeAcc != NULL ) timeAcc = new SlicedAcceptance( *(input.timeAcc) );

	if( histo != NULL )
	{

		f =  TFile::Open( fileName.c_str() );
		histo = (TH3D*) f->Get("histo");

		xaxis = histo->GetXaxis();
		xmin = xaxis->GetXmin();
		xmax = xaxis->GetXmax();
		nxbins = histo->GetNbinsX();
		deltax = (xmax-xmin)/nxbins;

		yaxis = histo->GetYaxis();
		ymin = yaxis->GetXmin();
		ymax = yaxis->GetXmax();
		nybins = histo->GetNbinsY();
		deltay = (ymax-ymin)/nybins;

		zaxis = histo->GetZaxis();
		zmin = zaxis->GetXmin();
		zmax = zaxis->GetXmax();
		nzbins = histo->GetNbinsZ();
		deltaz = (zmax-zmin)/nzbins;

		total_num_entries = histo->GetEntries();
	}
}


//Constructor
Bd2JpsiKstar_sWave_Fs::Bd2JpsiKstar_sWave_Fs(PDFConfigurator* configurator ) :
	cachedAzeroAzeroIntB(), cachedAparaAparaIntB(), cachedAperpAperpIntB(), cachedAparaAperpIntB(), cachedAzeroAparaIntB(), cachedAzeroAperpIntB(),
	cachedAsAsIntB(), cachedAparaAsIntB(), cachedAperpAsIntB(), cachedAzeroAsIntB(), AzeroAzeroB(), AparaAparaB(), AperpAperpB(), AsAsB(), ImAparaAperpB(),
	ReAzeroAparaB(), ImAzeroAperpB(), ReAparaAsB(), ImAperpAsB(), ReAzeroAsB(), cachedSinDeltaPerpPara(), cachedCosDeltaPara(), cachedSinDeltaPerp(),
	cachedCosDeltaParaS(), cachedSinDeltaPerpS(), cachedCosDeltaS(), cachedAzero(), cachedApara(), cachedAperp(), cachedAs(), timeAcc(NULL),
	// Physics parameters
	gammaName             ( configurator->getName("gamma" ))
	, Rzero_sqName ( configurator->getName("Rzero_sq"))
	, Rpara_sqName  ( configurator->getName("Rpara_sq" ))
	, Rperp_sqName  ( configurator->getName("Rperp_sq" ))
	, As_sqName     ( configurator->getName("As_sq" ))
	, delta_zeroName( configurator->getName("delta_zero" ))
	, delta_paraName( configurator->getName("delta_para" ))
	, delta_perpName( configurator->getName("delta_perp" ))
	, delta_sName   ( configurator->getName("delta_s" ))
	, CspName   ( configurator->getName("Csp" ))
	, angAccI1Name  ( configurator->getName("angAccI1" ))
	, angAccI2Name  ( configurator->getName("angAccI2" ))
	, angAccI3Name  ( configurator->getName("angAccI3" ))
	, angAccI4Name  ( configurator->getName("angAccI4" ))
	, angAccI5Name  ( configurator->getName("angAccI5" ))
	, angAccI6Name  ( configurator->getName("angAccI6" ))
	, angAccI7Name  ( configurator->getName("angAccI7" ))
	, angAccI8Name  ( configurator->getName("angAccI8" ))
	, angAccI9Name  ( configurator->getName("angAccI9" ))
	, angAccI10Name ( configurator->getName("angAccI10" ))
	, timeRes1Name  ( configurator->getName("timeResolution1") )
	, timeRes2Name  ( configurator->getName("timeResolution2" ))
	, timeRes1FractionName  ( configurator->getName("timeResolution1Fraction" ))
	, _useTimeAcceptance(false)
	, _plotAllComponents(false)

	// Observables (What we want to gain from the pdf after inserting physics parameter values)
	, normalisationCacheValid(false)
, evaluationCacheValid(false)
	, timeName      ( configurator->getName("time" ))
	, cosThetaName  ( configurator->getName("cosTheta" ))
	, phiName       ( configurator->getName("phi" ))
	, cosPsiName    ( configurator->getName("cosPsi" ))
	, KstarFlavourName  ( configurator->getName("KstarFlavour" ))

	, timeconstraintName( "time" )
, gamma(), Rzero_sq(), Rpara_sq(), Rperp_sq(), As_sq(), AzeroApara(), AzeroAperp(), AparaAperp(), AparaAs(), AperpAs(), AzeroAs(),
	delta_zero(), delta_para(), delta_perp(), delta_s(), omega(), timeRes(), timeRes1(), timeRes2(), timeRes1Frac(), angAccI1(), angAccI2(),
	angAccI3(), angAccI4(), angAccI5(), angAccI6(), angAccI7(), angAccI8(), angAccI9(), angAccI10(), Ap_sq(), Ap(), time(), cosTheta(), phi(),
	cosPsi(), KstarFlavour(), tlo(), thi(), useFlatAngularDistribution(true), _datapoint(NULL), f(NULL)
{
	componentIndex = 0;
	MakePrototypes();
	_useTimeAcceptance = configurator->isTrue( "UseTimeAcceptance" ) ;
	_plotAllComponents = configurator->isTrue( "PlotAllComponents" ) ;

	if( useTimeAcceptance() ) {
		timeAcc = new SlicedAcceptance( 0., 14.0, 0.0171 ) ;
		cout << "Bd2JpsiKstar_sWave_Fs:: Constructing timeAcc: Upper time acceptance beta=0.0171 [0 < t < 14] " << endl ;
	}
	else {
		timeAcc = new SlicedAcceptance( 0., 14. ) ;
		cout << "Bd2JpsiKstar_sWave_Fs:: Constructing timeAcc: DEFAULT FLAT [0 < t < 14]  " << endl ;
	}


	//AILSA
	//Find name of histogram needed to define 3-D angular distribution
	fileName = configurator->getConfigurationValue( "AngularDistributionHistogram" ) ;

	//Initialise depending upon whether configuration parameter was found
	if( fileName == "" )
	{
		cout << "   No AngularAcceptanceHisto found" << endl ;
		useFlatAngularDistribution = true ;
        histo = NULL ;
	}
	else
	{
		cout << "   AngularAcceptanceHisto requested: " << fileName << endl ;
		useFlatAngularDistribution = false ;



		ifstream input_file2;
		input_file2.open( fileName.c_str(), ifstream::in );
		input_file2.close();
		bool local_fail2 = input_file2.fail();
		if( !getenv("RAPIDFITROOT") && local_fail2 )
		{
			cerr << "\n" << endl;
			//cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cerr << "$RAPIDFITROOT NOT DEFINED, PLEASE DEFINE IT SO I CAN USE ACCEPTANCE DATA" << endl;
			//cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cerr << "\n" << endl;
			exit(-987);
		}


		string fullFileName;

		//File location
		//ifstream input_file2;
		//input_file2.open( fileName.c_str(), ifstream::in );
		//input_file2.close();
		//bool local_fail2 = input_file2.fail();

		if( getenv("RAPIDFITROOT") )
		{
			string path( getenv("RAPIDFITROOT") ) ;

			cout << "RAPIDFITROOT defined as: " << path << endl;
			fullFileName = path+"/pdfs/configdata/"+fileName ;

			input_file2.open( fullFileName.c_str(), ifstream::in );
			input_file2.close();
		}
		bool elsewhere_fail = input_file2.fail();

		if( elsewhere_fail && local_fail2 )
		{
			cerr << "\n\tFileName:\t" << fullFileName << "\t NOT FOUND PLEASE CHECK YOUR RAPIDFITROOT" << endl;
			cerr << "\t\tAlternativley make sure your XML points to the correct file, or that the file is in the current working directory\n" << endl;
			exit(-89);
		}

		if( fullFileName.empty() || !local_fail2 )
		{
			fullFileName = fileName;
		}

		fileName = fullFileName;
		f =  TFile::Open(fullFileName.c_str());
		histo = (TH3D*) f->Get("histo"); //(fileName.c_str())));

		xaxis = histo->GetXaxis();
		cout << "X axis Name: " << xaxis->GetName() << "\tTitle: " << xaxis->GetTitle() << endl;
		xmin = xaxis->GetXmin();
		xmax = xaxis->GetXmax();
		nxbins = histo->GetNbinsX();
		cout << "X axis Min: " << xmin << "\tMax: " << xmax << "\tBins: " << nxbins << endl;
		deltax = (xmax-xmin)/nxbins;

		yaxis = histo->GetYaxis();
		cout << "Y axis Name: " << yaxis->GetName() << "\tTitle: " << yaxis->GetTitle() << endl;
		ymin = yaxis->GetXmin();
		ymax = yaxis->GetXmax();
		nybins = histo->GetNbinsY();
		cout << "Y axis Min: " << ymin << "\tMax: " << ymax << "\tBins: " << nybins << endl;
		deltay = (ymax-ymin)/nybins;

		zaxis = histo->GetZaxis();
		cout << "Z axis Name: " << zaxis->GetName() << "\tTitle: " << zaxis->GetTitle() << endl;
		zmin = zaxis->GetXmin();
		zmax = zaxis->GetXmax();
		nzbins = histo->GetNbinsZ();
		cout << "Z axis Min: " << zmin << "\tMax: " << zmax << "\tBins: " << nzbins << endl;
		deltaz = (zmax-zmin)/nzbins;

		//method for Checking whether histogram is sensible

		total_num_entries = histo->GetEntries();
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
						zero_bins.push_back(1);
					}
					//cout << " Zero bins " << zero_bins.size() << endl;}
					else if (bin_content>0)                                        {
						sum += (int) bin_content;
					}
                }
            }
        }

	int average_bin_content = sum / total_num_bins;

	cout << "\n\n\t\t" << "****" << "For total number of bins " << total_num_bins << " there are " << zero_bins.size() << " bins containing zero events " << "****" << endl;
	cout <<  "\t\t\t" << "***" << "Average number of entries of non-zero bins: " << average_bin_content << "***" << endl;
	cout << endl;
	cout << endl;

	// Check.  This order works for both bases since phi is always the third one.
	if ((xmax-xmin) < 2. || (ymax-ymin) < 2. || (zmax-zmin) < 2.*Mathematics::Pi() )
	{
		cout << "In LongLivedBkg_3Dangular::LongLivedBkg_3Dangular: The full angular range is not used in this histogram - the PDF does not support this case" << endl;
		exit(1);
	}

	cout << "Finishing processing histo" << endl;
    }

}
// END AILSA


//Make the data point and parameter set
void Bd2JpsiKstar_sWave_Fs::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( cosThetaName );
	allObservables.push_back( phiName );
	allObservables.push_back( cosPsiName );
	allObservables.push_back( KstarFlavourName );
	// Need to think about additional parameters like
	// event-by-event propertime resolution and acceptance.
	// This will require event-by-event PDF normalisation,
	// but we are already doing this for tagging.

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName );
	parameterNames.push_back( Rpara_sqName );
	parameterNames.push_back( Rperp_sqName );
	parameterNames.push_back( As_sqName );
	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
	parameterNames.push_back( delta_sName );
	parameterNames.push_back( CspName );
	parameterNames.push_back( timeRes1Name );
	parameterNames.push_back( timeRes2Name );
	parameterNames.push_back( timeRes1FractionName );
	parameterNames.push_back( angAccI1Name );
	parameterNames.push_back( angAccI2Name );
	parameterNames.push_back( angAccI3Name );
	parameterNames.push_back( angAccI4Name );
	parameterNames.push_back( angAccI5Name );
	parameterNames.push_back( angAccI6Name );
	parameterNames.push_back( angAccI7Name );
	parameterNames.push_back( angAccI8Name );
	parameterNames.push_back( angAccI9Name );
	parameterNames.push_back( angAccI10Name );
	allParameters = ParameterSet(parameterNames);
}

//Destructor
Bd2JpsiKstar_sWave_Fs::~Bd2JpsiKstar_sWave_Fs()
{
	if( f != NULL ) f->Close();
	f = NULL;
	if( timeAcc != NULL ) delete timeAcc;
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bd2JpsiKstar_sWave_Fs::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	evaluationCacheValid = false;
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
	Rpara_sq   = allParameters.GetPhysicsParameter( Rpara_sqName )->GetValue();
	Rperp_sq   = allParameters.GetPhysicsParameter( Rperp_sqName )->GetValue();
	double Fs   = allParameters.GetPhysicsParameter( As_sqName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	delta_s = allParameters.GetPhysicsParameter( delta_sName )->GetValue();
	Csp = allParameters.GetPhysicsParameter( CspName )->GetValue();
	timeRes1 = allParameters.GetPhysicsParameter( timeRes1Name )->GetValue();
	timeRes2 = allParameters.GetPhysicsParameter( timeRes2Name )->GetValue();
	timeRes1Frac = allParameters.GetPhysicsParameter( timeRes1FractionName )->GetValue();
	angAccI1 = allParameters.GetPhysicsParameter( angAccI1Name )->GetValue();
	angAccI2 = allParameters.GetPhysicsParameter( angAccI2Name )->GetValue();
	angAccI3 = allParameters.GetPhysicsParameter( angAccI3Name )->GetValue();
	angAccI4 = allParameters.GetPhysicsParameter( angAccI4Name )->GetValue();
	angAccI5 = allParameters.GetPhysicsParameter( angAccI5Name )->GetValue();
	angAccI6 = allParameters.GetPhysicsParameter( angAccI6Name )->GetValue();
	angAccI7 = allParameters.GetPhysicsParameter( angAccI7Name )->GetValue();
	angAccI8 = allParameters.GetPhysicsParameter( angAccI8Name )->GetValue();
	angAccI9 = allParameters.GetPhysicsParameter( angAccI9Name )->GetValue();
	angAccI10 = allParameters.GetPhysicsParameter( angAccI10Name )->GetValue();

	//Rzero_sq = 1 - Rperp_sq - Rpara_sq;
	//Aperp_sq = Rperp_sq*(1-As_sq);
	//Apara_sq = Rpara_sq*(1-As_sq);
	//Azero_sq = Rzero_sq*(1-As_sq);

	Rzero_sq = 1 - Rperp_sq - Rpara_sq;
	Aperp_sq = Rperp_sq;
	Apara_sq = Rpara_sq;
	Azero_sq = Rzero_sq;
	As_sq = Fs/(1-Fs);

	AparaAperp = sqrt(Apara_sq)*sqrt(Aperp_sq);
	AzeroApara = sqrt(Azero_sq)*sqrt(Apara_sq);
	AzeroAperp = sqrt(Azero_sq)*sqrt(Aperp_sq);
	AparaAs    = sqrt(Apara_sq)*sqrt(As_sq);
	AperpAs	   = sqrt(Aperp_sq)*sqrt(As_sq);
	AzeroAs    = sqrt(Azero_sq)*sqrt(As_sq);

	return isOK;
}

//Return a list of parameters not to be integrated
vector<string> Bd2JpsiKstar_sWave_Fs::GetDoNotIntegrateList()
{
	vector<string> doNotIntList;
	doNotIntList.push_back( KstarFlavourName );
	//doNotIntList.push_back(timeName);
	return doNotIntList;
}

double Bd2JpsiKstar_sWave_Fs::q() const { return KstarFlavour;}

//Calculate the function value
double Bd2JpsiKstar_sWave_Fs::Evaluate(DataPoint * measurement)
{
	_datapoint = measurement;

	double returnValue;
	time = measurement->GetObservable( timeName )->GetValue();
	cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
	phi      = measurement->GetObservable( phiName )->GetValue();
	cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();
	KstarFlavour = measurement->GetObservable( KstarFlavourName )->GetValue();

	//cout << gamma << " " << Aperp_sq << " " << Azero_sq << endl;

	if(timeRes1Frac >= 0.9999)
	{
		// Set the member variable for time resolution to the first value and calculate
		timeRes = timeRes1;
		returnValue =  buildPDFnumerator();
	}
	else
	{
		// Set the member variable for time resolution to the first value and calculate
		timeRes = timeRes1;
		double val1 = buildPDFnumerator();
		// Set the member variable for time resolution to the second value and calculate
		timeRes = timeRes2;
		double val2 = buildPDFnumerator();
		//return timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
		returnValue = timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
	}


	if( componentIndex == 0 )
	{
		if( (returnValue <= 0.) || std::isnan(returnValue) )
		{
			cout << " Bd2JpsiKstar_sWave_Fs::Evaluate() returns <=0 or nan " << endl ;
			cout << " AT    " << Aperp_sq ;
			cout << " AP    " << Apara_sq ;
			cout << " A0    " << Azero_sq;
			cout << " As   " << As_sq;
			cout << "   Dperp    " << delta_perp;
			cout << "   Dpara    " << delta_para;
			cout << "   Ds     " << delta_s << endl;
			cout << "   gamma   " << gamma << endl;

			throw(10);
		}
	}

	return returnValue;
    
}


double Bd2JpsiKstar_sWave_Fs::buildPDFnumerator()
{
	// The angular functions f1->f6 as defined in roadmap Table 1.(same for Kstar)
	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10;

    Mathematics::getBs2JpsiPhiAngularFunctionsWithSwave( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, cosTheta, phi, cosPsi );
    //this->getAngularFunctionsTransversity( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, cosTheta, phi, cosPsi );

	// The time dependent amplitudes as defined in roadmap Eqns 48 -> 59  //No tagging so only need 2 (h± pg 72)
	// First for the B
	//double AzeroAzeroB, AparaAparaB, AperpAperpB, AsAsB;
	//double ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB;
	//double ReAparaAsB, ImAperpAsB, ReAzeroAsB;

	getTimeDependentAmplitudes( AzeroAzeroB, AparaAparaB, AperpAperpB
			, ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB
			, AsAsB, ReAparaAsB, ImAperpAsB, ReAzeroAsB
			);

	double v1(0.);
	if ( !_plotAllComponents )
	{
		switch (componentIndex)
		{
			case 1:
				v1  = f1 * AzeroAzeroB
					+ f2 * AparaAparaB
					+ f5 * ReAzeroAparaB;
				break;
			case 2:
				v1 = f3 * AperpAperpB;
				break;
			case 3:
				v1 = f7 * AsAsB;
				break;
			default:
				//q() tags the K* flavour - it changes the sign of f4, f6 and f9
				v1  = f1 * AzeroAzeroB
					+ f2 * AparaAparaB
					+ f3 * AperpAperpB
					+ f4 * ImAparaAperpB * q()
					+ f5 * ReAzeroAparaB
					+ f6 * ImAzeroAperpB * q()
					+ f7 * AsAsB
					+ Csp * f8 * ReAparaAsB
					+ Csp * f9 * ImAperpAsB * q()
					+ Csp * f10 * ReAzeroAsB;
				break;
		}
	}
	else
	{
		switch( componentIndex )
		{
			case 1:
				v1 = f1 * AzeroAzeroB;
				break;
			case 2:
				v1 = f2 * AparaAparaB;
				break;
			case 3:
				v1 = f3 * AperpAperpB;
				break;
			case 4:
				v1 = f4 * ImAparaAperpB * q();
				break;
			case 5:
				v1 = f5 * ReAzeroAparaB;
				break;
			case 6:
				v1 = f6 * ImAzeroAperpB * q();
				break;
			case 7:
				v1 = f7 * AsAsB;
				break;
			case 8:
				v1 = Csp * f8 * ReAparaAsB;
				break;
			case 9:
				v1 = Csp * f9 * ImAperpAsB * q();
				break;
			case 10:
				v1 = Csp * f10 * ReAzeroAsB;
				break;
			default:
				//q() tags the K* flavour - it changes the sign of f4, f6 and f9
				v1  = f1 * AzeroAzeroB
					+ f2 * AparaAparaB
					+ f3 * AperpAperpB
					+ f4 * ImAparaAperpB * q()
					+ f5 * ReAzeroAparaB
					+ f6 * ImAzeroAperpB * q()
					+ f7 * AsAsB
					+ Csp * f8 * ReAparaAsB
					+ Csp * f9 * ImAperpAsB * q()
					+ Csp * f10 * ReAzeroAsB;
				break;
		}
	}
    
	if( useTimeAcceptance() ) v1  = v1 * timeAcc->getValue(time);

    if( ! useFlatAngularDistribution ) v1  *=  angularFactor();
    
	return v1;
}


double Bd2JpsiKstar_sWave_Fs::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	_datapoint = measurement;

	double returnValue;
	IConstraint * timeBound = boundary->GetConstraint(timeconstraintName);

	time = measurement->GetObservable( timeName )->GetValue();
	KstarFlavour = measurement->GetObservable( KstarFlavourName )->GetValue();


	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		return -1.;
	}
	else
	{
		tlo = timeBound->GetMinimum();
		thi = timeBound->GetMaximum();
	}


	if(timeRes1Frac >= 0.9999)
	{
		// Set the member variable for time resolution to the first value and calculate
		timeRes = timeRes1;
		returnValue =  buildCompositePDFdenominator();

	}
	else
	{
		// Set the member variable for time resolution to the first value and calculate
		timeRes = timeRes1;
		double val1 = buildCompositePDFdenominator();
		// Set the member variable for time resolution to the second value and calculate
		timeRes = timeRes2;
		double val2 = buildCompositePDFdenominator();
		//return timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
		returnValue = timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;

	}


	if( componentIndex == 0 )
	{
		if( (returnValue <= 0.) || std::isnan(returnValue) )
		{
			cout << " Bd2JpsiKstar_sWave_Fs::Normalisation() returns <=0 or nan " << endl ;
			cout << " AT    " << Aperp_sq ;
			cout << " AP    " << Apara_sq ;
			cout << " A0    " << Azero_sq;
			cout << " As   " << As_sq;
			cout << "   Dperp    " << delta_perp;
			cout << "   Dpara    " << delta_para;
			cout << "   Ds     " << delta_s << endl;
			cout << "   gamma   " << gamma << endl;

			throw(1);
		}
	}

    return returnValue;
        
	//return -1 ;
}



//....................................................
// New method to calculate normalisation using a histogrammed "low-end" time acceptance function
// The acceptance function information is all contained in the timeAcceptance member object,

double Bd2JpsiKstar_sWave_Fs::buildCompositePDFdenominator( )
{
	double tlo_boundary = tlo ;
	double thi_boundary = thi ;
	double returnValue = 0;

	if( true /*useTimeAcceptance()*/ ) {                // Set to true because seleting false makes a single slice for 0 --> 14.
		//This loops over each time slice, does the normalisation between the limits, and accumulates
		for( unsigned int islice = 0; islice < timeAcc->numberOfSlices(); ++islice )
		{
			//Set the time integrals
			tlo = tlo_boundary > timeAcc->getSlice(islice)->tlow() ? tlo_boundary : timeAcc->getSlice(islice)->tlow() ;
			thi = thi_boundary < timeAcc->getSlice(islice)->thigh() ? thi_boundary : timeAcc->getSlice(islice)->thigh() ;
			if( thi> tlo ) returnValue+= this->buildPDFdenominator(  ) * timeAcc->getSlice(islice)->height() ;
		}
	}
	else {
		returnValue = this->buildPDFdenominator() ;
	}

	tlo = tlo_boundary;
	thi = thi_boundary ;
	return returnValue ;
}




double Bd2JpsiKstar_sWave_Fs::NormAnglesOnlyForAcceptanceWeights(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	_datapoint = measurement;
	(void) boundary;
	double returnValue;
	time = measurement->GetObservable( timeName )->GetValue();
	KstarFlavour = measurement->GetObservable( KstarFlavourName )->GetValue();
    
	//First job for any new set of parameters is to Cache the time integrals


	if(timeRes1Frac >= 0.9999)
	{
		// Set the member variable for time resolution to the first value and calculate
		timeRes = timeRes1;
		returnValue =  buildPDFdenominatorAngles();
	}
	else
	{
		// Set the member variable for time resolution to the first value and calculate
		timeRes = timeRes1;
		double val1 = buildPDFdenominatorAngles();
		// Set the member variable for time resolution to the second value and calculate
		timeRes = timeRes2;
		double val2 = buildPDFdenominatorAngles();
		//return timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
		returnValue = timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;

	}

	if( componentIndex == 0 )
	{
		if( (returnValue <= 0.) || std::isnan(returnValue) )
		{
			cout << " Bd2JpsiKstar_sWave_Fs::Normalisation() returns <=0 or nan " << endl ;
			cout << " AT    " << Aperp_sq ;
			cout << " AP    " << Apara_sq ;
			cout << " A0    " << Azero_sq;
			cout << " As   " << As_sq;
			cout << "   Dperp    " << delta_perp;
			cout << "   Dpara    " << delta_para;
			cout << "   Ds     " << delta_s << endl;
			cout << "   gamma   " << gamma << endl;

			throw(2);
		}
	}

	return returnValue;
}



double Bd2JpsiKstar_sWave_Fs::buildPDFdenominator()
{


	// The integrals of the time dependent amplitudes as defined in roadmap Eqns 48 -> 59
	getTimeAmplitudeIntegrals(  cachedAzeroAzeroIntB
			, cachedAparaAparaIntB
			, cachedAperpAperpIntB
			, cachedAparaAperpIntB
			, cachedAzeroAparaIntB
			, cachedAzeroAperpIntB
			, cachedAsAsIntB
			, cachedAparaAsIntB
			, cachedAperpAsIntB
			, cachedAzeroAsIntB
			);


	double v1= cachedAzeroAzeroIntB * angAccI1
		+ cachedAparaAparaIntB * angAccI2
		+ cachedAperpAperpIntB * angAccI3
		+ cachedAparaAperpIntB * angAccI4* q()
		+ cachedAzeroAparaIntB * angAccI5
		+ cachedAzeroAperpIntB * angAccI6 * q()
		+ cachedAsAsIntB * angAccI7
		+ Csp * cachedAparaAsIntB * angAccI8
		+ Csp * cachedAperpAsIntB * angAccI9 * q()
		+ Csp * cachedAzeroAsIntB * angAccI10;
	/*
	   if ( !_plotAllComponents )
	   {
	   switch( componentIndex )
	   {
	   case 1:
	   v1  = cachedAzeroAzeroIntB * angAccI1
	   + cachedAparaAparaIntB * angAccI2
	   + cachedAzeroAparaIntB * angAccI5;
	   break;
	   case 2:
	   v1 = cachedAperpAperpIntB * angAccI3;
	   break;
	   case 3:
	   v1 = cachedAsAsIntB * angAccI7;
	   break;
	   default:
	   v1 = cachedAzeroAzeroIntB * angAccI1
	   + cachedAparaAparaIntB * angAccI2
	   + cachedAperpAperpIntB * angAccI3
	   + cachedAparaAperpIntB * angAccI4* q()
	   + cachedAzeroAparaIntB * angAccI5
	   + cachedAzeroAperpIntB * angAccI6 * q()
	   + cachedAsAsIntB * angAccI7
	   + Csp * cachedAparaAsIntB * angAccI8
	   + Csp * cachedAperpAsIntB * angAccI9 * q()
	   + Csp * cachedAzeroAsIntB * angAccI10;
	   break;
	   }
	   }
	   else
	   {
	   switch ( componentIndex )
	   {
	   case 1:
	   v1 = cachedAzeroAzeroIntB * angAccI1;
	   break;
	   case 2:
	   v1 = cachedAparaAparaIntB * angAccI2;
	   break;
	   case 3:
	   v1 = cachedAperpAperpIntB * angAccI3;
	   break;
	   case 4:
	   v1 = cachedAparaAperpIntB * angAccI4* q();
	   break;
	   case 5:
	   v1 = cachedAzeroAparaIntB * angAccI5;
	   break;
	   case 6:
	   v1 = cachedAzeroAperpIntB * angAccI6 * q();
	   break;
	   case 7:
	   v1 = cachedAsAsIntB * angAccI7;
	   break;
	   case 8:
	   v1 = Csp * cachedAparaAsIntB * angAccI8;
	   break;
	   case 9:
	   v1 = Csp * cachedAperpAsIntB * angAccI9 * q();
	   break;
	   case 10:
	   v1 = Csp * cachedAzeroAsIntB * angAccI10;
	   break;
	   default:
	   v1 = cachedAzeroAzeroIntB * angAccI1
	   + cachedAparaAparaIntB * angAccI2
	   + cachedAperpAperpIntB * angAccI3
	   + cachedAparaAperpIntB * angAccI4* q()
	   + cachedAzeroAparaIntB * angAccI5
	   + cachedAzeroAperpIntB * angAccI6 * q()
	+ cachedAsAsIntB * angAccI7
		+ Csp * cachedAparaAsIntB * angAccI8
		+ Csp * cachedAperpAsIntB * angAccI9 * q()
		+ Csp * cachedAzeroAsIntB * angAccI10;
	break;
}
}
*/
return v1;
}

double Bd2JpsiKstar_sWave_Fs::buildPDFdenominatorAngles()  //test method
{

	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10;
    Mathematics::getBs2JpsiPhiAngularFunctionsWithSwave( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, cosTheta, phi, cosPsi );
    //this->getAngularFunctionsTransversity( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, cosTheta, phi, cosPsi );
    


	// The integrals of the time dependent amplitudes as defined in roadmap Eqns 48 -> 59
	getTimeDependentAmplitudes( AzeroAzeroB, AparaAparaB, AperpAperpB,
			ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB
			, AsAsB, ReAparaAsB, ImAperpAsB, ReAzeroAsB
			);


	double v1 =  AzeroAzeroB
		+ AparaAparaB
		+ AperpAperpB
		+ AsAsB
		;
	return v1;

}


void Bd2JpsiKstar_sWave_Fs::getTimeDependentAmplitudes(
		double & AzeroAzero
		, double & AparaApara
		, double & AperpAperp
		, double & ImAparaAperp
		, double & ReAzeroApara
		, double & ImAzeroAperp
		, double & AsAs
		, double & ReAparaAs
		, double & ImAperpAs
		, double & ReAzeroAs
		)
{
	// Quantities depending only on physics parameters can be cached
	if ( !evaluationCacheValid )
	{


		cachedAzero = sqrt( Azero_sq );
		cachedApara = sqrt( Apara_sq );
		cachedAperp = sqrt( Aperp_sq );
		cachedAs = sqrt (As_sq);

		cachedSinDeltaPerpPara	= sin( delta_perp - delta_para );
		cachedCosDeltaPara	= cos( delta_para );
		cachedSinDeltaPerp	= sin( delta_perp );
		cachedCosDeltaParaS	= cos( delta_para - delta_s );
		cachedSinDeltaPerpS 	= sin( delta_perp - delta_s );
		cachedCosDeltaS		= cos( delta_s);

		evaluationCacheValid = true;
	}



	// Now calculate the amplitudes

	double Exp = Mathematics::Exp(time, gamma, timeRes);

	AzeroAzero = Azero_sq * Exp;  // changed- see note 2009-015 eq 11-13
	AparaApara = Apara_sq * Exp;  //
	AperpAperp = Aperp_sq * Exp;  //
	AsAs = As_sq * Exp;

	ImAparaAperp = cachedApara*cachedAperp * cachedSinDeltaPerpPara * Exp;    //See http://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=33933 page14
	ReAzeroApara = cachedAzero*cachedApara * cachedCosDeltaPara * Exp;
	ImAzeroAperp = cachedAzero*cachedAperp * cachedSinDeltaPerp * Exp;
	ReAparaAs = cachedApara*cachedAs * cachedCosDeltaParaS * Exp; //AILSA_ NOT SURE
	ImAperpAs = cachedAperp*cachedAs * cachedSinDeltaPerpS * Exp; //AILSA
	ReAzeroAs = cachedAzero*cachedAs * cachedCosDeltaS * Exp; //AILSA


	//if ( std::isnan(ImAparaAperp)) cout << Azero_sq << " " << Apara_sq << " " << Aperp_sq << " " << Exp << endl;

	return;
}

void Bd2JpsiKstar_sWave_Fs::getTimeAmplitudeIntegrals(
		double & AzeroAzeroInt
		, double & AparaAparaInt
		, double & AperpAperpInt
		, double & AparaAperpInt
		, double & AzeroAparaInt
		, double & AzeroAperpInt
		, double & AsAsInt
		, double & AparaAsInt
		, double & AperpAsInt
		, double & AzeroAsInt
		)
{

	if ( !evaluationCacheValid )
	{


		cachedAzero = sqrt( Azero_sq );
		cachedApara = sqrt( Apara_sq );
		cachedAperp = sqrt( Aperp_sq );
		cachedAs = sqrt (As_sq);

		cachedSinDeltaPerpPara  = sin( delta_perp - delta_para );
		cachedCosDeltaPara      = cos( delta_para );
		cachedSinDeltaPerp      = sin( delta_perp );
		cachedCosDeltaParaS     = cos( delta_para - delta_s );
		cachedSinDeltaPerpS     = sin( delta_perp - delta_s );
		cachedCosDeltaS         = cos( delta_s);

		evaluationCacheValid = true;
	}




	double ExpInt = Mathematics::ExpInt(tlo, thi, gamma, timeRes);


	AzeroAzeroInt = Azero_sq * ExpInt;
	AparaAparaInt = Apara_sq  * ExpInt;
	AperpAperpInt = Aperp_sq * ExpInt;
	AsAsInt = As_sq * ExpInt;

	AparaAperpInt = AparaAperp * cachedSinDeltaPerpPara * ExpInt;
	AzeroAparaInt = AzeroApara * cachedCosDeltaPara * ExpInt;
	AzeroAperpInt = AzeroAperp * cachedSinDeltaPerp * ExpInt;
	AparaAsInt = AparaAs * cachedCosDeltaParaS * ExpInt;
	AperpAsInt = AperpAs * cachedSinDeltaPerpS * ExpInt;
	AzeroAsInt = AzeroAs * cachedCosDeltaS * ExpInt;

	return;
}

double Bd2JpsiKstar_sWave_Fs::angularFactor( )
{
	double returnValue=0.;

	int globalbin=-1;
	int xbin=-1, ybin=-1, zbin=-1;
	double num_entries_bin=-1.;

	//Observable* cos1Obs = cos1Obs = _datapoint->GetObservable( cosPsiName );
	Observable* cos1Obs = _datapoint->GetObservable( cosPsiName );

	if( cos1Obs->GetBkgBinNumber() != -1 )
	{
		return cos1Obs->GetBkgAcceptance();
	}
	else
	{
		Observable* cos2Obs = _datapoint->GetObservable( cosThetaName );
		Observable* phiObs = _datapoint->GetObservable( phiName );

		//Find global bin number for values of angles, find number of entries per bin, divide by volume per bin and normalise with total number of entries in the histogram
		xbin = xaxis->FindFixBin( cosPsi ); if( xbin > nxbins ) xbin = nxbins; if( xbin == 0 ) xbin = 1;
		ybin = yaxis->FindFixBin( cosTheta ); if( ybin > nybins ) ybin = nybins; if( ybin == 0 ) ybin = 1;
		zbin = zaxis->FindFixBin( phi ); if( zbin > nzbins ) zbin = nzbins; if( zbin == 0 ) zbin = 1;

		globalbin = histo->GetBin( xbin, ybin, zbin );
		num_entries_bin = (double)histo->GetBinContent(globalbin);

		//Angular factor normalized with phase space of histogram and total number of entries in the histogram
		returnValue = (double)num_entries_bin;// / (deltax * deltay * deltaz) / (double)total_num_entries;
		//returnValue = (double)num_entries_bin / histo->Integral();

		cos1Obs->SetBkgBinNumber( xbin ); cos1Obs->SetBkgAcceptance( returnValue );
		cos2Obs->SetBkgBinNumber( ybin ); cos2Obs->SetBkgAcceptance( returnValue );
		phiObs->SetBkgBinNumber( zbin ); phiObs->SetBkgAcceptance( returnValue );
	}

	return returnValue;
}

/*
//Angular distribution function
double Bd2JpsiKstar_sWave_Fs::angularFactor( )
{
double returnValue=0.;

int globalbin=-1;
int xbin=-1, ybin=-1, zbin=-1;
double num_entries_bin=-1.;

if( useFlatAngularDistribution ) {
returnValue = 1.0; /// 8.0 / TMath::Pi() ;
}
else {
//Find global bin number for values of angles, find number of entries per bin, divide by volume per bin and normalise with total number of entries in the histogram
xbin = xaxis->FindFixBin( cosPsi ); if( xbin > nxbins ) xbin = nxbins;
ybin = yaxis->FindFixBin( cosTheta ); if( ybin > nybins ) ybin = nybins;
zbin = zaxis->FindFixBin( phi ); if( zbin > nzbins ) zbin = nzbins;

globalbin = histo->GetBin( xbin, ybin, zbin );
num_entries_bin = histo->GetBinContent(globalbin);

//Angular factor normalized with phase space of histogram and total number of entries in the histogram
returnValue = num_entries_bin; /// (deltax * deltay * deltaz) / total_num_entries ;
}

return returnValue;
}*/

vector<string> Bd2JpsiKstar_sWave_Fs::PDFComponents()
{
	vector<string> this_component_list;
	this_component_list.push_back( "0" );
	if( !_plotAllComponents ) {
		this_component_list.push_back( "P-even" );
		this_component_list.push_back( "P-odd" );
		this_component_list.push_back( "S-wave" );
	}
	else {
		this_component_list.push_back( "f1" );
		this_component_list.push_back( "f2" );
		this_component_list.push_back( "f3" );
		this_component_list.push_back( "f4" );
		this_component_list.push_back( "f5" );
		this_component_list.push_back( "f6" );
		this_component_list.push_back( "f7" );
		this_component_list.push_back( "f8" );
		this_component_list.push_back( "f9" );
		this_component_list.push_back( "f10" );
	}
	return this_component_list;
}

string Bd2JpsiKstar_sWave_Fs::GetComponentName( ComponentRef* Component )
{
	return Component->getComponentName();
}

double Bd2JpsiKstar_sWave_Fs::EvaluateComponent(DataPoint * measurement, ComponentRef* Component)
{
	componentIndex = Component->getComponentNumber();
	if ( !_plotAllComponents )
	{
		if( componentIndex == -1 )
		{
			string ComponentName = Component->getComponentName();
			if( ComponentName.compare( "P-even" ) == 0 )
			{
				Component->setComponentNumber( 1 );
				componentIndex = 1;
			}
			else if( ComponentName.compare( "P-odd" ) == 0 )
			{
				Component->setComponentNumber( 2 );
				componentIndex = 2;
			}
			else if( ComponentName.compare( "S-wave" ) == 0 )
			{
				Component->setComponentNumber( 3 );
				componentIndex = 3;
			}
			else
			{
				Component->setComponentNumber( 0 );
				componentIndex = 0;
			}
		}
	}
	else
	{
		if( componentIndex == -1 )
		{
			string ComponentName = Component->getComponentName();
			if( ComponentName.compare( "f1" ) == 0 )
			{
				Component->setComponentNumber( 1 );
				componentIndex = 1;
			}
			else if( ComponentName.compare( "f2" ) == 0 )
			{
				Component->setComponentNumber( 2 );
				componentIndex = 2;
			}
			else if( ComponentName.compare( "f3" ) == 0 )
			{
				Component->setComponentNumber( 3 );
				componentIndex = 3;
			}
			else if( ComponentName.compare( "f4" ) == 0 )
			{
				Component->setComponentNumber( 4 );
				componentIndex = 4;
			}
			else if( ComponentName.compare( "f5" ) == 0 )
			{
				Component->setComponentNumber( 5 );
				componentIndex = 5;
			}
			else if( ComponentName.compare( "f6" ) == 0 )
			{
				Component->setComponentNumber( 6 );
				componentIndex = 6;
			}
			else if( ComponentName.compare( "f7" ) == 0 )
			{
				Component->setComponentNumber( 7 );
				componentIndex = 7;
			}
			else if( ComponentName.compare( "f8" ) == 0 )
			{
				Component->setComponentNumber( 8 );
				componentIndex = 8;
			}
			else if( ComponentName.compare( "f9" ) == 0 )
			{
				Component->setComponentNumber( 9 );
				componentIndex = 9;
			}
			else if( ComponentName.compare( "f10" ) == 0 )
			{
				Component->setComponentNumber( 10 );
				componentIndex = 10;
			}
			else
			{
				Component->setComponentNumber( 0 );
				componentIndex = 0;
			}
		}
	}
	return this->Evaluate( measurement );
}

void Bd2JpsiKstar_sWave_Fs::getAngularFunctionsHelicity( double &f1, double &f2, double &f3, double &f4, double &f5, double &f6, double &f7, double &f8, double &f9, double &f10, double _cosThetal, double _phi, double _cosThetaK )
{
   	vector<double> input ;
    input.push_back(_cosThetaK) ;
    input.push_back(_cosThetal) ;
    input.push_back( _phi ) ;
    f1 = Bs2JpsiPhi_Angular_Terms::HangleFactorA0A0( input ) ;
    f2 = Bs2JpsiPhi_Angular_Terms::HangleFactorAPAP( input ) ;
    f3 = Bs2JpsiPhi_Angular_Terms::HangleFactorATAT( input ) ;
    f4 = Bs2JpsiPhi_Angular_Terms::HangleFactorImAPAT( input ) ;
    f5 = Bs2JpsiPhi_Angular_Terms::HangleFactorReA0AP( input ) ;
    f6 = Bs2JpsiPhi_Angular_Terms::HangleFactorImA0AT( input ) ;
    f7 = Bs2JpsiPhi_Angular_Terms::HangleFactorASAS( input ) ;
    f8 = Bs2JpsiPhi_Angular_Terms::HangleFactorReASAP( input ) ;
    f9 = Bs2JpsiPhi_Angular_Terms::HangleFactorImASAT( input ) ;
    f10 = Bs2JpsiPhi_Angular_Terms::HangleFactorReASA0( input ) ;
    
    return ;
}

void Bd2JpsiKstar_sWave_Fs::getAngularFunctionsTransversity( double &f1, double &f2, double &f3, double &f4, double &f5, double &f6, double &f7, double &f8, double &f9, double &f10, double _cosTheta, double _phi, double _cosPsi )
{
   	vector<double> input ;
    input.push_back(_cosTheta) ;
    input.push_back(_cosPsi) ;
    input.push_back( _phi ) ;
    f1 = Bs2JpsiPhi_Angular_Terms::TangleFactorA0A0( input ) ;
    f2 = Bs2JpsiPhi_Angular_Terms::TangleFactorAPAP( input ) ;
    f3 = Bs2JpsiPhi_Angular_Terms::TangleFactorATAT( input ) ;
    f4 = Bs2JpsiPhi_Angular_Terms::TangleFactorImAPAT( input ) ;
    f5 = Bs2JpsiPhi_Angular_Terms::TangleFactorReA0AP( input ) ;
    f6 = Bs2JpsiPhi_Angular_Terms::TangleFactorImA0AT( input ) ;
    f7 = Bs2JpsiPhi_Angular_Terms::TangleFactorASAS( input ) ;
    f8 = Bs2JpsiPhi_Angular_Terms::TangleFactorReASAP( input ) ;
    f9 = Bs2JpsiPhi_Angular_Terms::TangleFactorImASAT( input ) ;
    f10 = Bs2JpsiPhi_Angular_Terms::TangleFactorReASA0( input ) ;
    
    return ;
}
