/** @class Bd2JpsiKstar_sWave Bd2JpsiKstar_sWave.cpp
 *
 *  RapidFit PDF for Bd2JpsiPhi with average angular acceptance input as a paramter
 *
 *  @author Ailsa Sparkes asparkes@cern.ch
 *  @date 2011-01-26
 */

#include "TMath.h"
#include <cmath>

#include "Bd2JpsiKstar_sWave.h"
#include "Mathematics.h"
#include "SlicedAcceptance.h"
#include <iostream>
#include <fstream>
#include "math.h"
#include "TFile.h"
#include "TH3D.h"
#include "TROOT.h"

PDF_CREATOR( Bd2JpsiKstar_sWave );
//Constructor
Bd2JpsiKstar_sWave::Bd2JpsiKstar_sWave(PDFConfigurator* configurator ) :
	cachedAzeroAzeroIntB(), cachedAparaAparaIntB(), cachedAperpAperpIntB(), cachedAparaAperpIntB(), cachedAzeroAparaIntB(), cachedAzeroAperpIntB(),
	cachedAsAsIntB(), cachedAparaAsIntB(), cachedAperpAsIntB(), cachedAzeroAsIntB(), AzeroAzeroB(), AparaAparaB(), AperpAperpB(), AsAsB(), ImAparaAperpB(),
	ReAzeroAparaB(), ImAzeroAperpB(), ReAparaAsB(), ImAperpAsB(), ReAzeroAsB(), cachedSinDeltaPerpPara(), cachedCosDeltaPara(), cachedSinDeltaPerp(),
	cachedCosDeltaParaS(), cachedSinDeltaPerpS(), cachedCosDeltaS(), cachedAzero(), cachedApara(), cachedAperp(), cachedAs(), timeAcc(NULL),
	// Physics parameters
	  gammaName     ( "gamma" )
	, Azero_sqName ( "Azero_sq")
	, Apara_sqName  ( "Apara_sq" )
	, Aperp_sqName  ( "Aperp_sq" )
	, As_sqName	( "As_sq" )
	, delta_zeroName( "delta_zero" )
	, delta_paraName( "delta_para" )
	, delta_perpName( "delta_perp" )
	, delta_sName	( "delta_s" )
	, angAccI1Name	( "angAccI1" )
	, angAccI2Name	( "angAccI2" )
	, angAccI3Name	( "angAccI3" )
	, angAccI4Name	( "angAccI4" )
	, angAccI5Name	( "angAccI5" )
	, angAccI6Name	( "angAccI6" )
	, angAccI7Name  ( "angAccI7" )
	, angAccI8Name  ( "angAccI8" )
	, angAccI9Name  ( "angAccI9" )
	, angAccI10Name ( "angAccI10" )
	, timeRes1Name	( "timeResolution1" )
	, timeRes2Name	( "timeResolution2" )
	, timeRes1FractionName	( "timeResolution1Fraction" )
	, _useTimeAcceptance(false)

	// Observables (What we want to gain from the pdf after inserting physics parameter values)
	, normalisationCacheValid(false)
	, evaluationCacheValid(false)
	, timeName	( "time" )
	, cosThetaName	( "cosTheta" )
	, phiName	( "phi" )
	, cosPsiName	( "cosPsi" )
	, KstarFlavourName  ( "KstarFlavour" )
	//, timeres	( "resolution" )

	, timeconstraintName( "time" )
	, gamma(), Azero_sq(), Apara_sq(), Aperp_sq(), As_sq(), AzeroApara(), AzeroAperp(), AparaAperp(), AparaAs(), AperpAs(), AzeroAs(),
	delta_zero(), delta_para(), delta_perp(), delta_s(), omega(), timeRes(), timeRes1(), timeRes2(), timeRes1Frac(), angAccI1(), angAccI2(),
	angAccI3(), angAccI4(), angAccI5(), angAccI6(), angAccI7(), angAccI8(), angAccI9(), angAccI10(), Ap_sq(), Ap(), time(), cosTheta(), phi(),
	cosPsi(), KstarFlavour(), tlo(), thi(), useFlatAngularDistribution(true)
{
	MakePrototypes();
	_useTimeAcceptance = configurator->isTrue( "UseTimeAcceptance" ) ;

       if( useTimeAcceptance() ) {
                        timeAcc = new SlicedAcceptance( 0., 14.0, 0.0171 ) ;
                        cout << "Bd2JpsiKstar_sWave:: Constructing timeAcc: Upper time acceptance beta=0.0171 [0 < t < 14] " << endl ;
                }
        else {
                        timeAcc = new SlicedAcceptance( 0., 14. ) ;
                        cout << "Bd2JpsiKstar_sWave:: Constructing timeAcc: DEFAULT FLAT [0 < t < 14]  " << endl ;
        }


//AILSA
	//Find name of histogram needed to define 3-D angular distribution
        string fileName = configurator->getConfigurationValue( "AngularDistributionHistogram" ) ;

        //Initialise depending upon whether configuration parameter was found
        if( fileName == "" )
        {
                cout << "   No AngularAcceptanceHisto found" << endl ;
                useFlatAngularDistribution = true ;
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
                        cerr << "\t\tAlternativley make sure your XML points to the correct file, or that the file is in the current working directory\n" <<
 endl;
                        exit(-89);
                }

                if( fullFileName.empty() || !local_fail2 )
                {
                        fullFileName = fileName;
                }
	
		TFile* f =  TFile::Open(fullFileName.c_str());
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
                if ((xmax-xmin) < 2. || (ymax-ymin) < 2. || (zmax-zmin) < 2.*TMath::Pi() )
                {
                        cout << "In LongLivedBkg_3Dangular::LongLivedBkg_3Dangular: The full angular range is not used in this histogram - the PDF does not support this case" << endl;
                        exit(1);
                }

                cout << "Finishing processing histo" << endl;
        }

}
// END AILSA


//Make the data point and parameter set
void Bd2JpsiKstar_sWave::MakePrototypes()
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
	parameterNames.push_back( Apara_sqName );
	parameterNames.push_back( Aperp_sqName );
	parameterNames.push_back( As_sqName );
	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
	parameterNames.push_back( delta_sName );
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
Bd2JpsiKstar_sWave::~Bd2JpsiKstar_sWave()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bd2JpsiKstar_sWave::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	evaluationCacheValid = false;
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
	Apara_sq   = allParameters.GetPhysicsParameter( Apara_sqName )->GetValue();
	Aperp_sq   = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
	As_sq   = allParameters.GetPhysicsParameter( As_sqName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	delta_s = allParameters.GetPhysicsParameter( delta_sName )->GetValue();
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

	Azero_sq = 1 - Aperp_sq - Apara_sq - As_sq ;
	AparaAperp = sqrt(Apara_sq)*sqrt(Aperp_sq);
	AzeroApara = sqrt(Azero_sq)*sqrt(Apara_sq);
	AzeroAperp = sqrt(Azero_sq)*sqrt(Aperp_sq);
	AparaAs    = sqrt(Apara_sq)*sqrt(As_sq);
	AperpAs	   = sqrt(Aperp_sq)*sqrt(As_sq);
	AzeroAs    = sqrt(Azero_sq)*sqrt(As_sq);

	return isOK;
}

//Return a list of parameters not to be integrated
vector<string> Bd2JpsiKstar_sWave::GetDoNotIntegrateList()
{
	vector<string> doNotIntList;
	doNotIntList.push_back(timeName);

	return doNotIntList;
}

double Bd2JpsiKstar_sWave::q() const { return KstarFlavour;}

//Calculate the function value
double Bd2JpsiKstar_sWave::Evaluate(DataPoint * measurement)
{

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

		 if( (returnValue <= 0.) || std::isnan(returnValue) ) {
                cout << " Bd2JpsiKstar_sWave::Evaluate() returns <=0 or nan " << endl ;
                cout << " AT    " << Aperp_sq ;
                cout << " AP    " << Apara_sq ;
                cout << " A0    " << Azero_sq;
                cout << " As   " << As_sq;
                cout << "   Dperp    " << delta_perp;
                cout << "   Dpara    " << delta_para;
                cout << "   Ds     " << delta_s << endl;
                cout << "   gamma   " << gamma << endl;

                throw 10 ;


}

 return returnValue;
}


double Bd2JpsiKstar_sWave::buildPDFnumerator()
{
	// The angular functions f1->f6 as defined in roadmap Table 1.(same for Kstar)
	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10;
	Mathematics::getBs2JpsiPhiAngularFunctionsWithSwave( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, cosTheta, phi, cosPsi );

	// The time dependent amplitudes as defined in roadmap Eqns 48 -> 59  //No tagging so only need 2 (h± pg 72)
	// First for the B
	double AzeroAzeroB, AparaAparaB, AperpAperpB, AsAsB;
	double ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB;
	double ReAparaAsB, ImAperpAsB, ReAzeroAsB;

	getTimeDependentAmplitudes( AzeroAzeroB, AparaAparaB, AperpAperpB
			, ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB
			, AsAsB, ReAparaAsB, ImAperpAsB, ReAzeroAsB
			);

	//q() tags the K* flavour - it changes the sign of f4, f6 and f9
	double v1 = f1 * AzeroAzeroB
		+ f2 * AparaAparaB
		+ f3 * AperpAperpB
		+ f4 * ImAparaAperpB * q()
		+ f5 * ReAzeroAparaB
		+ f6 * ImAzeroAperpB * q()
		+ f7 * AsAsB
		+ f8 * ReAparaAsB
		+ f9 * ImAperpAsB * q()
		+ f10 * ReAzeroAsB
		;
	if( useTimeAcceptance() ) v1  = v1 * timeAcc->getValue(time);
	
	v1  *=  angularFactor();
	return v1;
}


double Bd2JpsiKstar_sWave::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{


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


		 if( (returnValue <= 0.) || std::isnan(returnValue) ) {
                cout << " Bd2JpsiKstar_sWave::Normalisation() returns <=0 or nan " << endl ;
                cout << " AT    " << Aperp_sq ;
                cout << " AP    " << Apara_sq ;
                cout << " A0    " << Azero_sq;
                cout << " As   " << As_sq;
                cout << "   Dperp    " << delta_perp;
                cout << "   Dpara    " << delta_para;
                cout << "   Ds     " << delta_s << endl;
                cout << "   gamma   " << gamma << endl;

                exit(1) ;
}

return returnValue;
}



//....................................................
// New method to calculate normalisation using a histogrammed "low-end" time acceptance function
// The acceptance function information is all contained in the timeAcceptance member object,

double Bd2JpsiKstar_sWave::buildCompositePDFdenominator( )
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




double Bd2JpsiKstar_sWave::NormAnglesOnlyForAcceptanceWeights(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
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


                 if( (returnValue <= 0.) || std::isnan(returnValue) ) {
                cout << " Bd2JpsiKstar_sWave::Normalisation() returns <=0 or nan " << endl ;
                cout << " AT    " << Aperp_sq ;
                cout << " AP    " << Apara_sq ;
                cout << " A0    " << Azero_sq;
                cout << " As   " << As_sq;
                cout << "   Dperp    " << delta_perp;
                cout << "   Dpara    " << delta_para;
                cout << "   Ds     " << delta_s << endl;
                cout << "   gamma   " << gamma << endl;

                exit(1) ;
}

return returnValue;
}



double Bd2JpsiKstar_sWave::buildPDFdenominator()
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




	double v1 = cachedAzeroAzeroIntB * angAccI1
		+ cachedAparaAparaIntB * angAccI2
		+ cachedAperpAperpIntB * angAccI3
		+ cachedAparaAperpIntB * angAccI4* q()
		+ cachedAzeroAparaIntB * angAccI5
		+ cachedAzeroAperpIntB * angAccI6 * q()
		+ cachedAsAsIntB * angAccI7
		+ cachedAparaAsIntB * angAccI8
		+ cachedAperpAsIntB * angAccI9 * q()
		+ cachedAzeroAsIntB * angAccI10
		;
	return v1;




}

double Bd2JpsiKstar_sWave::buildPDFdenominatorAngles()  //test method
{

	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10;
        Mathematics::getBs2JpsiPhiAngularFunctionsWithSwave( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, cosTheta, phi, cosPsi );


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


void Bd2JpsiKstar_sWave::getTimeDependentAmplitudes(
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

void Bd2JpsiKstar_sWave::getTimeAmplitudeIntegrals(
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

//Angular distribution function
double Bd2JpsiKstar_sWave::angularFactor( )
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
}

