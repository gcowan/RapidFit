// $Id: DPHistoBackground.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class DPHistoBackground DPHistoBackground.cpp
 *
 *  PDF for  long lived background with 3D histogram angular description. A 3D root Histogram required provided by user and input in the xml file. The transversity angles in the order (x,y,z) --> (cosPsi, cosTheta, phi)
 *
 *  @author Ailsa Sparkes
 *  @date 2011-05-30
 */

#include "DPHistoBackground.h"
#include "Mathematics.h"
#include <iostream>
#include <fstream>
#include "math.h"
#include "TMath.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TH3D.h"
#include "TAxis.h"
#include "TH1.h"

PDF_CREATOR( DPHistoBackground );

//Constructor
DPHistoBackground::DPHistoBackground( PDFConfigurator* config ) :

	// Observables
	  massName		( config->getName("m23") )
	, cosTheta1Name		( config->getName("cosTheta1") )
	, phiName		( config->getName("phi") )
	, cosTheta2Name		( config->getName("cosTheta2") )
	, mass(), cos1(), cos2(), phi(), histo()
	, xaxis(), yaxis(), zaxis(), maxis()
	, nxbins(), nybins(), nzbins(), nmbins()
	, xmin(), xmax()
	, ymin(), ymax()
	, zmin(), zmax()
	, mmin(), mmax()
	, deltax(), deltay(), deltaz(), deltam()
	, total_num_entries(), useFlatAngularDistribution(true)
	, histogramFile()
{

	cout << "Constructing PDF: DPHistoBackground  " << endl ;

	//Make prototypes
	MakePrototypes();
	
	//Find name of histogram needed to define 3-D angular distribution
	string fileName = config->getConfigurationValue( "AngularDistributionHistogram" ) ;

	//Initialise depending upon whether configuration parameter was found
	if( fileName == "" )
	{
		cout << "   No AngularDistributionHistogram found: using flat background " << endl ;
		useFlatAngularDistribution = true ;
	}
	else
	{
		cout << "   AngularDistributionHistogram requested: " << fileName << endl ;
		useFlatAngularDistribution = false ;

		//File location
		ifstream input_file;
		input_file.open( fileName.c_str(), ifstream::in );
		input_file.close();

		bool local_fail = input_file.fail();

		if( !getenv("RAPIDFITROOT") && local_fail )
		{
			cerr << "\n" << endl;
			//cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cerr << "$RAPIDFITROOT NOT DEFINED, PLEASE DEFINE IT SO I CAN USE ACCEPTANCE DATA" << endl;
			//cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cerr << "\n" << endl;
			exit(-987);
		}

		if( getenv("RAPIDFITROOT") )
		{
			string path( getenv("RAPIDFITROOT") ) ;

			cout << "RAPIDFITROOT defined as: " << path << endl;

			fullFileName = path+"/pdfs/configdata/"+fileName ;

			input_file.open( fullFileName.c_str(), ifstream::in );
			input_file.close();
		}
		bool elsewhere_fail = input_file.fail();

		if( elsewhere_fail && local_fail )
		{
			cerr << "\n\tFileName:\t" << fullFileName << "\t NOT FOUND PLEASE CHECK YOUR RAPIDFITROOT" << endl;
			cerr << "\t\tAlternativley make sure your XML points to the correct file, or that the file is in the current working directory\n" << endl;
			exit(-89);
		}

		if( fullFileName.empty() || !local_fail )
		{
			fullFileName = fileName;
		}

		//Read in histo
		histogramFile = TFile::Open(fullFileName.c_str());
		histo = (THnSparse*)histogramFile->Get("histo_4var_SB"); //(fileName.c_str())));

		// cos k
		xaxis = histo->GetAxis(0);
		cout << "X axis Name: " << xaxis->GetName() << "\tTitle: " << xaxis->GetTitle() << endl;
		xmin = xaxis->GetXmin();
		xmax = xaxis->GetXmax();
		nxbins = xaxis->GetNbins();
		cout << "X axis Min: " << xmin << "\tMax: " << xmax << "\tBins: " << nxbins << endl;
		deltax = (xmax-xmin)/nxbins;

		// cos mu
		yaxis = histo->GetAxis(1);
		cout << "Y axis Name: " << yaxis->GetName() << "\tTitle: " << yaxis->GetTitle() << endl;
		ymin = yaxis->GetXmin();
		ymax = yaxis->GetXmax();
		nybins = yaxis->GetNbins();
		cout << "Y axis Min: " << ymin << "\tMax: " << ymax << "\tBins: " << nybins << endl;
		deltay = (ymax-ymin)/nybins;

		// delta phi
		zaxis = histo->GetAxis(2);
		cout << "Z axis Name: " << zaxis->GetName() << "\tTitle: " << zaxis->GetTitle() << endl;
		zmin = zaxis->GetXmin();
		zmax = zaxis->GetXmax();
		nzbins = zaxis->GetNbins();
		cout << "Z axis Min: " << zmin << "\tMax: " << zmax << "\tBins: " << nzbins << endl;
		deltaz = (zmax-zmin)/nzbins;

		// m Kpi
		maxis = histo->GetAxis(3);
		cout << "M axis Name: " << maxis->GetName() << "\tTitle: " << maxis->GetTitle() << endl;
		mmin = maxis->GetXmin();
		mmax = maxis->GetXmax();
		nmbins = maxis->GetNbins();
		cout << "M axis Min: " << mmin << "\tMax: " << mmax << "\tBins: " << nmbins << endl;
		deltam = (mmax-mmin)/nmbins;

		//method for Checking whether histogram is sensible

		total_num_entries = histo->GetEntries();
		int total_num_bins = nxbins * nybins * nzbins * nmbins;
		int sum = 0;

		vector<int> zero_bins;

		//loop over each bin in histogram and print out how many zero bins there are
		int idx[4] = {0,0,0,0};
		for (int i=1; i < nxbins+1; ++i)
		{
			for (int j=1; j < nybins+1; ++j)
			{
				for (int k=1; k < nzbins+1; ++k)
				{
					for (int l=1; l < nmbins+1; ++l)
					{
						idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
						double bin_content = histo->GetBinContent(idx);
						//cout << "Bin content: " << bin_content << endl;
						if(bin_content<=0)
						{
							zero_bins.push_back(1);
						}
						//cout << " Zero bins " << zero_bins.size() << endl;}
						else if (bin_content>0)
						{
							sum += (int) bin_content;
						}
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
		if ((xmax-xmin) < 2. || (ymax-ymin) < 2. || (zmax-zmin) < 2.*3.14159 || (mmax-mmin) < 1.35 )
		{
			cout << "In DPHistoBackground::DPHistoBackground: The full angular range is not used in this histogram - the PDF does not support this case" << endl;
			exit(1);
		}

		cout << "Finishing processing histo" << endl;
	}
}

// Copy
DPHistoBackground::DPHistoBackground( const DPHistoBackground& input ) : BasePDF( (BasePDF) input ),
xaxis(), yaxis(), zaxis(), maxis(),
nxbins(input.nxbins), nybins(input.nybins), nzbins(input.nzbins), nmbins(input.nmbins),
xmin(input.xmin), xmax(input.xmax), 
ymin(input.ymin), ymax(input.ymax), 
zmin(input.zmin), zmax(input.zmax), 
mmin(input.zmin), mmax(input.zmax), 
deltax(input.deltax), deltay(input.deltay), deltaz(input.deltaz), deltam(input.deltam),
total_num_entries(input.total_num_entries), useFlatAngularDistribution(input.useFlatAngularDistribution),
massName(input.massName), cosTheta1Name(input.cosTheta1Name), phiName(input.phiName), cosTheta2Name(input.cosTheta2Name), 
mass(input.mass), cos1(input.cos1), cos2(input.cos2), phi(input.phi),
histogramFile(), histo()
,fullFileName(input.fullFileName)
{
	if ( !useFlatAngularDistribution ) {
                histogramFile = TFile::Open(fullFileName.c_str());
                histo = (THnSparse*)histogramFile->Get("histo_4var_SB");
		xaxis = histo->GetAxis(0);
		yaxis = histo->GetAxis(1);
		zaxis = histo->GetAxis(2);
		maxis = histo->GetAxis(3);
	}
}

//Destructor
DPHistoBackground::~DPHistoBackground()
{
        if( !useFlatAngularDistribution ) {
		delete histogramFile;
		delete histo;
	}
}

//Make the data point and parameter set
void DPHistoBackground::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( massName );
	allObservables.push_back( cosTheta1Name );
	allObservables.push_back( phiName );
	allObservables.push_back( cosTheta2Name );

	//Make the parameter set
	vector<string> parameterNames;
	allParameters = ParameterSet(parameterNames);
}

bool DPHistoBackground::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	return isOK;
}

//Main method to build the PDF return value
double DPHistoBackground::Evaluate(DataPoint * measurement)
{
	// Observable
	mass = measurement->GetObservable( massName )->GetValue();
	cos1 = measurement->GetObservable( cosTheta1Name )->GetValue();
	phi  = measurement->GetObservable( phiName )->GetValue();
	cos2 = measurement->GetObservable( cosTheta2Name )->GetValue();
	return this->angleMassFactor();
}


// Normlisation
double DPHistoBackground::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void)boundary;
	return buildPDFdenominator();
}

double DPHistoBackground::buildPDFdenominator()
{
	return 1.;
}


//Angular distribution function
double DPHistoBackground::angleMassFactor( )
{
	double returnValue = 0.;

	int globalbin = -1;
	int xbin = -1, ybin = -1, zbin = -1, mbin = -1;
	double num_entries_bin = -1.;

	if( useFlatAngularDistribution ) {
		returnValue = 1.0 / 8.0 / TMath::Pi(); // probably need some mass factor in here also
	}
	else {
		//Find global bin number for values of angles, find number of entries per bin, divide by volume per bin and normalise with total number of entries in the histogram
		xbin = xaxis->FindFixBin( cos2 ); if( xbin > nxbins ) xbin = nxbins;
		ybin = yaxis->FindFixBin( cos1 ); if( ybin > nybins ) ybin = nybins;
		zbin = zaxis->FindFixBin( phi  ); if( zbin > nzbins ) zbin = nzbins;
		mbin = maxis->FindFixBin( mass ); if( mbin > nmbins ) mbin = nmbins;
		
		int idx[4] = { xbin, ybin, zbin, mbin};
		globalbin = (int)histo->GetBin( idx );
		num_entries_bin = histo->GetBinContent(globalbin);

		//Angular factor normalized with phase space of histogram and total number of entries in the histogram
		returnValue = num_entries_bin / (deltax * deltay * deltaz * deltam) / total_num_entries ;
	}
	return returnValue;
}
