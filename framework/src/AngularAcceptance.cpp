/**
  @class AngularAcceptance

  A class for holding angular acceptance information

  @author Pete Clarke
  @data 2012-03-29
 */


#include "AngularAcceptance.h"
#include "Mathematics.h"
#include "StringProcessing.h"

#include "TChain.h"

using namespace::std;

AngularAcceptance::AngularAcceptance( const AngularAcceptance& input ) :
	_af1(input._af1), _af2(input._af2), _af3(input._af3), _af4(input._af4), _af5(input._af5), _af6(input._af6)
	, _af7(input._af7), _af8(input._af8), _af9(input._af9), _af10(input._af10), useFlatAngularAcceptance(input.useFlatAngularAcceptance)
	, histo(input.histo) , xaxis(input.xaxis), yaxis(input.yaxis), zaxis(input.zaxis), nxbins(input.nxbins), nybins(input.nybins), nzbins(input.nzbins)
	, xmin(input.xmin), xmax(input.xmax), ymin(input.ymin), ymax(input.ymax), zmin(input.zmin), zmax(input.zmax), deltax(input.deltax), deltay(input.deltay), deltaz(input.deltaz)
	, total_num_entries(input.total_num_entries), average_bin_content(input.average_bin_content)
{
	TString Histo_Name="histo_";
	TString XAxis_Name="Xaxis_";
	TString YAxis_Name="Yaxis_";
	TString ZAxis_Name="Zaxis_";
	size_t uniqueNum = reinterpret_cast<size_t>(this);
	Histo_Name+=uniqueNum;
	XAxis_Name+=uniqueNum;
	YAxis_Name+=uniqueNum;
	ZAxis_Name+=uniqueNum;
	if( input.histo != NULL )
	{
		histo = (TH3D*)input.histo->Clone(Histo_Name);
		xaxis = (TAxis*)histo->GetXaxis()->Clone(XAxis_Name);
		yaxis = (TAxis*)histo->GetYaxis()->Clone(YAxis_Name);
		zaxis = (TAxis*)histo->GetZaxis()->Clone(ZAxis_Name);
	}
}

AngularAcceptance::~AngularAcceptance()
{
	if( histo != NULL ) delete histo;
	//	Objects claim to be looked after by the Histogram so cannot be 'doubly deleted'
	//	This sounds suspicious but we can't delete them here due to segfaults
	//if( xaxis != NULL ) delete xaxis;
	//if( yaxis != NULL ) delete yaxis;
	//if( zaxis != NULL ) delete zaxis;
}

//............................................
// Constructor for accpetance from a file
AngularAcceptance::AngularAcceptance( string fileName, bool useHelicityBasis, bool quiet ) :
	_af1(1), _af2(1), _af3(1), _af4(0), _af5(0), _af6(0), _af7(1), _af8(0), _af9(0), _af10(0), useFlatAngularAcceptance(false)
	, histo(), xaxis(), yaxis(), zaxis(), nxbins(), nybins(), nzbins(), xmin(), xmax(), ymin(), ymax(), zmin(), zmax(), deltax(), deltay(), deltaz(), total_num_entries(), average_bin_content()
{

	//Initialise depending upon whether configuration parameter was found
	if( fileName == "" )
	{
		useFlatAngularAcceptance = true ;
	}
	else
	{
		useFlatAngularAcceptance = false ;

		// Get the acceptance histogram

		string fullFileName = StringProcessing::FindFileName( fileName, quiet );//this->openFile( fileName, quiet ) ;

		TFile* f =  TFile::Open(fullFileName.c_str());

		//cout << " AngularAcceptance::AngularAcceptance fileName: " <<  fullFileName << endl;

		if( useHelicityBasis ) {
			histo = (TH3D*) f->Get("helacc"); //(fileName.c_str())));
			if( !quiet ) cout << " AngularAcceptance::  Using heleicity basis" << endl ;
		}
		else {
			histo = (TH3D*) f->Get("tracc"); //(fileName.c_str())));
			if( !quiet ) cout << " AngularAcceptance::  Using transversity basis" << endl ;
		}

		if( histo == NULL ) histo = (TH3D*) f->Get("acc");

		if( histo == NULL )
		{
			cerr << "Cannot Open a Valid NTuple" << endl;
			exit(0);
		}
		histo->SetDirectory(0);
		size_t uniqueNum = reinterpret_cast<size_t>(this);
		TString Histo_Name="Histo_";Histo_Name+=uniqueNum;
		this->processHistogram( quiet );


		// Get the 10 angular factors

		TChain* decayTree ;
		decayTree = new TChain("tree");
		decayTree->Add(fullFileName.c_str());

		vector<double> *pvect = new vector<double>() ;
		decayTree->SetBranchAddress("weights", &pvect ) ;
		decayTree->GetEntry(0);

		_af1=(*pvect)[0], _af2=(*pvect)[1], _af3=(*pvect)[2], _af4=(*pvect)[3], _af5=(*pvect)[4],
			_af6=(*pvect)[5], _af7=(*pvect)[6], _af8=(*pvect)[7], _af9=(*pvect)[8], _af10=(*pvect)[9] ;

		//for( int ii=0; ii <10; ii++) {
		//	cout << "AcceptanceWeight "<<ii+1<< " = "  << (*pvect)[ii] << endl ;
		//}
		f->Close();
	}

}



//............................................
// Return numerator for evaluate
double AngularAcceptance::getValue( double cosPsi, double cosTheta, double phi ) const
{
	if( useFlatAngularAcceptance ) return 1. ;

	double returnValue=0.;

	int globalbin=-1;
	int xbin=-1, ybin=-1, zbin=-1;
	double num_entries_bin=-1.;

	//Find global bin number for values of angles, find number of entries per bin, divide by volume per bin and normalise with total number of entries in the histogram
	xbin = xaxis->FindFixBin( cosPsi ); if( xbin > nxbins ) xbin = nxbins;
	ybin = yaxis->FindFixBin( cosTheta ); if( ybin > nybins ) ybin = nybins;
	zbin = zaxis->FindFixBin( phi ); if( zbin > nzbins ) zbin = nzbins;

	globalbin = histo->GetBin( xbin, ybin, zbin );
	num_entries_bin = histo->GetBinContent(globalbin);

	returnValue = num_entries_bin; /// (deltax * deltay * deltaz) / total_num_entries ;

	return returnValue / average_bin_content;
}

double AngularAcceptance::getValue( Observable* cosPsi, Observable* cosTheta, Observable* phi ) const
{
	if( useFlatAngularAcceptance ) return 1.;

	int psi_num = cosPsi->GetBinNumber();
	//int theta_num = cosTheta->GetBinNumber();
	//int phi_num = phi->GetBinNumber();

	if( psi_num > -1  ) //	&& theta_num > -1) && phi_num > -1 )	These are implicitly +ve if psi is when using this class
	{
		//	By definition doesn't matter which Observable we interrogate here,
		//	The user can be intentionally stupid and break this, but we trust that they won't
		return cosPsi->GetAcceptance();
	}
	else
	{
		//	This has to be here to protect ROOT from breaking everything
		//	GetBinContent is NOT a const function!!!
		//	It will break the copy of the histogram in memory if you request an object out of scope
		int xbin=-1, ybin=-1, zbin=-1;
		xbin = xaxis->FindFixBin( cosPsi->GetValue() ); if( xbin > nxbins ) xbin = nxbins;
		ybin = yaxis->FindFixBin( cosTheta->GetValue() ); if( ybin > nybins ) ybin = nybins;
		zbin = zaxis->FindFixBin( phi->GetValue() ); if( zbin > nzbins ) zbin = nzbins;

		cosPsi->SetBinNumber( xbin );
		cosTheta->SetBinNumber( ybin );
		phi->SetBinNumber( zbin );

		int globalbin = histo->GetBin( xbin, ybin, zbin );
		double num_entries_bin = histo->GetBinContent(globalbin);

		double acc = (double)num_entries_bin / average_bin_content;

		cosPsi->SetAcceptance( acc );
		cosTheta->SetAcceptance( acc );
		phi->SetAcceptance( acc );

		return acc;
	}

	return -1.;
}

//............................................
// Open the input file containing the acceptance
string AngularAcceptance::openFile( string fileName, bool quiet )
{
	ifstream input_file2;
	input_file2.open( fileName.c_str(), ifstream::in );
	input_file2.close();
	bool local_fail2 = input_file2.fail();

	if( !quiet )
	{
		cout << "Looking For: " << fileName << endl;
		if( !local_fail2 ) cout << "Found Locally!" << endl;
	}

	string fileName_pwd = "pdfs/configdata/";
	fileName_pwd.append( fileName );
	input_file2.open( fileName_pwd.c_str(), ifstream::in );
	input_file2.close();
	bool pwd_fail = input_file2.fail();

	if( !getenv("RAPIDFITROOT") && ( local_fail2 && pwd_fail ) )
	{
		cerr << "\n" << endl;
		//cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cerr << "$RAPIDFITROOT NOT DEFINED, PLEASE DEFINE IT SO I CAN FIND ACCEPTANCE FILE" << endl;
		//cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cerr << "\n" << endl;
		exit(-987);
	}

	string fullFileName;

	if( (local_fail2 && !pwd_fail) || !getenv("RAPIDFITROOT") ) fullFileName=fileName_pwd;

	if( getenv("RAPIDFITROOT") )
	{
		string path( getenv("RAPIDFITROOT") ) ;

		//cout << "RAPIDFITROOT defined as: " << path << endl;
		fullFileName = path+"/pdfs/configdata/"+fileName ;

		input_file2.open( fullFileName.c_str(), ifstream::in );
		input_file2.close();
	}
	bool elsewhere_fail = input_file2.fail();

	if( elsewhere_fail && ( local_fail2 && pwd_fail ) )
	{
		cerr << "AngularAcceptance::AngularAcceptance: \n\tFileName:\t" << fullFileName << "\t NOT FOUND PLEASE CHECK YOUR RAPIDFITROOT" << endl;
		cerr << "\t\tAlternativley make sure your XML points to the correct file, or that the file is in the current working directory\n" <<
			endl;
		exit(-89);
	}


	if( fullFileName.empty() || !local_fail2 )
	{
		fullFileName = fileName;
	}

	return fullFileName;
}



//............................................
// Open the input file containing the acceptance
void AngularAcceptance::processHistogram( bool quiet )
{
	TString XAxis_Name="XAxis_";
	TString YAxis_Name="YAxis_";
	TString ZAxis_Name="ZAxis_";
	size_t uniqueNum = reinterpret_cast<size_t>(this);
	XAxis_Name+=uniqueNum;
	YAxis_Name+=uniqueNum;
	ZAxis_Name+=uniqueNum;

	xaxis = histo->GetXaxis();
	xaxis->SetName(XAxis_Name);
	xmin = xaxis->GetXmin();
	xmax = xaxis->GetXmax();
	nxbins = histo->GetNbinsX();
	deltax = (xmax-xmin)/nxbins;
	if( !quiet ) cout << " X axis Name: " << xaxis->GetName() << "\tTitle: " << xaxis->GetTitle() << "\t\t" << "X axis Min: " << xmin << "\tMax: " << xmax << "\tBins: " << nxbins << endl;

	yaxis = histo->GetYaxis();
	yaxis->SetName(YAxis_Name);
	ymin = yaxis->GetXmin();
	ymax = yaxis->GetXmax();
	nybins = histo->GetNbinsY();
	deltay = (ymax-ymin)/nybins;
	if( !quiet ) cout << " Y axis Name: " << yaxis->GetName() << "\tTitle: " << yaxis->GetTitle() << "\t\t" << "Y axis Min: " << ymin << "\tMax: " << ymax << "\tBins: " << nybins << endl;

	zaxis = histo->GetZaxis();
	zaxis->SetName(ZAxis_Name);
	zmin = zaxis->GetXmin();
	zmax = zaxis->GetXmax();
	nzbins = histo->GetNbinsZ();
	deltaz = (zmax-zmin)/nzbins;
	if( !quiet ) cout << " Z axis Name: " << zaxis->GetName() << "\tTitle: " << zaxis->GetTitle() << "\t\t" << "Z axis Min: " << zmin << "\tMax: " << zmax << "\tBins: " << nzbins << "\t\t\t";

	//method for Checking whether histogram is sensible

	total_num_entries = histo->GetEntries();
	int total_num_bins = nxbins * nybins * nzbins;
	double sum = 0;

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
				else if (bin_content>0)
				{
					sum += bin_content;
				}
			}
		}
	}

	average_bin_content = sum / total_num_bins;

	if( !quiet ) cout << "Total Bins: " << total_num_bins << "\tEmpty Bins: " << zero_bins.size() << "\tAvg Bin Content: " << average_bin_content << endl;

	//cout << "For total number of bins " << total_num_bins << " there are " << zero_bins.size() << " bins containing zero events "  << endl;
	//cout << "Average number of entries of non-zero bins: " << average_bin_content <<  endl;

	// Check.  This order works for both bases since phi is always the third one.
	if( (xmax-xmin) < 2. || (ymax-ymin) < 2. || (zmax-zmin) < (2.*TMath::Pi()-0.01) )
	{
		cout << endl << "In AngularAcceptance::processHistogram: The full angular range is not used in this histogram - the PDF does not support this case" << endl;
		exit(1);
	}

	//cout << "Finishing processing angular acceptance histo" << endl;
}


