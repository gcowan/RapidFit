/**
  @class SlicedAcceptance

  A class for holding a sliced propertime acceptance

  @author Pete Clarke
  @data 2011-06-07
  */


#include "SlicedAcceptance.h"
#include "StringProcessing.h"
#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"

#include <cmath>
#include <stdio.h>
#include <vector>
#include <string>

using namespace::std;

AcceptanceSlice::AcceptanceSlice( const AcceptanceSlice& input ) :
	_tlow(input._tlow), _thigh(input._thigh), _height(input._height)
{
}

//............................................
// Constructor for flat acceptance
SlicedAcceptance::SlicedAcceptance( double tl, double th, bool quiet ) :
	slices(), nullSlice( new AcceptanceSlice(0.,0.,0.) ), tlow(tl), thigh(th), beta(0), _sortedSlices(false), maxminset(false), t_min(0.), t_max(0.)
{
	(void)quiet;
	//Reality checks
	if( tl > th )
	{
		cout << "SlicedAcceptance::SlicedAcceptance :  tl > th  - exiting" << endl;
		exit(1);
	}

	//Create single slice with acceptance=1.0
	slices.push_back( new AcceptanceSlice( tlow, thigh, 1.0 ) );

	//....done.....

	_sortedSlices = true;
}

SlicedAcceptance::SlicedAcceptance( const SlicedAcceptance& input ) :
	slices(), nullSlice( new AcceptanceSlice(0.,0.,0.) ), tlow( input.tlow ), thigh( input.thigh ), beta( input.beta ), _sortedSlices( input._sortedSlices ), maxminset(input.maxminset), t_min(input.t_min), t_max(input.t_max), _hasChecked(input._hasChecked), _storedDecision(input._storedDecision)
{
	for( unsigned int i=0; i< input.slices.size(); ++i )
	{
		if( (input.slices[i]) != NULL ) slices.push_back( new AcceptanceSlice( *(input.slices[i]) ) );
		else slices.push_back( NULL );
	}
}

SlicedAcceptance::~SlicedAcceptance()
{
	if( nullSlice != NULL ) delete nullSlice;
	while( !slices.empty() )
	{
		if( slices.back() != NULL ) delete slices.back();
		slices.pop_back();
	}
}

//............................................
// Constructor for simple upper time acceptance only
SlicedAcceptance::SlicedAcceptance( double tl, double th, double b, bool quiet ) :
	slices(), nullSlice(new AcceptanceSlice(0.,0.,0.)), tlow(tl), thigh(th), beta(b), _sortedSlices(false), maxminset(false), t_min(0.), t_max(0.), _hasChecked(false), _storedDecision(false)
{
	//Reality checks
	if( tlow > thigh )
	{
		cout << "SlicedAcceptance::SlicedAcceptance :  tl > th  - exiting" << endl;
		exit(1);
	}
	if( beta < 0. )
	{
		cout << "SlicedAcceptance::SlicedAcceptance :  beta < 0  - exiting" << endl;
		exit(1);
	}

	// The upper time acceptance formula is (1 - beta*time) so determine highest and lowest acceptance
	double accUp = 1.0 - beta*tlow;
	double accLow = 1.0 - beta*thigh;


	//Create N more, equispaced
	int N = 15 ;
	double dh = (accUp-accLow) / (double)N;
	double dt = (thigh - tlow) / (double)N;

	//First slice which goes underneath the whole slope at the accLow level +dh/2
	slices.push_back( new AcceptanceSlice( tlow, thigh, (accLow+(dh/2.0)) ) );

	for( int is=1; is <= (N-1) ; is++ )    //The N-1 is important as the N one is included as the +dh/2 in first slice
	{
		double this_th = tlow + is*dt;
		slices.push_back( new AcceptanceSlice( tlow, this_th, dh ) );
		//cout << tlow << "   " <<  this_th << "  " << dh << endl;
	}

	//....done.....
	_sortedSlices = this->isSorted();

	if( _sortedSlices )
	{
		if( !quiet ) cout << "Sliced Acceptance is using sorted horizontal slices" << endl;
	}
	else
	{
		if( !quiet ) cout << "Sliced Acceptance is NOT using sorted horizontal slices" << endl;
	}
}


//............................................
// Constructor for simple 2010 version of lower time acceptance only
SlicedAcceptance::SlicedAcceptance( string s, bool quiet ) :
	slices(), nullSlice(new AcceptanceSlice(0.,0.,0.)), tlow(), thigh(), beta(), _sortedSlices(false), maxminset(false), t_min(0.), t_max(0.), _hasChecked(false), _storedDecision(false)
{
	(void)s;
	int N = 31;

	double ts[31] = {0, 0.05,  0.1, 0.15, 0.2, 0.25, 0.3, 0.34, 0.38, 0.42, 0.46, 0.5, 0.54, 0.58, 0.62, 0.66,
		0.7,  0.75,  0.8, 0.85 , 0.9, 0.95, 1.0, 1.05, 1.1,  1.2, 1.3, 1.4,1.6 , 1.8, 2};

	double ac[31] = {0.000122069 , 0.00436829 , 0.022288 , 0.0633692 , 0.132094 , 0.225401 , 0.321979 , 0.409769 ,
		0.494438 , 0.570185 , 0.637997 , 0.695043 , 0.744576 , 0.784799 , 0.816805 , 0.844606 , 0.869917 ,
		0.891761 , 0.910512 , 0.926418 , 0.937439 , 0.947105 , 0.953599 , 0.960726 , 0.968278 , 0.976149 ,
		0.980624 , 0.98685 , 0.991765 , 0.995575 , 1 };


	for( int is=0; is < N ; is++ )
	{
		double this_tlow = ts[is];
		double this_thigh = 14.0;
		double height = is > 0 ? ac[is] - ac[is-1] : ac[is];
		slices.push_back( new AcceptanceSlice( this_tlow, this_thigh, height ) );
	}

	//....done.....
	_sortedSlices = this->isSorted();

	if( _sortedSlices )
	{
		if( !quiet ) cout << "Sliced Acceptance is using sorted horizontal slices" << endl;
	}
	else
	{
		if( !quiet ) cout << "Sliced Acceptance is NOT using sorted horizontal slices" << endl;
	}
}

//............................................
// Constructor for accpetance from a file
SlicedAcceptance::SlicedAcceptance( string type, string fileName, bool quiet ) :
	slices(), nullSlice(new AcceptanceSlice(0.,0.,0.)), tlow(), thigh(), beta(), _sortedSlices(false), maxminset(false), t_min(0.), t_max(0.), _hasChecked(false), _storedDecision(false)
{
	(void)type;
	if( type != "File" ) {   }//do nothing for now

	string fullFileName = StringProcessing::FindFileName( fileName, quiet );

	if( !quiet ) cout << "Opening: " << fullFileName << endl;

	int number_of_lines = 0;
	std::string line;
	std::ifstream myfile( fullFileName.c_str() );

	while (std::getline(myfile, line))	++number_of_lines;

	myfile.close();


	ifstream in;
	in.open(fullFileName.c_str());
	if( in.fail() )
	{
		cout << "SlicedAcceptance::SlicedAcceptance : failed to open acceptance file  '  " << fullFileName  << "  '  " << endl;
		exit(1);
	}

	vector<double> lowEdge; vector<double> hiEdge; vector<double> binContent;
	int ok = true;
	unsigned int this_line=0;
	while (ok)
	{
		lowEdge.push_back(stream(in));
		hiEdge.push_back(stream(in));
		binContent.push_back(stream(in) );
		++this_line;
		if( in.eof()  == true || in.good() == false || this_line == (unsigned int)number_of_lines ) ok =false;
	}
	in.close();

	for( unsigned int is=0; is < lowEdge.size() ; ++is )
	{
		double this_tlow = lowEdge[is];
		double this_thigh = hiEdge[is];
		double height = binContent[is];
		//cout << " Adding slice " << tlow << " /  " << thigh << " /  " << height << endl ;
		slices.push_back( new AcceptanceSlice( this_tlow, this_thigh, height ) );
		//cout << this_tlow << "  " << this_thigh << "  " << height << endl;
	}
	//exit(0);
	if( !quiet ) cout << "Time Acc Slices: " << lowEdge.size() << endl;
	if( lowEdge.size() == 1 )
	{
		cout << "SlicedAcceptance: SERIOUS ERROR" << endl;
		exit(0);
	}
	//....done.....


	_sortedSlices = this->isSorted();

	if( _sortedSlices )
	{
		if( !quiet ) cout << "Sliced Acceptance is using sorted horizontal slices" << endl;
	}
	else
	{
		if( !quiet ) cout << "Sliced Acceptance is NOT using sorted horizontal slices" << endl;
	}
}

//............................................
// Constructor for accpetance from a ROOT Tfile
SlicedAcceptance::SlicedAcceptance( string type, string fileName,string histName, bool fluctuate, bool quiet ) :
	slices(), nullSlice(new AcceptanceSlice(0.,0.,0.)), tlow(), thigh(), beta(), _sortedSlices(false), maxminset(false), t_min(0.), t_max(0.), _hasChecked(false), _storedDecision(false)
{
	if(!quiet) cout << "Root file being used for acceptance" << endl;
	(void)type;
	if( type != "RootFile" ) {   }//do nothing for now

	string fullFileName = StringProcessing::FindFileName( fileName, quiet );

	if( !quiet ) cout << "Opening: " << fullFileName << endl;

	TFile* file = TFile::Open(TString(fullFileName));
 	if(!quiet) cout << "File " << fullFileName << " opened!" << endl;
	TH1F* histo = (TH1F*)file->Get(TString(histName));
	if(!quiet) cout << "Histo " << histName << " opened!" << endl;
	histo->Draw();
	if(!quiet) cout << "Histo " << histName << " drawn!" << endl;
//	histo->Sumw2();

	if(fluctuate){
	cout << "WARNING! You have fluctuated the acceptance." << endl;
	cout << "WARNING! This is for systematic studies only. " << endl;
	cout << "WARNING! Projections and pull fits will have a different fluctuated acceptance to the PDF you fit with." << endl;
	cout << "WARNING! ONLY USE FluctuateAcceptance:True FOR SYSTEMATIC STUDIES" << endl;
	TRandom3 * rng = new TRandom3(0);	
	//Randomly fluctuate bin contents within error: 
	for (int l = 1; l <= histo->GetNbinsX(); ++l){
	if(!quiet)	cout << "Bin content and error before: " << histo->GetBinContent(l) << "+/-" << histo->GetBinError(l);
		histo->SetBinContent(l,rng->Gaus(histo->GetBinContent(l),histo->GetBinError(l)));
	if(!quiet)	cout <<  " and after: " << histo->GetBinContent(l) << "+/-" << histo->GetBinError(l) << endl;
	}
	delete rng;
	}

	histo->Scale(1./(histo->GetBinContent(histo->GetMaximumBin())));
	histo->SetMinimum(0);


	double maxend = histo->GetBinLowEdge(histo->GetNbinsX()) + histo->GetBinWidth(histo->GetNbinsX());
	double height;
	double start;
	double end  = histo->GetBinLowEdge(histo->GetNbinsX()) + histo->GetBinWidth(histo->GetNbinsX());
	double dheight;
	for (int l = 1; l <= histo->GetNbinsX(); ++l){
		height = histo->GetBinContent(l);
		dheight = height;
		for (int n = l; n>0; n--){
			if(histo->GetBinContent(n)<height){
				dheight = height - histo->GetBinContent(n);
				cout << l << " " << n << " " << dheight << endl;
				break;
			}
		}
		start = histo->GetBinLowEdge(l);
		end = maxend;
		for (int m = l; m <= histo->GetNbinsX(); ++m){
			double thisbinheight = histo->GetBinContent(m);
			if(thisbinheight<height){
				end = histo->GetBinLowEdge(m);
				break;
			}
		}
		slices.push_back( new AcceptanceSlice( start, end, dheight ) );
	if(!quiet)	cout << start << "	" << end << "      " << dheight << endl;
	}
	histo->Delete();
//	delete histo;
	file->Close();
	delete file;

	if( !quiet ) cout << "Time Acc Slices: " << slices.size() << endl;
	if( slices.size() == 1 )
	{
		cout << "SlicedAcceptance: SERIOUS ERROR" << endl;
		exit(0);
	}
	//....done.....


	_sortedSlices = this->isSorted();

	if( _sortedSlices )
	{
		if( !quiet ) cout << "Sliced Acceptance is using sorted horizontal slices" << endl;
	}
	else
	{
		if( !quiet ) cout << "Sliced Acceptance is NOT using sorted horizontal slices" << endl;
	}
}

//............................................
// Return numerator for evaluate
double SlicedAcceptance::getValue( const double t ) const
{
	returnValue = 0;
	for( unsigned int is = 0; is < slices.size() ; ++is )
	{
		if( (t>=slices[is]->tlow()) && (t<slices[is]->thigh()) ) returnValue += slices[is]->height();
		if( _sortedSlices ) if( t<slices[is]->thigh() ) break;
	}
	if( fabs( t - slices.front()->tlow() ) < 1E-6 ) returnValue += slices.front()->height();
	if( fabs( t - slices.back()->thigh() ) < 1E-6 ) returnValue += slices.back()->height();
	return returnValue;
}

double SlicedAcceptance::getValue( const Observable* time, const double timeOffset ) const
{
	//return this->getValue( time->GetValue() - timeOffset );

	if( time->GetBinNumber() != -1 )
	{
		//#pragma GCC diagnostic ignored "-Wfloat-equal"
		if( abs(time->GetOffSet() - timeOffset) <= numeric_limits<float>::epsilon() ) return time->GetAcceptance();
		//#pragma GCC diagnostic pop
		else
		{
			time->SetBinNumber( -1 );
			time->SetOffSet( -999. );
		}
	}

	double t = time->GetValue() - timeOffset;
	returnValue = 0;
	finalBin = -1;

	if( time->GetBinNumber() < 0 )
	{
		if( t < this->GetMin() )
		{
			if( time->GetValue() > this->GetMin() )
			{
				cout << "TIME OFFSET PUSHING VALUE BELOW ACCEPTANCE HISTOGRAM!!!" << endl;
			}
			else
			{
				cout << "TIME BELOW ACCEPTANCE HISTO!!!" << endl;
			}

			cout << " time: " << time->GetValue() << " offset: " << timeOffset << " min: " << this->GetMin() << endl;
			this->Print();
			throw(-987643);
		}
		if( t > this->GetMax() )
		{
			if( time->GetValue() > this->GetMax() )
			{
				cout << "TIME OFFSET PUSHING VALUE ABOVE ACCEPTANCE HISTOGRAM!!!" << endl;
			}
			else
			{
				cout << "TIME ABOVE ACCEPTANCE HISTO!!!" << endl;
			}
			cout << " time: " << time->GetValue() << " offset: " << timeOffset << " max: " << this->GetMax() << endl;

			this->Print();
			throw(-987643);
		}
	}

	if( time->GetBinNumber() >= 0 )
	{
		finalBin = time->GetBinNumber();
	}
	else
	{
		finalBin = (int)this->findSliceNum( time, timeOffset );
	}

	//cout << finalBin << endl;

	if( _sortedSlices )
	{
		returnValue += slices[(unsigned)finalBin]->height();
	}
	else
	{
		for( unsigned int is=0; is < slices.size(); ++is )
		{
			if( (t >= slices[is]->tlow() ) && ( t < slices[is]->thigh() ) )
			{
				returnValue += slices[is]->height();
			}
		}

		if( abs(t - slices.back()->thigh()) <= numeric_limits<float>::epsilon() )
		{
			returnValue += slices.back()->height();
		}
	}

	if( fabs( time->GetValue() - slices[0]->tlow() ) < 1E-5 ) returnValue = slices[0]->height();

	time->SetBinNumber( finalBin );
	time->SetAcceptance( returnValue );
	time->SetOffSet( timeOffset );
	return returnValue;
}

//............................................
// Return the number of slices
unsigned int SlicedAcceptance::numberOfSlices() const
{
	return (unsigned)slices.size();
}

//............................................
// Return a slice
AcceptanceSlice * SlicedAcceptance::getSlice( const unsigned int s ) const
{
	if( s>=(unsigned)slices.size() ) return nullSlice;
	return slices[s];
}

//............................................
//Helpers	Only Used in Constructor
double SlicedAcceptance::stream(ifstream& thisStream)
{
	double tmpVal;
	thisStream >> tmpVal;
	return tmpVal;
}

unsigned int SlicedAcceptance::findSliceNum( const Observable* time, const double timeOffset ) const
{
	double t = time->GetValue() - timeOffset;
	if( time->GetBinNumber() >= 0 ) return (unsigned) time->GetBinNumber();
	unsigned int thisNum=0;
	for( unsigned int is=0; is < slices.size(); ++is )
	{
		if( (t >= slices[is]->tlow() ) && ( t < slices[is]->thigh() ) )
		{
			break;
		}
		else
		{
			continue;
		}
		++thisNum;
	}
	if( abs(t - slices.back()->thigh()) <= numeric_limits<float>::epsilon() ) --thisNum;
	time->SetBinNumber( (int)thisNum );
	return thisNum;
}

bool SlicedAcceptance::isSorted() const
{
	if( _hasChecked ) return _storedDecision;
	else
	{
		bool isReallySorted = true;
		double last_thigh = -999.;
		for( unsigned int i=0; i< slices.size(); ++i )
		{
			if( isReallySorted )
			{
				if( last_thigh <= slices[i]->tlow() )
				{
					isReallySorted = true;
					last_thigh = slices[i]->thigh();
				}
				else
				{
					isReallySorted = false;
				}
				continue;
			}
			else
			{
				break;
			}
		}
		_storedDecision = isReallySorted;
		_hasChecked = true;
		return isReallySorted;
	}
}

bool SlicedAcceptance::GetIsSorted() const
{
	return _sortedSlices;
}

double SlicedAcceptance::GetMax() const
{
	if( !maxminset ) this->FindMaxMin();
	return t_max;
}

double SlicedAcceptance::GetMin() const
{
	if( !maxminset ) this->FindMaxMin();
	return t_min;
}

void SlicedAcceptance::FindMaxMin() const
{
	double temp_max = slices[0]->tlow();
	double temp_min = slices[0]->thigh();
	if( _sortedSlices )
	{
		temp_min = (*slices.begin())->tlow();
		temp_max = (*slices.back()).thigh();
	}
	else
	{
		for( unsigned int i=0; i< slices.size(); ++i )
		{
			if( slices[i]->tlow() < temp_min ) temp_min = slices[i]->tlow();
			if( slices[i]->thigh() > temp_max ) temp_max = slices[i]->thigh();
		}
	}
	t_min = temp_min;
	t_max = temp_max;
	maxminset = true;
}

void SlicedAcceptance::Print() const
{
	cout << "SlicedAcceptance:" << endl;
	for( unsigned int i=0; i< slices.size(); ++i )
	{
		cout << "Slice: " << i << endl;
		slices[i]->Print();
		cout << endl;
	}
}

