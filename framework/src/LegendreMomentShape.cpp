#include "LegendreMomentShape.h"
#include "DPHelpers.hh"
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "TBranch.h"
double LegendreMomentShape::Moment(const int l, const int i, const int k, const int j, const double mKK_mapped, const double phi, const double ctheta_1, const double ctheta_2)
{
	if(j < k)
		return 0;
	if(std::abs(mKK_mapped) > 1)
		return 0; 
	double Q_l  = gsl_sf_legendre_Pl     (l, mKK_mapped );
	double P_i  = gsl_sf_legendre_Pl     (i, ctheta_2   );
	double Y_jk = gsl_sf_legendre_sphPlm (j, k, ctheta_1);
	if(k != 0)
		Y_jk *= sqrt(2) * cos(k * phi);
	return Q_l * P_i * Y_jk;
}
LegendreMomentShape::LegendreMomentShape() : init(true), copied(false)
{
}
LegendreMomentShape::LegendreMomentShape(std::string filename) : init(true), copied(false)
{
	Open(filename);
}
LegendreMomentShape::LegendreMomentShape(const LegendreMomentShape& copy) :
	  mKK_min(copy.mKK_min)
	, mKK_max(copy.mKK_max)
	, coeffs(copy.coeffs)
	, init(copy.init)
	, copied(true)
{
}
LegendreMomentShape::~LegendreMomentShape()
{
}
void LegendreMomentShape::Open(const std::string filename)
{
	TFile* file;
	if(!filename.empty())
	{
		std::cout << "Opening " << filename << std::endl;
		file = TFile::Open(filename.c_str());
	}
	else return;
	if(file->IsZombie())
	{
		std::cerr << "No file found. Defaulting to uniform shape." << std::endl;
		delete file;
		return;
	}
	TTree* tree = (TTree*)file->Get("LegendreMomentsTree");
	if(tree == nullptr) throw std::runtime_error("LegendreMomentsTree not found");
	tree->SetBranchAddress("mKK_min",&mKK_min);
	tree->SetBranchAddress("mKK_max",&mKK_max);
	std::string limbranchtitle = tree->GetBranch("c")->GetTitle();
	// Read the index maxima from the name of the branch
	size_t found = 0;
	for(int* maximum: {&l_max, &i_max, &k_max, &j_max})
	{
		found = limbranchtitle.find('[',found+1);
		limbranchtitle.find(']',found);
		*maximum = atoi(limbranchtitle.substr(found+1,1).c_str());
	}
	double**** c = newcoefficients();
		// Set up the 4D array and prepare to read from the tree
	char branchtitle[10]; // the letter "c" + four 2-digit numbers + 1 for luck
	// Seriously I don't expect *any* 2-digit numbers
	for ( int l = 0; l < l_max; l++ )
		for ( int i = 0; i < i_max; i++ )
			for ( int k = 0; k < k_max; k++ )
				for ( int j = 0; j < j_max; j++ )
				{
					sprintf(branchtitle,"c%d%d%d%d",l,i,k,j);
					tree->SetBranchAddress(branchtitle,&c[l][i][k][j]);
				}
	tree->GetEntry(0);
	storecoefficients(c);
	deletecoefficients(c);
	if(coeffs.size() == 0)
	{
		std::cerr << "No coefficients found. Defaulting to uniform shape." << std::endl;
		return;
	}
	printcoefficients();
	delete tree;
	delete file;
	init = false;
}
void LegendreMomentShape::Save(const std::string filename)
{
	TTree* outputTree = new TTree("LegendreMomentsTree","");
	char branchtitle[20];
	double**** c = newcoefficients();
	sprintf(branchtitle,"c[%d][%d][%d][%d]/D",l_max,i_max,k_max,j_max);
	outputTree->Branch("c",c,branchtitle);
	outputTree->Branch("mKK_min",&mKK_min,"mKK_min/D");
	outputTree->Branch("mKK_max",&mKK_max,"mKK_max/D");
	// Set up the branches
	for ( int l = 0; l < l_max; l++ )
		for ( int i = 0; i < i_max; i++ )
			for ( int k = 0; k < k_max; k++ )
				for ( int j = 0; j < j_max; j++ )
				{
					char branchname[5];
					sprintf(branchname,"c%d%d%d%d",l,i,k,j);
					outputTree->Branch(branchname,&c[l][i][k][j],((std::string)branchname+"/D").c_str());
				}
	// Pass the non-zero coefficients
	for(auto coeff: coeffs)
		c[coeff.l][coeff.i][coeff.k][coeff.j] = coeff.val;
	outputTree->Fill();
	// Save the tree to a file
	outputTree->SaveAs(filename.c_str());
	deletecoefficients(c);
}
void LegendreMomentShape::Generate(IDataSet* dataSet, const PhaseSpaceBoundary* boundary, const std::string mKKname, const std::string phiname, const std::string ctheta_1name, const std::string ctheta_2name)
{
	double**** c    = newcoefficients();
	double**** c_sq = newcoefficients(); // Used in caclulating the error
	mKK_min = boundary->GetConstraint(mKKname)->GetMinimum();
	mKK_max = boundary->GetConstraint(mKKname)->GetMaximum();
	const double mK  = 0.493677; // TODO: read these from config somehow
	const double mBs  = 5.36677;
	const double mPhi= 1.019461;
	int numEvents = dataSet->GetDataNumber();
	// Calculate the coefficients by summing over the dataset
	std::cout << "Sum over " << numEvents << " events" << std::endl;
	for (int e = 0; e < numEvents; e++)
	{
		// Retrieve the data point
		DataPoint * event = dataSet->GetDataPoint(e);
		double phi        = event->GetObservable(phiname)->GetValue();
		double ctheta_1   = event->GetObservable(ctheta_1name)->GetValue();
		double ctheta_2   = event->GetObservable(ctheta_2name)->GetValue();
		double mKK        = event->GetObservable(mKKname)->GetValue();
		double mKK_mapped = (mKK - mKK_min)/(mKK_max-mKK_min)*2.+ (-1);
		// Calculate phase space element
		double p1_st = DPHelpers::daughterMomentum(mKK, mK, mK);
		double p3    = DPHelpers::daughterMomentum(mBs,mKK,mPhi);
		double val = p1_st*p3;
		for ( int l = 0; l < l_max; l++ )
			for ( int i = 0; i < i_max; i++ )
				for ( int k = 0; k < k_max; k++ )
					for ( int j = 0; j < j_max; j++ )
					{
						double coeff = ((2*l + 1)/2.)*((2*i + 1)/2.)*Moment(l, i, k, j, mKK_mapped, phi, ctheta_1, ctheta_2)/val;
						c[l][i][k][j] += coeff;
						c_sq[l][i][k][j] += coeff * coeff;
					}
	}
	// Accept or reject the coefficients
	double threshold = 4; // TODO: read from config
	std::cout << "Keeping coefficients more significant than " << threshold << "σ" << std::endl;
	for ( int l = 0; l < l_max; l++ )
		for ( int i = 0; i < i_max; i++ )
			for ( int k = 0; k < k_max; k++ )
				for ( int j = 0; j < j_max; j++ )
				{
					if(std::abs(c[l][i][k][j]) < 1e-12)
					{
						c[l][i][k][j] = 0.;
						continue;
					}
					double error = sqrt(1./numEvents/numEvents * ( c_sq[l][i][k][j] - c[l][i][k][j]*c[l][i][k][j]/numEvents) );
					double signif = std::abs(c[l][i][k][j]/numEvents)/error;
					if ( signif > threshold )
					{
						printf("c[%d][%d][%d][%d] = %f;// ± %f with significance %fσ\n", l, i, k, j, c[l][i][k][j]/numEvents, error, signif );
						c[l][i][k][j] /= numEvents;
					}
					else
						c[l][i][k][j] = 0.;
				}
	storecoefficients(c);
	deletecoefficients(c);
	deletecoefficients(c_sq);
	init = false;
}
double LegendreMomentShape::Evaluate(const std::array<double,4>& datapoint) const
{
	return Evaluate(datapoint[0],datapoint[1],datapoint[2],datapoint[3]);
}
double LegendreMomentShape::Evaluate(const double mKK, const double phi, const double ctheta_1, const double ctheta_2) const
{
	if(init) return 1;
	double result = 0;
	double mKK_mapped = (mKK - mKK_min) / (mKK_max - mKK_min)*2 - 1;
	if(std::abs(mKK_mapped) > 1) return 0; // I could print a warning here, but it gets tedious when you just want a mass projection with sensibly-sized bins that includes the threshold
	for(auto coeff : coeffs)
		result += coeff.val*Moment(coeff.l, coeff.i, coeff.k, coeff.j, mKK_mapped, phi, ctheta_1, ctheta_2);
	return result;
}
double**** LegendreMomentShape::newcoefficients() const
{
	double**** c = new double***[l_max];
	for ( int l = 0; l < l_max; l++ )
	{
		c[l] = new double**[i_max];
		for ( int i = 0; i < i_max; i++ )
		{
			c[l][i] = new double*[k_max];
			for ( int k = 0; k < k_max; k++ )
			{
				c[l][i][k] = new double[j_max];
				for ( int j = 0; j < j_max; j++ )
					c[l][i][k][j] = 0;
			}
		}
	}
	return c;
}
void LegendreMomentShape::deletecoefficients(double**** c) const
{
	for ( int l = 0; l < l_max; l++ )
	{
		for ( int i = 0; i < i_max; i++ )
		{
			for ( int k = 0; k < k_max; k++ )
				delete c[l][i][k];
			delete c[l][i];
		}
		delete c[l];
	}
	delete c;
}
void LegendreMomentShape::storecoefficients(double**** c)
{
	std::cout << "Storing coefficients" << std::endl;
	// Create std::vector of non-zero coefficients after reading from tree
	for ( int l = 0; l < l_max; l++ )
		for ( int i = 0; i < i_max; i++ )
			for ( int k = 0; k < k_max; k++ )
				for ( int j = 0; j < j_max; j++ )
				{
					if(std::abs(c[l][i][k][j]) < 1e-12)
						continue;
					coefficient coeff;
					coeff.l = l;
					coeff.i = i;
					coeff.k = k;
					coeff.j = j;
					coeff.val = c[l][i][k][j];
					coeffs.push_back(coeff);
				}
}
void LegendreMomentShape::printcoefficients() const
{
	for(auto coeff : coeffs)
		coeff.print();
}
