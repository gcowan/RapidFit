// Self
#include "Bs2PhiKKSignal.h"
// Std Libraries
#include <iostream>
#include <stdexcept>
#include <complex>
// ROOT Libraries
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
// RapidFit
#include "StringProcessing.h"
#include "PDFConfigurator.h"
#include "DPHelpers.hh"
// GSL
#include <gsl/gsl_randist.h>

PDF_CREATOR( Bs2PhiKKSignal )
/*****************************************************************************/
// Constructor
Bs2PhiKKSignal::Bs2PhiKKSignal(PDFConfigurator* config) :
	// Dependent variable names
	 mKKName(config->getName("mKK"))
	,phiName(config->getName("phi"))
	,ctheta_1Name(config->getName("ctheta_1"))
	,ctheta_2Name(config->getName("ctheta_2"))
	,acceptance_moments((std::string)config->getConfigurationValue("CoefficientsFile") != "")
	,convolve(config->isTrue("convolve"))
	,outofrange(false)
{
	std::cout << "\nBuilding Bs → ϕ K+ K− signal PDF\n\n";
	std::string phiname = config->getConfigurationValue("phiname");
	thraccscale = Bs2PhiKKComponent::PhysPar(config,"thraccscale");
	phimass = Bs2PhiKKComponent::PhysPar(config,phiname+"_mass");
	dGsGs = Bs2PhiKKComponent::PhysPar(config,"dGsGs");
	if(convolve)
		for(std::string name: {"sigma1","sigma2","frac","nsteps","nsigma"})
			mKKrespars[name] = std::stod(config->getConfigurationValue("mKKres_"+name));
	std::vector<std::string> reslist = StringProcessing::SplitString(config->getConfigurationValue("resonances"), ' ');
	std::cout << "┏━━━━━━━━━━━━━━━┯━━━━━━━┯━━━━━━━━━━━━━━━┓\n";
	std::cout << "┃ Component\t│ Spin\t│ Lineshape\t┃\n";
	std::cout << "┠───────────────┼───────┼───────────────┨\n";
	for(const auto& name: reslist)
	{
		if(name=="") continue;
		components[name] = ParseComponent(config,phiname,name);
		componentnames.push_back(name);
	}
	if(components.size() > 1) componentnames.push_back("interference");
	std::cout << "┗━━━━━━━━━━━━━━━┷━━━━━━━┷━━━━━━━━━━━━━━━┛" << std::endl;
	if(acceptance_moments) acc_m = std::unique_ptr<LegendreMomentShape>(new LegendreMomentShape(config->getConfigurationValue("CoefficientsFile")));
	Initialise();
	MakePrototypes();
}
/*****************************************************************************/
// Copy constructor
Bs2PhiKKSignal::Bs2PhiKKSignal(const Bs2PhiKKSignal& copy) : BasePDF( (BasePDF) copy)
	// Dependent variable names
	,mKKName(copy.mKKName)
	,phiName(copy.phiName)
	,ctheta_1Name(copy.ctheta_1Name)
	,ctheta_2Name(copy.ctheta_2Name)
	// Width splitting
	,dGsGs(copy.dGsGs)
	// Phi mass
	,phimass(copy.phimass)
	// threshold acceptance scale
	,thraccscale(copy.thraccscale)
	// mass resolution parameters
	,mKKrespars(copy.mKKrespars)
	// PDF components
	,components(copy.components)
	,componentnames(copy.componentnames)
	// Plotting components
	// Options
	,acceptance_moments(copy.acceptance_moments)
	,convolve(copy.convolve)
	// Status
	,outofrange(copy.outofrange)
{
	if(acceptance_moments) acc_m = std::unique_ptr<LegendreMomentShape>(new LegendreMomentShape(*copy.acc_m));
	Initialise();
}
/*****************************************************************************/
// Destructor
Bs2PhiKKSignal::~Bs2PhiKKSignal()
{
}
/*****************************************************************************/
// Code common to the constructors
void Bs2PhiKKSignal::Initialise()
{
	// Enable numerical normalisation and disable caching
	this->SetNumericalNormalisation( true );
	this->TurnCachingOff();
}
/*****************************************************************************/
// Build a component object from a passed option
Bs2PhiKKComponent Bs2PhiKKSignal::ParseComponent(PDFConfigurator* config, std::string phiname, std::string option) const
{
	// Syntax: <resonance name>(<spin>,<lineshape>)
	// - the list is space-delimited: no extra spaces, please!
	// - resonance name must be alphanumeric
	// - spin must be a single digit 0, 1 or 2, or you have to implement the code for higher spins
	// - lineshape name must be 2 capital letters
	// - see Bs2PhiKKComponent.cpp for implemented lineshapes (BW, FT, NR... but more can be added)
	// Example: "phi1020(1,BW)"
	std::string KKname = option.substr(0,option.find('('));
	int JKK = atoi(option.substr(option.find('(')+1,1).c_str());
	std::string lineshape = option.substr(option.find(',')+1,2);
	std::cout << "┃ " << KKname << "\t│    " << JKK << "\t│ " << lineshape << " shape \t┃\n";
	return Bs2PhiKKComponent(config, phiname, KKname, JKK, lineshape);
}
/*****************************************************************************/
// Make the data point and parameter set
void Bs2PhiKKSignal::MakePrototypes()
{
	// Make the DataPoint prototype
	// The ordering here matters. It has to be the same as the XML file, apparently.
	allObservables.push_back(mKKName);
	allObservables.push_back(phiName);
	allObservables.push_back(ctheta_1Name);
	allObservables.push_back(ctheta_2Name);
	// Make the parameter set
	std::vector<std::string> parameterNames;
	// Resonance parameters
	parameterNames.push_back(dGsGs.name);
	parameterNames.push_back(phimass.name);
	parameterNames.push_back(thraccscale.name);
	for(const auto& comp: components)
		for(std::string par: comp.second.GetPhysicsParameters())
			parameterNames.push_back(par);
	std::sort(parameterNames.begin(),parameterNames.end());
	parameterNames.erase(std::unique(parameterNames.begin(),parameterNames.end()),parameterNames.end());
	allParameters = ParameterSet(parameterNames);
}
/*****************************************************************************/
// Set the physics parameters
bool Bs2PhiKKSignal::SetPhysicsParameters(ParameterSet* NewParameterSet)
{
	outofrange = false;
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	dGsGs.Update(&allParameters);
	phimass.Update(&allParameters);
	thraccscale.Update(&allParameters);
	for(auto& comp: components)
	{
		try
		{
			comp.second.SetPhysicsParameters(&allParameters);
		}
		catch(...)
		{
			outofrange = true;
		}
	}
	return isOK;
}
/*****************************************************************************/
// List of components
std::vector<std::string> Bs2PhiKKSignal::PDFComponents()
{
	// Avoid redundant plotting for single-component PDFs
	if(components.size() == 1) return {};
	return componentnames;
}
/*****************************************************************************/
// Evaluate a single component
double Bs2PhiKKSignal::EvaluateComponent(DataPoint* measurement, ComponentRef* component)
{
	std::string compName = component->getComponentName();
	if(compName=="0") // should quicken things up slightly?
		return Evaluate(measurement);
	const Bs2PhiKKComponent::datapoint_t datapoint = ReadDataPoint(measurement);
	// Evaluate the PDF at this point
	double MatrixElementSquared;
	if(compName=="interference")
		MatrixElementSquared = convolve? Convolve(&Bs2PhiKKSignal::InterferenceMsq,datapoint,"") : InterferenceMsq(datapoint);
	else
		MatrixElementSquared = convolve? Convolve(&Bs2PhiKKSignal::ComponentMsq,datapoint,compName) : ComponentMsq(datapoint,compName);
	return Evaluate_Base(MatrixElementSquared, datapoint);
}
// Evaluate the entire PDF
double Bs2PhiKKSignal::Evaluate(DataPoint* measurement)
{
	if(outofrange)
		return 1e-100;
	const Bs2PhiKKComponent::datapoint_t datapoint = ReadDataPoint(measurement);
	double MatrixElementSquared = convolve? Convolve(&Bs2PhiKKSignal::TotalMsq,datapoint,"") : TotalMsq(datapoint);
	return Evaluate_Base(MatrixElementSquared, datapoint);
}
// The stuff common to both Evaluate() and EvaluateComponent()
double Bs2PhiKKSignal::Evaluate_Base(const double MatrixElementSquared, const Bs2PhiKKComponent::datapoint_t& datapoint) const
{
	return MatrixElementSquared * p1stp3(datapoint[0]) * Acceptance(datapoint);
}
/*****************************************************************************/
Bs2PhiKKComponent::datapoint_t Bs2PhiKKSignal::ReadDataPoint(DataPoint* measurement) const
{
	// Get values from the datapoint
	double mKK      = measurement->GetObservable(mKKName     )->GetValue();
	double phi      = measurement->GetObservable(phiName     )->GetValue();
	double ctheta_1 = measurement->GetObservable(ctheta_1Name)->GetValue();
	double ctheta_2 = measurement->GetObservable(ctheta_2Name)->GetValue();
	phi+=M_PI;
	return {mKK, phi, ctheta_1, ctheta_2};
}
/*Calculate matrix elements***************************************************/
// Total |M|²: coherent sum of all amplitudes
double Bs2PhiKKSignal::TotalMsq(const Bs2PhiKKComponent::datapoint_t& datapoint, const std::string& dummy) const
{
	(void)dummy;
	Bs2PhiKKComponent::amplitude_t TotalAmp = {std::complex<double>(0, 0), std::complex<double>(0, 0)};
	for(const auto& comp : components)
	{
		Bs2PhiKKComponent::amplitude_t CompAmp = comp.second.Amplitude(datapoint);
		if(std::isnan(CompAmp[0].real()) || std::isnan(CompAmp[0].imag())){ std::cerr << comp.first << " amplitude evaluates to " << CompAmp[0] << std::endl; std::exit(1);}
		TotalAmp[false] += CompAmp[false];
		TotalAmp[true] += CompAmp[true];
	}
	return TimeIntegratedMsq(TotalAmp);
}
// Single-component |M|²
double Bs2PhiKKSignal::ComponentMsq(const Bs2PhiKKComponent::datapoint_t& datapoint, const std::string& compName) const
{
	return TimeIntegratedMsq(components.at(compName).Amplitude(datapoint, compName));
}
// Interference |M|²: difference between incoherent and coherent sums of components
double Bs2PhiKKSignal::InterferenceMsq(const Bs2PhiKKComponent::datapoint_t& datapoint, const std::string& dummy) const
{
	(void)dummy;
	double MatrixElementSquared = TotalMsq(datapoint);
		for(const auto& comp : components)
			MatrixElementSquared -= TimeIntegratedMsq(comp.second.Amplitude(datapoint));
	return MatrixElementSquared;
}
// Take the amplitudes of the B and Bbar decay and return the time-integrated |M|²
double Bs2PhiKKSignal::TimeIntegratedMsq(const Bs2PhiKKComponent::amplitude_t& Amp) const
{
	double GH = (2 - dGsGs.value); // Actually Γ/2ΓH but who cares about an overall factor Γ?
	double GL = (2 + dGsGs.value); // Actually Γ/2ΓL
	double termone = (std::norm(Amp[false]) + std::norm(Amp[true])) * (1./GL + 1./GH);
	double termtwo = 2*std::real(Amp[true] * std::conj(Amp[false])) * (1./GL - 1./GH);
	return termone + termtwo;
}
// Convolution of the matrix element function with a Gaussian... the slow integral way
double Bs2PhiKKSignal::Convolve(MsqFunc_t EvaluateMsq, const Bs2PhiKKComponent::datapoint_t& datapoint, const std::string& compName) const
{
	const double res1 = mKKrespars.at("sigma1");
	const double res2 = mKKrespars.at("sigma2");
	const double frac = mKKrespars.at("frac");
	const double nsigma = mKKrespars.at("nsigma");
	const int nsteps = mKKrespars.at("nsteps");
	const double resolution = frac*res1+(1-frac)*res2;
	// Can't do this if we're too close to threshold: it starts returning nan
	if(datapoint[0] - nsigma*resolution > 2*Bs2PhiKKComponent::mK)
	{
		double Msq_conv = 0.;
		const double stepsize = 2.*nsigma*resolution/nsteps;
		// Integrate over range −nσ to +nσ
		for(double x = -nsigma*resolution; x < nsigma*resolution; x += stepsize)
			Msq_conv += (frac*gsl_ran_gaussian_pdf(x,res1)+(1-frac)*gsl_ran_gaussian_pdf(x,res2)) * (this->*EvaluateMsq)({datapoint[0]-x,datapoint[1],datapoint[2],datapoint[3]},compName) * stepsize; // FML
		return Msq_conv;
	}
	else
		return (this->*EvaluateMsq)(datapoint,compName);
}
/*Stuff that factors out of the time integral*********************************/
double Bs2PhiKKSignal::Acceptance(const Bs2PhiKKComponent::datapoint_t& datapoint) const
{
	double acceptance;
	if(acceptance_moments)
	{
		acceptance = acc_m->Evaluate(datapoint);
		acceptance *= std::erf(thraccscale.value*(datapoint[0]-2*Bs2PhiKKComponent::mK));
//		acceptance *= std::tanh(thraccscale.value*(datapoint[0]-2*Bs2PhiKKComponent::mK));
//		acceptance *= std::atan(thraccscale.value*(datapoint[0]-2*Bs2PhiKKComponent::mK))*2.0/M_PI;
	}
	if(std::isnan(acceptance)) std::cerr << "Acceptance evaluates to nan" << std::endl;
	return acceptance;
}
double Bs2PhiKKSignal::p1stp3(const double& mKK) const
{
	const double mK   = Bs2PhiKKComponent::mK;
	const double mBs  = Bs2PhiKKComponent::mBs;
	const double mPhi = phimass.value;
	if(mKK < 2*mK || mKK > mBs-mPhi) return 0;
	double pR = DPHelpers::daughterMomentum(mKK, mK,  mK);
	double pB = DPHelpers::daughterMomentum(mBs, mKK, mPhi);
	double pRpB = pR * pB;
	if(std::isnan(pRpB)) std::cerr << "p1stp3 evaluates to nan" << std::endl;
	return pRpB;
}
/*****************************************************************************/
double Bs2PhiKKSignal::Normalisation(PhaseSpaceBoundary* boundary)
{
	(void)boundary;
	return -1;
}

