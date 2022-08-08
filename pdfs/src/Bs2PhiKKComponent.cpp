// Self
#include "Bs2PhiKKComponent.h"
// Std Libraries
#include <iostream>
#include <stdexcept>
// RapidFit Dalitz Plot Libraries
#include "DPBWResonanceShape.hh"
#include "DPFlatteShape.hh"
#include "DPNonresonant.hh"
#include "DPHelpers.hh"
#include "DPWignerFunctionJ0.hh"
#include "DPWignerFunctionJ1.hh"
#include "DPWignerFunctionJ2.hh"

double Bs2PhiKKComponent::mBs  = 5.36677;
double Bs2PhiKKComponent::mK   = 0.493677;
double Bs2PhiKKComponent::mpi  = 0.139570;

// Constructor
Bs2PhiKKComponent::Bs2PhiKKComponent(PDFConfigurator* config, std::string _phiname, std::string KKname, int _JKK, std::string _lineshape) :
	  phimass(PhysPar(config,_phiname+"_mass"))
	, fraction(PhysPar(config,KKname+"_fraction"))
	, JKK(_JKK)
	, lineshape(_lineshape)
{
	// Barrier factors
	if(lineshape != "NR")
	{
		BsBFradius = PhysPar(config,"BsBFradius");
		KKBFradius = PhysPar(config,"KKBFradius");
		Bsbarrier = DPBarrierFactor(  0,1.0,0);
		KKbarrier = DPBarrierFactor((unsigned int)JKK,3.0,0);
	}
	// KK resonance parameters
	// Breit Wigner
	if(lineshape=="BW")
	{
		KKpars.push_back(PhysPar(config,KKname+"_mass"));
		KKpars.push_back(PhysPar(config,KKname+"_width"));
	}
	// Flatte
	else if(lineshape=="FT")
	{
		KKpars.push_back(PhysPar(config,KKname+"_mass")); 
		KKpars.push_back(PhysPar(config,KKname+"_gpipi"));
		KKpars.push_back(PhysPar(config,KKname+"_Rg"));
	}
	else if(lineshape!="NR")
	{
		lineshape = "NR";
		std::cerr << "Bs2PhiKKComponent WARNING: unknown lineshape '" << lineshape << "'. Treating this component non-resonant." << std::endl;
	}
	// Helicity amplitude observables
	int lambda_max = std::min(1, _JKK); // Maximum helicity
	std::vector<int> helicities;
	if(lineshape!="NR")
		for(int lambda = -lambda_max; lambda <= lambda_max; lambda++)
			helicities.push_back(lambda);
	long unsigned int n = helicities.size();
	switch(n)
	{
		case 0:
			break;
		case 1:
			phases.push_back(PhysPar(config,KKname+"_deltazero"));
			break;
		case 3:
			magsqs.push_back(PhysPar(config,KKname+"_Aperpsq"));
			magsqs.push_back(PhysPar(config,KKname+"_Azerosq"));
			phases.push_back(PhysPar(config,KKname+"_deltaperp"));
			phases.push_back(PhysPar(config,KKname+"_deltazero"));
			phases.push_back(PhysPar(config,KKname+"_deltapara"));
			break;
		default:
			std::cerr << "Bs2PhiKKComponent can't handle this many helicities: " << n << std::endl;
			std::exit(-1);
			break;
	}
	// Make the helicity amplitude vector
	for(const auto& lambda: helicities)
		Ahel[lambda] = std::polar<double>(sqrt(1. / (double)n), 0);
	Initialise();
}
// Copy constructor
Bs2PhiKKComponent::Bs2PhiKKComponent(const Bs2PhiKKComponent& other) :
	// Floatable parameters
	  fraction(other.fraction)
	, Ahel(other.Ahel)
	, magsqs(other.magsqs)
	, phases(other.phases)
	, phimass(other.phimass)
	, KKpars(other.KKpars)
	, BsBFradius(other.BsBFradius)
	, KKBFradius(other.KKBFradius)
	// Fixed parameters
	, JKK(other.JKK)
	// DP objects
	, Bsbarrier(other.Bsbarrier)
	, KKbarrier(other.KKbarrier)
	// Options
	, lineshape(other.lineshape)
{
	Initialise();
}
// Copy by assignment
Bs2PhiKKComponent& Bs2PhiKKComponent::operator=(const Bs2PhiKKComponent& other)
{
	// Floatable parameters
	fraction = other.fraction;
	Ahel = other.Ahel;
	magsqs = other.magsqs;
	phases = other.phases;
	phimass = other.phimass;
	KKpars = other.KKpars;
	BsBFradius = other.BsBFradius;
	KKBFradius = other.KKBFradius;
	// Fixed parameters
	JKK = other.JKK;
	// DP objects
	Bsbarrier = other.Bsbarrier;
	KKbarrier = other.KKbarrier;
	// Options
	lineshape = other.lineshape;
	Initialise();
	return *this;
}
Bs2PhiKKComponent::~Bs2PhiKKComponent()
{
}
void Bs2PhiKKComponent::Initialise()
{
	// Breit Wigner
	if(lineshape=="BW")
		KKLineShape = std::unique_ptr<DPBWResonanceShape>(new DPBWResonanceShape(KKpars[0].value, KKpars[2].value, JKK, mK, mK, KKBFradius.value));
	// Flatte
	else if(lineshape=="FT")
		KKLineShape = std::unique_ptr<DPFlatteShape>(new DPFlatteShape(KKpars[0].value, KKpars[1].value, mpi, mpi, KKpars[1].value*KKpars[2].value, mK, mK));
	else
		KKLineShape = std::unique_ptr<DPNonresonant>(new DPNonresonant());
	// Build the barrier factor and Wigner function objects
	wignerPhi = std::unique_ptr<DPWignerFunctionJ1>(new DPWignerFunctionJ1());
	switch (JKK) // I hate this but I'd rather it just worked...
	{
		case 0:
			wignerKK  = std::unique_ptr<DPWignerFunctionJ0>(new DPWignerFunctionJ0());
			break;
		case 1:
			wignerKK  = std::unique_ptr<DPWignerFunctionJ1>(new DPWignerFunctionJ1());
			break;
		case 2:
			wignerKK  = std::unique_ptr<DPWignerFunctionJ2>(new DPWignerFunctionJ2());
			break;
		default:
			wignerKK  = std::unique_ptr<DPWignerFunctionGeneral>(new DPWignerFunctionGeneral(JKK)); // This should only happen for the rho_3 (1690)
			break;
	}
	UpdateAmplitudes();
	UpdateLineshape();
}
// Angular part of the amplitude
std::complex<double> Bs2PhiKKComponent::F(const int lambda, const double Phi, const double ctheta_1, const double ctheta_2) const
{
	const double d_phi = wignerPhi->function(ctheta_1, lambda, 0);
	const double d_KK = wignerKK->function(ctheta_2, lambda, 0);
	if(std::isnan(d_phi)) std::cerr << "\tWigner function for the phi evaluates to nan" << std::endl;
	if(std::isnan(d_phi)) std::cerr << "\tWigner function for the KK resonance evaluates to nan" << std::endl;
	return d_phi * d_KK * std::polar<double>(1, lambda*Phi);
}
std::complex<double> Bs2PhiKKComponent::AngularPart(const double phi, const double ctheta_1, const double ctheta_2) const
{
	std::complex<double> angularPart(0, 0);
	if(Ahel.empty()) return std::complex<double>(1, 0); // Must be non-resonant
	for(const auto& A : Ahel)
	{
		int lambda = A.first;
		if(std::isnan(A.second.real()) || std::isnan(A.second.imag())) std::cerr << "\tA(" << lambda << ") is " << A.second << std::endl;
		angularPart += A.second * F(lambda, phi, ctheta_1, ctheta_2);
	}
	return angularPart;
}
// Orbital and barrier factor
double Bs2PhiKKComponent::OFBF(const double mKK) const
{
	if(mKK < 2*mK) return 0;
	if(lineshape=="NR")
		return 1;
	// Orbital factor
	// Masses
	double mphi = phimass.value;
	double mKK0 = KKpars[0].value;
	// Momenta
	double pBs  = DPHelpers::daughterMomentum(mBs,  mphi, mKK   );
	double pKK  = DPHelpers::daughterMomentum(mKK,  mK,   mK    );
	double orbitalFactor = //std::pow(pBs/mBs,   0)* // == 1 so don't bother
	                       std::pow(pKK/mKK0, JKK);
	// Barrier factors
	double barrierFactor = Bsbarrier.barrier(pBs)*
	                       KKbarrier.barrier(pKK);
	if(std::isnan(orbitalFactor)) std::cerr << "\tOrbital factor evaluates to nan" << std::endl;
	if(std::isnan(barrierFactor)) std::cerr << "\tBarrier factor evaluates to nan" << std::endl;
	return orbitalFactor*barrierFactor;
}
// The full amplitude.
Bs2PhiKKComponent::amplitude_t Bs2PhiKKComponent::Amplitude(const datapoint_t& datapoint) const
{
	double mKK = datapoint[0];
	double phi = datapoint[1];
	double ctheta_1 = datapoint[2];
	double ctheta_2 = datapoint[3];
	amplitude_t angularPart = {AngularPart(phi, ctheta_1, ctheta_2), AngularPart(-phi, -ctheta_1, -ctheta_2)};
	std::complex<double> massPart = KKLineShape->massShape(mKK);
	if(std::isnan(fraction.value)) std::cerr << "\tFraction is nan" << std::endl;
	if(std::isnan(massPart.real()) || std::isnan(massPart.imag())) std::cerr << "\tLineshape evaluates to " << massPart << std::endl;
	massPart *= fraction.value * OFBF(mKK);
	return {massPart*angularPart[false], massPart*angularPart[true]};
}
// The full amplitude with an option.
Bs2PhiKKComponent::amplitude_t Bs2PhiKKComponent::Amplitude(const datapoint_t& datapoint, const std::string option) const
{
	if(Ahel.empty() || option == "" || option.find("odd") == std::string::npos || option.find("even") == std::string::npos ) return Amplitude(datapoint);
	double mKK = datapoint[0];
	double phi = datapoint[1];
	double ctheta_1 = datapoint[2];
	double ctheta_2 = datapoint[3];
	// Angular part
	amplitude_t angularPart = {std::complex<double>(0, 0), std::complex<double>(0, 0)};
	std::complex<double> Aperp = std::polar(sqrt(magsqs[0].value),phases[0].value);
	std::complex<double> Apara = std::polar(sqrt(1. - magsqs[0].value - magsqs[1].value),phases[2].value);
	// Temporary helicity amplitudes
	std::vector<std::complex<double>> HelAmp;
	// CP-odd component
	if(option.find("odd") != std::string::npos)
		HelAmp = {-Aperp/sqrt(2.), std::complex<double>(0, 0), Aperp/sqrt(2.)};
	// CP-even component
	else
		HelAmp = {Apara/sqrt(2.), Ahel.at(0), Apara/sqrt(2.)};
	for(const auto& A : Ahel)
	{
		int lambda = A.first;
		angularPart[false] += HelAmp[(unsigned int)lambda+1] * F(lambda, phi, ctheta_1, ctheta_2);
		angularPart[true] += HelAmp[(unsigned int)lambda+1] * F(lambda, -phi, -ctheta_1, -ctheta_2);
	}
	std::complex<double> massPart = KKLineShape->massShape(mKK);
	if(std::isnan(fraction.value)) std::cerr << "\tFraction is nan" << std::endl;
	if(std::isnan(massPart.real()) || std::isnan(massPart.imag())) std::cerr << "\tLineshape evaluates to " << massPart << std::endl;
	massPart *= fraction.value * OFBF(mKK);
	return {massPart*angularPart[false], massPart*angularPart[true]};
}
// Update everything from the parameter set
void Bs2PhiKKComponent::SetPhysicsParameters(ParameterSet* fitpars)
{
	// Update the parameters objects first
	fraction.Update(fitpars);
	phimass.Update(fitpars);
	if(lineshape != "NR")
	{
		KKBFradius.Update(fitpars);
		BsBFradius.Update(fitpars);
	}
	for(auto& par: magsqs) par.Update(fitpars);
	for(auto& par: phases) par.Update(fitpars);
	for(auto& par: KKpars) par.Update(fitpars);
	UpdateBarriers();
	UpdateAmplitudes();
	UpdateLineshape();
}
void Bs2PhiKKComponent::UpdateAmplitudes()
{
	// Update the helicity amplitudes
	if(Ahel.empty())
		return;
	else if(Ahel.size() == 1)
		Ahel[0] = std::polar<double>(1.0,phases[0].value);
	else if(Ahel.size() == 3)
	{
		std::complex<double> Aperp = std::polar(sqrt(magsqs[0].value),phases[0].value);
		double Apara_mag = sqrt(1. - magsqs[0].value - magsqs[1].value);
		if(std::isnan(Apara_mag))
		{
			Apara_mag = 0;
			throw std::out_of_range("sum of |A|² > 1");
		}
		std::complex<double> Apara = std::polar(Apara_mag,phases[2].value);
		Ahel[-1] = (Apara - Aperp)/sqrt(2.); // A− = (A‖ − A⊥)/sqrt(2)
		Ahel[ 0] = std::polar(sqrt(magsqs[1].value),phases[1].value);;
		Ahel[+1] = (Apara + Aperp)/sqrt(2.); // A+ = (A‖ + A⊥)/sqrt(2)
	}
	else
	{
		std::cerr << "Bs2PhiKKComponent can't handle this many helicities: " << Ahel.size() << std::endl;
		std::exit(-1);
	}
}
void Bs2PhiKKComponent::UpdateLineshape()
{
	// Update the resonance line shape
	vector<double> respars;
	if(lineshape == "BW")
	{
		respars.push_back(KKpars[0].value); // mass
		respars.push_back(KKpars[1].value); // width
		respars.push_back(KKBFradius.value); // barrier factor radius
	}
	else if(lineshape == "FT")
	{
		respars.push_back(KKpars[0].value); // mass
		respars.push_back(KKpars[1].value); // gpipi
		respars.push_back(KKpars[1].value*KKpars[2].value); // gKK = gpipi*Rg
	}
	KKLineShape->setParameters(respars);
}
void Bs2PhiKKComponent::UpdateBarriers()
{
	if(lineshape == "NR") return;
	double mphi = phimass.value;
	double Mres = KKpars[0].value;
	double m_min  = mK + mK;
	double m_max  = mBs - mphi;
	double m0_eff = m_min + (m_max - m_min) * (1 + std::tanh((Mres - (m_min + m_max) / 2) / (m_max - m_min))) / 2;
	double pBs0 = DPHelpers::daughterMomentum(mBs,  mphi, m0_eff);
	double pKK0 = DPHelpers::daughterMomentum(Mres, mK,   mK    );
	Bsbarrier.setparameters(BsBFradius.value,pBs0);
	KKbarrier.setparameters(KKBFradius.value,pKK0);
}
vector<ObservableRef> Bs2PhiKKComponent::GetPhysicsParameters() const
{
	vector<ObservableRef> parameters;
	for(const auto& set: {{fraction,phimass},magsqs,phases,KKpars})
		for(const auto& par: set)
			parameters.push_back(par.name);
	// Add barrier factors if this is a resonant component
	if(lineshape != "NR")
		for(const auto& par: {BsBFradius,KKBFradius})
			parameters.push_back(par.name);
	return parameters;
}

