#include "FitFractionCalculator.h"
#include <memory>
#include "TTree.h"

FitFractionCalculator::FitFractionCalculator(IPDF& pdf, PhaseSpaceBoundary& boundary) : storedBoundary(boundary)
{
	RapidFitIntegrator testIntegrator( &pdf, true, true );
	vector<std::string> doNotIntegrate = pdf.GetDoNotIntegrateList();
	for(const auto& combination: storedBoundary.GetDiscreteCombinations())
	{
		int i = GetIndex(combination);
		double TotalIntegral = testIntegrator.Integral(combination, &storedBoundary);
		for(const auto& ComponentName: pdf.PDFComponents())
		{
			ComponentRef thisRef( ComponentName, "dummyObservable" );
			ComponentFitFractions[i][ComponentName] = testIntegrator.NumericallyIntegrateDataPoint(combination, &storedBoundary, doNotIntegrate, &thisRef);
			ComponentFitFractions[i][ComponentName] /= TotalIntegral;
		}
	}
}
std::map<std::string, double> FitFractionCalculator::GetFitFractions(DataPoint* combination) const
{
	return ComponentFitFractions.at(GetIndex(combination));
}
double FitFractionCalculator::GetComponentFitFraction(DataPoint* combination, std::string ComponentName) const
{
	return ComponentFitFractions.at(GetIndex(combination)).at(ComponentName);
}
void FitFractionCalculator::Print() const
{
	for(const auto& combination: ComponentFitFractions)
	{
		std::cout << "Fit fractions in combination " << combination.first << "\n";
		for(const auto& fraction: combination.second)
			std::cout << fraction.first << "\t" << fraction.second << "\n";
		std::cout << std::endl;
	}
}
void FitFractionCalculator::WriteToFile(std::string filename)
{
	auto OutputFile = std::unique_ptr<TFile>(TFile::Open(filename.c_str(),"UPDATE"));
	for(auto& combination: ComponentFitFractions)
	{
		std::string treename = "combination_"+std::to_string(combination.first);
		auto OutputTree = std::unique_ptr<TTree>(dynamic_cast<TTree*>(OutputFile->Get(treename.c_str())));
		bool createtree = OutputTree == nullptr;
		// If the tree already exists, update it with a new entry
		// If it doesn't exist, create it
		if(createtree) OutputTree = std::unique_ptr<TTree>(new TTree(treename.c_str(),"Fit fractions"));
		for(auto& fraction: combination.second)
		{
			std::string branchname = fraction.first;
			// Remove all punctuation because ROOT interprets the branch name as a TFormula when drawings
			branchname.erase(std::remove_if(branchname.begin(), branchname.end(), &ispunct), branchname.end());
			// Avoid purely-numerical branch names for the same reason as above
			branchname = "component_" + branchname;
			if(createtree)
				OutputTree->Branch(branchname.c_str(),&fraction.second);
			else
				OutputTree->SetBranchAddress(branchname.c_str(),&fraction.second);
		}
		OutputTree->Fill();
		OutputTree->Write();
	}
	OutputFile->Close();
}
int FitFractionCalculator::GetIndex(DataPoint* combination) const
{
	return storedBoundary.GetDiscreteIndex(combination);
}

