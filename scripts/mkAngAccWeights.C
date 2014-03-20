
{
	vector<double> weights;
	weights.push_back(  1.0000 );
	weights.push_back(  1.0514 );
	weights.push_back(  1.0532 );
	weights.push_back( -0.0006 );
	weights.push_back(  0.0001 );
	weights.push_back( -0.0003 );
	weights.push_back(  1.0155 );
	weights.push_back(  0.0007 );
	weights.push_back( -0.0004 );
	weights.push_back(  0.0008 );
	/*weights.push_back(  1.0      );
	weights.push_back(  1.0472   );
	weights.push_back(  1.0489   );
	weights.push_back( -0.0008   );
	weights.push_back( -0.0002   );
	weights.push_back( -0.0003   );
	weights.push_back(  1.0130   );
	weights.push_back(  0.0007   );
	weights.push_back( -0.0002   );
	weights.push_back( -0.0020   );*/

	TFile* theseWeights = new TFile( "newAngAccWeights.root", "RECREATE" );

	TTree* Thistree = new TTree( "tree", "tree containing acceptance weights and histo" );

	Thistree->Branch("weights", "std::vector<double>", &weights);

	Thistree->Fill();
	Thistree->Write("",TObject::kOverwrite);
	theseWeights->Write("",TObject::kOverwrite);
}
