void RapidFit_TOYS(char* XMLStr="XML_UB_B.xml", char* repeat_num="", char* pull_plots="PullPlots.root", char* SEED="12345" )
{
	gSystem->Load("libFoam");
	gSystem->Load("libRooFit");
 	gSystem->Load("libMinuit");
	gSystem->Load("libMinuit2");
	gSystem->Load("libRapidRun");
	//set up main, eg command line opts
	TList* args = new TList();
	args->Add(new TObjString("RapidFit"));//program name
	args->Add(new TObjString( "-f" ));
	args->Add(new TObjString( XMLStr ));
	args->Add(new TObjString( "-repeats" ));
	args->Add(new TObjString( repeat_num ));
	args->Add(new TObjString( "--doPulls" ));
	args->Add(new TObjString( pull_plots ));
	args->Add(new TObjString( "--SetSeed" ));
	args->Add(new TObjString( SEED ));

	cout << "Constructing Fit" <<endl;
	cout << XMLStr << "\t" << repeat_num << "\t" << pull_plots << "\t" << SEED << endl;
	gDirectory->ls();

	RapidRun* fitting = new RapidRun( args );
	int result = fitting->run();
	gApplication->Terminate();
}
