void RapidFit_FC( char* XML_Str="XML_UB_B.xml", char* param_string_1="Phi_s,0,3.14159,50", char* param_string_2="deltaGamma,-0.7,0.7,50", char* repeat_events="10", char* SEED="12345" )
{
	gSystem->Load("libFoam");
	gSystem->Load("libRooFit");
 	gSystem->Load("libMinuit");
	gSystem->Load("libMinuit2");
	gSystem->Load("libRapidRun");
	//set up main, eg command line opts
	TList* args = new TList();
	args->Add(new TObjString("RapidFit"));//program name
	args->Add(new TObjString("-f"));
	args->Add(new TObjString( XML_Str ));
	args->Add(new TObjString("--defineContour"));
	args->Add(new TObjString( param_string_1 ));
	args->Add(new TObjString( param_string_2 ));
	args->Add(new TObjString("--doFCscan"));
	args->Add(new TObjString("-repeats"));
	args->Add(new TObjString( repeat_events ));
	args->Add(new TObjString("--useUUID"));
	args->Add(new TObjString("--SetSeed"));
	args->Add(new TObjString( SEED ));

	cout << "Constructing Fit" <<endl;
	cout << param_string_1 << "\t" << param_string_2 << "\t" << repeat_events << "\t" << SEED << endl;
	gDirectory->ls();

	RapidRun* fitting = new RapidRun( args );
	int result = fitting->run();
	gApplication->Terminate();
}
