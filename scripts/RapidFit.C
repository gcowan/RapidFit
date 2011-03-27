void RapidFit(char* XMLStr="XML_UB_B.xml", char* param_string_1="", char* param_string_2="" )
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
	args->Add(new TObjString( XMLStr ));
	args->Add(new TObjString("--defineContour"));
	args->Add(new TObjString( param_string_1 ));
	args->Add(new TObjString( param_string_2 ));
	args->Add(new TObjString("--doLLcontour"));

	cout << "Constructing Fit" <<endl;
	cout << param_string_1 << "\t" << param_string_2 << endl;
	gDirectory->ls();

	RapidRun* fitting = new RapidRun( args );
	int result = fitting->run();
	gApplication->Terminate();
}
