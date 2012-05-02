void test()
{
	TString XMLStr("Production_Pass3_UandB_v4_newAF.xml");

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

	RapidRun* fitting = new RapidRun( args );
	int result = fitting->run();
	gApplication->Terminate();

}
