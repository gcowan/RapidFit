
#include <string>

void runRapidScripts(std::string filename = "theFile.root", std::string treename = "DecayTree", std::string trigger = "ALL",  std::string trailer = "_Selected.root"){

  //  gSystem->Load("./ksfun_h.so");
  gSystem->Load("./AddRapidFitInfo_C.so");
  gSystem->Load("./ApplySelection_C.so");
  gSystem->Load("./PlotRapidVar_C.so");


  // add the rapid fit info
  std::string rapidFile = AddRapidFitInfo(filename,treename); 
  
  // make the selection
  std::string selectedFile = ApplySelection(rapidFile,treename, trigger,trailer);

  // make plots
  PlotRapidVar(selectedFile,treename);

}
