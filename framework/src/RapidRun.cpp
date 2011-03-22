#include <iostream>
using std::cout;
using std::endl;

#include "main.h"
#include "RapidRun.h"
#include "XMLConfigReader.h"
#include "TList.h"
#include "TObject.h"
#include "TObjString.h"
#include "TString.h"

#include <cstring>

ClassImp(RapidRun)//needed for Cint

RapidRun::RapidRun( TList* _args ):
    args(_args)
{
}

int RapidRun::run(){
  cout << TString("HELLO") <<endl;
  //convert the TList args into a main compatible list
  const Int_t argc = args->GetSize();
  //memory is dynamically allocated
  char** argv = new char*[argc];

  TObjLink *lnk = args->FirstLink();
  Int_t count = 0;
  while (lnk) {
    const char* str = ((TObjString*)lnk->GetObject())->String();
    argv[count] = new char[strlen(str)];
    strcpy(argv[count],str);
    lnk = lnk->Next();
    count++;
  }
  cout << count << "\t" << argv[0] << "\t" << argv[1] << "\t" << argv[2] <<endl;
 //                XMLConfigReader* xmlFile=NULL;
 //                       xmlFile = new XMLConfigReader("/home/robert/PHYS/Tagged-Pass9-UB-XML1/Tagged-Pass9-UB-XML1/DataFile_Tagged_cFit_Working-forMainFit-mistagObservable-Unbiased-XML1.xml"); 

  Int_t result = RapidFit(count, argv);
  //clean up
  //for(int i = 0; i < argc; i++){
  //  delete argv[i];
  //}
  //delete argv;

  return result;
}

