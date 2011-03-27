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
  cout << "RapidRun Launching RapidFit!" <<endl;

  Int_t result = RapidFit(count, argv);

  return result;
}

