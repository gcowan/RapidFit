
#include "TList.h"
#include "TObject.h"
#include "TClass.h"

#include "Template_Functions.h"

#include <iostream>

using namespace::std;

void PrintPrimatives( TList* thisList )
{
        //TList* thisList = thisObject->GetListOfPrimitives();

        //cout << "TObject: " << thisObject->GetName() << "\t of type: " << thisObject->GetTypeName() << "\tat: " << thisObject << endl;
        cout << "I have knowledge of: " << endl;

        for( unsigned int i=0; i< (unsigned) thisList->GetSize(); ++i )
        {
                TObject* objPointer = thisList->At( i );
                TString Name = objPointer->GetName();
                TString Type = objPointer->IsA()->GetImplFileName();
                cout << "Name: " << Name << "\tType: " << Type << "\tAt: " << objPointer << endl;
        }

        cout << endl;
}

