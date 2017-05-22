//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "ClassifierRepresentation.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)
Classifier *CARTBaggingRep::CreateClassifier()
{
  return 0;
}
//---------------------------------------------------------------------------
Classifier *IGPBaggingRep::CreateClassifier()
{
  return 0;
}
//---------------------------------------------------------------------------
ListaDeRepresentaciones::ListaDeRepresentaciones()
{
  ClassifierFactoryRepresentation *cr = new IGPBaggingRep();
  Representaciones.push_back(cr);
  cr = new CARTBaggingRep();
  Representaciones.push_back(cr);
}
