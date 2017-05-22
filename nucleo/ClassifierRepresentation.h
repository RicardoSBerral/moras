//---------------------------------------------------------------------------

#ifndef ClassifierRepresentationH
#define ClassifierRepresentationH
#include "Classifier.h"
#include <stdlib.h>
//---------------------------------------------------------------------------
/*class Propiedad
{
  protected:
    string nombre;
  public:
    Propiedad(string nombre){this->nombre = nombre;}
    virtual void DarValor(string valor)=0; //throws valor no válido
    virtual string LeerValor()=0;
    virtual string Nombre() {return nombre;}
};
class PropiedadTextoCorto : public Propiedad
{
  protected:
    string valor;
  public:
    PropiedadTextoCorto(string kk){}
    virtual void DarValor(string valor) {this->valor = valor;}
    virtual string LeerValor() {return valor;}
};
class PropiedadEntero : public Propiedad
{
  protected:
    int valor;
  public:
    PropiedadEntero(){}
    virtual void DarValor(string valor) {this->valor = itoa(valor,0,10);}
    virtual string LeerValor() {return atoi(valor);}
};*/
class ClassifierFactoryRepresentation
{
  public:
    ClassifierRepresentation(){}
    virtual void Properties(){}
    virtual Classifier *CreateClassifier()=0;
    virtual AnsiString GetName()=0;
};
class RepresentacionDeConjunto : public ClassifierFactoryRepresentation
{
  int numClass;
  public:
    AnsiString GetName() {return "Conjunto de clasificadores";}
    Classifier *CreateClassifier() {return 0;}
};
class CARTBaggingRep : public RepresentacionDeConjunto
{
  public:
    AnsiString GetName() {return "Bagging CART";}
    Classifier *CreateClassifier();
};
class IGPBaggingRep : public RepresentacionDeConjunto
{
  public:
    AnsiString GetName() {return "Bagging IGP";}
    Classifier *CreateClassifier();
};
class ListaDeRepresentaciones
{
  public:
    vector<ClassifierFactoryRepresentation*> Representaciones;
    ListaDeRepresentaciones();
};
#endif
