//---------------------------------------------------------------------------


#include "Lector.h"
#include "Matriz.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <exception>
#include <stdexcept>
#include <iterator>
#include <typeinfo>
#include <cctype>
#include <functional>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

string Mandato::VERDADERO("true");
string Mandato::FALSO("false");
int Mandato::NivelSalida = 1;
EventoSalida Mandato::Salida = 0;
vector<Mandato*> Mandato::Mandatos;
//---------------------------------------------------------------------------
Mandato::Mandato()
{
  Mandatos.push_back(this);
}
Mandato::~Mandato()
{
  vector<Mandato*>::iterator pos = find(Mandatos.begin(), Mandatos.end(), this);
  if (pos!=Mandatos.end()) Mandatos.erase(pos);
//  Variable::EliminarVars();
}
//---------------------------------------------------------------------------
Tipo* Mandato::DameTipo()
{
  if (EsCadena())   return Tipo::TCadena();
  if (EsNumero())   return Tipo::TNumerico();
  if (EsEntero())   return Tipo::TEntero();
  if (EsBooleano()) return Tipo::TBooleano();
  if (EsDatos())    return Tipo::TDatos();

  return Tipo::TNinguno();
}
bool Mandato::EsTipo(Tipo *tipo)
{
  return this->DameTipo()->EsTipo(tipo);//typeid(*tipo)==typeid(*this);
}
//---------------------------------------------------------------------------
    /**
      * Busca cod dentro de exp y devuelve la posiciOn
      * siempre que cod estE fuera de parEntesis y comillas.
      * si no devuelve -1;
      */
//---------------------------------------------------------------------------
int Mandato::BuscaFueraDeParentesis(string exp, string cod)
{
/*  unsigned I;
  unsigned encontrado=0;
  for (unsigned i=0;i<exp.size();i++) {
    if (exp.at(i)=='(') {
      I=exp.find(")", i+1);
      if (I>i) i=I;
      else throw logic_error(ACad((int)i));
      encontrado = 0;
    }
    else if (exp.at(i)=='\"') {
      I=exp.find("\"", i+1);
      if (I>i) i=I;
      else throw logic_error(ACad((int)i));
      encontrado = 0;
    }
    else if (exp.at(i)==cod.at(encontrado)) encontrado++;
    else encontrado = 0;
    if (encontrado==cod.length())
      return i-encontrado+1;
  }
  return -1;*/
  string::size_type I;
  int nparen=0;
  int encontrado=0;
  for (int i=0;i<(int)exp.length();i++) {
    if (exp[i]=='(') nparen++;
    else if (exp[i]==')') nparen--;
    else if (exp[i]=='\"') {
      I=exp.find("\"", i+1);
      if ((int)I>i) i=I;
      else throw logic_error(ACad((int)i));
      encontrado = 0;
    }
    else if (exp[i]==cod[encontrado]) encontrado++;
    else encontrado = 0;
    if (encontrado==(int)cod.length()) {
      if (nparen==0) return i-encontrado+1;
      else encontrado = 0;
    }
  }
  return -1;
}
//---------------------------------------------------------------------------
    /**
      * Busca cod dentro de exp y devuelve la posiciOn
      * siempre que cod estE fuera de parEntesis y comillas.
      * si no devuelve -1;
      */
//---------------------------------------------------------------------------
int Mandato::BuscaFueraDeParentesisInv(string exp, string cod)
{
  string::size_type I;
  int nparen=0;
  int encontrado=cod.length()-1;
  for (int i=exp.length()-1;i>=0;i--) {
    if (exp[i]==')') nparen++;
    else if (exp[i]=='(') nparen--;
    else if (exp[i]=='\"') {
      I=exp.find_last_of("\"", i-1);
      if ((int)I<i) i=I;
      else throw logic_error(ACad((int)i));
      encontrado = cod.length()-1;
    }
    else if (exp[i]==cod[encontrado]) encontrado--;
    else encontrado = cod.length()-1;
    if (encontrado==-1) {
      if (nparen==0) return i;
      else encontrado = cod.length()-1;
    }
  }
  return -1;
  /*string::size_type I;
  int encontrado=cod.length()-1;
  for (int i=exp.size()-1;i>=0;i--) {
    if (exp.at(i)==')') {
      I=exp.find_last_of("(", i-1);
      if (I<i) i=I;
      else throw logic_error(ACad((int)i));
      encontrado = cod.length()-1;
    }
    else if (exp.at(i)=='\"') {
      I=exp.find_last_of("\"", i-1);
      if (I<i) i=I;
      else throw logic_error(ACad((int)i));
      encontrado = cod.length()-1;
    }
    else if (exp.at(i)==cod.at(encontrado)) encontrado--;
    else encontrado = cod.length()-1;
    if (encontrado==-1)
      return i;
  }
  return -1;*/
}
//---------------------------------------------------------------------------
    /*
     * Covierte un nUmero a cadena
     */
//---------------------------------------------------------------------------
string Mandato::ACad(double num)
{
  char str[25];
  gcvt(num, 15, str);
  string cad(str);
  return cad;
}
string Mandato::ACad(int num)
{
  char str[25];
  sprintf(str, "%d", num);//itoa(num, str, 10);
  string cad(str);
  return cad;
}
//---------------------------------------------------------------------------
    /*
     * Quita espacios del final e inicio de una expresiOn
     */
//---------------------------------------------------------------------------
void Mandato::Trim(string *exp)
{
  unsigned ini, fin, i, i2;
  if (exp->size()==0) return;
  char *cad = new char[exp->size()+1];
  for(ini=0;ini<exp->size();ini++) if ((*exp)[ini]!=' ') break;
  for(fin=exp->size()-1;fin>ini;fin--) if ((*exp)[fin]!=' ') break;
  for(i=ini, i2=0;i<=fin;i++, i2++) cad[i2] = (*exp)[i];
  cad[i2] = '\0';
  (*exp) = cad;
  delete []cad;
}
//---------------------------------------------------------------------------
    /*
     * Quita espacios y parEntesis que engloben a una exesiOn
     */
//---------------------------------------------------------------------------
void Mandato::TrimExtendida(string *exp)
{
  unsigned ini=0, fin=0, nparen=0, i;

  /*for(int ikk=0;ikk<exp.size();ikk++)printf( "(ANTES:%d-%d" , exp->size(),
		  (int)(*exp)[exp->size()-1]);*/ 
  if ((*exp)[exp->size()-1]==13){ 
	  exp->erase(exp->end()-1);
  /*for(int ikk=0;ikk<exp.size();ikk++)printf( "DESPU:%d-%d)\n" ,exp->size(),
		  (int)(*exp)[exp->size()-1]);*/ 
  }
  Trim(exp);
  char *cad = new char[exp->size()+1];
  strcpy(cad, exp->c_str());
  ini = fin = 0;
  if (cad[0]=='+') ini++;
  if (cad[ini]=='(') {
    for (i=ini;i<strlen(cad);i++) {
      if (cad[i]=='(') nparen++;
      else if (cad[i]==')') nparen--;
      if (nparen==0) break;
    }
    i++;
    if (nparen==0 && i==strlen(cad)) {
      ini++;
      fin++;
    }
    else if (nparen>0 && i==strlen(cad)) {
      throw logic_error("Parentesis no cerrado");
    }
  }
  if (fin>0) {
    (*exp) = cad;
    (*exp) = exp->substr(ini, exp->size() - fin - 1);
    TrimExtendida(exp);
  }
  Trim(exp);
  delete []cad;
}

//---------------------------------------------------------------------------
Mandato* Mandato::valueOf(string exp)
{
  Mandato* mand = 0;
#ifdef _DEBUG
  printf("MND: %s", exp.c_str()); fflush(stdout);
#endif
  TrimExtendida(&exp);
  if (exp[0]=='/' && exp[1]=='/') return 0; //Comentario
  if (exp.size()==0) return 0;
  int i;
  if (0 < (i=BuscaFueraDeParentesis(exp, "=")) ) {
    string kk = exp.substr(0, i);
    TrimExtendida(&kk);
    mand = new MandatoAsignacion(valueOf(exp.substr(i+1)), kk);
  }
  else if (0 < (i=BuscaFueraDeParentesisInv(exp, "+")) ) {
    Mandato *s1 = valueOf(exp.substr(0, i));
    Mandato *s2 = valueOf(exp.substr(i+1));
    if (s1->EsCadena()) mand = new Concatenacion(s1, s2);
    else                mand = new Suma(s1, s2);
  }
  else if (0 == (i=BuscaFueraDeParentesisInv(exp, "-")) )
    mand = new Producto(new Numero(-1.0), valueOf(exp.substr(1)));
  else if (i >0)
    mand = new Resta(valueOf(exp.substr(0, i)), valueOf(exp.substr(i+1)));
  else if (0 < (i=BuscaFueraDeParentesisInv(exp, "*")) )
    mand = new Producto(valueOf(exp.substr(0, i)), valueOf(exp.substr(i+1)));
  else if (0 < (i=BuscaFueraDeParentesisInv(exp, "/")) )
    mand = new Division(valueOf(exp.substr(0, i)), valueOf(exp.substr(i+1)));
  else if (0 != (mand=Variable::VariablePorNombre(exp)) ) {
  }
  else if (0 != (mand=ifns().CrearFuncionDeCadena(&exp)) ) {
  }
  else if (exp[0]=='\"' && exp[exp.size()-1]=='\"') {
    char *kk = new char[exp.size()];
    exp = exp.substr(1, exp.size()-2);
    sprintf(kk, "%s", exp.c_str());
    mand = new Cadena(kk);
    delete []kk;
  }
  else
    mand = new Numero((char*)exp.c_str());

  return mand;
}
//---------------------------------------------------------------------------
void Mandato::salida(char* sal)
{
  if (Salida) Salida(sal);
}
//---------------------------------------------------------------------------
void Mandato::EliminarMandatos()
{
  for(unsigned i=0;i<Mandatos.size();i++)
    delete Mandatos[i];
  Mandatos.clear();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Tipo* Tipo::_TCadena    = new TipoCadena();
Tipo* Tipo::_TNumerico  = new TipoNumerico();
Tipo* Tipo::_TEntero    = new TipoEntero();
Tipo* Tipo::_TBooleano  = new TipoBooleano();
Tipo* Tipo::_TDatos      = new TipoDatos();
Tipo* Tipo::_TMatriz      = new TipoMatriz();
Tipo* Tipo::_TNinguno      = new TipoNinguno();
Tipo* Tipo::_TCualquiera    = new TipoCualquiera();
Tipo* Tipo::_TFicheroEntrada = new TipoCadenaFicheroLectura();
Tipo* Tipo::_TFicheroSalida  = new TipoCadenaFicheroEscritura();
//---------------------------------------------------------------------------
bool Tipo::EsTipo(Tipo *tipo)
{
  return tipo->ConvertibleA(this);//typeid(*tipo)==typeid(*this);
}
bool Tipo::EsTipoBasico(Tipo *tipo)
{
  return tipo==TCadena() || tipo==TNumerico() || tipo==TEntero() || tipo==TBooleano();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Variable *TipoCadena::NuevaVariable(std::string nombre, Mandato *inic)
{
  VariableCadena *var = new VariableCadena(nombre);
  if (inic) var->Set(inic);
  return var;
}
Variable *TipoNumerico::NuevaVariable(std::string nombre, Mandato *inic)
{
  VariableNumero *var = new VariableNumero(nombre, 0.0);
  if (inic) var->Set(inic);
  return var;
}
Variable *TipoEntero::NuevaVariable(std::string nombre, Mandato *inic)
{
  VariableEntero *var = new VariableEntero(nombre, 0);
  if (inic) var->Set(inic);
  return var;
}
Variable *TipoBooleano::NuevaVariable(std::string nombre, Mandato *inic)
{
  VariableBooleano *var = new VariableBooleano(nombre, false);
  if (inic) var->Set(inic);
  return var;
}
Variable *TipoDatos::NuevaVariable(std::string nombre, Mandato *inic)
{
  Tipo *t = inic && inic->DameTipo()->EsTipo(this) ? inic->DameTipo() : 0;
  if (!t) t = this;
  VariableDatos *var = new VariableDatos(nombre, 0, t);
  if (inic) var->Set(inic);
  return var;
}
Variable *TipoMatriz::NuevaVariable(std::string nombre, Mandato *inic)
{
  VariableMatriz *var = new VariableMatriz(nombre);
  if (inic) var->Set(inic);
  return var;
}
/*Variable *TipoCompuesto::NuevaVariable(std::string nombre, Mandato *inic=0)
{
  VariableCompuesta *var = new VariableCompuesta(nombre);
  if (inic) { //Ver si inic es de tipo compuesto ...
  }
  return var;
}*/
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
map<std::string, Variable*> Variable::Variables;
map<std::string, std::vector<Variable *> * > Variable::esps;
int Variable::VarsAnon = 0;
//---------------------------------------------------------------------------
Variable::Variable(string nombre) : Mandato()
{
  if (nombre=="") {//Creamos una variable anonima
    do {
      char kk[24];
      memset(kk, 0, 24);
      sprintf(kk, "Var%d", VarsAnon++);
      nombre=kk;
    }
    while (Variables[nombre]!=0);
    Variables[nombre] = this;
    ParteNombre(nombre, FEspacioDeNombres, this->FNombre);
  }
  else {
    if (!NombreValido(nombre)) {
      throw invalid_argument(string("Nombre de variable (") + nombre +
                                                               ") no valido");
    }

    ParteNombre(nombre, FEspacioDeNombres, this->FNombre);
  
    if (!VariablePorNombre(FNombre, FEspacioDeNombres)) {
      Variables[nombre] = this;
      if (!edn(FEspacioDeNombres)) { 
//        cout << "Nuevo espacio <" << FEspacioDeNombres <<  ">" << endl;
        esps[FEspacioDeNombres] = new vector<Variable*>;
      }
      edn(FEspacioDeNombres)->push_back(this);
    }
  }

}
void Variable::print()
{
  map<string, vector<Variable *> * >::iterator ii;
  
  cout << "*************************************" << endl;
  for(ii=esps.begin();ii!=esps.end();++ii) {
    cout << "<" << ii->first << ">[" << ii->second->size() << "]" << endl;
    vector<Variable *>::iterator vv = ii->second->begin();
    for(;vv!=ii->second->end();++vv)  {
      cout << "\t" << (*vv)->NombreCorto() << " = " << (*vv)->ComoCadena();
      cout  << "\t[" << (*vv)->DameTipo()->nombre() << "]" << endl;
    }
  }
}
Variable::Variable() : Mandato()
{
  char kk[24];
  sprintf(kk, "Var%d", VarsAnon);
  this->FNombre = kk;
//  if (!VariablePorNombre(kk)) Variables.push_back(this);
  VarsAnon++;
}
Variable::~Variable()
{
  EliminarVarDeListas(NombreCompleto());
  /*for (unsigned i=0;i<Variables.size();i++) {
    if(Variables[i]==this) {
      Variables[i] = Variables[Variables.size()-1];
      Variables.pop_back();
    }
  }*/
}
vector<Variable *> * Variable::edn(std::string Namespace)
{
  return esps[Namespace];
}
string Variable::UneNombre(const string &nombre, const string &EspacioDeNombres)
{
  if (EspacioDeNombres.size()==0) return nombre;

  return EspacioDeNombres + "." + nombre;
}
void Variable::ParteNombre(const string NombreLargo, string &EspacioDeNombres,
                                                           string &NombreCorto)
{
  int i = NombreLargo.find_last_of(".");
  NombreCorto = NombreLargo.substr(i<0?0:i+1);
  EspacioDeNombres = NombreLargo.substr(0, i<0?0:i);
}

Variable* Variable::VariablePorNombre(string nom, string Namespace, bool force,
 Mandato *inic)
{
  string nl = UneNombre(nom, Namespace);

  if (!NombreValido(nl)) return 0;

  Variable *v = Variables[nl];

//  cout << "Buscando " << nl <<  ": ";
//  cout << (v ? "ENCONTRADA" : "") << endl;

  /*for(unsigned i=0;i<Variables.size();i++) {
    if (Variables[i]->FNombre.compare(nom)==0) {
      return Variables[i];
    }
  }*/

  if (!v && force) {
    v = inic ? inic->DameTipo()->NuevaVariable(nom, inic) : 
                                                new VariableNumero(nom, 0.0);
  }

  return v; 
}
int Variable::NumVariables(string Namespace)
{
  if (edn(Namespace)) return edn(Namespace)->size();
  return 0;
}
Variable* Variable::VariablePorNumero(int i, string Namespace)
{
  if (!edn(Namespace)) return 0;
  return (*edn(Namespace))[i];

}
void Variable::EliminarVarDeListas(string nombre)
{
  Variable *v = VariablePorNombre(nombre);

  if (v) {
    //Eliminamos de la lista general
    Variables[nombre] = 0;

    //Eliminamos de la lista de espacios de nombres
    string esp, nom;
    ParteNombre(nombre, esp, nom);
    vector<Variable*>::iterator it;
    it = find(esps[esp]->begin(), esps[esp]->end(), v);
    if (it!=esps[esp]->end()) esps[esp]->erase(it);
  }
}
void Variable::EliminarVars(string Namespace)
{
  if (edn(Namespace)) {
    while(esps[Namespace]->begin()!=esps[Namespace]->end()) {
      delete *esps[Namespace]->begin();
    }
  }
}
bool Variable::NombreValido0(string nom)
{
  if (nom.size()==0 || (!isalpha(nom[0]) && nom[0]!='_') ) return false;

  for(unsigned i=1;i<nom.size();i++) {
    if (!isalnum(nom[i]) && nom[i]!='_') return false;
  }

  return true; 
}
bool Variable::NombreValido(string nom, string Namespace)
{
  string nl;

  nl = UneNombre(nom, Namespace);

  int i, pos = 0;
  do {
    i = nl.find(".", pos);
    string g = nl.substr(pos, i<0 ? nl.size() : i-pos);
    if (!NombreValido0(g)) return false;
    pos = i+1;
  }
  while (i>0);

  return true; 
}
vector<Variable*> Variable::VariablesPorTipo(Tipo *tipo, string Namespace)
{
  vector<Variable*> ret;

  if (edn(Namespace)) {
    vector<Variable*> *Variables = edn(Namespace);
    vector<Variable*>::iterator i = Variables->begin();

    for( ;i!=Variables->end();i++)
      if ((*i)->DameTipo()->EsTipo(tipo))
        ret.push_back(*i);
  }

  return ret;
}
//---------------------------------------------------------------------------
Tipo* VariableMutante::DameTipo()
{
  return FVar ? FVar->DameTipo() : Tipo::TCualquiera();
}
bool VariableMutante::EsCadena()
{
  return FVar ? FVar->EsCadena() : false;
}
bool VariableMutante::EsNumero()
{
  return FVar ? FVar->EsNumero() : false;
}
bool VariableMutante::EsEntero()
{
  return FVar ? FVar->EsEntero() : false;
}
bool VariableMutante::EsBooleano()
{
  return FVar ? FVar->EsBooleano() : false;
}
bool VariableMutante::EsDatos()
{
  return FVar ? FVar->EsDatos() : false;
}
bool VariableMutante::EsVectorNumeros()
{
  return FVar ? FVar->EsVectorNumeros() : false;
}
bool VariableMutante::EsMatrizNumeros()
{
  return FVar ? FVar->EsMatrizNumeros() : false;
}
//---------------------------------------------------------------------------
string VariableCadena::ComoCadena()
{
  return FValor;
}
double VariableCadena::ComoNumero()
{
  return atof(FValor.c_str());
}
int VariableCadena::ComoEntero()
{
  return atoi(FValor.c_str());
}
bool VariableCadena::ComoBooleano()
{
  return FValor.compare(VERDADERO)==0;
}
void* VariableCadena::ComoDatos()
{
  return (void*)FValor.c_str();
}
void VariableCadena::SetCadena(string valor)
{
  FValor = valor;
}
void VariableCadena::SetNumero(double valor)
{
  char str[25];
  gcvt(valor, 5, str);
  FValor = str;
}
void VariableCadena::SetEntero(int valor)
{
  char str[25];
  sprintf(str, "%d", valor);//itoa(valor, str, 10);
  FValor = str;
}
void VariableCadena::SetBooleano(bool valor)
{
  if (valor) FValor = Mandato::VERDADERO;
  else FValor = Mandato::FALSO;
}
void VariableCadena::SetDatos(void* valor)
{
  FValor = (char*)valor;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
VariableBooleano *VariableBooleano::VarTrue  = GenVarTrue();
VariableBooleano *VariableBooleano::VarFalse = GenVarFalse();
VariableBooleano *VariableBooleano::GenVarTrue()
{
  return new VariableBooleano(Mandato::VERDADERO, true);
}
VariableBooleano *VariableBooleano::GenVarFalse()
{
  return new VariableBooleano(Mandato::FALSO, false);
}
string VariableBooleano::ComoCadena()
{
  return ComoBooleano() ? Mandato::VERDADERO : Mandato::FALSO;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
std::string VariableMatriz::ComoCadena()
{
  std::ostringstream buf;
  if (GetMatriz()) GetMatriz()->saveToStream(buf);
  return buf.str();
}
void VariableMatriz::SetCadena(std::string valor)
{
//Algo habra' que hacer
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
MandatoAsignacion::MandatoAsignacion(Mandato *mand, string NomVar) : Mandato()
{
  FMandato = mand;

  if (FMandato->DameTipo()==Tipo::TNinguno()) 
    throw MandatoError("Mandato no asignable, no devuelve ningún valor");

  FVar=Variable::VariablePorNombre(NomVar);
  if (FVar && FVar->DameTipo()==FMandato->DameTipo()) 
    ;
  else {
    if (FVar) delete FVar;//Si existe se cambia de tipo
    FVar = FMandato->DameTipo()->NuevaVariable("", FMandato);
    FVar =  new VariableMutante(NomVar, FVar);//La encapsulo
  }

/*if (FMandato->EsCadena())
    FVar = new VariableCadena(NomVar, mand->ComoCadena());
  else if (FMandato->EsNumero())
    FVar = new VariableNumero(NomVar, mand->ComoNumero());
  else if (FMandato->EsDatos())
    FVar = new VariableDatos(NomVar, mand->ComoDatos());
  else
    FVar = new VariableCadena(NomVar);*/
}
MandatoAsignacion::MandatoAsignacion(Mandato *mand, Variable *Var) : Mandato()
{
  FMandato = mand;
  FVar = Var;
}
Mandato* MandatoAsignacion::Ejecutar()
{
  FMandato->Ejecutar();

  //Mutamos el tipo de la variable si FMandato ha cambiado
  if (FVar->DameTipo()!=FMandato->DameTipo()) {
    VariableMutante *vm = dynamic_cast<VariableMutante*>(FVar);
    Variable *v = vm ? vm->GetVar() : FVar;
    string NomVar = v->NombreCompleto();
    delete v;
    v = FMandato->DameTipo()->NuevaVariable("", FMandato);
    if (vm) vm->SetVar(v);
    else FVar = v;
  }

  FVar->Set(FMandato);

/*  if (FVar->EsCadena()) FVar->ComoCadena(FMandato->ComoCadena());
  else if (FVar->EsNumero()) FVar->ComoNumero(FMandato->ComoNumero());
  else if (FVar->EsDatos()) FVar->ComoDatos(FMandato->ComoDatos());*/

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

string Cadena::ComoCadena()
{
  return cad;
}
//---------------------------------------------------------------------------
double Cadena::ComoNumero()
{
  return atof(cad.c_str());
}
//---------------------------------------------------------------------------
int Cadena::ComoEntero()
{
  return atoi(cad.c_str());
}
//---------------------------------------------------------------------------
bool Cadena::ComoBooleano()
{
  return cad.compare(VERDADERO)==0;
}
//---------------------------------------------------------------------------
void* Cadena::ComoDatos()
{
  return (void*)cad.c_str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
bool Numero::ComoBooleano()
{
  if (!ejec) Ejecutar();
  return num!=0.0;
}
//---------------------------------------------------------------------------
string Numero::ComoCadena()
{
  if (!ejec) Ejecutar();
  return ACad(num);
}
//---------------------------------------------------------------------------
void* Numero::ComoDatos()
{
  if (!ejec) Ejecutar();
  return (void*)&num;
}
//---------------------------------------------------------------------------
int Numero::ComoEntero()
{
  if (!ejec) Ejecutar();
  return num>0.0 ? (int)(num + 0.5) : (int)(num - 0.5);
}
//---------------------------------------------------------------------------
double Numero::ComoNumero()
{
  if (!ejec) Ejecutar();
  return num;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Parametro::checkValue(Mandato *mnd)
{
  if (!DameTipo()->EsTipo(Tipo::TNinguno()) && mnd) {
    if (!mnd->EsTipo(DameTipo())) throw TipoIncorrecto("Tipo incorrecto");
  }
  else if (!mnd && !Opcional()) throw ParametroNoOpcional("");
}
Mandato* Parametro::Ejecutar()
{
  if (mnd) mnd->Ejecutar();
  checkValue(mnd);
  return this;
}
void Parametro::PonPropiedades(bool o, Tipo *t, string n, string d)
{
  mnd = 0;
  opcional = o;
  tipo = t;
  nombre = n;
  descripcion = d;
  poromision = "";
}
//---------------------------------------------------------------------------
void ParametroNumero::checkValue(Mandato *mnd)
{
   Parametro::checkValue(mnd);
  if (menor<=mayor && (mnd->ComoNumero()>mayor || mnd->ComoNumero()<menor))
    throw RangoNoValido(string("Parametro \"") + nombre + "\" fuera de rango");
}
//---------------------------------------------------------------------------
void ParametroEntero::checkValue(Mandato *mnd)
{
  if (!mnd && !Opcional()) throw ParametroNoOpcional("");
  if (!mnd) return;
  if (!mnd->EsTipo(Tipo::TNumerico()))
    throw TipoIncorrecto("Parametro no numerico");
  if (menor<=mayor && (mnd->ComoEntero()>mayor || mnd->ComoEntero()<menor))
    throw RangoNoValido(string("Parametro \"") + nombre + "\" fuera de rango");
}
//---------------------------------------------------------------------------
void ParametroCadena::checkValue(Mandato *mnd)
{
  if (!mnd && !Opcional()) throw ParametroNoOpcional("");
  if (!mnd) return;
  if (CadenasValidas.size()>0 && CadenasValidas.end()==
            find(CadenasValidas.begin(), CadenasValidas.end(), ComoCadena()))
    throw RangoNoValido("Texto no valido");
}
string ParametroCadena::ComoCadena()
{
  string cad = mnd->ComoCadena();

  if (CadenasValidas.size()==0) return cad;

  if (CadenasValidas.end()!=
           find(CadenasValidas.begin(), CadenasValidas.end(), cad)) return cad;

  int n = (int)mnd->ComoNumero();
  if (n>0 && n<=(int)CadenasValidas.size()) 
    return CadenasValidas[n-1];

  return cad;
}
int ParametroCadena::ComoEntero()
{
  if (!mnd) return 0;
  if (CadenasValidas.size()==0) return mnd->ComoEntero();
  string cad = ComoCadena();
  for (unsigned i=0;i<CadenasValidas.size();i++) {
    if (CadenasValidas[i]==cad) return i+1;
  }
  return 0;
}
//---------------------------------------------------------------------------
void ParametroFicheroEntrada::checkValue(Mandato *mnd)
{
  if (mnd->ComoCadena()=="" && Opcional()) return;
  FILE *f = fopen(mnd->ComoCadena().c_str(), "r");
  if (!f) throw NombreFicheroNoValido( mnd->ComoCadena() );
  fclose(f);
}
//---------------------------------------------------------------------------
void ParametroFicheroSalida::checkValue(Mandato *mnd)
{
  if (mnd->ComoCadena()=="" && Opcional()) return;
  FILE *f = fopen(mnd->ComoCadena().c_str(), "a");
  if (!f) throw NombreFicheroNoValido( mnd->ComoCadena() );
  fclose(f);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Funcion::Funcion(vector<Mandato *> params) : Mandato()
{
  Parametros = params;//new Mandato*[NumParams()];
//  for(unsigned i=Parametros.size();i<NumParams();i++)
//    Parametros.push_back(0);
//  for(unsigned i=0;i<NumParams() && i<numperams();i++)
//    Parametros[i] = params[i];
}                             
Funcion::Funcion() : Mandato()
{
//  Parametros = 0;
//  for(unsigned i=0;i<NumParams();i++)
//    Parametros.push_back(0);
}
Funcion::~Funcion()
{
//  if (Parametros) delete []Parametros;
}
/*vector<Mandato *> *Funcion::GetParametros()
{
  for(unsigned i=Parametros.size();i<NumParams();i++)
    Parametros.push_back(0);

  return &Parametros;
}*/
int Funcion::NumParamsAsignados()
{
  int pa = 0;
  
  for(unsigned i=0;i<Parametros.size();i++)
    if (ParamInfo(i)->Asignado()) pa++;

  return pa;
}
Parametro* Funcion::ParamInfo(int idx)
{
  if (int(Parametros.size())<=idx) {
    for(int i=Parametros.size();i<idx+1;i++)
      Parametros.push_back(new Parametro());
  }

  return dynamic_cast<Parametro*>(Parametros[idx]);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
ListaInterfazFunciones ListaInterfazFunciones::_ListaInterfazFunciones;
ListaInterfazFunciones &ifns()
{
  return ListaInterfazFunciones::lif();
}
//---------------------------------------------------------------------------
void ListaInterfazFunciones::AnadirInterfazFuncion(IntefazFuncion *inf, 
                                                             std::string cat) 
{
  iFunciones.push_back(inf);
  FunsCat.push_back(cat);
  InfoFunciones.push_back(inf->CrearFuncion());
  if ( Categorias.end()==find(Categorias.begin(), Categorias.end(), cat)) 
    Categorias.push_back(cat);
}
Funcion *ListaInterfazFunciones::CrearFuncionDeCadena(string *cadena)
{
  Mandato::Trim(cadena);
  string::size_type I = cadena->find("(");
  if (I==string::npos) return 0;
  string fn = cadena->substr(0, I);  Mandato::Trim(&fn);
  string pars = cadena->substr(I);   Mandato::Trim(&pars);

  for(unsigned i=0;i<iFunciones.size();i++) {
    if (iFunciones[i]->NombreFuncion().compare(fn)==0) {
      Funcion * func = iFunciones[i]->CrearFuncion();
      CrearParametrosDeCadena(&pars, func);
      return func;
    }
  }
  return 0;
}
void ListaInterfazFunciones::CrearParametrosDeCadena(string *cadena,
                                                                  Funcion *f)
{
  if ((*cadena)[0]!='(') throw "Error on CrearParametrosDeCadena";
  Mandato::TrimExtendida(cadena);
  int pos;
//  vector<Mandato*> *pars = f->GetParametros();
//  pars->clear();
  Parametro *par;

  int i=0;
  while( 0<=(pos = Mandato::BuscaFueraDeParentesis(*cadena, ",")) ) {
    string kk = cadena->substr(0, pos);
    par = f->ParamInfo(i);
    par->SetMandato(Mandato::valueOf(kk));
    (*cadena) = cadena->substr(pos+1);
    i++;
  }
  f->ParamInfo(i)->SetMandato(Mandato::valueOf(*cadena));
}
Funcion *ListaInterfazFunciones::InfoFuncionPorNombre(std::string nom)
{
  for(unsigned i=0;i<iFunciones.size();i++) {
    if (iFunciones[i]->NombreFuncion().compare(nom)==0) {
      return InfoFuncion(i);
    }
  }

  return 0;
}
string ListaInterfazFunciones::DameCategoria(int i)
{
  if (i<0 || i>=NumCategorias()) return "";
  return Categorias[i];
}
vector<int> ListaInterfazFunciones::IndicesFuncionesPorCat(string cat)
{
  vector<int> ret;

  for(unsigned i = 0; i<iFunciones.size(); i++) {
    if (FunsCat[i]==cat) ret.push_back(i);
  }

  return ret;
}
vector<IntefazFuncion*> ListaInterfazFunciones::InterfacesPorCat(string cat)
{
  vector<IntefazFuncion*> ret;

  for(unsigned i = 0; i<iFunciones.size(); i++) {
    if (FunsCat[i]==cat) ret.push_back(iFunciones[i]);
  }

  return ret;
}
vector<IntefazFuncion*> ListaInterfazFunciones::InterfacesPorTipo(Tipo *tipo)
{
  vector<IntefazFuncion*> ret;

  for(unsigned i = 0; i<iFunciones.size(); i++) {
    if (InfoFunciones[i]->TipoRetorno()->EsTipo(tipo))
      ret.push_back(iFunciones[i]);
  }

  return ret;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Seno::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Seno());
}
//---------------------------------------------------------------------------
Mandato* Seno::Ejecutar()
{
  Param(0)->Ejecutar();
  res = sin(Param(0)->ComoNumero());
  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Coseno::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Coseno());
}
//---------------------------------------------------------------------------
Mandato* Coseno::Ejecutar()
{
  Param(0)->Ejecutar();
  res = cos(Param(0)->ComoNumero());
  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

