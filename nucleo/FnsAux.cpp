//---------------------------------------------------------------------------


#include "FnsAux.h"
#include "util20.h"
#include "Utils.h"
#include "Matriz.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <typeinfo>

#if defined(__STRICT_ANSI__) or defined(WIN32)
  #include <time.h>
#else
  #include <sys/times.h>
  #include <unistd.h>
#endif


//---------------------------------------------------------------------------
using namespace std;
//---------------------------------------------------------------------------
string &replace_all(string s_find, string s_repl, string &s)
{
  size_t ini = 0;
  int lon_f = s_find.size();
  int lon_r = s_repl.size();

  while( (ini=s.find(s_find, ini)) != string::npos) {
    s.replace(ini, lon_f, s_repl);
    ini += lon_r;
  };

  return s;
}
//---------------------------------------------------------------------------
void DarDeAltaFnsAux()
{
  DameNombreFichero::DarDeAlta();
  NombreProc::DarDeAlta();
  Append::DarDeAlta();
  Mean::DarDeAlta();
  StartRand::DarDeAlta();
  Call::DarDeAlta();
  Time::DarDeAlta();
  Repeat::DarDeAlta();
  If::DarDeAlta();
  Print::DarDeAlta();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void DameNombreFichero::CreateParams()
{
  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("NewFileName");
  CadenasValidas.push_back("PIDFileName");
  CadenasValidas.push_back("ProcNameFileName");
  CadenasValidas.push_back("tmpfile");
  Parametro *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(false, 0, "Tipo de nombre", "");
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
void DameNombreFichero::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new DameNombreFichero(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* DameNombreFichero::Ejecutar()
{
  Param(0)->Ejecutar();

  string prev = "";
  string ext = "";
  string dir = ".";

  if ( NumParamsAsignados() > 1 ) {
    prev = Param(1)->Ejecutar()->ComoCadena();
  }
  if ( NumParamsAsignados() > 2 ) {
    ext = Param(2)->Ejecutar()->ComoCadena();
  }
  if ( NumParamsAsignados() > 3 ) {
    dir = Param(3)->Ejecutar()->ComoCadena();
  }

  int what = Param(0)->ComoEntero();
  if (what==1) { //NewFileName
    txt = GenerateNewFileName(prev, ext, dir);
  }
  else if (what==2) { //PIDFileName
    txt = GeneratePIDFileName(prev, ext, dir);
  }
  else if (what==3) { //ProcNameFileName
    txt = GenerateProcNameFileName(prev, ext, dir);
  }
  else if (what==4) { //tmpfile
    txt = tmpnam(0);
  }
  
  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void NombreProc::CreateParams()
{
  ParametroCadena *p2 = new ParametroCadena();
  p2->PonPropiedades(false, 0, "Nombre del proceso", 
     "Se usa para dar nombres a ficheros de salida y para lo que haga falta");
  Parametros.push_back(p2);

}
//---------------------------------------------------------------------------
void NombreProc::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new NombreProc(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* NombreProc::Ejecutar()
{
  Param(0)->Ejecutar();

  SetProcName(Param(0)->ComoCadena());

  return this;
}
string NombreProc::ComoCadena()
{
  return GetProcName();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Append::CreateParams()
{
  ParametroCadena *p2 = new ParametroCadena();
  p2->PonPropiedades(false, 0, "Texto a añadir", 
                    "Texto que se añadirá al final del fichero especificado");
  Parametros.push_back(p2);

  Parametro *pfs = new ParametroFicheroSalida();
  pfs->PonPropiedades(false, 0, "Fichero de salida",
                                  "Nombre del fichero para guardar el texto");
  Parametros.push_back(pfs);
}
//---------------------------------------------------------------------------
void Append::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Append(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* Append::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  FILE *f=fopen(Param(1)->ComoCadena().c_str(), "at");
  if (!f) return this;

  string txt = Param(0)->ComoCadena();

  if (txt[txt.length()-1]!='\n') 
    txt = txt + "\n";

  fprintf(f, "%s", txt.c_str());

  fclose(f);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Mean::CreateParams()
{
  string fil = "Todos |*.csv;*.txt;";
  Parametro *pfe = new ParametroFicheroEntrada(fil);
  pfe->PonPropiedades(false, 0, "Fichero con una matriz",
                       "Nombre del fichero que contiene una matriz numérica");
  Parametros.push_back(pfe);
}
//--------------------------------------------------------------------------
void Mean::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Mean(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* Mean::Ejecutar()
{
  Param(0)->Ejecutar();

  char f[1024];
  strcpy(f, Param(0)->ComoCadena().c_str());
  Matriz *aux = Matriz::leerDeFichero(f);
  if (!aux) return this;

  if (means) delete means;
  means = aux->mean(true, false);
  delete aux;

  return this;
}
//---------------------------------------------------------------------------
string Mean::ComoCadena()
{
  string cad = "";
  string sep = "";
  for(int i=0;i<means->columnas();i++) {
    cad = cad + sep + Mandato::ACad((*means)[0][i]);
    sep = "\t";
  }
  return cad;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void StartRand::CreateParams()
{
  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("Normal");
  CadenasValidas.push_back("Normal guardando en fichero");
  CadenasValidas.push_back("Recuperardo los numeros aleatorios del fichero");
  ParametroCadena *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(false, 0, "Numeros aleatorios",
                               "Metodo de generacion de números aleatorios");
  Parametros.push_back(p2);

  Parametro *p = new ParametroEntero();
  p->PonPropiedades(true, 0, "Semilla de inic.", 
                                 "Solo valido para los dos primeros metodos");
  p->SetPorOmision("0");
  Parametros.push_back(p);
}
//--------------------------------------------------------------------------
void StartRand::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new StartRand(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* StartRand::Ejecutar()
{
  int modo, semilla;
  FILE *f;
  
  Param(0)->Ejecutar();
  modo = Param(0)->ComoEntero();


  semilla=NumParamsAsignados()>1?Param(1)->Ejecutar()->ComoEntero():time(0);

printf("--------\nPeti=%d\n----------\n", TDebugRand::Peticiones);
  switch (modo) {
    case 0:
      if (semilla>=0) srand(semilla);
      TDebugRand::SetModo(moNormal);
      TDebugRand::Peticiones = 0;
      break;
    case 1:
      if (semilla>=0) srand(semilla);
      f=fopen("semillas.txt", "wt");
      fprintf(f, "%d\n", semilla);
      fclose(f);
      TDebugRand::SetModo(moDeRandYGuardaEnFichero);
      TDebugRand::Start();
      TDebugRand::Peticiones = 0;
      break;
    case 2:
      TDebugRand::SetModo(moDeFichero);
      TDebugRand::Start();
      TDebugRand::Peticiones = 0;
      break;
    default:
      break;
  }

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Print::CreateParams()
{
  ParametroCadena *p2 = new ParametroCadena();
  p2->PonPropiedades(false, 0, "Texto a imprimir"); 
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
Mandato* Print::Ejecutar()
{
  txt = Param(0)->Ejecutar()->ComoCadena();
  printf("%s\n", txt.c_str());
  return this;
}
//---------------------------------------------------------------------------
void Print::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Print());
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Call::CreateParams()
{
  string fil = string("Script MKO|*.mKo");
  Parametro *pfe = new ParametroFicheroEntrada(fil);
  pfe->PonPropiedades(false, 0, "Fichero script",
                               "Nombre del fichero con el script a ejecutar");
  Parametros.push_back(pfe);
}
//---------------------------------------------------------------------------
void Call::LoadScript()
{
  vector<string> params;
  string nf;
  string kk;
  char fnd[24];
  
  //Leemos los parametros del script
  nf = Param(0)->Ejecutar()->ComoCadena();
  for (int i=1;i<NumParamsAsignados();i++) {
    params.push_back(Param(i)->Ejecutar()->ComoCadena());
  }

  //Leemos y compilamos el script
  ifstream ifs(nf.c_str(), ios_base::in);
  while(getline(ifs, kk)) {
    int ip=NumParamsAsignados() - 2;
    for (int i = NumParamsAsignados()-1; i >= 1; i--) {
      sprintf(fnd, "$%d", ip--);
      replace_all (fnd, params[i-1], kk);
    }
    txt.push_back(kk);
    
    if (0<=(ip=kk.find("//"))) {//Comentario
      kk = kk.substr(0, ip);
    }
    script.push_back(Mandato::valueOf(kk));
    if (script[script.size()-1]!=0 && 
                              typeid(*script[script.size()-1])==typeid(*this))
      ((Call*)script[script.size()-1])->LoadScript(); 
  };
  ifs.close();

}
//---------------------------------------------------------------------------
bool Call::ScriptLoaded()
{
	return script.size()>0;
}
//---------------------------------------------------------------------------
Mandato* Call::Ejecutar()
{

  if(!ScriptLoaded()) LoadScript();
  //Ejecutamos el script
  for(unsigned i=0;i<script.size();i++) {
    cout << i <<  ": " << txt[i] << endl;
    if (script[i]) {
	    script[i]->Ejecutar();
	    cout << script[i]->ComoCadena() << endl;
    }
  }

  //Liberamos memoria
/*  for(unsigned i=0;i<script.size();i++)
    delete script[i];
  script.clear();
  */
  return this;
}
//---------------------------------------------------------------------------
void Call::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Call(), "File");
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Time::CreateParams()
{
  Parametro *p = new Parametro(false, Tipo::TNinguno(), "Mandato a ejecutar");
  Parametros.push_back(p);

  p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(true, 0, "Numero de veces que se ejecuta", 
                     "Mas ejecuciones producen mayor precision en la medida");
  p->SetPorOmision("1");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* Time::Ejecutar()
{
  int rep = 1;

  if (NumParamsAsignados()>1) {
    rep = Param(1)->Ejecutar()->ComoEntero();
  }

#if defined(__STRICT_ANSI__) or defined(WIN32)
  double antes = clock();
#else
  double antes = times(0);
#endif
  double antes2 = time(0);
  for (int i=0;i<rep;i++) {
    Param(0)->Ejecutar();
//    cout << i << ": " << Param(0)->ComoCadena() << endl;
  }
#if defined(__STRICT_ANSI__) or defined(WIN32)
  double despues = clock();
  double div = rep*CLOCKS_PER_SEC;
#else
  double despues = times(0);
  double div = rep*sysconf(_SC_CLK_TCK);
#endif
  double despues2 = time(0);
  printf("TIEMPO ABSOLUTO: %g seg\n", (despues2-antes2)/rep);
  
  secs = (despues-antes)/div;

  return this;
}
//---------------------------------------------------------------------------
void Time::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Time(), "Tools");
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Repeat::CreateParams()
{
  Parametro *p = new Parametro(false, Tipo::TNinguno(), "Mandato a ejecutar");
  Parametros.push_back(p);

  p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Numero de vueltas del bucle");
  p->SetPorOmision("1");
  Parametros.push_back(p);

  p = new Parametro(true, Tipo::TNinguno(), "Ejecutar al final del bucle");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* Repeat::Ejecutar()
{
  int rep = Param(1)->Ejecutar()->ComoEntero();
  Mandato *mnt=0;

  if (NumParamsAsignados()>2)
    mnt = Param(2);

  for (int i=0;i<rep;i++) {
    Param(0)->Ejecutar();
    Mandato *output = Variable::VariablePorNombre("OUTPUT");
    if (!output || output->ComoBooleano())
      cout << i << ": " << Param(0)->ComoCadena() << endl;
    if (mnt) mnt->Ejecutar();
  }

  return this;
}
//---------------------------------------------------------------------------
void Repeat::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Repeat(), "Tools");
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void If::CreateParams()
{
  Parametro *p = new Parametro(false, Tipo::TNinguno(), "Condicion");
  Parametros.push_back(p);

  p = new Parametro(false, Tipo::TNinguno(), "A ejecutar si verdadero");
  Parametros.push_back(p);

  p = new Parametro(false, Tipo::TNinguno(), "A ejecutar si falso");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* If::Ejecutar()
{
  bool cond = Param(0)->Ejecutar()->ComoBooleano();

  mand = cond ? Param(1)->Ejecutar() : Param(2)->Ejecutar();

  return this;
}
//---------------------------------------------------------------------------
void If::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new If());
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

