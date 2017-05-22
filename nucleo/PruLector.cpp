//---------------------------------------------------------------------------

#include <iostream>
#include "Lector.h"
#include "FnsOrdClas.h"
#include "FnsClasificacion.h"
#include <fstream>
#include <string.h>
#include "FnsAux.h"
#include "Arcing.h"

#ifndef VERSION
  #define VERSION "Unknown"
#endif

#ifndef FLAGS
  #define FLAGS "Unknown"
#endif

//---------------------------------------------------------------------------
using namespace std;
//---------------------------------------------------------------------------
void imprimir(vector<string> lista);
void ImprimirMandatos();
void ImprimirVars(string Namespace);
void RemplazarMandatos(string &kk);
void AyudaMandato(string txt);
Mandato *Run(string kk);
//---------------------------------------------------------------------------
class Shell : public IntefazFuncion, Funcion
{
  protected:
    int ret;

  public:
    Shell(std::vector<Mandato *> mands) : Funcion(mands){CreateParams();ret=0;}
    Shell() : Funcion() {CreateParams();ret=0;}

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 1;}
    virtual Funcion* CrearFuncion() {return new Shell();}
    virtual std::string NombreFuncion() {return "Shell";}

    virtual Mandato* Ejecutar();

    virtual bool EsEntero(){return true;}

    virtual std::string ComoCadena() { return Mandato::ACad(ret);}
    virtual double ComoNumero(){return ret;}
    virtual int ComoEntero(){return ret;}
    virtual bool ComoBooleano(){return ret==0;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Shell::CreateParams()
{
  Parametro *pc = new ParametroCadena();
  Parametros.push_back(pc);
}
//--------------------------------------------------------------------------
void Shell::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Shell(), "");
}
//---------------------------------------------------------------------------
Mandato* Shell::Ejecutar()
{
  Param(0)->Ejecutar();

  ret = system(Param(0)->ComoCadena().c_str());

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    bool seguir=true;
    vector<string> mnds;

    if (argc>1) {
      if ( strcmp(argv[1], "-v") == 0 ) {
        printf("Executable compiled:\n\t\t%s\n", VERSION);
        printf("Compilation flags:\n\t\t%s\n", FLAGS);
        exit(0);
      }
    }

    mnds.push_back("Call(\"./test.mKo\")");
    
    DarDeAltaFnsOrdClas();
    DarDeAltaFnsClasificacion();
    DarDeAltaFnsAux();
    Shell::DarDeAlta();

    Mandato *OUTPUT = Mandato::valueOf("OUTPUT=true");
	    
    string kk;
    FILE *f;

    try {
      if (argc>1) {
        string kk = "Call(\"";
        kk = kk + argv[1] + "\"";
        for(int i=2;i<argc;i++) kk = kk + ", \"" + argv[i] + "\"";
        kk = kk + ")";
        printf("%s\n", kk.c_str());
        Mandato *n = Mandato::valueOf(kk);
        if (n!=0) n->Ejecutar();
        exit(0);
      }
    
      f=fopen("last.mKo", "wt");
      fclose(f);

      printf("Escriba una expresiOn:\n");
      while (seguir) {
        printf("> ");
        if ( !getline(cin, kk) ) return 0;
        if (kk.size()==0) continue;
        if (kk=="#remove_all") {
          Mandato::EliminarMandatos();
          continue;
        }
        else if (kk.compare(0,2,"##")==0 || kk.compare(0,4,"help")==0) {
          if (kk=="##" || kk=="help") //Listamos los mandatos dados de alta
            ImprimirMandatos();
          else if (kk[2]=='!') //Listamos los mandatos dados de alta
            ImprimirVars(kk.substr(3));
          else 
            AyudaMandato(kk.substr(kk.compare(0,2,"##")==0?3:5));
          continue;
        }
        else if (kk[0]=='#') {
          if (kk.size()==1) {
            imprimir(mnds);
            continue;
          }
          else {
            unsigned ii = atoi(kk.c_str()+1);
            if (ii>mnds.size()) continue;
            kk = mnds[mnds.size()-ii];
            cout << " " << kk << endl;
          }
        }
        if (kk[0]=='!') {
          system(kk.c_str()+1);
        }
        else if (kk == "fin" || kk == "exit") seguir = false;
        else {
          try {
            RemplazarMandatos(kk);
#ifdef _DEBUG
            cout << "Llamando a valueOf(" << kk << ")" << endl;
            fflush(stdout);
#endif
            Mandato *n = Mandato::valueOf(kk);
            if (n!=0) {
              //Solo se anade a la lista de mandatos si es distinto
              //al ultimo mandato
              if (*(mnds.end()-1)!=kk) {
                mnds.push_back(kk);
                FILE *f=fopen("last.mKo", "at");
                fprintf(f, "%s\n", kk.c_str());
                fclose(f);
              }
              n->Ejecutar();
              if (OUTPUT->ComoBooleano()) cout << n->ComoCadena() << endl;
            }
            else {
              cout << "ERROR en: " << kk;
            }
          }
          catch(char *msg){
            cout << "ERROR 0: " << msg << endl;
          }
          catch(MandatoError &err){
            cout << "ERROR a: " << err.que();
          }
          catch(exception &exp){
            cout << "ERROR b: " << exp.what();
          }
        }
      }
    }
    catch(char *msg){
      cout << "ERROR c: " << msg << endl;
    }
    catch(exception &exp){
      cout << "ERROR d: " << exp.what();
    }
  // catch(...){
    //     }
    
    return 0;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void imprimir(vector<string> lista)
{
  for (unsigned i=0;i<lista.size();i++)
    cout << (lista.size()-i) << ": " << lista[i] << endl;
}
//---------------------------------------------------------------------------
void AyudaMandato(string txt)
{
  Mandato::Trim(&txt);

  Funcion *f = ifns().InfoFuncionPorNombre(txt);
  if (!f) {
    cout << txt << ": Nombre de funcion no encotrada" << endl;
    return;
  }

  cout << txt << "[" << f->TipoRetorno()->nombre()<< "]. Parametros: " << endl;
  for(unsigned k=0;k<f->NumParams();k++) {
    Parametro *p = f->ParamInfo(k);
    cout << "\t-" << p->Nombre() << "[" << p->DameTipo()->nombre()<< "]: ";
    cout << p->Descripcion() << (p->Opcional() ? " [Opcional]":"") << endl;
    ParametroCadena *pc = dynamic_cast<ParametroCadena *>(p);
    if (pc && pc->NumCadenasValidas()>0) {
      cout << "\t  " << "Posibles valores:" << endl;
      for(int j=0;j<pc->NumCadenasValidas();j++) 
        cout << "\t\t" << (j+1) << " = " << pc->CadenaValida(j) << endl;
    }
  }


/*  for(int i=0;i<ifns().size();i++){
    cout << "##" << (i+1) << ": ";
    cout << ifns().iFuncion(i)->NombreFuncion() 
      << "\t\t"*/ /*<< ifns().InfoFuncion(i)->NumParams()*/ /*<< endl;
    }*/
}
//---------------------------------------------------------------------------
void ImprimirMandatos()
{
  for(int i=0;i<ifns().NumCategorias();i++){
    string cat = ifns().DameCategoria(i);
    vector<int> inds = ifns().IndicesFuncionesPorCat(cat);
    cout << "Categoria: " << cat << endl;
    for(unsigned j=0;j<inds.size();j++) {
//      Funcion *f = ifns().InfoFuncion(inds[j]);
      IntefazFuncion *fi = ifns().iFuncion(inds[j]);
//      cout << "\t" << f->TipoRetorno()->nombre() << " ";
      cout << "  " << fi->NombreFuncion() << "(";
/*      string sep = ", ";
      for(int k=0;k<f->NumParams();k++) {
        if (k==f->NumParams()-1) sep = "";
        cout << f->ParamInfo(k)->DameTipo()->nombre() << sep;
      }*/
      cout << "...)" << endl;
    }
  }


/*  for(int i=0;i<ifns().size();i++){
    cout << "##" << (i+1) << ": ";
    cout << ifns().iFuncion(i)->NombreFuncion() 
      << "\t\t"*/ /*<< ifns().InfoFuncion(i)->NumParams()*/ /*<< endl;
    }*/
}
//---------------------------------------------------------------------------
void RemplazarMandatos(string &kk)
{
  char fnd[24];
  for(int i=0;i<ifns().size();i++){
    sprintf(fnd, "##%d", i+1);
    replace_all(fnd, ifns().iFuncion(i)->NombreFuncion(), kk);
  }
}
//---------------------------------------------------------------------------
void ImprimirVars(string Namespace)
{
  Variable::print();
/*  for(int i=0;i<Variable::NumVariables(Namespace);i++){
    Variable *var = Variable::VariablePorNumero(i, Namespace);
    cout << var->NombreCorto() <<"\t["<< var->DameTipo()->nombre()
                                     << "]: \t"  << var->ComoCadena() << endl;
  }*/
}
//---------------------------------------------------------------------------
Mandato *Run(string kk)
{
  return 0;
}
//---------------------------------------------------------------------------

