//---------------------------------------------------------------------------


#include "FnsClasificacion.h"
#include "data20.h"
#include "node20.h"

#include "stats.h"
#undef Min
#undef Max

#include "Arcing.h"
#include "GPTree.h"
#include "NNet.h"
#include "Matriz.h"
#include "C45Tree.h"
#include "UtilsGraf.h"
#include "Utils.h"
#include "SpaceCubes.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <iterator>
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
using namespace std;
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
TipoData* TipoData::_TData = new TipoData();
TipoClasificador* TipoClasificador::_TClasificador = new TipoClasificador();
TipoEnsemble* TipoEnsemble::_TEnsemble = new TipoEnsemble();
//---------------------------------------------------------------------------
Variable *TipoEnsemble::NuevaVariable(std::string nombre, Mandato *inic)
{
  VariableEnsemble *var = new VariableEnsemble(nombre);
  if (inic) var->Set(inic);
  return var;
}
//---------------------------------------------------------------------------
VariableEnsemble::~VariableEnsemble()
{
  clear();
}
void VariableEnsemble::clear()
{
  for(int i=0;i<NumElementos();i++) {
    delete clsfs[i];
  }
  clsfs.clear();
  pos.clear();
}
void VariableEnsemble::init(Ensemble *ens)
{
  clear();
  if (!ens) return;

  Variable *v;
  VariableDatos *c = new VariableDatos();

  for (int i=0;i<ens->Count();i++) {
    std::ostringstream buf;
    buf << NombreCorto() << ".";
    buf << TipoClasificador::TClasificador()->nombre() << "_" << i;
    c->SetDatos(ens->GetClassifier(i));
    v = TipoClasificador::TClasificador()-> NuevaVariable(buf.str(), c);
    pos[ens->GetClassifier(i)] = i;
    clsfs.push_back(v);
  } 

  delete c;
}
int VariableEnsemble::NumElementos()
{
  if (!FValor) return 0;
  return clsfs.size();
}
Variable* VariableEnsemble::Elemento(int i)
{
  return clsfs[pos[((Ensemble*)FValor)->GetClassifier(i)]];
}
//---------------------------------------------------------------------------
void DarDeAltaFnsClasificacion()
{
  CargaClasificador::DarDeAlta();
  CargaDatos::DarDeAlta();
  CambiaColumnaDependiente::DarDeAlta();
  CopiarColumnaDependiente::DarDeAlta();
  FiltraDatos::DarDeAlta();
  Clasificar::DarDeAlta();
  ClasErr::DarDeAlta();
  ErrSecClas::DarDeAlta();
  MatrizClasif::DarDeAlta();
  MatrizConfusion::DarDeAlta();
  Margen::DarDeAlta();
  Certeza::DarDeAlta();
  SumaMargen::DarDeAlta();
  EstimErrorVal::DarDeAlta();
  EstimErrorTest::DarDeAlta();
  DiversityMeasures::DarDeAlta();
  EnsembleMeasures::DarDeAlta();
  SetPropertyValue::DarDeAlta();
  SetNTrain::DarDeAlta();
  CreaConjuntosTrainAleat::DarDeAlta();
  GuardaDatos::DarDeAlta();
  GuardaClasificador::DarDeAlta();
  InfoClasificador::DarDeAlta();
  ConstruyeBaggingCART::DarDeAlta();
  ConstruyeBaggingC45::DarDeAlta();
  ConstruyeBaggingNNet::DarDeAlta();
  ConstruyeRandomForest::DarDeAlta();
//  ConstruyeBaggingComMod::DarDeAlta();
//  ConstruyeBoostingWithBaggingInfo::DarDeAlta();
  ConstruyeBoosting::DarDeAlta();
  ConstruyeBoostingC45::DarDeAlta();
  ConstruyeClassSwitching::DarDeAlta();
  ConstruyeCART::DarDeAlta();
  ConstruyeC45::DarDeAlta();
  ConstruyeNNet::DarDeAlta();
  MapaDeClasif2D::DarDeAlta();
  TestAll::DarDeAlta();
}
//---------------------------------------------------------------------------
void DoCreaConjuntosTrainAleat(Data *data, int numdatos,
                                      int numconjuntos, string salida, int opts)
{
  bool stratify = opts==1 ? true : false;
  bool eqdistrib = opts==2 ? true : false;
  char salts[300], saltr[300];

  //Barajamos un poco
  data->Scramble(0, data->GetNTotal()-1);
  data->Scramble(0, data->GetNTotal()-1);
  data->Scramble(0, data->GetNTotal()-1);
  data->Scramble(0, data->GetNTotal()-1);

  //Si numdatos<0 hacemos numconjuntos-folds para CV
  if (numdatos<=0) {
    data->SetNTrain(data->GetNTotal());
    data->SetNGrow(data->GetNTrain());
    data->ResetCV(numconjuntos);
    for (int i=0;i<numconjuntos;i++) {
      sprintf(saltr, "%strain_cv%d.cre", salida.c_str(), i);
      sprintf(salts, "%stest_cv%d.cre", salida.c_str(), i);
      data->SetCV(i+1);
      data->SaveToFile(saltr, 0, data->GetNumData()-1);
      data->SaveToFile(salts, data->GetNumData(), data->GetNTotal()-1);
    }
  }
  else {
    data->SetNTrain(numdatos);
    for (int i=0;i<numconjuntos;i++) {
      sprintf(saltr, "%strain_%d.cre", salida.c_str(), i);
      sprintf(salts, "%stest_%d.cre", salida.c_str(), i);
      data->Scramble(0, data->GetTotalData()-1);
      if (stratify) data->SelectNewTrainingSet();
      else if (eqdistrib) data->SelectNewTrainingSetEqDis();
      data->SaveToFile(saltr, 0, numdatos-1);
      data->SaveToFile(salts, numdatos, data->GetTotalData()-1);
    }
  }
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void CargaClasificador::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new CargaClasificador(), "File");
}
//---------------------------------------------------------------------------
void CargaClasificador::CreateParams()
{
  string fil = string("Todos|*.*|") + "Formato rapido (*.clf)|*.clf";
  Parametro *pfe = new ParametroFicheroEntrada(fil);
  pfe->PonPropiedades(false, 0, "Fichero", 
                                    "Nombre del fichero con el clasificador");
  Parametros.push_back(pfe);

  Parametro *p = new ParametroData();
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Tipo* CargaClasificador::DameTipo()
{
  if (dynamic_cast<Ensemble*>(Clasif))
    return TipoEnsemble::TEnsemble(); 

  return TipoClasificador::TClasificador();
}
//---------------------------------------------------------------------------
Mandato* CargaClasificador::Ejecutar()
{
  Param(0)->Ejecutar();

  //Leemos el clasificador
  ifstream ifs(Param(0)->ComoCadena().c_str(), ios_base::in | ios_base::binary);
  if (Clasif) delete Clasif;
  Clasif = Classifier::LeerClasificadorGeneral(ifs);
  ifs.close();

  if (Clasif && ParamInfo(1)->Asignado()) {
    Param(1)->Ejecutar();
    Clasif->SetData((Data*)Param(1)->ComoDatos());
//    Clasif->ResetOriginalTrainingData();
  }


  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
CargaDatos::~CargaDatos()
{
  if(nd) delete nd;
}
//---------------------------------------------------------------------------
void CargaDatos::CreateParams()
{
  string fil = string("Todos|*.cre;*.asc;*.names;|") + "Propio (*.cre)|*.cre|" +
                              "C45 (*.names)|*.names|" + "Raetsch(*.asc)|*.asc";
  Parametro *pfe = new ParametroFicheroEntrada(fil);
  pfe->PonPropiedades(false, 0, "Fichero", "Nombre del fichero con la base de datos");
  Parametros.push_back(pfe);
}
//---------------------------------------------------------------------------
void CargaDatos::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new CargaDatos(), "File");
}
//---------------------------------------------------------------------------
Mandato* CargaDatos::Ejecutar()
{
  Param(0)->Ejecutar();

  if (nd) delete nd;

  nd = Data::DataFromFile((char*)Param(0)->ComoCadena().c_str());
  nd->SetNTrain(nd->GetNTotal());

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
CambiaColumnaDependiente::~CambiaColumnaDependiente()
{
  if(nd) delete nd;
}
//---------------------------------------------------------------------------
void CambiaColumnaDependiente::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                             "Nombre de la base de datos con los ejemplos"));
  
  Parametro *p = new ParametroEntero();
  p->PonPropiedades(false, 0, "Column index for dependent variable",
                              "starts in 0");
  Parametros.push_back(p);

}
//---------------------------------------------------------------------------
void CambiaColumnaDependiente::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new CambiaColumnaDependiente(), "File");
}
//---------------------------------------------------------------------------
Mandato* CambiaColumnaDependiente::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nueva = Param(1)->ComoEntero();

  nd->ChangeDependentColumn(nueva);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
CopiarColumnaDependiente::~CopiarColumnaDependiente()
{
  if(nd) delete nd;
}
//---------------------------------------------------------------------------
void CopiarColumnaDependiente::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                             "Nombre de la base de datos con los ejemplos"));

}
//---------------------------------------------------------------------------
void CopiarColumnaDependiente::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new CopiarColumnaDependiente(), "File");
}
//---------------------------------------------------------------------------
Mandato* CopiarColumnaDependiente::Ejecutar()
{
  Param(0)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  nd->CopiarColumnaDependiente();

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
FiltraDatos::~FiltraDatos()
{
  if(nd) delete nd;
}
//---------------------------------------------------------------------------
void FiltraDatos::CreateParams()
{
  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("Dispersar");
  CadenasValidas.push_back("NaiveBayes");
  CadenasValidas.push_back("Clonar");
  CadenasValidas.push_back("NN");
  CadenasValidas.push_back("Margen");
  CadenasValidas.push_back("OOB");
  CadenasValidas.push_back("Wrong");
  CadenasValidas.push_back("Ok");
  CadenasValidas.push_back("Border");
  CadenasValidas.push_back("ClassNoise");
  Parametro *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(false, 0, "Filtro", "Filtro a aplicar");
  Parametros.push_back(p2);

  p2 = new ParametroData();
  p2->PonPropiedades(false, 0, "Datos a filtros", 
                                    "Datos sobre los que se aplica el filtro");
  Parametros.push_back(p2);

  p2 = new ParametroCadena();
  p2->PonPropiedades(true, 0, "Parametro", "Parametro del filtro (si aplica)");
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
void FiltraDatos::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new FiltraDatos(), "Tools");
}
//---------------------------------------------------------------------------
Matriz* ClasificationMatrix(Ensemble *ens, NomData *data, bool ExcluirTrain)
{
  int nDatos = data->GetNTotal();
  int nClasf = ens->GetClassifiersToUse();
  int M = nClasf;
  int N = nDatos;

  Matriz *Cmn = new Matriz(M+1, N);
  for(int m=0;m<M;m++) {
    Classifier *c = ens->GetClassifier(m);
    for(int n=0;n<N;n++) {
      int clase = c->Classify(n);
      bool Val = ExcluirTrain &&
                    c->UsedInOriginalTrainingData(data->GetDatIniPos(n));
      (*Cmn)[m][n] = Val ? -1 : clase;
    }
  }
  for(int n=0;n<N;n++) {
    vector<double> votos(data->NumClass, 0.0);
    for(int m=0;m<M;m++) {
      if ((*Cmn)[m][n]>=0) {
        votos[(int)(*Cmn)[m][n]]++;
      }
    }
    (*Cmn)[M][n] = Classifier::WhichClass(votos);
  }

  return Cmn;
}
Mandato* FiltraDatos::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  NomData *data = (NomData*)Param(1)->ComoDatos();

  if (nd) delete nd;

  int what = Param(0)->ComoEntero();
  if (what==1) {
    int mg = Param(2)->Ejecutar()->ComoEntero();
    nd = data->GenerateSinteticData(mg, 0, data->GetNTotal()-1);
  }
  else if (what==2) {
    int ntot = Param(2)->Ejecutar()->ComoEntero();
    nd = data->SinteticDataNaiveBayes(ntot, 0, data->GetNTotal()-1);
  }
  else if (what==3) {
    nd = data->Clone();
  }
  else if (what==4) {
    int ntot = Param(2)->Ejecutar()->ComoEntero();
    bool ByClass = NumParamsAsignados()>3 ? 
                                   Param(3)->Ejecutar()->ComoBooleano() : true;
    nd = data->SinteticDataNN(ntot, 0, data->GetNTotal()-1, ByClass);
  }
  else if (what==5) {
    Classifier *cls = (Classifier*)Param(2)->Ejecutar()->ComoDatos();
    int ncls = ((NomData*)data)->NumClass;
    nd = data->Clone(0, data->GetNTotal()-1, ncls);
    for(int i=0;i<data->GetNTotal();i++) {
      vector<double> votos = cls->Distribution(i);
      int clase = nd->GetDatClass(i);
      nd->SetDatWeight(i, votos[clase]);
      for(int j=1;j<ncls;j++) {
        clase++;
        clase %= ((NomData*)data)->NumClass;
        nd->SetDatClass(j*data->GetNTotal() +i, clase);
        nd->SetDatWeight(j*data->GetNTotal()+i, votos[clase]);
      }
    } 
  }
  else if (what==6) {
    Ensemble *ens = (Ensemble*)Param(2)->Ejecutar()->ComoDatos();

    int ncls = ((NomData*)data)->NumClass;
    //Se ordena por orden inicial
    data->SortOn(data->GetNumVar()+data->GetIniPosIndex());
    Matriz *clases = ClasificationMatrix(ens, data, true);

    nd = data->Clone(0, data->GetNTotal()-1, ncls);
    for(int i=0;i<data->GetNTotal();i++) {
      for(int j=0;j<ncls;j++) {
        nd->SetDatClass(j*data->GetNTotal() +i, j);
        nd->SetDatWeight(j*data->GetNTotal()+i, 0.0);
      }
      if ((*clases)[ens->GetClassifiersToUse()][i]!=data->GetDatClass(i)) continue;
      for(int j=0;j<ens->GetClassifiersToUse();j++) {
        if ((*clases)[j][i]>=0) {
          int clase = (int)(*clases)[j][i];
          double w = 1.0 + nd->GetDatWeight(clase*data->GetNTotal()+i);
          nd->SetDatWeight(clase*data->GetNTotal()+i, w);
        }
      }
    } 
/*    FILE *f=fopen("oob.csv", "a");
    for(int j=0;j<ncls;j++) {
      for(int i=0;i<data->GetNTotal();i++) {
        fprintf(f, "%g, ", nd->GetDatWeight(j*data->GetNTotal()+i));
      }
      fprintf(f, "\n");
    }
    fclose(f);*/
    delete clases;
    //solo multiplos enteros de los datos originales
/*    int T = ens->GetClassifiersToUse();
    int factor = (int)(0.5+(double)ntot/T);
    ntot = T*factor;
    d->redim(ntot);
    vector<int> from(data->GetNTotal());
    vector<int> to(data->GetNTotal());
    for(int i=0;i<T;i++) {
      Classifier *cls = ens->GetClassifier(i);
      vector<bool> set(data->GetNTotal(), false);
      for(int j=0;j<cls->OriginalTrainingDataCount();j++) {
        set[cls->OriginalTrainingDataPos(j)] = true;
      }
      int n_oob = 0;
      for(unsigned j=0;j<set.size();j++) {
        if (!set[j]) {
          from[n_oob++] = j;
        }
      }
      random_sample_n(from->first(), from->first()+n_oob, to, factor);
      
      for(int j=1;j<ncls;j++) {
        clase++;
        clase %= ((NomData*)data)->NumClass;
        nd->SetDatClass(j*data->GetNTotal() +i, clase);
        nd->SetDatWeight(j*data->GetNTotal()+i, votos[clase]);
      }
    } */
  }
  else if (what==7 || what==8) {
    Classifier *cls = (Classifier*)Param(2)->Ejecutar()->ComoDatos();

    int n_ok=0, n_nook=0;
    double w_ok, w_nook;
    if (what==7) {
      w_ok = 1.0;
      w_nook = 0.0;
    }
    else {
      w_ok = 0.0;
      w_nook = 1.0;
    }
    int nvar = data->GetNumVar();
    for(int i=0;i<data->GetNTotal();i++) {
      double w;
      if (cls->Classify(i)==data->GetDatClass(i)) {
        w = w_ok;
        n_ok++;
      }
      else {
        w = w_nook;
        n_nook++;
      }
      data->SetValueVar(i, nvar, w);
    }
    data->SortOn(nvar);
    nd = data->Clone(0, what==7 ? n_nook-1 : n_ok-1);
  }
  else if (what==9) {
    int knn = NumParamsAsignados()>2 ? Param(2)->Ejecutar()->ComoEntero() : 5;
    nd = data->BorderData(0, data->GetNTotal()-1, knn);
    printf("No. datos: %d\n", nd->GetNTotal());
  }
  else if (what==10) {//ClassNoise
    double prob = Param(2)->Ejecutar()->ComoNumero();
    nd = ((NomData*)data)->ClassNoise(prob);
  }

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Clasificar::~Clasificar()
{
  if(mat) delete mat;
}
//---------------------------------------------------------------------------
void Clasificar::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Clasificar(), "Classifiers");
}
//---------------------------------------------------------------------------
void Clasificar::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
Mandato* Clasificar::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Classifier *c = (Classifier *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  if (mat) delete mat;
  mat = new Matriz(1, dat->GetNTotal());

  c->SetData(dat);

  for(int i=0;i<dat->GetNTotal();i++)
    (*mat)[0][i] = c->Classify(i);

  return this;
}
//---------------------------------------------------------------------------
string Clasificar::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ClasErr::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ClasErr(), "Classifiers");
}
//---------------------------------------------------------------------------
void ClasErr::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
Mandato* ClasErr::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Classifier *c = (Classifier *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  c->SetData(dat);
  err = c->Error(0, dat->GetNTotal()-1);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
ErrSecClas::~ErrSecClas()
{
  if(err) delete err;
}
void ErrSecClas::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ErrSecClas(), "Classifiers");
}
//---------------------------------------------------------------------------
void ErrSecClas::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
  Parametro *pfs = new ParametroFicheroSalida();
  pfs->PonPropiedades(true, 0, "Fichero clasificaciones",
                "Nombre del fichero para guardar las clasificacion individual");
  Parametros.push_back(pfs);
}
//---------------------------------------------------------------------------
Mandato* ErrSecClas::Ejecutar()
{
  int ini, fin;
  string nom_fich;
  vector<double> *inderrs = new vector<double>();

  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();

  Data *dat = (Data*)Param(1)->ComoDatos();
  ens->SetData(dat);
  ini = 0;
  fin = dat->GetNTotal()-1;

  vector<int> *clas = 0;
  if (NumParamsAsignados()==3) {
    Param(2)->Ejecutar();
    nom_fich = Param(2)->ComoCadena();
    clas = new vector<int>();
  }

  vector<double> errores;
  ens->SecuencialError(ini, fin, &errores, inderrs, clas);

  if (err) delete err;
  err = new Matriz(1, (int)errores.size());
  for (unsigned i=0; i<errores.size(); i++) {
    (*err)[0][i] = errores[i];
  }

/*  //Generamos variables para su acceso externo
  double min = 1.0;
  int pos;
  for (unsigned i=0; i<err->size(); i++) {
    if (min>(*err)[i]) {
      min = (*err)[i];
      pos = i;
    }
  }

  Variable *mmin, *mpos;
  Mandato::valueOf("min=1");
  Mandato::valueOf("pos=1");
  mmin=Variable::VariablePorNombre("min");
  mpos=Variable::VariablePorNombre("pos");
  mmin->ComoNumero(min);
  mpos->ComoNumero(pos);*/

/*  FILE *f=fopen("inderrs.txt", "a");
  for (unsigned i=0; i<inderrs->size(); i++) {
    fprintf(f, "%g\t", (*inderrs)[i]);
  }
  fprintf(f, "\n");
  fclose(f);
*/
  if (clas) {
    FILE *f=fopen(nom_fich.c_str(), "a");
    for (unsigned i=0; i<clas->size(); i++) {
      fprintf(f, "%d\t", (*clas)[i]);
    }
    fprintf(f, "\n");
    fclose(f);
    delete clas;
  }

  delete inderrs;

  return this;
}
//---------------------------------------------------------------------------
string ErrSecClas::ComoCadena()
{
  string cad="", sep="";
  for(int i=0;i<err->columnas();i++) {
    cad = cad + sep + Mandato::ACad((*err)[0][i]);
    sep = "\t";
  }
  return cad;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
MatrizClasif::~MatrizClasif()
{
  if (mat) delete mat;
}
void MatrizClasif::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new MatrizClasif(), "Classifiers");
}
//---------------------------------------------------------------------------
void MatrizClasif::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
  Parametro *p = new ParametroBooleano();
  p->PonPropiedades(true, 0, "Secuencial/individual",
                             "Da la clase/error acumulado o no del conjunto");
  Parametros.push_back(p);
  p = new ParametroBooleano();
  p->PonPropiedades(true, 0, "Clase/error",
                             "Da la clase/error acumulado o no del conjunto");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* MatrizClasif::Ejecutar()
{
  int ini, fin;

  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data*)Param(1)->ComoDatos();
  ens->SetData(dat);
  ini = 0;
  fin = dat->GetNTotal()-1;

  bool Acum = NumParamsAsignados()>2 ?  Param(2)->Ejecutar()->ComoBooleano() : 
                                                                          false;
  bool Erro = NumParamsAsignados()>3 ?  Param(3)->Ejecutar()->ComoBooleano() : 
                                                                          false;

  if (mat) delete mat;
  mat = new Matriz(ens->GetClassifiersToUse(), dat->GetNTotal());
  if (Acum) {
    vector<int> classes;
    for(int j=0;j<dat->GetNTotal();j++) {
      ens->SecuencialClassify(j, &classes);
      for(int i=0;i<ens->GetClassifiersToUse();i++){
        if (Erro) {
          (*mat)[i][j] = classes[i]==dat->GetDatClass(j) ? 1 : 0;
        }
        else {
          (*mat)[i][j] = classes[i];
        }
      }
    }
  }
  else {
    int **M = ens->MatClasif0(dat);
    int **M2 = Erro ? M + ens->GetClassifiersToUse() : M;
    for(int i=0;i<ens->GetClassifiersToUse();i++){
      for(int j=0;j<dat->GetNTotal();j++)
        (*mat)[i][j] = M2[i][j];
      delete []M[i];
      delete []M[i+ens->GetClassifiersToUse()];
    }
    delete []M;
  }

  return this;
}
//---------------------------------------------------------------------------
string MatrizClasif::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
MatrizConfusion::~MatrizConfusion()
{
  if (mat) delete mat;
}
void MatrizConfusion::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new MatrizConfusion(), "Classifiers");
}
//---------------------------------------------------------------------------
void MatrizConfusion::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
Mandato* MatrizConfusion::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Classifier *c = (Classifier *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  if (mat) delete mat;
  mat = new Matriz(((NomData*)dat)->NumClass);

  c->SetData(dat);

  for(int i=0;i<dat->GetNTotal();i++)
    (*mat)[dat->GetDatClass(i)][c->Classify(i)] += 1.0;

  return this;
}
//---------------------------------------------------------------------------
string MatrizConfusion::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Margen::~Margen()
{
  if (mat) delete mat;
}
void Margen::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Margen(), "Classifiers");
}
//---------------------------------------------------------------------------
void Margen::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
  Parametro *p = new ParametroEntero(2, 2100000000);
  p->PonPropiedades(true, 0, "Divisiones", "");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* Margen::Ejecutar()
{
  int ini, fin;

  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data*)Param(1)->ComoDatos();
  ens->SetData(dat);
  ini = 0;
  fin = dat->GetNTotal()-1;

  int div=-1;
  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    div = Param(2)->ComoEntero();
  }

//  std::vector<double> margen = ens->Margen(ini, fin, div);

  if (mat) delete mat;
/*  mat = new Matriz(2, (int)margen.size());
  (*mat)[0][0] = (*mat)[1][0] = margen[0];
  for (unsigned i=1; i<margen.size(); i++) {
    (*mat)[0][i] = margen[i];
    (*mat)[1][i] = (*mat)[0][i-1] + margen[i];
  }*/
  MarginDistribution md(ens, dat, true);
  MarginDistribution md2(ens, dat, false);
  mat = new Matriz(2, div+1);
  double x = -1.0;
  for (int i=0;i<=div;i++) {
    (*mat)[0][i] = md.fm(x);
    (*mat)[1][i] = md2.fm(x);
    x+=2.0/div;
  }

  return this;
}
//---------------------------------------------------------------------------
string Margen::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Certeza::~Certeza()
{
  if (mat) delete mat;
}
void Certeza::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Certeza(), "Classifiers");
}
//---------------------------------------------------------------------------
void Certeza::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
Mandato* Certeza::Ejecutar()
{
  int ini, fin;

  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();

  if (NumParamsAsignados()==3) {
    Param(2)->Ejecutar();
    ini = Param(1)->ComoEntero();
    fin = Param(2)->ComoEntero();
  }
  else {
    Data *dat = (Data*)Param(1)->ComoDatos();
    ens->SetData(dat);
    ini = 0;
    fin = dat->GetNTotal()-1;
  }


  Data *dat = ens->GetData();
  dat->SortOn(dat->GetNumVar() + Data::IniPosIndex, ini, fin);
  if (mat) delete mat;
  mat = new Matriz(2, fin-ini+1);
  for(int i=ini; i<=fin;i++) {
    int clase;
    double cer = ens->ClassificationCertainty(i, clase/*, ExcluirTrain*/);
    if (clase!=dat->GetDatClass(i)) cer = -cer;
    (*mat)[0][i-ini] = cer;
  }

  return this;
}
//---------------------------------------------------------------------------
string Certeza::ComoCadena()
{
  string cad = "";
  string sep = "";
  for(int i=0;i<mat->columnas();i++) {
    cad = cad + sep + Mandato::ACad((*mat)[0][i]);
    sep = "\t";
  }
  return cad;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
SumaMargen::~SumaMargen()
{
  if (mat) delete mat;
}
void SumaMargen::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new SumaMargen(), "Classifiers");
}
//---------------------------------------------------------------------------
void SumaMargen::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
  Parametro *p = new ParametroNumero(-1.0, 1.0);
//  p->SetPorOmision("0.0");
  p->PonPropiedades(true, 0, "Divisiones", "");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* SumaMargen::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data*)Param(1)->ComoDatos();
  ens->SetData(dat);

  double f = 0.0;
  if (ParamInfo(2)->Asignado()) f=Param(2)->Ejecutar()->ComoNumero();

  int nClasf = ens->GetClassifiersToUse();
  int nDatos = dat->GetNTotal();
  int nClases = ((NomData*)dat)->NumClass;
  vector<double> w = ens->GetWeights();

  Matriz *mrg_sum = new Matriz(3, nClasf);
  Matriz cTrain(nDatos, nClases + 1, 0.0);
  for(int i=0;i<nClasf;i++) {
    Classifier *c = ens->GetClassifier(i);
    (*mrg_sum)[0][i] = 0.0;
    (*mrg_sum)[1][i] = 0.0;
    (*mrg_sum)[2][i] = 0.0;
    double num = 0.0;
    double num10 = 0.0;
    double num02 = 0.0;
    for(int j=0;j<nDatos;j++) {
      int clase_ij = c->Classify(j);
      cTrain[j][clase_ij] += w[i];
      //int clase = Classifier::WhichClass(cVal[j], nClases);
      double mrg_ij = 2.0*cTrain[j][dat->GetDatClass(j)]/(i+1)-1.0;
      if (mrg_ij>f) {
        (*mrg_sum)[0][i] += 1.0-mrg_ij;
        num++;
      }
      if (mrg_ij>0.2) {
        (*mrg_sum)[1][i] += 1.0-mrg_ij;
        num02++;
      }
      (*mrg_sum)[2][i] += 1.0-mrg_ij;
      num10++;
    }
    (*mrg_sum)[0][i] /= num;
    (*mrg_sum)[1][i] /= num02;
    (*mrg_sum)[2][i] /= num10;
  }

  if (mat) delete mat;
  mat = mrg_sum;

  return this;
}
//---------------------------------------------------------------------------
string SumaMargen::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
EstimErrorVal::~EstimErrorVal()
{
  if (mat) delete mat;
}
void EstimErrorVal::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new EstimErrorVal(), "Classifiers");
}
//---------------------------------------------------------------------------
void EstimErrorVal::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
Mandato* EstimErrorVal::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data*)Param(1)->ComoDatos();
  ens->SetData(dat);

  int nClasf = ens->GetClassifiersToUse();
  int nDatos = dat->GetNTotal();
  int nClases = ((NomData*)dat)->NumClass;

  if (mat) delete mat;
  mat = new Matriz(3, nClasf);
  Matriz cVal(nDatos, nClases + 1, 0.0);
  Matriz cTrain(nDatos, nClases + 1, 0.0);
  int winners[1024];
  for(int i=0;i<nClasf;i++) {
    Classifier *c = ens->GetClassifier(i);
    double err = 0.0, denom = 0.0, draws = 0.0;
    for(int j=0;j<nDatos;j++) {
      int clase_ij = c->Classify(j);
      cTrain[j][clase_ij] += 1.0;
      if (!c->UsedInOriginalTrainingData(dat->GetDatIniPos(j))) {
        cVal[j][clase_ij] += 1.0;
        cVal[j][nClases] += 1.0;
      }
      int clase;
      if (cVal[j][nClases]==0.0) { //Si nadie se ha pronunciado pasamos de el
        continue;
      } 
      clase = Classifier::WhichClass(cVal[j], nClases, winners);
      if (winners[0]>1) {
        //Si hay empate se.. 
        int k;
        for(k=0;k<winners[0];k++) {
          if (winners[k+1]==dat->GetDatClass(j))
            break;
        }
        if (k==winners[0]) {//La clase correcta no estA entre las elegidas
          err++;
          denom++;
        }
        else {//Si que estA
          draws += 1.0;
        }
      }
      else if (clase!=dat->GetDatClass(j)) {
        err++;
      }
      if (winners[0]==1) denom++;
    }

    (*mat)[0][i] = err/denom;
    (*mat)[1][i] = (err+draws)/(denom+draws);
    (*mat)[2][i] = err/(denom+draws);
  }

  return this;
}
//---------------------------------------------------------------------------
string EstimErrorVal::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
EstimErrorTest::~EstimErrorTest()
{
  if (mat) delete mat;
}
void EstimErrorTest::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new EstimErrorTest(), "Classifiers");
}
//---------------------------------------------------------------------------
void EstimErrorTest::CreateParams()
{
  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("GammaBound");
  CadenasValidas.push_back("DeltaBound");
  CadenasValidas.push_back("MarginError");
  Parametro *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(false, 0, "Bound", "Limite a aplicar");
  Parametros.push_back(p2);

  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
/*Matriz *CalculateGammaBound(Ensemble *ens, Data *dat, double gamma)
{
  ens->SetData(dat);
  int nClasf = ens->GetClassifiersToUse();
  int n = dat->GetNTotal();
  double bound = pow(n, -1.0+gamma/2.0);
  double delta, delta_incr, Pn;

  Matriz *m = new Matriz(1, nClasf);
  for(int i=0;i<nClasf;i++) {
    ens->SetClassifiersToUse(i+1);
    vector<double> mrg = ens->Margen(0, n-1);
    delta_incr = 2.0/(mrg.size()-1.0);
    delta = -1.0;
    Pn = mrg[0];
    for(int j=1;j<(int)mrg.size();j++) {
      if (pow(delta+delta_incr, gamma)*(Pn+mrg[j])>bound) {
        break;
      }
      delta += delta_incr; 
      Pn += mrg[j];
    }
    if (bound/Pn>=pow(delta+delta_incr, gamma))
      (*m)[0][i] = bound*(1.0/pow(delta+delta_incr, gamma));
    else
      (*m)[0][i] = Pn;//bound*(1.0/pow(delta, gamma));
  }

  return m;
}*/
Matriz *CalculateGammaBound(Ensemble *ens, Data *dat, double gamma)
{
  ens->SetData(dat);
  int nClasf = ens->GetClassifiersToUse();
  int n = dat->GetNTotal();
  double bound = pow(n, -1.0+gamma/2.0);
  double delta;
delta= 1.0 - gamma/2.0;
printf("%g, %g, %g\n", bound, pow(1000.0, -0.5), 1.0/pow(n, delta));

  Matriz *m = new Matriz(2, nClasf);
  for(int i=0;i<nClasf;i++) {
    ens->SetClassifiersToUse(i+1);
    MarginDistribution md1(ens, dat, true);
    MarginDistribution md2(ens, dat, false);
    // 
    double x0, x1;
    for(int j=0;j<2;j++) {
      MarginDistribution &md = j==0 ? md1 : md2;
      int ite=0;
      x0 = 0.0, x1 = 1.0;
if(j==0) {
      while(true) {
        ite++;
        if (x1-x0<0.000001) {
          delta = x0;
printf("Aqui: %g\n", x0);
          break;
        }
        delta = (x1+x0)/2.0;
        if (pow(delta, gamma)*md.fm(delta)<bound) {
          if (pow(delta, gamma)*md.fm(delta)+0.00001>bound) break;
          x0 = delta;
        }
        else {
          x1 = delta;
        }
      }
}
else {
double min = 5.1234;
for(int k=0;k<100;k++) {
delta=0.01+k/100.0;
double kl=md2.fm(delta)+bound*(1.0/pow(delta, gamma));
if (kl<min)
min=kl;
}
delta=min;
}
      (*m)[j][i] = j==0 ? bound*(1.0/pow(delta, gamma)) : delta;
    }
  }

  return m;
}
Matriz* ErrorMargen(Ensemble *ens, Data *dat, double omega)
{
  ens->SetData(dat);

  int nClasf = ens->GetClassifiersToUse();
  int nDatos = dat->GetNTotal();
  int nClases = ((NomData*)dat)->NumClass;
  vector<double> w = ens->GetWeights();

  Matriz *mrg_sum = new Matriz(1, nClasf);
  Matriz cTrain(nDatos, nClases + 1, 0.0);
  for(int i=0;i<nClasf;i++) {
    Classifier *c = ens->GetClassifier(i);
    (*mrg_sum)[0][i] = 0.0;
//    (*mrg_sum)[1][i] = -100.0;
//    (*mrg_sum)[2][i] = 0.0;
    for(int j=0;j<nDatos;j++) {
      int clase_ij = c->Classify(j);
      cTrain[j][clase_ij] += w[i];
      cTrain[j][nClases] += w[i];
      double maxmalo = 0.0;
      int clase_real = dat->GetDatClass(j);
      for(int k=0;k<nClases;k++) {
        if (k==clase_real) continue;
        if (cTrain[j][k]>maxmalo) 
          maxmalo = cTrain[j][k];
      }
      double mrg_ij = (cTrain[j][clase_real] - maxmalo)/cTrain[j][nClases];
      if (mrg_ij>=0.0 && mrg_ij<0.05)
        (*mrg_sum)[0][i] += 1.0;
/*      if (mrg_ij<-0.05)
        (*mrg_sum)[1][i] += 1.0;
      if (cTrain[j][clase_real] == maxmalo)
        (*mrg_sum)[0][i] += 1.0;*/
/*      if (mrg_ij<0.05)
        (*mrg_sum)[2][i] += 1.0;
      if (mrg_ij<0.10)
        (*mrg_sum)[3][i] += 1.0;
      if (mrg_ij<0.15)
        (*mrg_sum)[4][i] += 1.0;*/
    }
 //   (*mrg_sum)[0][i] /= nDatos;
 //   (*mrg_sum)[1][i] /= nDatos;
//    (*mrg_sum)[2][i] /= nDatos;
  }

  return mrg_sum;
}
Matriz *CalculateDeltaBound(Ensemble *ens, Data *dat)
{
  return 0;
}
Mandato* EstimErrorTest::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();

  string what = Param(0)->ComoCadena();
  Ensemble *ens = (Ensemble *)Param(1)->ComoDatos();
  Data *dat = (Data*)Param(2)->ComoDatos();

  if (mat) delete mat;
  mat = 0;
  if (what=="GammaBound") 
    mat = CalculateGammaBound(ens, dat, Param(3)->Ejecutar()->ComoNumero());
  else if (what=="DeltaBound")
    mat = CalculateDeltaBound(ens, dat);
  else if (what=="MarginError")
    mat = ErrorMargen(ens, dat, 0.0);

  return this;
}
//---------------------------------------------------------------------------
string EstimErrorTest::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
DiversityMeasures::~DiversityMeasures()
{
  if(div) delete div;
}
void DiversityMeasures::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new DiversityMeasures(), "Classifiers");
}
//---------------------------------------------------------------------------
void DiversityMeasures::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());

  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("Todos");
  CadenasValidas.push_back("DiferenciasError");
  CadenasValidas.push_back("Secuencial Kappa");
  CadenasValidas.push_back("N00-N01-N10-N11");
  CadenasValidas.push_back("Secuencial Q[Conj(n),Conj(n+1)]");
  CadenasValidas.push_back("Secuencial Kappa [Conj(n),1-Conj(n)]");
  CadenasValidas.push_back("MelvilleMooney");
  CadenasValidas.push_back("KW");
  ParametroCadena *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(false, 0, "Estadistico", "Medida estadistica a utilizar");
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
Mandato* DiversityMeasures::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();

  Ensemble *ens = (Ensemble*)Param(0)->ComoDatos();
  Data* dt = (Data*)Param(1)->ComoDatos();
  string what = Param(2)->ComoCadena();
  int what_i = Param(2)->ComoEntero();
    
  int nClases = ((NomData*)dt)->NumClass;
  int nDatos = dt->GetNTotal();
  int nClasf = ens->GetClassifiersToUse();

  int **Mcd = ens->MatClasif0(dt);

  int **Res = new int*[nClasf];
  for(int i=0;i<nClasf;i++) 
    Res[i] = new int[nDatos];

  int *clase_ok = new int[nDatos];
  for(int i=0;i<nDatos;i++) 
    clase_ok[i] = dt->GetDatClass(i);

  for(int j=0;j<nDatos;j++) {
    vector<double> votes;
    votes.resize(nClases, 0.0);
    for(int i=0;i<nClasf;i++) {
      votes[Mcd[i][j]] += 1.0;
      Mcd[i+nClasf][j] = Classifier::WhichClass(votes);
    }
    for(int i=0;i<nClasf;i++) {
      Res[i][j] = Classifier::WhichClass(votes);
      votes[Mcd[i][j]] -= 1.0;
    }
  }

  if (div) delete div;
  div = new Matriz(10, nClasf);
  int **CInd = Mcd;
  int **CEns = Mcd+nClasf;
  int f=0, c;
  for(c=0;c<nClasf-2;c++) {
    f=0;
    if (what=="Todos" || what=="DiferenciasError") {
      int err = 0;
      for(int k=0;k<nDatos;k++)
        if (CEns[c][k]!=CEns[c+2][k]) err++;

      (*div)[f++][c] = (double)err/nDatos;
    }
    if (what=="Todos" || what=="Secuencial Kappa")
      (*div)[f++][c] = ens->KappaStatistic(CEns[c], CEns[c+2], nDatos, nClases);
    if (what=="Todos" || what=="N00-N01-N10-N11"){
      int N00, N01, N10, N11;
      N00 = N01 = N10 = N11 = 0;
      for(int k=0;k<nDatos;k++) {
        if (clase_ok[k]==CEns[c][k] && clase_ok[k]==CEns[c+2][k])      N11++;
        else if (clase_ok[k]!=CEns[c][k] && clase_ok[k]!=CEns[c+2][k]) N00++;
        else if (clase_ok[k]!=CEns[c][k])                              N01++;
        else                                                           N10++;
      }

      (*div)[f++][c] = (double)N00/nDatos;
      (*div)[f++][c] = (double)N01/nDatos;
      (*div)[f++][c] = (double)N10/nDatos;
      (*div)[f++][c] = (double)N11/nDatos;
    }
    if (what=="Todos" || what_i==5)
      (*div)[f++][c] = QStatistic(CEns[c], CEns[c+2], clase_ok, nDatos);
    if (what=="Todos" || what_i==6)
      (*div)[f++][c] = ens->KappaStatistic(CEns[c], Res[c+2], nDatos, nClases);
  }

  if (what=="Todos" || what=="MelvilleMooney") {
    int err = 0;
    for(c=0;c<nClasf;c++) {
      f=0;
      for(int k=0;k<nDatos;k++){
        if (CEns[c][k]!=CInd[c][k]) err++;
      }
      (*div)[f++][c] = (double)err/(nDatos*(c+1));
    } 
  }
  else if (what=="Todos" || what=="KW") {
    double *n_ok = new double[nDatos];
    for(int k=0;k<nDatos;k++) n_ok[k] = 0.0;

    for(c=0;c<nClasf;c++) {
      f=0;
      double a = 0.0;
      for(int k=0;k<nDatos;k++){
        if (CInd[c][k]==clase_ok[k]) n_ok[k] += 1.0;
        a += n_ok[k]*((c+1)-n_ok[k]);
      }
      (*div)[f++][c] = a/(nDatos*(c+1)*(c+1));
    } 

    delete []n_ok;
  }

  div->redim(f, c);

  for(int i=0;i<nClasf;i++) {
    delete []Mcd[i];
    delete []Mcd[i+nClasf];
    delete []Res[i];
  }
  delete []Mcd;
  delete []Res;

  return this;
}
//---------------------------------------------------------------------------
string DiversityMeasures::ComoCadena()
{
  std::ostringstream buf;
  div->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
EnsembleMeasures::~EnsembleMeasures()
{
  if(mat) delete mat;
}
void EnsembleMeasures::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new EnsembleMeasures(), "Classifiers");
}
//---------------------------------------------------------------------------
void EnsembleMeasures::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());

  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("Bartlett");
  CadenasValidas.push_back("AddErrs");
  ParametroCadena *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(false, 0, "Estadistico", "Medida estadistica a utilizar");
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
Mandato* EnsembleMeasures::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();

  Ensemble *ens = (Ensemble*)Param(0)->ComoDatos();
  Data* dt = (Data*)Param(1)->ComoDatos();
  string what = Param(2)->ComoCadena();
    
  int nClases = ((NomData*)dt)->NumClass;
  int nDatos = dt->GetNTotal();
  int nClasf = ens->GetClassifiersToUse();

  int **Mcd = ens->MatClasif0(dt);
  if (mat) delete mat;
  mat = new Matriz(1, nClasf);
  if (what=="Bartlett") {
    Matriz mrgs(nDatos, nClases, 0.0);
FILE *f=fopen("kk", "w");
    for(int i=0;i<nClasf;i++) {
      for(int j=0;j<nDatos;j++) {
        mrgs[j][Mcd[i][j]]++;
        double m_j = Classifier::Margin(mrgs[j], nClases, dt->GetDatClass(j));
fprintf(f, "%g\t", m_j);
        (*mat)[0][i] += exp(-nDatos*m_j*m_j/8);
      }
fprintf(f, "\n");
    }
fclose(f);
  }
  else if (what=="AddErrs") {
    Matriz mrgs(nDatos, nClases, 0.0);
    for(int i=0;i<nClasf;i++) {
      for(int j=0;j<nDatos;j++) {
        mrgs[j][Mcd[i][j]]++;
        double errs = 0.0;
        for(int k=0;k<nClases;k++) {
          if (k==dt->GetDatClass(j)) continue;
          errs += mrgs[j][k];
        }
        double adderrs = c45::AddErrs(errs+mrgs[j][dt->GetDatClass(j)], errs);
/*        vector<double> distr(nClases, 0.0);
        for(int k=0;k<nClases;k++) {
          distr[k] = mrgs[j][k];
          if (k!=dt->GetDatClass(j)) distr[k] += adderrs/(nClases-1.0);
        }
        if (Classifier::WhichClass(distr)!=dt->GetDatClass(j))
          (*mat)[0][i] += 1.0;*/
        (*mat)[0][i] += adderrs;
      }
      (*mat)[0][i] /= nDatos;
    }
  }


  for(int i=0;i<nClasf;i++) {
    delete []Mcd[i];
    delete []Mcd[i+nClasf];
  }
  delete []Mcd;

  return this;
}
//---------------------------------------------------------------------------
string EnsembleMeasures::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/*class ProgresoEnsemble : public FuncionDeProgreso
{
    protected:
        double llamadasPorTic;
        int longDeBarra;
        int iLlamadas;
        int iTics;

    public:
        ProgresoEnsemble(int longDeBarra=50){
            this->longDeBarra = longDeBarra;
        }
        virtual ~ProgresoEnsemble(){;}

    public:
        virtual void Start(int numeroDeLlamadas) {
            FuncionDeProgreso::Start(numeroDeLlamadas);
        }
        virtual void End() {
        }
        virtual bool operator()(int percent=-1) {
            return true;
         }
};*/
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void SetPropertyValue::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new SetPropertyValue(), "Classifiers");
}
//---------------------------------------------------------------------------
void SetPropertyValue::CreateParams()
{
  Parametros.push_back(new Parametro(false, Tipo::TDatos(), "Objeto", "Elemento a cambiar"));
  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("UseWeights");
  CadenasValidas.push_back("ClassifiersToUse");
  CadenasValidas.push_back("UsePrunedTrees");
  CadenasValidas.push_back("K");
  CadenasValidas.push_back("SetDefaultClassifiers");
  CadenasValidas.push_back("SetGlobalOpt");
  ParametroCadena *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(false, 0, "Propiedad", "Propiedad a modificar");
  Parametros.push_back(p2);
  p2 = new ParametroCadena();
  p2->PonPropiedades(false, 0, "Valor", "Nuevo valor");
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
Mandato* SetPropertyValue::Ejecutar()
{
  Param(0)->Ejecutar(); //Objecto
  Param(1)->Ejecutar(); //Propiedad
  Param(2)->Ejecutar(); //Valor
  
  string prop = Param(1)->ComoCadena();
  
  if (prop.compare("UseWeights")==0) {
    Ensemble *ens = (Ensemble*)Param(0)->ComoDatos();
    ens->SetUseWeights(Param(2)->ComoBooleano());
  }
  else if (prop.compare("ClassifiersToUse")==0) {
    Ensemble *ens = (Ensemble*)Param(0)->ComoDatos();
    ens->SetClassifiersToUse(Param(2)->ComoEntero());
  }
  else if (prop.compare("UsePrunedTrees")==0) {
    C45Tree::UsePrunedTrees = Param(2)->ComoEntero();
//    CARTTree::UsePrunedTrees = Param(2)->ComoEntero();
  }
  else if (prop.compare("K")==0) {
    CARTTree::K = Param(2)->ComoEntero();
  }
  else if (prop.compare("SetDefaultClassifiers")==0) {
    Boosting *boo = (Boosting*)Param(0)->ComoDatos();
    boo->SetDefaultClassifiers();
  }
  else if (prop.compare("SetGlobalOpt")==0) {
    int opt = (Param(0)->ComoCadena().c_str())[0];
    char val[1024];
    strcpy(val, Param(2)->ComoCadena().c_str());
    C45Tree::SetGlobalOpt(opt, val);
  }
  else {
    printf("Unknown property: %s", prop.c_str());
  }

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void SetNTrain::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new SetNTrain());
}
//---------------------------------------------------------------------------
void SetNTrain::CreateParams()
{
  Parametros.push_back(new ParametroData());
  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Numero de ejemplos", "Ejemplos a utilizar como train");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* SetNTrain::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data *dt = (Data*)Param(0)->ComoDatos();
  int ndatos = Param(1)->ComoEntero();

  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    int fin = Param(2)->ComoEntero();
    dt->SetNTrain(ndatos, fin);
  }
  else {
    dt->SetNTrain(ndatos);
  }

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void CreaConjuntosTrainAleat::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new CreaConjuntosTrainAleat(), "Classifiers");
}
//---------------------------------------------------------------------------
void CreaConjuntosTrainAleat::CreateParams()
{
  Parametros.push_back(new ParametroData());
  Parametro *p = new ParametroEntero(-1, 2100000000);
  p->PonPropiedades(false, 0, "Tama�o train", "Ejemplos a utilizar como train");
  Parametros.push_back(p);
  p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Numero de particiones", "Particiones a realizar");
  Parametros.push_back(p);
  ParametroCadena *pc = new ParametroCadena();
  pc->PonPropiedades(false, 0, "Nombre salida", "Nombre base de los ficheros de salida");
  Parametros.push_back(pc);
}
//---------------------------------------------------------------------------
Mandato* CreaConjuntosTrainAleat::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();
  Param(3)->Ejecutar();

  Data *data = (Data*)Param(0)->ComoDatos();
  int ndatos = Param(1)->ComoEntero();
  int nparti =  Param(2)->ComoEntero();
  string pre = Param(3)->ComoCadena();

  int opts = 1;
  if (NumParamsAsignados()>4) {
    Param(4)->Ejecutar();
    opts = Param(4)->ComoEntero();
  }

  DoCreaConjuntosTrainAleat( data, ndatos, nparti, pre, opts);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void GuardaDatos::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new GuardaDatos(), "File");
}
//---------------------------------------------------------------------------
void GuardaDatos::CreateParams()
{
  Parametros.push_back(new ParametroData());
  Parametro *pfs = new ParametroFicheroSalida();
  pfs->PonPropiedades(false, 0, "Guardar con nombre", "Nombre del fichero");
  Parametros.push_back(pfs);
}
//---------------------------------------------------------------------------
Mandato* GuardaDatos::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data *data = (Data*)Param(0)->ComoDatos();
  string filename = Param(1)->ComoCadena();

  int ndatos = data->GetNTotal();
  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    ndatos = Param(2)->ComoEntero();
  }

  data->SaveToFile((char*)filename.c_str(), 0, ndatos-1);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void GuardaClasificador::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new GuardaClasificador(), "File");
}
//---------------------------------------------------------------------------
void GuardaClasificador::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());

  string fil = string("Cualquiera|*.*|") + "Formato rapido (*.clf)|*.clf";
  Parametro *pfs = new ParametroFicheroSalida(fil);
  pfs->PonPropiedades(false, 0, "Guardar con nombre", "Nombre del fichero");
  Parametros.push_back(pfs);
}
//---------------------------------------------------------------------------
Mandato* GuardaClasificador::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Classifier *c = (Classifier*)Param(0)->ComoDatos();
  string f = Param(1)->ComoCadena();
  ofstream ofs(f.c_str(), ios_base::out | ios_base::binary);
  int version = f.rfind(".clf")==f.length()-4 ? 1 : 0;
  c->Guardar(ofs, version);
  ofs.close();

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void InfoClasificador::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new InfoClasificador(), "Classifiers");
}
//---------------------------------------------------------------------------
void InfoClasificador::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());
  Parametro *p = new ParametroEntero(0, 2100000000);
  p->PonPropiedades(false, 0, "Nivel de informacion", "");
  Parametros.push_back(p);
}//---------------------------------------------------------------------------
Mandato* InfoClasificador::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Classifier *c = (Classifier*)Param(0)->ComoDatos();
  int p = Param(1)->ComoEntero();

  kk = c->Info(p);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ConstruyeBaggingCART::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ConstruyeBaggingCART(), "Classifiers");
}
//---------------------------------------------------------------------------
void ConstruyeBaggingCART::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                             "Nombre de la base de datos con los ejemplos"));

  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Numero de arboles CART", 
                                        "Numero que compondran el conjunto");
  Parametros.push_back(p);

  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("Stumps");
  CadenasValidas.push_back("Pruned");
  CadenasValidas.push_back("Unpruned");
  ParametroCadena *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(true, 0, "Podar arboles", "Por omision: El arbol se poda");
  Parametros.push_back(p2);

  p = new ParametroNumero(0.0, 1.0);
  p->PonPropiedades(true, 0, "Tamano de muestreo [0.0;1.0]", "");
  Parametros.push_back(p);

//  p = new ParametroBooleano(true, "Podar arboles",
//                                           "Por omision: El arbol se poda");
//  Parametros.push_back(p);

  p = new ParametroBooleano(true, "Multisplits",
                            "Por omision: se hacen corte unidimensionales");
  Parametros.push_back(p);

  p = new ParametroBooleano(true, "SE 0 rule",
                                             "Por omision: se hace SE 1.0");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* ConstruyeBaggingCART::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();

  int Podar = 2;
  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    Podar = Param(2)->ComoEntero();
  }

  double TamanoMuestreo = 1.0;
  if (NumParamsAsignados()>3) {
    Param(3)->Ejecutar();
    TamanoMuestreo = Param(3)->ComoNumero();
  }
    
  bool MultiSplits = false;
  if (NumParamsAsignados()>4) {
    Param(4)->Ejecutar();
    MultiSplits = Param(4)->ComoBooleano();
  }

  bool SE_0 = false;
  if (NumParamsAsignados()>5) {
    Param(5)->Ejecutar();
    SE_0 = Param(5)->ComoBooleano();
  }

  nd->SetNTrain(nd->GetNTotal());

  if (Clasif) delete Clasif;
  Bagging *b = new Bagging(0.5 + TamanoMuestreo*nd->GetNTrain());

  for(int i=0;i<nArb;i++) {
    if (Podar==1) {
      b->AddClassifier(new DecisionStump());
    }
    else if (Podar==2) {
      b->AddClassifier(new CARTTree(0, 1, true, MultiSplits, SE_0));
    }
    else if (Podar==3) {
      b->AddClassifier(new CARTTree(0, 1, false, MultiSplits, SE_0));
    }
  }

  ProgresoConsola pc;
  b->Build(nd, &pc);

  Clasif = b;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
ConstruyeBaggingC45::~ConstruyeBaggingC45()
{
  //printf("%x\n", Clasif);
  if (Clasif) delete Clasif; 
}
//---------------------------------------------------------------------------
void ConstruyeBaggingC45::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ConstruyeBaggingC45(), "Classifiers");
}
//---------------------------------------------------------------------------
void ConstruyeBaggingC45::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                               "Nombre de la base de datos con los ejemplos"));
  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Numero de arboles C45", "Numero que compondran el conjunto");
  Parametros.push_back(p);
  p = new ParametroBooleano(true, "Podar arboles",
                                            "Por omision: El arbol se poda");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* ConstruyeBaggingC45::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();

  if (Clasif) delete Clasif;
  Bagging *b = new Bagging();

  bool Podar = true;
  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    Podar = Param(2)->ComoBooleano();
  }

  int ndatos = nd->GetNTotal();
  if (NumParamsAsignados()>3) {
    Param(3)->Ejecutar();
    ndatos = Param(3)->ComoEntero();
  }

  nd->SetNTrain(ndatos);
  for(int i=0;i<nArb;i++) {
    b->AddClassifier(new C45Tree(Podar));
  }

  ProgresoConsola pc;
  b->Build(nd, &pc);

  Clasif = b;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
ConstruyeBaggingNNet::~ConstruyeBaggingNNet()
{
  //printf("%x\n", Clasif);
  if (Clasif) delete Clasif; 
}
//---------------------------------------------------------------------------
void ConstruyeBaggingNNet::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ConstruyeBaggingNNet(), "Classifiers");
}
//---------------------------------------------------------------------------
void ConstruyeBaggingNNet::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                               "Nombre de la base de datos con los ejemplos"));
  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Numero de redes neuronales", 
                                          "Numero que compondran el conjunto");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* ConstruyeBaggingNNet::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nNN = Param(1)->ComoEntero();

  if (Clasif) delete Clasif;
  Bagging *b = new Bagging();

  int ndatos = nd->GetNTotal();
  nd->SetNTrain(ndatos);
  for(int i=0;i<nNN;i++) {
    b->AddClassifier(new NNet());
  }

  ProgresoConsola pc;
  b->Build(nd, &pc);

  Clasif = b;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ConstruyeRandomForest::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ConstruyeRandomForest(), "Classifiers");
}
//---------------------------------------------------------------------------
void ConstruyeRandomForest::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                             "Nombre de la base de datos con los ejemplos"));

  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Numero de arboles aleatorios", 
                                        "Numero que compondran el conjunto");
  Parametros.push_back(p);

  p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(true, 0, "Numero de atributos", 
                  "Se seleccionan aleatoriomente en cada nodo, Por omision: 1");
  Parametros.push_back(p);

  p = new ParametroBooleano(true, "Multisplits",
                            "Por omision: se hacen corte unidimensionales");
  Parametros.push_back(p);

  p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(true, 0, "Numero de instancias por hoja", "Por omision: 1");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* ConstruyeRandomForest::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();

  if (Clasif) delete Clasif;
  Bagging *b = new Bagging();

  int NumAtributos = 1;
  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    NumAtributos = Param(2)->ComoEntero();
  }

  bool MultiSplits = false;
  if (NumParamsAsignados()>3) {
    Param(3)->Ejecutar();
    MultiSplits = Param(3)->ComoBooleano();
  }

  int EjemplosPorHoja = 1;
  if (NumParamsAsignados()>4) {
    Param(4)->Ejecutar();
    EjemplosPorHoja = Param(4)->ComoEntero();
  }
    
  nd->SetNTrain(nd->GetNTotal());
  for(int i=0;i<nArb;i++) {
    b->AddClassifier(new RandomForestTree(0, EjemplosPorHoja, NumAtributos, 
                                                                MultiSplits));
  }

  ProgresoConsola pc;
  b->Build(nd, &pc);

  Clasif = b;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/*void ConstruyeBaggingComMod::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ConstruyeBaggingComMod());
}
//---------------------------------------------------------------------------
Mandato* ConstruyeBaggingComMod::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();

  if (Clasif) delete Clasif;
  CompetentModelBagging *b = new CompetentModelBagging();

  int ndatos = nd->GetNTotal();
  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    ndatos = Param(2)->ComoEntero();
  }

  nd->SetNTrain(ndatos);
  for(int i=0;i<nArb;i++) {
    b->AddClassifier(new CARTTree());
  }
  b->Build(nd);

  Clasif = b;

  return this;
}                          */
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ConstruyeBoosting::DarDeAlta()
{
  ifns().AnadirInterfazFuncion( new ConstruyeBoosting(), "Classifiers");
}
//---------------------------------------------------------------------------
void ConstruyeBoosting::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                               "Nombre de la base de datos con los ejemplos"));

  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Numero de arboles CART", "Numero de arboles que compondran el conjunto");
  Parametros.push_back(p);

//  p = new ParametroBooleano(true, "Podar arboles", "Por omision: El arbol se poda");
  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("Stumps");
  CadenasValidas.push_back("Pruned");
  CadenasValidas.push_back("Unpruned");
  ParametroCadena *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(true, 0, "Podar arboles", "Por omision: El arbol se poda");
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
Mandato* ConstruyeBoosting::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();

  if (Clasif) delete Clasif;
  Boosting *b = new Boosting();

  int Podar = 2;
  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    Podar = Param(2)->ComoEntero();
  }

  int ndatos = nd->GetNTotal();
  if (NumParamsAsignados()>3) {
    Param(3)->Ejecutar();
    ndatos = Param(3)->ComoEntero();
  }

  nd->SetNTrain(ndatos);
  for(int i=0;i<nArb;i++) {
    if (Podar==1) {
      b->AddClassifier(new DecisionStump());
    }
    else if (Podar==2) {
      b->AddClassifier(new CARTTree(0, 1, true));
    }
    else if (Podar==3) {
      b->AddClassifier(new CARTTree(0, 1, false));
    }
  }

  ProgresoConsola pc;
  nd->ResetWeights();
  nd->SortOn(nd->GetNumVar() + Data::IniPosIndex, 0, nd->GetNTrain()-1);
  b->Build(nd, &pc);
  nd->ResetWeights();

  Clasif = b;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ConstruyeBoostingC45::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ConstruyeBoostingC45(), "Classifiers");
}
//---------------------------------------------------------------------------
void ConstruyeBoostingC45::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                               "Nombre de la base de datos con los ejemplos"));
  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Numero de arboles C45", "Numero que compondran el conjunto");
  Parametros.push_back(p);
  p = new ParametroBooleano(true, "Podar arboles",
                                            "Por omision: El arbol se poda");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* ConstruyeBoostingC45::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();

  if (Clasif) delete Clasif;
  Boosting *b = new Boosting();

  bool Podar = true;
  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    Podar = Param(2)->ComoBooleano();
  }

  int ndatos = nd->GetNTotal();
  if (NumParamsAsignados()>3) {
    Param(3)->Ejecutar();
    ndatos = Param(3)->ComoEntero();
  }

  nd->SetNTrain(ndatos);
  for(int i=0;i<nArb;i++) {
    b->AddClassifier(new C45Tree(Podar));
  }

  ProgresoConsola pc;
  nd->ResetWeights();
  nd->SortOn(nd->GetNumVar() + Data::IniPosIndex, 0, nd->GetNTrain()-1);
  b->Build(nd, &pc);
  nd->ResetWeights();

  Clasif = b;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ConstruyeClassSwitching::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ConstruyeClassSwitching(), "Classifiers");
}
//---------------------------------------------------------------------------
void ConstruyeClassSwitching::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                               "Nombre de la base de datos con los ejemplos"));

  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Numero de arboles (C45 o CART)", 
                                          "Numero que compondran el conjunto");
  Parametros.push_back(p);

  p = new ParametroNumero(0.0, 1.0);
  p->PonPropiedades(false, 0, "Parametro de switching [0.0;0.5]", "");
  Parametros.push_back(p);

  p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Tipo de arboles (C45=1, CART=-1)", "");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* ConstruyeClassSwitching::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();
  double err = Param(2)->ComoNumero();

  int nClasf = 1;
  if (NumParamsAsignados()>3) {
    Param(3)->Ejecutar();
    nClasf = Param(3)->ComoEntero();
  }
  bool AddC45Trees = nClasf > 0 ? true : false;
  nClasf = nClasf > 0 ? nClasf : -nClasf;

  bool UseDummyCls = false;
  if (NumParamsAsignados()>4) {
    Param(4)->Ejecutar();
    UseDummyCls = Param(4)->ComoBooleano();
  }

  double f_incr = 0.0;
  if (NumParamsAsignados()>5) {
    Param(5)->Ejecutar();
     f_incr = Param(5)->ComoNumero();
  }

  Data *test = 0;
  if (NumParamsAsignados()>6) {
    Param(6)->Ejecutar();
    test = (Data*)Param(6)->ComoDatos();
  }

  if (Clasif) delete Clasif;
  IndepClassifiers *b = new IndepClassifiers(err, UseDummyCls, f_incr);

  bool Podar = false;
  int ndatos = nd->GetNTotal();
  nd->SetNTrain(ndatos);
  if (AddC45Trees) {
    C45Tree::SetGlobalOpt('m', "1");
    C45Tree::SetGlobalOpt('d', "0");
  }
  for(int i=0;i<nArb;i++) {
//    Classifier *c = new NNet();
    Classifier *c = AddC45Trees ? (Classifier *) new C45Tree(Podar) : 
                                  (Classifier *) new CARTTree(0, 1, Podar);
    b->AddClassifier(c);
  }

  ProgresoEnsemble *pe = test ? new ProgresoEnsemble(b, test) : 0;
  ProgresoConsola pc(pe);

  b->Build(nd, &pc);
  if (AddC45Trees) {
    //Reset Default value without checking what was there before
    C45Tree::SetGlobalOpt('m', "2");
    C45Tree::SetGlobalOpt('d', "1");
  }

  Clasif = b;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/*void ConstruyeBoostingWithBaggingInfo::DarDeAlta()
{
  ifns().AnadirInterfazFuncion( new ConstruyeBoostingWithBaggingInfo(),
                                                              "Classifiers");
}
//---------------------------------------------------------------------------
void ConstruyeBoostingWithBaggingInfo::CreateParams()
{
}
//---------------------------------------------------------------------------
Mandato* ConstruyeBoostingWithBaggingInfo::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  double umbral = NumParamsAsignados()>2 ? 
	  Param(2)->Ejecutar()->ComoNumero() : 0.90;
  Ensemble *ens =(Ensemble*)( NumParamsAsignados()>3 ? 
	  Param(3)->Ejecutar()->ComoDatos() : 0);

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();

  if (Clasif) delete Clasif;
  Clasif = 0;

  nd->SetNTrain(nd->GetNTotal());

  Data *reduced = ReduceData(nd, ens, umbral);
  if (reduced->GetNTrain()<nd->GetNTrain()) {
    BoostingWithBaggingInfo *b = new BoostingWithBaggingInfo(umbral, ens);

    bool Podar = true;
    if (NumParamsAsignados()>4) {
      Param(4)->Ejecutar();
      Podar = Param(4)->ComoBooleano();
    }

    for(int i=0;i<nArb;i++) {
      b->AddClassifier(new CARTTree());//C45Tree(Podar));
    }

    reduced->SetNTrain(reduced->GetNTotal());
    ProgresoConsola pc;
    reduced->ResetWeights();
    reduced->SortOn(reduced->GetNumVar() + Data::IniPosIndex, 0, 
                                                    reduced->GetNTrain()-1);
    b->Build(reduced, &pc);
    reduced->ResetWeights();

    Clasif = b;
  }

  return this;
}                   */
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ConstruyeCART::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                               "Nombre de la base de datos con los ejemplos"));
  Parametro *p = new ParametroBooleano(true, "Podar",
                                            "Por omision: El arbol se poda");
  Parametros.push_back(p);
  p = new ParametroBooleano(true, "Multisplits",
                            "Por omision: se hacen corte unidimensionales");
  Parametros.push_back(p);
  p = new ParametroBooleano(true, "SE 0 rule",
                                             "Por omision: se hace SE 1.0");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
void ConstruyeCART::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ConstruyeCART(), "Classifiers");
}
//---------------------------------------------------------------------------
Mandato* ConstruyeCART::Ejecutar()
{
  Param(0)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();

  if (Clasif) delete Clasif;

  bool Podar = true;
  if (NumParamsAsignados()>1) {
    Param(1)->Ejecutar();
    Podar = Param(1)->ComoBooleano();
  }

  bool MultiSplits = false;
  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    MultiSplits = Param(2)->ComoBooleano();
  }

  bool SE_0 = false;
  if (NumParamsAsignados()>3) {
    Param(3)->Ejecutar();
    SE_0 = Param(3)->ComoBooleano();
  }

  int ndatos = nd->GetNTotal();
  if (NumParamsAsignados()>4) {
    Param(4)->Ejecutar();
    ndatos = Param(4)->ComoEntero();
  }

  nd->SetNTrain(ndatos);
  Clasif = new CARTTree(0, 1, Podar, MultiSplits, SE_0);
  Clasif->Build(nd);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ConstruyeC45::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ConstruyeC45(), "Classifiers");
}
//---------------------------------------------------------------------------
void ConstruyeC45::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                               "Nombre de la base de datos con los ejemplos"));
  Parametro *p = new ParametroBooleano(true, "Podar",
                                            "Por omision: El arbol se poda");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* ConstruyeC45::Ejecutar()
{
  Param(0)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();

  if (c45tree) delete c45tree;

  bool podar = true;
  if (NumParamsAsignados()>1) {
    Param(1)->Ejecutar();
    podar = Param(1)->ComoBooleano();
  }

  int ndatos = nd->GetNTotal();
  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    ndatos = Param(2)->ComoEntero();
  }
  nd->SetNTrain(ndatos);

  c45tree = new C45Tree(podar);
  c45tree->Build(nd);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ConstruyeNNet::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ConstruyeNNet(), "Classifiers");
}
//---------------------------------------------------------------------------
void ConstruyeNNet::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                               "Nombre de la base de datos con los ejemplos"));
}
//---------------------------------------------------------------------------
Mandato* ConstruyeNNet::Ejecutar()
{
  Param(0)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();

  if (nnet) delete nnet;

  int ndatos = nd->GetNTotal();
  nd->SetNTrain(ndatos);

  nnet = new NNet();
  nnet->Build(nd);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/*void EnsembleReport::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new EnsembleReport());
}
//---------------------------------------------------------------------------
Mandato* EnsembleReport::Ejecutar()
{ */
/*  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *train   = (Data*)Param(1)->ComoDatos();
  Data *test    = (Data*)Param(2)->ComoDatos();
*/
  /*return this;
} */
//---------------------------------------------------------------------------
/*string EnsembleReport::ComoCadena()
{
  string cad="";
//  for(unsigned i=0;i<certezas.size();i++)
//    cad = cad + "\t" + Mandato::ACad(certezas[i]);
  return cad;
}                              */
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void MapaDeClasif2D::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new MapaDeClasif2D(), "Tools");
}
//---------------------------------------------------------------------------
void MapaDeClasif2D::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());

  Parametro *p = new ParametroEntero(0, 2100000000);
  p->PonPropiedades(false, 0, "Atrib. 1", "");
  Parametros.push_back(p);
  p = new ParametroEntero(0, 2100000000);
  p->PonPropiedades(false, 0, "Atrib. 2", "");
  Parametros.push_back(p);

  p = new ParametroNumero();
  p->PonPropiedades(false, 0, "Minimo (Atrib. 1)", "");
  Parametros.push_back(p);
  p = new ParametroNumero();
  p->PonPropiedades(false, 0, "Minimo (Atrib. 2)", "");
  Parametros.push_back(p);
  p = new ParametroNumero();
  p->PonPropiedades(false, 0, "Maximo (Atrib. 1)", "");
  Parametros.push_back(p);
  p = new ParametroNumero();
  p->PonPropiedades(false, 0, "Maximo (Atrib. 2)", "");
  Parametros.push_back(p);

  p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Ancho imagen de salida", "");
  Parametros.push_back(p);
  p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Alto imagen de salida", "");
  Parametros.push_back(p);
  p = new ParametroCadena();
  p->PonPropiedades(false, 0, "Nombre imagen de salida", "");
  Parametros.push_back(p);

}
//---------------------------------------------------------------------------
Mandato* MapaDeClasif2D::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();
  Param(3)->Ejecutar();
  Param(4)->Ejecutar();
  Param(5)->Ejecutar();
  Param(6)->Ejecutar();
  Param(7)->Ejecutar();
  Param(8)->Ejecutar();
  Param(9)->Ejecutar();

  int x1 = Param(1)->ComoEntero();
  int x2 = Param(2)->ComoEntero();
  double min1 = Param(3)->ComoNumero();
  double min2 = Param(4)->ComoNumero();
  double max1 = Param(5)->ComoNumero();
  double max2 = Param(6)->ComoNumero();
  int Ancho = Param(7)->ComoEntero();
  int Alto  = Param(8)->ComoEntero();
  string nom  = Param(9)->ComoCadena();
//  double div1 = (max1-min1)/Ancho;
//  double div2 = (max2-min2)/Alto;
  Classifier* c = (Classifier*)Param(0)->ComoDatos();
  vector<double> dato;
  dato.push_back(min1);
  dato.push_back(min2);

  MargenFunction *fm = new MargenFunction(c);
  MargenClasfFunction *fmc = new MargenClasfFunction(c);

  DoMapaDeClasif2D(c, x1, x2, min1, min2, max1, max2, Ancho, Alto, dato,
                                                  (nom + "clasf.bmp").c_str());
/*  DoMapaDeAlturas2D(fm, x1, x2, min1, min2, max1, max2, Ancho, Alto, dato,
                                                 (nom + "margen.bmp").c_str());
  DoMapaDeAlturas2D(fmc, x1, x2, min1, min2, max1, max2, Ancho, Alto, dato,
                                            (nom + "margenclasf.bmp").c_str());
*/
  delete fm;
  delete fmc;
  
  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void TestAll::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new TestAll(), "Tools");
}
//---------------------------------------------------------------------------
void TestAll::CreateParams()
{
//
}
//---------------------------------------------------------------------------
void savematriz(int **Mcd, int nClasf, int nDatos, char *nom)
{
    FILE *fmat=fopen(nom, "wt");
    for(int i=0;i<nClasf;i++) {
      for(int j=0;j<nDatos;j++) {
        fprintf(fmat, "%d\t", Mcd[i][j]);
      }
      fprintf(fmat, "\n");
    }
    for(int i=0;i<nClasf;i++) {
      for(int j=0;j<nDatos;j++) {
        fprintf(fmat, "%d\t", Mcd[i+nClasf][j]);
      }
      fprintf(fmat, "\n");
    }
    fclose(fmat);
}
void deletematriz(int **Mcd, int nClasf)
{
    for(int i=0;i<nClasf;i++) {
      delete []Mcd[i];
      delete []Mcd[i+nClasf];
    }
    delete []Mcd;
}
Mandato* TestAll::Ejecutar()
{
  Param(0)->Ejecutar();

  string what = Param(0)->ComoCadena();
  
  if (what.compare("ComModGraf")==0) { //             ---  ComModGraf  ---
    class FunVotes : public Function {
      CompetentModelBagging *k;
      int t;
      public:
        FunVotes(CompetentModelBagging *_k, int _t):Function(2){k = _k; t = _t;}
        virtual double f(double *xs){
        double *x = new double[3];
        double val=0.0;
        vector<int> A, B;
        Data *dat = k->GetData();
        x[0] = dat->GetValueVar(0,0); x[1] = dat->GetValueVar(0,1);
        dat->SetValueVar(0,0,xs[0]);dat->SetValueVar(0,1,xs[1]);
        k->GetVotes(0, A, B);
        if (t==1) val=A[0]+A[1];
        else if (t==2) val = ((double)(A[0]>A[1] ? A[0] : A[1]))/(A[0]+A[1]);
        else if (t==3) val = ((double)(B[0]>B[1] ? B[0] : B[1]))/(B[0]+B[1]);
        dat->SetValueVar(0,0,x[0]);dat->SetValueVar(0,1,x[1]);
        delete []x;
        return val; 
              }
              virtual double max(){return t==1 ? 200.0 : 1.0;}
              virtual double min(){return t==1 ? 0.0 : 0.5;}
          };
          int x1 = 0;
          int x2 = 1;
          double min1 = -2.5;
          double min2 = -2.5;
          double max1 = 2.5;
          double max2 = 2.5;
          int Ancho = 160;
          int Alto  = 160;
          string nom  = "kk";

          Param(1)->Ejecutar();
          CompetentModelBagging* c = 
                                (CompetentModelBagging*)Param(1)->ComoDatos();
          vector<double> dato;
          dato.push_back(min1);
          dato.push_back(min2);
          
          FunVotes *f = new FunVotes(c, 1);
//    DoMapaDeAlturas2D(f, x1, x2, min1, min2, max1, max2, Ancho, Alto, dato,
  //                                              (nom + "numclasf.bmp").c_str());
          f = new FunVotes(c, 2);
          DoMapaDeAlturas2D(f, x1, x2, min1, min2, max1, max2, Ancho, Alto, 
                                        dato, (nom + "margenRefs.bmp").c_str());
//    f = new FunVotes(c, 3);
//    DoMapaDeAlturas2D(f, x1, x2, min1, min2, max1, max2, Ancho, Alto, dato,
  //                                          (nom + "margenNoRefs.bmp").c_str());
  
    delete f;
  }
  else if (what.compare("ExpandData")==0) { //          ---  ExpandData  ---
    Param(1)->Ejecutar();
    Param(2)->Ejecutar();
    Data* dt = (Data*)Param(1)->ComoDatos();
    Classifier *cls = (Classifier*)Param(2)->ComoDatos();
    algo = ExpandData(dt, cls);
  }
  else if (what.compare("ExpandData2")==0) { //        ---  ExpandData2  ---
    Param(1)->Ejecutar();
    Param(2)->Ejecutar();
    Data* dt = (Data*)Param(1)->ComoDatos();
    Classifier *cls = (Classifier*)Param(2)->ComoDatos();
    algo = ExpandData2(dt, cls);
  }
  else if (what.compare("ConjuntoIGP")==0) {//---  ConjuntoIGP  ---
    Param(1)->Ejecutar();
    Param(2)->Ejecutar();
    Param(3)->Ejecutar();

    Data* nd = (Data*)Param(1)->ComoDatos();
    int nArb = Param(2)->ComoEntero();
    Data* ts = (Data*)Param(3)->ComoDatos();

    GroupBagging *b = new GroupBagging();

    for(int i=0;i<nArb;i++) {
      Classifier *c = new GPTree(); 
      b->AddClassifier(c);
    }

    ProgresoConsola pc;
    b->Build(nd, &pc);
    
    vector<double> errores;
    char nf[256];

    b->SetData(nd);
    b->SecuencialError(0, nd->GetNTrain()-1, &errores);
    sprintf(nf, "train_igp_%d.txt", nd->GetNTrain());
    FILE *f = fopen(nf, "a");
    FILE *ft = fopen("tam_arbol.txt", "a");
    for (unsigned i=0; i<errores.size(); i++) {
      fprintf(f, "%g\t", errores[i]);
      fprintf(ft, "%d\t", ((GPTree*)b->GetClassifier(i))->GetTree()->CountNodes());
      fprintf(ft, "%d\t", ((GPTree*)b->GetClassifier(i))->GetTree()->Tmax);
    }
    fprintf(f, "\n");
    fprintf(ft, "\n");
    fclose(f);
    fclose(ft);

    b->SetData(ts);
    b->SecuencialError(0, ts->GetNTotal()-1, &errores);
    sprintf(nf, "test_igp_%d.txt", nd->GetNTrain());
    f = fopen(nf, "a");
    for (unsigned i=0; i<errores.size(); i++) {
      fprintf(f, "%g\t", errores[i]);
    }
    fprintf(f, "\n");
    fclose(f);

  }
  else if (what.compare("IndepClassifiers")==0) {//---  IndepClassifiers  ---
    Param(1)->Ejecutar();
    Param(2)->Ejecutar();
    Param(3)->Ejecutar();

    int nClasf = 1;
    if (NumParamsAsignados()>4) {
      Param(4)->Ejecutar();
      nClasf = Param(4)->ComoEntero();
    }
    bool AddC45Trees = nClasf > 0 ? true : false;
    nClasf = nClasf > 0 ? nClasf : -nClasf;

    bool UseDummyCls = false;
    if (NumParamsAsignados()>5) {
      Param(5)->Ejecutar();
      UseDummyCls = Param(5)->ComoBooleano();
    }

    double f_incr = 0.0;
    if (NumParamsAsignados()>6) {
      Param(6)->Ejecutar();
       f_incr = Param(6)->ComoNumero();
    }

    Data *test = 0;
    if (NumParamsAsignados()>7) {
      Param(7)->Ejecutar();
      test = (Data*)Param(7)->ComoDatos();
    }

    Data* nd = (Data*)Param(1)->ComoDatos();
    int nArb = Param(2)->ComoEntero();
    double err = Param(3)->ComoNumero();

    if (algo) delete (IndepClassifiers*)algo;
    IndepClassifiers *b = new IndepClassifiers(err, UseDummyCls, f_incr);

    bool Podar = false;
    int ndatos = nd->GetNTotal();
    nd->SetNTrain(ndatos);
    if (AddC45Trees) {
      C45Tree::SetGlobalOpt('m', "1");
      C45Tree::SetGlobalOpt('d', "0");
    }
    for(int i=0;i<nArb;i++) {
//      Boosting *bo = new Boosting();
  //    for(int i=0;i<nClasf;i++) {
    //    if (AddC45Trees) bo->AddClassifier(new C45Tree(Podar)); 
      //  else             bo->AddClassifier(new CARTTree(0, 1, Podar));
//      }
  //    bo->SetForceContinuationByReweighting(false);
    //  b->AddClassifier(bo);
      Classifier *c = AddC45Trees ? (Classifier *) new C45Tree(Podar) : 
                                      (Classifier *) new CARTTree(0, 1, Podar);
      b->AddClassifier(c);
    }

    ProgresoEnsemble *pe = test ? new ProgresoEnsemble(b, test) : 0;
    ProgresoConsola pc(pe);

    b->Build(nd, &pc);
    if (AddC45Trees) {
      //Reset Default value without checking what was there before
      C45Tree::SetGlobalOpt('m', "2");
      C45Tree::SetGlobalOpt('d', "1");
    }

    algo = b;
  }
  else if (what.compare("Multiclassing")==0) { // ---  Multiclassing  ---
    Param(1)->Ejecutar();
    Param(2)->Ejecutar();
    Param(3)->Ejecutar();
    Param(4)->Ejecutar();

    Data* dtr = (Data*)Param(1)->ComoDatos();
    int nArb = Param(2)->ComoEntero();
    double p = Param(3)->ComoNumero();
    Data* dts = (Data*)Param(4)->ComoDatos();

    //Duplicamos el numero de clases
    NomData *nd = (NomData*)dtr->Clone(0, dtr->GetNTotal()-1);
    nd->SetNTrain(nd->GetNTotal());
    int nC = nd->NumClass;
    int icls = nd->GetNumVarNom();
    for(int i=0;i<nC;i++) {
       string nomC = string("~~~") + nd->GetTerm(i, icls);
       nd->AddClass(nomC);
       ((NomData*)dtr)->AddClass(nomC);
       ((NomData*)dts)->AddClass(nomC);
    }
    int kk = nd->GetNumVar()-1;
    nd->SortOn(kk, 0, nd->GetNTotal()-1);
    for(int i=0;i<nd->GetNTotal();i++) {
      int nuevaC = (int)(0.5 + nd->GetDatClass(i) + (Flip() ? nC : 0));
      nd->SetValueVar(i, kk, nuevaC);
    }

    if (algo) delete (IndepClassifiers*)algo;
    IndepClassifiers *b = new IndepClassifiers(p);

    int ndatos = nd->GetNTotal();
    nd->SetNTrain(ndatos);
    C45Tree::SetGlobalOpt('m', "1");
    C45Tree::SetGlobalOpt('d', "0");
    for(int i=0;i<nArb;i++) 
      b->AddClassifier(new C45Tree(false));

    ProgresoConsola pc;
    b->Build(nd, &pc);
    C45Tree::SetGlobalOpt('m', "2");
    C45Tree::SetGlobalOpt('d', "1");

    int err = 0;
    b->SetData(dtr);
    int **Mcd = b->MatClasif0(dtr);
    savematriz(Mcd, nArb, dtr->GetNTotal(), "mtr.csv");
    deletematriz(Mcd, nArb);
    for(int i=0;i<dtr->GetNTotal();i++) {
      int c = b->Classify(i);
      if (c%nC!=dtr->GetDatClass(i))
        err++;
    }
    FILE *f = fopen("train.txt","a");
    fprintf(f, "%.10g\n", (double)err/dtr->GetNTrain());
    fclose(f);

    err = 0;
    b->SetData(dts);
    Mcd = b->MatClasif0(dts);
    savematriz(Mcd, nArb, dts->GetNTotal(), "mts.csv");
    deletematriz(Mcd, nArb);
    for(int i=0;i<dts->GetNTotal();i++) {
      int c = b->Classify(i);
      if (c%nC!=dts->GetDatClass(i))
        err++;
    }
    f = fopen("test.txt","a");
    fprintf(f, "%.10g\n", (double)err/dts->GetNTrain());
    fclose(f);

    f = fopen("info.txt","a");
    fprintf(f, "%s\n", b->Info(1).c_str());
    fclose(f);

    b->SetData(nd);

    algo = b;
  }
  else if (what.compare("Flipping")==0) {          //---  Flipping  ---
    Param(1)->Ejecutar();
    Param(2)->Ejecutar();
    Param(3)->Ejecutar();

    int nClasf = 1;
    if (NumParamsAsignados()>4) {
      Param(4)->Ejecutar();
      nClasf = Param(4)->ComoEntero();
    }
    bool AddC45Trees = nClasf > 0 ? true : false;
    nClasf = nClasf > 0 ? nClasf : -nClasf;

    Data* nd = (Data*)Param(1)->ComoDatos();
    int nArb = Param(2)->ComoEntero();
    double err = Param(3)->ComoNumero();

    if (algo) delete (IndepClassifiers*)algo;
    IndepClassifiers *b = new IndepClassifiers(err);

    bool Podar = false;
    int ndatos = nd->GetNTotal();
    nd->SetNTrain(ndatos);
    if (AddC45Trees) {
      C45Tree::SetGlobalOpt('m', "1");
      C45Tree::SetGlobalOpt('d', "0");
    }
    for(int i=0;i<nArb;i++) {
    //  Boosting *bo = new Boosting();
      //for(int i=0;i<nClasf;i++) {
  //      if (AddC45Trees) bo->AddClassifier(new C45Tree(Podar)); 
//        else             bo->AddClassifier(new CARTTree(0, 1, Podar));
//      }
  //    bo->SetForceContinuationByReweighting(false);
//      b->AddClassifier(bo);
      Classifier *c = AddC45Trees ? (Classifier *) new C45Tree(Podar) : 
                                      (Classifier *) new CARTTree(0, 1, Podar);
      b->AddClassifier(c);
    }

    b->SetAsBreiman(true);
    ProgresoConsola pc;
    b->Build(nd, &pc);
    if (AddC45Trees) {
      //Reset Default value without checking what was there before
      C45Tree::SetGlobalOpt('m', "2");
      C45Tree::SetGlobalOpt('d', "1");
    }

    algo = b;
  }
  else if (what.compare("KappaStatistic")==0) { //---  KappaStatistic  ---
    Param(1)->Ejecutar();
    Param(2)->Ejecutar();
    Param(3)->Ejecutar();

    Ensemble *ens = (Ensemble*)Param(1)->ComoDatos();
    Data* dt = (Data*)Param(2)->ComoDatos();
    string nf = Param(3)->ComoCadena();
    
  int nClases = ((NomData*)dt)->NumClass;
  int nDatos = dt->GetNTotal();
  int nClasf = ens->Count();

    int **Mcd = ens->MatClasif0(dt);
  vector<int> e;
  for(int i=0;i<nClasf;i++) {
    e.push_back(0);
    for(int j=0;j<nDatos;j++) {
      if (Mcd[i][j]!=dt->GetDatClass(j)) e[i]++;
    }
  }
  double k, q;
  FILE *fke = fopen(nf.c_str(), "wt");
  FILE *fkp = fopen("kappa.csv", "wt");
  for(int i=nClasf-1;i>=0;i--) {
    for(int j=nClasf-1;j>i;j--)
      fprintf(fkp, "\t");
       
    for(int j=i-1;j>=0;j--) {
      k = Classifier::KappaStatistic(Mcd[i], Mcd[j], nDatos, nClases);
      q = QStatistic(Mcd[i+nClasf], Mcd[j+nClasf], nDatos);
      fprintf(fkp, "%g\t", k);
      fprintf(fke, "%g\t", k);
      fprintf(fke, "%g\t", q);
      fprintf(fke, "%g\t", (double)e[i]/nDatos);
      fprintf(fke, "%g\t", (double)e[j]/nDatos);
      fprintf(fke, "%g\t", ((double)e[i]+e[j])/(2.0*nDatos));
      fprintf(fke, "%d\t", i);
      fprintf(fke, "%d\n", j); 
    }
    fprintf(fkp, "\n");
  }
  fclose(fke);
  fclose(fkp);

    FILE *fmat=fopen("mat0.txt", "wt");
    for(int i=0;i<nClasf;i++) {
      for(int j=0;j<nDatos;j++) {
        fprintf(fmat, "%d\t", Mcd[i][j]);
      }
      delete []Mcd[i];
      fprintf(fmat, "\n");
    }
    for(int i=0;i<nClasf;i++) {
      for(int j=0;j<nDatos;j++) {
        fprintf(fmat, "%d\t", Mcd[i+nClasf][j]);
      }
      delete []Mcd[i+nClasf];
      fprintf(fmat, "\n");
    }
    delete []Mcd;
    fclose(fmat);
  }
  else if (what.compare("Divisiones")==0) { //---  Divisiones  ---
    Param(1)->Ejecutar();
    Param(2)->Ejecutar();
    
    Ensemble *ens = (Ensemble*)Param(1)->ComoDatos();
    Data* dt = (Data*)Param(2)->ComoDatos();

    for(int i=0;i<dt->GetNumVar()-1;i++){
    }
  }
  else if (what.compare("AllTreeTest")==0) { //---  AllTreeTest  ---

    Param(1)->Ejecutar();
    Ensemble *e = (Ensemble*)Param(1)->ComoDatos();
    vector<set<double> > splits(e->GetData()->GetNumVar());

    for(int i=0;i<e->GetClassifiersToUse();i++) {
      Tree *t = ((DecisionTree*)e->GetClassifier(i))->GetTree();
      int K = ((DecisionTree*)e->GetClassifier(i))->GetK();
      Node *cur = t->GetRoot();
      int numLeafs = 0;
      while(cur->child) cur = cur->child;
      while(cur) {
        if (!cur->IsLeaf(K)) {
          splits[cur->att].insert(cur->fSplit);
        }
        else {
          numLeafs++;
        }
        cur = cur->nextDown(K);
      }
      cout << numLeafs << "*";
    }
    cout << endl;
    for (unsigned i=0;i<splits.size();i++) {
      cout << "Atts: " << i << " - Elements: " << splits[i].size() << endl;
      copy(splits[i].begin(), splits[i].end(), 
                                          ostream_iterator<double>(cout, " "));
      cout << endl << "----------------" << endl;
    }
  }
  else if (what.compare("SpaceCubes")==0) { //---  SpaceCubes  ---
    Param(1)->Ejecutar();
    Param(2)->Ejecutar();

    Ensemble *e = (Ensemble*)Param(1)->ComoDatos();
    Data *d = (Data*)Param(2)->ComoDatos();

/*  for(int i=0;i<e->Count();i++) {
    Tree *t = ((CARTTree*)e->GetClassifier(i))->GetTree();
    e->SetClassifierWeight(i, 1.0/t->GetRoot()->T);
  }

  e->OrdenarClasificadoresPorPeso();*/

    CARTTree* t =  (CARTTree*)e->GetClassifier(0);
    vector<double> mins, maxs;
    for(int ik=0;ik<t->GetData()->GetNumVar()-1;ik++) {
      mins.push_back(t->GetData()->Min(ik));
      maxs.push_back(t->GetData()->Max(ik));
    }
    SpaceCubes sc(t->GetTree(), mins, maxs);
//M    printf("Acumulad + tree      = resultad -> after merging\n");
//printf("%s", sc.Print().c_str());
//SpaceCubes sc2(sc.ToTree(), mins, maxs);
//sc.ToTree();
//printf("%s", sc2.Print().c_str());
    for(int i=1;i<e->GetClassifiersToUse();i++) {
      int ac, tr, re;
      ac = sc.size();
      t = (CARTTree*)e->GetClassifier(i);
      SpaceCubes spt(t->GetTree(), mins, maxs);
      tr = spt.size();
      sc = sc & spt; 
//sc.ToTree();
      re = sc.size();
      while (sc.Merge()) ;
      printf("  %6d +    %5d =   %6d ->  %7d\n", ac, tr, re, sc.size());
    }
//    printf("%s", sc.Print().c_str());
/*SpaceCubes sc2(sc.ToTree(), mins, maxs);
printf("\nSIZE: %d", sc2.size());
while (sc2.Merge()) ;
printf("\nSIZE: %d", sc2.size());*/

/*M    for(int i=0;i<sc.size();i++) {
      HyperCube *hc = sc.GetHyperCube(i);
      hc->tag = 255.99*hc->tag/e->Count();
    }
    FILE *ff=fopen("too", "w");
    fprintf(ff, "%s", sc.Print().c_str());
    fclose(ff);
*/
    for(int i=0;i<sc.size();i++) {
      HyperCube *hc = sc.GetHyperCube(i);
      hc->tag = int(hc->tag);
    }
    printf("#cubos: %d\n", sc.size());

    while (sc.Merge()>0) 
      printf("#cubos: %d\n", sc.size());

//    printf("%s", sc.Print().c_str());
    algo = 0;//sc.ToData();
SpaceCubes sc2(sc.ToTree(), mins, maxs);
printf("\nSIZE: %d", sc2.size());
while (sc2.Merge()) ;
printf("\nSIZE: %d", sc2.size());
    algo = sc2.ToTree();
char buff[32000];
    ((Tree*)algo)->Write(d, -1, buff);
printf("%s\n", buff);
  }
  else if (what.compare("MinimumCostTrees")==0) { //---  MinimumCostTrees  ---
    Param(1)->Ejecutar();

    CARTTree *cart = (CARTTree*)Param(1)->ComoDatos();
    map<Node*, PruneInfo*> trees = cart->GetTree()->MinimumCostTrees();

    Node *cur = cart->GetTree()->GetRoot();
    while(cur->child) cur = cur->child;
    while (cur) {
      PruneInfo *pi = trees[cur];
//      cout << "---------------------" << cur->Rleaf;
//      cout << "(" << cur->RsubTree << ") ----\n";
      char buf[1000];
      cart->GetData()->WriteNode(cur, buf, -1);
      cout << "-----------------------------------" << "\n";
      cout << buf;
      cout << "l:";
      for(unsigned i=0;i<pi->min_cost_left.size();i++) 
        cout << "\t" << pi->min_cost_left[i];
      cout << "\nr:";
      for(unsigned i=0;i<pi->min_cost_right.size();i++) 
        cout << "\t" << pi->min_cost_right[i];
      cout << "\ne:";
      for(unsigned i=0;i<pi->min_cost.size();i++) 
        cout << "\t" << pi->min_cost[i];
      cout << "\n";
      cur = cur->nextDown(-1);
    }
    cart->GetTree()->SelectMinimumCostTrees(trees);
  } else if (what.compare("ValidacionMCT")==0) { // Validacion cruzada para MCT
    Param(1)->Ejecutar();
    Data *data = (Data*)Param(2)->ComoDatos(); // Seguro que se puede sacar de Param(1)...
	  
    CARTTree *cart = (CARTTree*)Param(1)->ComoDatos();
    map<Node*, PruneInfo*> trees = cart->GetTree()->MinimumCostTrees();	  
	  
	cart->GetTree()->SelectMinimumCostTrees(trees);
	  	
	  cout << "Llego" << endl;
	  
	int Nmin = 10; // valor de esto?
	bool SE_0 = false; 
	cart->GetTree()->PruneMCT(data, Nmin, SE_0);
	  
	// Faltan las funciones de test.
  }
  
  if (what.compare("Twonormcilla")==0) { //---  Twonormcilla  --- 
/*    double d[2];
    d[0] = d[1] = 1.0;
    Gaussiana *g11 = new Gaussiana(2, 0, d);
    d[0] = d[1] = -1.0;
    Gaussiana *g12 = new Gaussiana(2, 0, d);
    GaussianaMultiple *clase1= new GaussianaMultiple();
    clase1->AnadirGaussiana(g11, 0.5);
    clase1->AnadirGaussiana(g12, 0.5);

    d[1] = 1.0;
    Gaussiana *clase2 = new Gaussiana(2, 0, d);

    vector<double> APrioriProb(2, 0.5);
    vector<ProbFunction*> ProbFuns;
    ProbFuns.push_back(clase2);
    ProbFuns.push_back(clase1);
    BayesClassifier *clasf = new BayesClassifier(APrioriProb, ProbFuns);
    vector<double> dato;
    dato.push_back(-2);
    dato.push_back(-2);
    DoMapaDeClasif2D(clasf, 0, 1, -2, -2, 2, 2, 600, 600, dato,
                                                      "twonormcilla_clasf.bmp");
*/
    Gaussiana *g11 = new Gaussiana(2);
    double **d = new double*[2];
    d[0] = new double[2];
    d[1] = new double[2];
    d[0][0] = d[1][1] = 2.0;
    d[1][0] = d[0][1] = 0.0;
    Gaussiana *g12 = new Gaussiana(2, d);
    delete []d[0];
    delete []d[1];
    delete []d;

    vector<double> APrioriProb(2, 0.5);
    vector<ProbFunction*> ProbFuns;
    ProbFuns.push_back(g11);
    ProbFuns.push_back(g12);
    BayesClassifier *clasf = new BayesClassifier(APrioriProb, ProbFuns);
    vector<double> dato;
    dato.push_back(-3);
    dato.push_back(-3);
    DoMapaDeClasif2D(clasf, 0, 1, -2, -2, 2, 2, 600, 600, dato,
                                                      "circle_gauss_clasf.bmp");
  }
  else {
  }
  
  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
