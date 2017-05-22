//---------------------------------------------------------------------------


#include "FnsClasificacion.h"
#include "data20.h"
#include "node20.h"

#include <typeinfo>

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
#include "Evaluations.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <iterator>
#include <string.h>
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
  LoadClassifier::DarDeAlta();
  LoadDataset::DarDeAlta();
  GenerateWrongLabels::DarDeAlta();
  Corrected::DarDeAlta();
  LoadLabels::DarDeAlta();
  FilterData::DarDeAlta();
  SetSplitCriterium::DarDeAlta();
  Classify::DarDeAlta();
  Error::DarDeAlta();
  ResubstituteNodeStats::DarDeAlta();
  EvaluateClassifier::DarDeAlta();
  Mrse::DarDeAlta();
  SequentialError::DarDeAlta();
  ClassificationMatrix::DarDeAlta();
  ConfusionMatrix::DarDeAlta();
  Margin::DarDeAlta();
  Certainty::DarDeAlta();
  MarginOfInstances::DarDeAlta();
  MarginSum::DarDeAlta();
  EstimErrorVal::DarDeAlta();
  EstimErrorTest::DarDeAlta();
  DiversityMeasures::DarDeAlta();
  EnsembleMeasures::DarDeAlta();
  SetPropertyValue::DarDeAlta();
  SetNTrain::DarDeAlta();
  GeneratePartitions::DarDeAlta();
  SaveDataset::DarDeAlta();
  SaveClassifier::DarDeAlta();
  ClassifierInfo::DarDeAlta();
  BuildIGPEnsemble::DarDeAlta();
  BuildBaggingCART::DarDeAlta();
  BuildBaggingC45::DarDeAlta();
  BuildBaggingNNet::DarDeAlta();
  BuildRandomForest::DarDeAlta();
//  BuildBaggingComMod::DarDeAlta();
//  BuildBoostingWithBaggingInfo::DarDeAlta();
  BuildBoosting::DarDeAlta();
  BuildBoostingC45::DarDeAlta();
  BuildClassSwitching::DarDeAlta();
  BuildCART::DarDeAlta();
  BuildC45::DarDeAlta();
  BuildNNet::DarDeAlta();
  ClassificationMap2D::DarDeAlta();
  TestAll::DarDeAlta();
  Clear::DarDeAlta();
  GValue::DarDeAlta();
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
      data->Scramble(0, numdatos-1);
      data->SaveToFile(saltr, 0, numdatos-1);
      data->SaveToFile(salts, numdatos, data->GetTotalData()-1);
    }
  }
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void LoadClassifier::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new LoadClassifier(), "File");
}
//---------------------------------------------------------------------------
void LoadClassifier::CreateParams()
{
  string fil = string("All|*.*|") + "Fast format (*.clf)|*.clf";
  Parametro *pfe = new ParametroFicheroEntrada(fil);
  pfe->PonPropiedades(false, 0, "File Name", 
                                    "Name of the file with the classifier");
  Parametros.push_back(pfe);

  Parametro *p = new ParametroData();
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Tipo* LoadClassifier::DameTipo()
{
  if (dynamic_cast<Ensemble*>(Clasif))
    return TipoEnsemble::TEnsemble(); 

  return TipoClasificador::TClasificador();
}
//---------------------------------------------------------------------------
Mandato* LoadClassifier::Ejecutar()
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
LoadDataset::~LoadDataset()
{
  if(nd) delete nd;
}
//---------------------------------------------------------------------------
void LoadDataset::CreateParams()
{
  string fil = string("All|*.cre;*.asc;*.names;|") + "Self (*.cre)|*.cre|" +
                              "C45 (*.names)|*.names|" + "Raetsch(*.asc)|*.asc";
  Parametro *pfe = new ParametroFicheroEntrada(fil);
  pfe->PonPropiedades(false, 0, "File", "File name of the dataset");
  Parametros.push_back(pfe);

  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("BasicInstance");
  CadenasValidas.push_back("UInt8Instance");
  CadenasValidas.push_back("Int32Instance");
  CadenasValidas.push_back("FloatInstance");
  CadenasValidas.push_back("ShortInstance");
  CadenasValidas.push_back("ComboInstance");
  CadenasValidas.push_back("MultilabelInstance");
  CadenasValidas.push_back("MultiobjectiveInstance");
  Parametro *p2 = new ParametroCadena(&CadenasValidas);
  //p2->PonPropiedades(false, 0, "Tipo de dato donde cargar las instancias", "Por omision es double");
  p2->PonPropiedades(false, 0, "Instance data type, BasicInstance by default (double)", 
                     "This field allows us to use less memory by using less precision");
  Parametros.push_back(p2);

  Parametro *p = new ParametroEntero();
  p->PonPropiedades(false, 0, "Column index with the multi-instance information (for MI datasets)"
                              "If negative, number of classes for multi-label datasets",
      "starts in 0");
  Parametros.push_back(p);

}
//---------------------------------------------------------------------------
void LoadDataset::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new LoadDataset(), "File");
}
//---------------------------------------------------------------------------
Mandato* LoadDataset::Ejecutar()
{
  Param(0)->Ejecutar();

  if (nd) delete nd;

  Instance *instanceBuilder = 0;
  if (NumParamsAsignados() > 1) {
    Param(1)->Ejecutar();
    instanceBuilder = Instance::NewInstanceByName(Param(1)->ComoCadena());
    cout << "Instance type is " << Param(1)->ComoCadena();
  }

  int attribute_index = -1;
  if (NumParamsAsignados() > 2) {
    int val = Param(2)->ComoEntero();
    if ( val < 0 ) {
        MultilabelInstance::SetNumMultilabels(-val);
        MultiobjectiveInstance::SetNumMultiobjectives(-val);
    }
    else {
      attribute_index = val;
    }
  }


  if ( attribute_index >= 0 ) {
    Param(2)->Ejecutar();
    nd = new MINomData((char*)Param(0)->ComoCadena().c_str(), attribute_index, instanceBuilder);
    nd->SetNTrain(nd->GetNTotal());
  }
  else {
    nd = Data::DataFromFile((char*)Param(0)->ComoCadena().c_str(), instanceBuilder);
    nd->SetNTrain(nd->GetNTotal());
  }

  if (instanceBuilder) delete instanceBuilder;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
GenerateWrongLabels::~GenerateWrongLabels()
{
  if(mat) delete mat;
}
//---------------------------------------------------------------------------
void GenerateWrongLabels::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new GenerateWrongLabels(), "Tools");
}
//---------------------------------------------------------------------------
void GenerateWrongLabels::CreateParams()
{
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
Mandato* GenerateWrongLabels::Ejecutar()
{
  Param(0)->Ejecutar();

  Data *data = (Data*)Param(0)->ComoDatos();

  int ndatos = data->GetNTotal();
  int nclasses = ((NomData*)data)->NumClass;

  if (mat) delete mat;
  mat = new Matriz(1, 2*ndatos);

  data->Scramble(0, ndatos-1);

  for(int i=0;i<ndatos;i++) {
    int new_class = (int) ((double)(nclasses-1)*rand()/(RAND_MAX+1.0));
    if (new_class>=data->GetDatClass(i)) new_class++;
    (*mat)[0][i*2] = data->GetDatIniPos(i);
    (*mat)[0][i*2+1] = new_class;
  }

  return this;
}
//---------------------------------------------------------------------------
string GenerateWrongLabels::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Corrected::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Corrected(), "File");
}
//---------------------------------------------------------------------------
void Corrected::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());
  Parametros.push_back(new ParametroData());

  Parametro *pfs = new ParametroFicheroEntrada();
  pfs->PonPropiedades(false, 0, "labels file name", "File name");
  Parametros.push_back(pfs);

  Parametro *p = new ParametroEntero();
  p->PonPropiedades(false, 0, "Row index for the row with the labels", "it starts in 0");
  Parametros.push_back(p);

  p = new ParametroNumero(0.0, 1.0);
  p->PonPropiedades(false, 0, "Percentage of instances that get their label changed", "");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* Corrected::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();
  Param(3)->Ejecutar();
  Param(4)->Ejecutar();

  Classifier *c = (Classifier*)Param(0)->ComoDatos();
  Data *data = (Data*)Param(1)->ComoDatos();
  string filename = Param(2)->ComoCadena();
  int i_fila = Param(3)->ComoEntero();
  int n_ejemplos = (int)(Param(4)->ComoNumero()*data->GetNTotal() + 0.5);

  FILE *f = fopen(filename.c_str(), "r");
  if (!f) return this;
  char linea[0x100000];
  int i = 0;
  while(!feof(f)) {
    fgets(linea, 0x100000-1, f);
    if (feof(f)) return this;
    if (i == i_fila) {
      break;
    }
    i++;
  }
  fclose(f);

  data->SortOn(data->GetNumVar()+data->GetIniPosIndex());
  int i_ej = atoi(strtok(linea, " \t\n"));
  int clase = atoi(strtok(0, " \t\n"));
  int ok = 0;
  for(int i = 0 ; i < n_ejemplos - 1 ; i++ )
  {
    if (c->Classify(i_ej) == data->GetDatClass(i_ej)) ok++;
    i_ej = atoi(strtok(0, " \t\n"));
    clase = atoi(strtok(0, " \t\n"));
  }
  if (c->Classify(i_ej) == data->GetDatClass(i_ej)) ok++;

  corr = (double)ok/n_ejemplos;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void LoadLabels::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new LoadLabels(), "File");
}
//---------------------------------------------------------------------------
void LoadLabels::CreateParams()
{
  Parametros.push_back(new ParametroData());

  Parametro *pfs = new ParametroFicheroEntrada();
  pfs->PonPropiedades(false, 0, "File name with the new labels", "File name");
  Parametros.push_back(pfs);

  Parametro *p = new ParametroEntero();
  p->PonPropiedades(false, 0, "Indice de la fila con las etiquetas", "indice empezando en 0");
  Parametros.push_back(p);

  p = new ParametroNumero(0.0, 1.0);
  p->PonPropiedades(false, 0, "Porcentaje de ejemplos a los que se cambia la etiqueta", "");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* LoadLabels::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();
  Param(3)->Ejecutar();

  Data *data = (Data*)Param(0)->ComoDatos();
  string filename = Param(1)->ComoCadena();
  int i_fila = Param(2)->ComoEntero();
  int n_ejemplos = (int)(Param(3)->ComoNumero()*data->GetNTotal() + 0.5);
cout << n_ejemplos << endl;
  FILE *f = fopen(filename.c_str(), "r");
  if (!f) return this;
  char linea[0x100000];
  int i = 0;
  while(!feof(f)) {
    fgets(linea, 0x100000-1, f);
    if (feof(f)) return this;
    if (i == i_fila) {
      break;
    }
    i++;
  }
  fclose(f);

  data->SortOn(data->GetNumVar()+data->GetIniPosIndex());
  int i_ej = atoi(strtok(linea, " \t\n"));
  int clase = atoi(strtok(0, " \t\n"));
  for(int i = 0 ; i < n_ejemplos - 1 ; i++ )
  {
    data->SetDatClass(i_ej, clase);
    i_ej = atoi(strtok(0, " \t\n"));
    clase = atoi(strtok(0, " \t\n"));
  }
  data->SetDatClass(i_ej, clase);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
FilterData::~FilterData()
{
  if(nd) delete nd;
}
//---------------------------------------------------------------------------
void FilterData::CreateParams()
{
  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("Dispersar");
  CadenasValidas.push_back("NaiveBayes");
  CadenasValidas.push_back("Clonar");
  CadenasValidas.push_back("NN");
  CadenasValidas.push_back("Margin");
  CadenasValidas.push_back("OOB");
  CadenasValidas.push_back("Wrong");
  CadenasValidas.push_back("Ok");
  CadenasValidas.push_back("Border");
  CadenasValidas.push_back("ClassNoise");
  CadenasValidas.push_back("Bootstrap");
  CadenasValidas.push_back("Ponderate");
  CadenasValidas.push_back("DeleteFeature");
  CadenasValidas.push_back("MultiInstanceData");
  CadenasValidas.push_back("Merge");
  CadenasValidas.push_back("RandomizeAttributeValues");
  CadenasValidas.push_back("NormalizeVectors");
  CadenasValidas.push_back("Paste");
  CadenasValidas.push_back("NormalizeMIWeights");
  CadenasValidas.push_back("SetClass");
  CadenasValidas.push_back("RemoveInstancesOfClass");
  CadenasValidas.push_back("GenerateBackgroundDistributionData_Breiman");
  Parametro *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(false, 0, "Filtro", "Filtro a aplicar");
  Parametros.push_back(p2);

  p2 = new ParametroData();
  p2->PonPropiedades(false, 0, "Datos a filtros", 
                                    "Datos sobre los que se aplica el filtro");
  Parametros.push_back(p2);

  p2 = new Parametro();
  p2->PonPropiedades(true, Tipo::TCualquiera(), "Parametro", 
                                           "Parametro del filtro (si aplica)");
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
void FilterData::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new FilterData(), "Tools");
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
Mandato* FilterData::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  NomData *data = (NomData*)Param(1)->ComoDatos();

  if (nd) {
    delete nd;
    nd = 0;
  }

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
    int ii, ff;

    ii = 0;
    ff = data->GetNTotal() - 1;

    if ( NumParamsAsignados() > 3 ) {
      ii = Param(2)->Ejecutar()->ComoEntero();
      ff = Param(3)->Ejecutar()->ComoEntero();
    }
    else if ( NumParamsAsignados() > 2 ) {
      ff = Param(2)->Ejecutar()->ComoEntero() - 1;
    }

    data->SortOn(data->GetNumVar() + Data::IniPosIndex, 0, data->GetNTotal() - 1 );

    nd = data->Clone(ii, ff);
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
  else if (what==11) {//Bootstrap
//    double prob = Param(2)->Ejecutar()->ComoNumero();
    nd = data->GetBootstrapSample();
  }
  else if (what==12) {//Ponderate
    nd = data->Clone();
    nd->EqualizeClassWeight();
  }
  else if (what==13) {//DeleteFeature
    int icol = Param(2)->Ejecutar()->ComoEntero();
    nd = data->Clone();
    nd->DeleteColumn(icol);
  }
  else if (what==14) {//MultiInstanceData
    int icol = Param(2)->Ejecutar()->ComoEntero();
    nd = NomData::GenerateMIDataSet((NomData*)data, icol);
cout << ((NomData*)nd)->NumClass << endl;
  }
  else if (what==15) {//Merge
    Data *d2 = (Data*)Param(2)->Ejecutar()->ComoDatos();
    data->AddData(d2);
    nd = data;
cout << ((NomData*)nd)->NumClass << endl;
cout << ((MINomData*)nd)->GetMITrain() << endl;
  }
  else if (what==16) {//RandomizeAttributeValues
    int icol = Param(2)->Ejecutar()->ComoEntero();
    nd = data->Clone();
    nd->PermuteAttributeValues(icol);
  }
  else if (what==17) {//NormalizeVectors
    nd = data->Clone();
cout << nd->GetNumVarOrd() << endl;
    for(int i = 0; i < nd->GetNTotal(); i++) {
      double w = 0.0;
      for(int j = 0; j < nd->GetNumVarOrd(); j++) {
        w += fabs(nd->GetValueVar(i, j));
      }
      for(int j = 0; j < nd->GetNumVarOrd(); j++) {
        if (w != 0.0) {
          nd->SetValueVar(i, j, nd->GetValueVar(i,j) / w );
        }
      }
    }
  }
  else if (what==18) {//Paste
    //Data *d2[128];
    Data *d2 = (Data*)Param(2)->Ejecutar()->ComoDatos();
    if (d2->GetNTotal() != data->GetNTotal()) {
      cout << "Different Number of instances: Nothing done." << endl;
    }
    else {
      bool free_input = NumParamsAsignados() > 3 ? Param(3)->Ejecutar()->ComoBooleano() : false;
      //idatum = 1;
      //for(int i = 4; i < NumParamsAsignados(); i++) {
      //  d2[idatum++] = (Data*)Param(i)->Ejecutar()->ComoDatos();
      //}
      if (data->GetNumVarOrd()+1 != data->GetNumVar() || 
            d2->GetNumVarOrd()+1 !=   d2->GetNumVar()  ) {
        cout << "Warning: only pasting ordinal attributes" << endl;
      }
      data->OriginalOrder();
      d2->OriginalOrder();
      MINomData *mind = dynamic_cast<MINomData*>(data);
      int pasos[5];
      for(int i = 0; i < 5; i++) pasos[i] = (data->GetNTotal()*(i+1))/5;
      if (mind) {
         nd = new MINomData(pasos[0], data->GetNumVarOrd() +
                                            d2->GetNumVarOrd() + 1);
      }
      else {
         nd = new NomData(pasos[0], data->GetNumVarOrd() +
                                            d2->GetNumVarOrd() + 1);
      }
      for(int i = 0; i < ((NomData*)data)->NumClass; i++) {
        ((NomData*)nd)->AddClass( data->GetTerm(i, data->GetNumVarNom()) );
      }
      ((NomData*)nd)->init();
      int kk = data->GetNumVarOrd();
      int kk2 = d2->GetNumVarOrd();
      int kk3 = data->GetNumVar();

      int i_first_wrong = -1;
      int n_wrong = 0;
      for (int i = 0, io = data->GetNTotal() - 1, ip = 0; i < pasos[4]; i++, io--) {
        if (i == nd->GetNTotal() ) {
          if (free_input) {
            data->redim(io+1);
            d2->redim(io+1);
            cout << "reducing " << io+1 << endl;
          }
          ip++;
          nd->redim(pasos[ip]);
        }
        for (int j = 0; j < kk; j++) {
          nd->SetValueVar(i, j, data->GetValueVar(io, j));
        }
        for (int j = 0; j < kk2; j++) {
          nd->SetValueVar(i, j + kk, d2->GetValueVar(io, j));
        }
        for (int j = 0; j < 7; j++) {
          nd->SetValueVar(i, 1 + j + kk + kk2, data->GetValueVar(io, kk3 + j));
        }
        nd->SetValueVar(i, kk2 + kk, data->GetValueVar(io, data->GetNumVar()-1));
        if (data->GetDatClass(io) != d2->GetDatClass(io)) { 
          if ( n_wrong == 0 ) i_first_wrong = i+1;
          n_wrong++;
        }
      }
      if ( n_wrong > 0 ) {
        cout << "Warning: classes of " << n_wrong << " pasting vectors do ";
        cout << "no coincide. First is: " << i_first_wrong << endl;
      }
      nd->OriginalOrder();
      nd->SetNTrain(pasos[4]);
            cout << "to " << pasos[4] << endl;
      if (free_input) {
        delete data;
        delete d2;
      }
    }
  }
  else if (what==19) {// NormalizeMIWeights
    MINomData *mind = dynamic_cast<MINomData*>(data);
    if (!mind) {
      cout << "Data is not Multi-instance data: Nothing done" << endl;
      return this;
    }
    mind->SortByGroup();
    int igroup = mind->GetDatGroup(0);
    for(int i = 0, ii = 0; i < mind->GetNTotal(); ) {
      while (ii < mind->GetNTotal() && mind->GetDatGroup(ii) == igroup) {
        ii++;
      }
      igroup = ii < mind->GetNTotal() ? mind->GetDatGroup(ii) : 0;
      double mi_weight = 1.0*mind->GetNTotal()/((ii-i)*mind->GetMITotal());
      for ( ; i < ii; i++) {
        mind->SetDatWeight(i, mi_weight);
      }
    }
    cout << endl;
  }
  else if (what==20) {//SetClass
    if ( NumParamsAsignados() > 3 ) {
      Data *d2 = (Data*)Param(2)->Ejecutar()->ComoDatos();
      int attribute_idx = Param(3)->Ejecutar()->ComoEntero();
      data->SubstituteClassLabels(d2, attribute_idx);
    }
    else {
      int class_idx = Param(2)->Ejecutar()->ComoEntero();
      int N = data->GetNTotal();
      for(int i = 0; i < N; i++) {
        data->SetDatClass(i, class_idx);
      }
    }
  }
  else if (what==21) {//RemoveInstancesOfClass
    int class_index = Param(2)->Ejecutar()->ComoEntero();
    data->RemoveInstancesOfClass(class_index);
  }
  else if (what==22) {//GenerateBackgroundDistributionData_Breiman
    int N = NumParamsAsignados() > 2 ? Param(2)->Ejecutar()->ComoEntero() : data->GetNTotal();
    nd = data->Clone();
    nd->redim(N);
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < data->GetNumVar()-1; j++) {
        int value_index = (int) ((double)(data->GetNTotal())*rand()/(RAND_MAX+1.0));
        nd->SetValueVar(i, j, data->GetValueVar(value_index, j));
      }
    }
  }

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void SetSplitCriterium::CreateParams()
{
  Parametro *p2 = new ParametroData();
  p2->PonPropiedades(false, 0, "Dataset in which the new split criterium is set", 
                                    "Datos sobre los que se aplica el split");
  Parametros.push_back(p2);

  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("GiniCriterium");
  CadenasValidas.push_back("MSECriterium");
  CadenasValidas.push_back("MultipleMSECriterium");
  CadenasValidas.push_back("ComboCriterium");
  CadenasValidas.push_back("MultilabelGiniCriterium");
  CadenasValidas.push_back("BoostedGiniCriterium");
  p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(false, 0, "Split name", "Split to be applied");
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
void SetSplitCriterium::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new SetSplitCriterium(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* SetSplitCriterium::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  nd = (NomData*)Param(0)->ComoDatos();
  nd->SetSplitCriterium(Param(1)->ComoCadena());

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Classify::~Classify()
{
  if(mat) delete mat;
}
//---------------------------------------------------------------------------
void Classify::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Classify(), "Classifiers");
}
//---------------------------------------------------------------------------
void Classify::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());
  Parametros.push_back(new ParametroData());
  Parametros.push_back(new ParametroBooleano());
  Parametros.push_back(new ParametroBooleano());
  Parametros.push_back(new ParametroBooleano());
}
//---------------------------------------------------------------------------
Mandato* Classify::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Classifier *c = (Classifier *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  c->SetData(dat);
  dat->OriginalOrder();

  if (mat) delete mat;

  bool WriteClass = NumParamsAsignados()>=4 ? Param(3)->Ejecutar()->ComoBooleano() : false;

  if (NumParamsAsignados()>=3 && Param(2)->Ejecutar()->ComoBooleano()) {
    Ensemble *e = (Ensemble*)c;

    bool OOB = NumParamsAsignados()>=5 ? Param(4)->Ejecutar()->ComoBooleano() : false;

    mat = new Matriz(dat->GetNTotal(), e->GetClassifiersToUse() + (WriteClass ? 1 : 0));
    for(int i = 0; i < e->GetClassifiersToUse(); i++) {
      c = e->GetClassifier(i);
      for(int j = 0; j < dat->GetNTotal(); j++) {
        (*mat)[j][i] = OOB && c->UsedInOriginalTrainingData(j) ? -1 : c->Classify(j);
      }
    }
    if (WriteClass) {
      for(int j = 0; j < dat->GetNTotal(); j++) {
        (*mat)[j][e->GetClassifiersToUse()] = dat->GetDatClass(j);
      }
    }
  }
  else {
    mat = new Matriz(dat->GetNTotal(), 1 + (WriteClass ? 1 : 0));

    for(int i=0;i<dat->GetNTotal();i++) {
      (*mat)[i][0] = c->Classify(i);
      if (WriteClass) {
        (*mat)[i][1] = dat->GetDatClass(i);
      }
    }
  }

  return this;
}
//---------------------------------------------------------------------------
string Classify::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Mrse::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Mrse(), "Classifiers");
}
//---------------------------------------------------------------------------
void Mrse::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
Mandato* Mrse::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Classifier *c = (Classifier *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  c->SetData(dat);
  dat->OriginalOrder();

  tot = 0.0;
  for(int i=0;i<dat->GetNTotal();i++) {
    double val = c->Average(i);
    tot += (dat->GetDatClass(i)-val) * (dat->GetDatClass(i)-val);
    
  }

  tot = sqrt(tot/dat->GetNTotal());

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Error::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Error(), "Classifiers");
}
//---------------------------------------------------------------------------
void Error::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
Mandato* Error::Ejecutar()
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
void ResubstituteNodeStats::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ResubstituteNodeStats(), "Classifiers");
}
//---------------------------------------------------------------------------
void ResubstituteNodeStats::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());
  Parametros.push_back(new ParametroData());
  Parametros.push_back(new ParametroBooleano());
  Parametros.push_back(new ParametroBooleano());
}
//---------------------------------------------------------------------------
Mandato* ResubstituteNodeStats::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  bool WithOOB = NumParamsAsignados() > 2 ? Param(2)->Ejecutar()->ComoBooleano() : false;
  bool AddInfo = NumParamsAsignados() > 3 ? Param(3)->Ejecutar()->ComoBooleano() : false;

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  ens->SetData(dat);
  int ninst = dat->GetNTotal();
  for (int i = 0; i < ens->GetClassifiersToUse(); i++) {
    DecisionTree *dt = (DecisionTree*)ens->GetClassifier(i);
    if (WithOOB) {
      ninst = 0;
      for (int j = 0; j < dat->GetNTotal(); j++) {
        int iori = dat->GetDatIniPos(j);
        dat->GetInstance(i)(dat->GetNumVar(), dt->UsedInOriginalTrainingData(iori) ? 1.0 : 0.0);
        if (!dt->UsedInOriginalTrainingData(iori)) ninst++;
      }
      dat->SortOn(dat->GetNumVar());
    }
    dt->GetTree()->RecomputeNodeStats(dat, 0, ninst - 1, AddInfo);
  }


  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void EvaluateClassifier::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new EvaluateClassifier(), "Classifiers");
}
//---------------------------------------------------------------------------
void EvaluateClassifier::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());
  Parametros.push_back(new ParametroData());
  Parametros.push_back(new ParametroCadena());
  Parametros.push_back(new ParametroCadena());
}
//---------------------------------------------------------------------------
Mandato* EvaluateClassifier::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Classifier *c = (Classifier *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  string EvaluationName = "";
  if (NumParamsAsignados() > 2) {
    Param(2)->Ejecutar();
    EvaluationName = Param(2)->ComoCadena();
  }

  string Params = "";
  if (NumParamsAsignados() > 3) {
    Param(3)->Ejecutar();
    Params = Param(3)->ComoCadena();
  }

  Evaluation *eval = c->Evaluate(dat, EvaluationName, Params);
  output = eval->ToString();
  delete eval;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
SequentialError::~SequentialError()
{
  if(err) delete err;
}
void SequentialError::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new SequentialError(), "Classifiers");
}
//---------------------------------------------------------------------------
void SequentialError::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
  Parametro *pfs = new ParametroFicheroSalida();
  pfs->PonPropiedades(true, 0, "Output file",
                "File name where the individual classifications are stored");
  Parametros.push_back(pfs);
}
//---------------------------------------------------------------------------
Mandato* SequentialError::Ejecutar()
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

  

  if (clas) {
//    FILE *f=fopen(nom_fich.c_str(), "a");
//    for (unsigned i=0; i<inderrs->size(); i++) {
//      fprintf(f, "%g\t", (*inderrs)[i]);
//    }
//    fprintf(f, "\n");
//    fclose(f);
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
string SequentialError::ComoCadena()
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
ClassificationMatrix::~ClassificationMatrix()
{
  if (mat) delete mat;
}
void ClassificationMatrix::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ClassificationMatrix(), "Classifiers");
}
//---------------------------------------------------------------------------
void ClassificationMatrix::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
  Parametro *p = new ParametroBooleano();
  p->PonPropiedades(true, 0, "Sequential/individual",
                             "Da la clase/error acumulado o no del conjunto");
  Parametros.push_back(p);
  p = new ParametroBooleano();
  p->PonPropiedades(true, 0, "Clase/error",
                             "Da la clase/error acumulado o no del conjunto");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* ClassificationMatrix::Ejecutar()
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
string ClassificationMatrix::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
ConfusionMatrix::~ConfusionMatrix()
{
  if (mat) delete mat;
}
void ConfusionMatrix::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ConfusionMatrix(), "Classifiers");
}
//---------------------------------------------------------------------------
void ConfusionMatrix::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
Mandato* ConfusionMatrix::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Classifier *c = (Classifier *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  if (mat) delete mat;
  mat = new Matriz(((NomData*)dat)->NumClass);

  c->SetData(dat);

  c->Error(0, dat->GetNTotal() - 1);
 
  string n = GenerateProcNameFileName("predicted_class",".csv");
  FILE *f=fopen(n.c_str(),"a");
  for(int i = 0; i < (int)c->rc.size(); i++) {
    (*mat)[c->rc[i]][c->pc[i]] += 1.0;
    fprintf(f, "%d ", c->pc[i]);
  }
  fprintf(f, "\n");
  fclose(f);

//  for(int i=0;i<dat->GetNTotal();i++)
  //  (*mat)[dat->GetDatClass(i)][c->Classify(i)] += 1.0;

  return this;
}
//---------------------------------------------------------------------------
string ConfusionMatrix::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Margin::~Margin()
{
  if (mat) delete mat;
}
void Margin::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Margin(), "Classifiers");
}
//---------------------------------------------------------------------------
void Margin::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
  Parametro *p = new ParametroEntero(2, 2100000000);
  p->PonPropiedades(true, 0, "Divisiones", "");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* Margin::Ejecutar()
{
  int ini, fin;

  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data*)Param(1)->ComoDatos();
  ens->SetData(dat);
  ini = 0;
  fin = dat->GetNTotal()-1;

  int div = ens->GetClassifiersToUse();
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
//  MarginDistribution md(ens, dat, true);
  MarginDistribution md2(ens, dat, false);
  mat = new Matriz(1, div+1);
  double x = -1.0;
  for (int i=0;i<=div;i++) {
//    (*mat)[0][i] = md.fm(x);
    (*mat)[0][i] = md2.fm(x);
    x+=2.0/div;
  }

  bool acum = true;  
  if (!acum) {
    for (int i=div-1;i>=0;i--) {
      (*mat)[0][i+1] -= (*mat)[0][i];
    }
  }

  return this;
}
//---------------------------------------------------------------------------
string Margin::ComoCadena()
{
  std::ostringstream buf;
  mat->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
MarginOfInstances::~MarginOfInstances()
{
  if (mat) delete mat;
}
void MarginOfInstances::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new MarginOfInstances(), "Classifiers");
}
//---------------------------------------------------------------------------
void MarginOfInstances::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
  Parametros.push_back(new ParametroEntero());
}
//---------------------------------------------------------------------------
Mandato* MarginOfInstances::Ejecutar()
{

  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data*)Param(1)->ComoDatos();
  ens->SetData(dat);
  int ndatos = dat->GetNTotal();

  int clase = -2;
  if ( NumParamsAsignados() > 2 ) {
    Param(2)->Ejecutar();
    clase = Param(2)->ComoEntero();
  }

  dat->OriginalOrder();

  if (mat) delete mat;
  mat = new Matriz(ndatos, 1);

  for(int i = 0; i < ndatos; i++) {
    (*mat)[i][0] = Classifier::Margin(ens->Distribution(i), clase == -2 ? dat->GetDatClass(i) : clase);
  }

  return this;
}
//---------------------------------------------------------------------------
string MarginOfInstances::ComoCadena()
{
  string cad = "";
  string sep = "";
  for(int i = 0; i < mat->filas(); i++) {
    cad = cad + sep + Mandato::ACad((*mat)[i][0]);
    sep = "\t";
  }
  return cad;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Certainty::~Certainty()
{
  if (mat) delete mat;
}
void Certainty::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Certainty(), "Classifiers");
}
//---------------------------------------------------------------------------
void Certainty::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
Mandato* Certainty::Ejecutar()
{
  int ini, fin;

  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();

  if ( NumParamsAsignados() > 3 ) {
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
  mat = new Matriz(1, fin-ini+1);
  for(int i=ini; i<=fin;i++) {
    int clase;
    double cer = ens->ClassificationCertainty(i, clase/*, ExcluirTrain*/);
    if (clase!=dat->GetDatClass(i)) cer = 1.0 - cer;
    (*mat)[0][i-ini] = cer;
  }

  return this;
}
//---------------------------------------------------------------------------
string Certainty::ComoCadena()
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
MarginSum::~MarginSum()
{
  if (mat) delete mat;
}
void MarginSum::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new MarginSum(), "Classifiers");
}
//---------------------------------------------------------------------------
void MarginSum::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
  Parametro *p = new ParametroNumero(-1.0, 1.0);
//  p->SetPorOmision("0.0");
  p->PonPropiedades(true, 0, "Divisiones", "");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* MarginSum::Ejecutar()
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
string MarginSum::ComoCadena()
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

FILE *ff=fopen("dp.csv", "w");
for(int j=0;j<nDatos;j++) {
int clase = dat->GetDatClass(j); 
fprintf(ff, "%g\n", cVal[j][clase]/cVal[j][nClases]);
}
fclose(ff);

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
  CadenasValidas.push_back("MaxTreeDepth");
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
  else if (prop.compare("MaxTreeDepth")==0) {
    Node::max_depth = Param(2)->ComoEntero();
  }
  else if (prop.compare("K")==0) {
    CARTTree::K = Param(2)->ComoEntero();
    if (CARTTree::K == -2) Node::max_depth = 1;
    else if (Node::max_depth == 1) Node::max_depth = -100;
  }
  else if (prop.compare("SetDefaultClassifiers")==0) {
    Boosting *boo = (Boosting*)Param(0)->ComoDatos();
    boo->SetDefaultClassifiers();
    int max = Param(2)->ComoEntero();
    if (boo->GetClassifiersToUse() > max)
      boo->SetClassifiersToUse(max);
  }
  else if (prop.compare("SetGlobalOpt")==0) {
    int opt = (Param(2)->ComoCadena().c_str())[0];
    char val[1024];
    strcpy(val, Param(2)->ComoCadena().c_str()+1);
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
  p->PonPropiedades(false, 0, "Number of instances for train", "Ejemplos a utilizar como train");
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
void GeneratePartitions::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new GeneratePartitions(), "Classifiers");
}
//---------------------------------------------------------------------------
void GeneratePartitions::CreateParams()
{
  Parametros.push_back(new ParametroData());
  Parametro *p = new ParametroNumero(-1, 2100000000);
  p->PonPropiedades(false, 0, "Train set size", "if 0 this performs a n-fold cross validation, where n is given by the next parameter"
                              "if in range (0, 1) then uses this value as a % de the total number of instances");
  Parametros.push_back(p);
  p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Numbar of partitions", "Particiones a realizar");
  Parametros.push_back(p);
  ParametroCadena *pc = new ParametroCadena();
  pc->PonPropiedades(false, 0, "File prefix", "Prefix for the output files to be generated");
  Parametros.push_back(pc);
}
//---------------------------------------------------------------------------
Mandato* GeneratePartitions::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();
  Param(3)->Ejecutar();

  Data *data = (Data*)Param(0)->ComoDatos();
  double ndatos = Param(1)->ComoNumero();
  int nparti =  Param(2)->ComoEntero();
  string pre = Param(3)->ComoCadena();

  int opts = 1;
  if (NumParamsAsignados()>4) {
    Param(4)->Ejecutar();
    opts = Param(4)->ComoEntero();
  }

  if ( ndatos > 0 && ndatos < 1 )  {
    ndatos = (int)(0.5 + ndatos*data->GetNTotal());
  }
  DoCreaConjuntosTrainAleat( data, ndatos, nparti, pre, opts);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void SaveDataset::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new SaveDataset(), "File");
}
//---------------------------------------------------------------------------
void SaveDataset::CreateParams()
{
  Parametros.push_back(new ParametroData());
  Parametro *pfs = new ParametroFicheroSalida();
  pfs->PonPropiedades(false, 0, "Output file name", "Nombre del fichero");
  Parametros.push_back(pfs);
}
//---------------------------------------------------------------------------
Mandato* SaveDataset::Ejecutar()
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
void SaveClassifier::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new SaveClassifier(), "File");
}
//---------------------------------------------------------------------------
void SaveClassifier::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());

  string fil = string("Any|*.*|") + "Fast format (*.clf)|*.clf";
  Parametro *pfs = new ParametroFicheroSalida(fil);
  pfs->PonPropiedades(false, 0, "Output file name", "Nombre del fichero");
  Parametros.push_back(pfs);
}
//---------------------------------------------------------------------------
Mandato* SaveClassifier::Ejecutar()
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
void ClassifierInfo::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ClassifierInfo(), "Classifiers");
}
//---------------------------------------------------------------------------
void ClassifierInfo::CreateParams()
{
  Parametros.push_back(new ParametroClasificador());
  Parametro *p = new ParametroEntero(0, 2100000000);
  p->PonPropiedades(false, 0, "Info level", "");
  Parametros.push_back(p);
}//---------------------------------------------------------------------------
Mandato* ClassifierInfo::Ejecutar()
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
BuildIGPEnsemble::~BuildIGPEnsemble()
{
  //printf("%x\n", Clasif);
  if (Clasif) delete Clasif; 
}
//---------------------------------------------------------------------------
void BuildIGPEnsemble::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new BuildIGPEnsemble(), "Classifiers");
}
//---------------------------------------------------------------------------
void BuildIGPEnsemble::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Ejemplos",
                               "Nombre de la base de datos con los ejemplos"));
  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Numero de arboles IGP", "Numero que compondran el conjunto");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* BuildIGPEnsemble::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();

  if (Clasif) delete Clasif;
  GroupBagging *b = new GroupBagging();

  int ndatos = nd->GetNTotal();
  nd->SetNTrain(ndatos);
  for(int i=0;i<nArb;i++) {
    b->AddClassifier(new GPTree());
  }

  ProgresoConsola pc;
  b->Build(nd, &pc);

  Clasif = b;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void BuildBaggingCART::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new BuildBaggingCART(), "Classifiers");
}
//---------------------------------------------------------------------------
void BuildBaggingCART::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Training dataset",
                             "Nombre de la base de datos con los ejemplos"));

  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Number of CART trees to create", 
                                        "Numero que compondran el conjunto");
  Parametros.push_back(p);

  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("Stumps");
  CadenasValidas.push_back("Pruned");
  CadenasValidas.push_back("Unpruned");
  ParametroCadena *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(true, 0, "Prune trees?", "By default: Pruned");
  Parametros.push_back(p2);

  p = new ParametroNumero(0.0, 10000.0);
  p->PonPropiedades(true, 0, "Bootstrap size [0.0;...)", "By default 1 (=100%)");
  Parametros.push_back(p);

  p = new ParametroBooleano(true, "With Replacement?", 
                             "Bootstrap with or without replacement, by default with replacement");
  Parametros.push_back(p);

  p = new ParametroBooleano(false, "Equalize classes", 
                 "By default no equalization of the number of instances of each class is performed");
  Parametros.push_back(p);

//  p = new ParametroBooleano(true, "Podar arboles",
//                                           "Por omision: El arbol se poda");
//  Parametros.push_back(p);

  p = new ParametroBooleano(true, "Multisplits?",
                            "By default: unidimensional splits");
  Parametros.push_back(p);

  p = new ParametroBooleano(true, "SE 0 rule?",
                                             "By default: SE 1.0");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* BuildBaggingCART::Ejecutar()
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

  bool WithReplacement = true;
  if (NumParamsAsignados()>4) {
    Param(4)->Ejecutar();
    WithReplacement = Param(4)->ComoBooleano();
  }

  bool EqualizeClasses = false;
  if (NumParamsAsignados()>5) {
    Param(5)->Ejecutar();
    EqualizeClasses = Param(5)->ComoBooleano();
  }

  bool MultiSplits = false;
  if (NumParamsAsignados()>6) {
    Param(6)->Ejecutar();
    MultiSplits = Param(6)->ComoBooleano();
  }

  bool SE_0 = false;
  if (NumParamsAsignados()>7) {
    Param(7)->Ejecutar();
    SE_0 = Param(7)->ComoBooleano();
  }

  nd->SetNTrain(nd->GetNTotal());

  if (Clasif) delete Clasif;

  MINomData *mind = dynamic_cast<MINomData*>(nd);
  Ensemble *b;
  if (mind) {
    b = new MIBagging();
  }
  else {
    b = new Bagging((int)(0.5 + TamanoMuestreo*nd->GetNTrain()), 
                                            WithReplacement, EqualizeClasses);
  }

  for(int i=0;i<nArb;i++) {
    if (Podar==1) {
      b->AddClassifier(new DecisionStump());
    }
    else if (Podar==2) {
      b->AddClassifier(new CARTTree(0, 1, true, MultiSplits, SE_0));
    }
    else if (Podar==3) {
      //b->AddClassifier(new CARTTree(0, 1, false, MultiSplits, SE_0));
      b->AddClassifier(new CARTTree(0, 1, false, MultiSplits, SE_0));
    }
  }

  ProgresoConsola pc("tree");
  b->Build(nd, &pc);

  Clasif = b;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
BuildBaggingC45::~BuildBaggingC45()
{
  //printf("%x\n", Clasif);
  if (Clasif) delete Clasif; 
}
//---------------------------------------------------------------------------
void BuildBaggingC45::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new BuildBaggingC45(), "Classifiers");
}
//---------------------------------------------------------------------------
void BuildBaggingC45::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Dataset variable",
                               "Nombre de la base de datos con los ejemplos"));
  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Number of C45 trees", "Numero que compondran el conjunto");
  Parametros.push_back(p);
  p = new ParametroBooleano(true, "Pruned?",
                                            "By default: true");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* BuildBaggingC45::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();

  if (Clasif) delete Clasif;

  Ensemble *b;
  MINomData *mind = dynamic_cast<MINomData*>(nd);

  if (mind) {
    b = new MIBagging();
  }
  else {
    b = new Bagging();
  }

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

  ProgresoConsola pc("tree");
  b->Build(nd, &pc);

  Clasif = b;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
BuildBaggingNNet::~BuildBaggingNNet()
{
  //printf("%x\n", Clasif);
  if (Clasif) delete Clasif; 
}
//---------------------------------------------------------------------------
void BuildBaggingNNet::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new BuildBaggingNNet(), "Classifiers");
}
//---------------------------------------------------------------------------
void BuildBaggingNNet::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Dataset variable",
                               "Nombre de la base de datos con los ejemplos"));
  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Number of Neural Networks of the ensemble", 
                                          "Numero que compondran el conjunto");
  Parametros.push_back(p);
  
  p = new ParametroEntero(1, 32000);
  p->PonPropiedades(false, 0, "Number of nodes in the hidden layer", 
                                          "Pues eso");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* BuildBaggingNNet::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nNN = Param(1)->ComoEntero();
  int nnodos = Param(2)->ComoEntero();

  int epocas = 1000;
  if (NumParamsAsignados()>3) {
    Param(3)->Ejecutar();
    epocas = Param(3)->ComoEntero();
  }

  if (Clasif) delete Clasif;
  Bagging *b = new Bagging();

  int ndatos = nd->GetNTotal();
  nd->SetNTrain(ndatos);
  for(int i=0;i<nNN;i++) {
    b->AddClassifier(new NNet(nnodos, epocas));
  }

  ProgresoConsola pc("nn");
  b->Build(nd, &pc);

  Clasif = b;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void BuildRandomForest::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new BuildRandomForest(), "Classifiers");
}
//---------------------------------------------------------------------------
void BuildRandomForest::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Dataset variable",
                             "Nombre de la base de datos con los ejemplos"));

  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Number of ranom trees in the ensemble", 
                                        "Numero que compondran el conjunto");
  Parametros.push_back(p);

  p = new ParametroEntero(0, 2100000000);
  p->PonPropiedades(true, 0, "Number of random attributes in the splits", 
                  "Randomly selected in each node: By default: log2(number_of_attributes)");
  Parametros.push_back(p);

  p = new ParametroBooleano(true, "Multisplits?", "By default: false");
  Parametros.push_back(p);

  p = new ParametroEntero(-210000000, 2100000000);
  p->PonPropiedades(true, 0, "if positive then 'Minimum number of instances "
                             "per leaf'; if negative 'minimum number of inst"
                             "ances in node to continue spplitting'", "By default: 1");
  Parametros.push_back(p);

  p = new ParametroNumero(0.0, 10000.0);
  p->PonPropiedades(true, 0, "Bootstrap sample size [0.0;...)", "By default: 1 (=100%)");
  Parametros.push_back(p);

  p = new ParametroBooleano(false, "Equalize classes?", 
                 "By default: false (no equalization of class proportion is performed)");
  Parametros.push_back(p);

  p = new ParametroEntero(0, 2100000000);
  p->PonPropiedades(false, 0, "OutputLevel [0]", 
        "Only for MIBagging, put 4 to get votes for all instances and nodes to a file");
  Parametros.push_back(p);

  p = new ParametroEntero(0, 2100000000);
  p->PonPropiedades(false, 0, "Attributes from where to randomize [all]", 
        "Attribute range from where to select the random atributes in the splits. From 0 to the indicated value-1");
  Parametros.push_back(p);

}
//---------------------------------------------------------------------------
Mandato* BuildRandomForest::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();

  if (Clasif) delete Clasif;

  Ensemble *b;
  MINomData *mind = dynamic_cast<MINomData*>(nd);

  int NumAtributos = (int)(1 + log(nd->GetNumVar()-1)/log(2.0));
  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    if (Param(2)->ComoEntero() > 0) {
      NumAtributos = Param(2)->ComoEntero();
    }
  }
cout << "NAttribs " << NumAtributos << endl;
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
    
  double TamanoMuestreo = 1.0;
  if (NumParamsAsignados()>5) {
    Param(5)->Ejecutar();
    TamanoMuestreo = Param(5)->ComoNumero();
  }

  bool EqualizeClasses = false;
  if (NumParamsAsignados()>6) {
    Param(6)->Ejecutar();
    EqualizeClasses = Param(6)->ComoBooleano();
  }

  int OutputLevel = 0;
  if (NumParamsAsignados()>7) {
    Param(7)->Ejecutar();
    OutputLevel = Param(7)->ComoEntero();
  }

  int AttsToRandomize = -1;
  if (NumParamsAsignados()>8) {
    Param(8)->Ejecutar();
    AttsToRandomize = Param(8)->ComoEntero();
  }

  if (mind) {
    cout << "MI ensemble" << endl;
cout << ((MINomData*)nd)->GetMITrain() << endl;
    b = new MIBagging(((int)(0.5 + TamanoMuestreo*((MINomData*)nd)->GetMITrain())));
    ((MIBagging*)b)->SetOutputLevel(OutputLevel);
  }
  else {
    b = new Bagging(((int)(0.5 + TamanoMuestreo*nd->GetNTrain())), true, EqualizeClasses);
  }

  nd->SetNTrain(nd->GetNTotal());
  for(int i=0;i<nArb;i++) {
    b->AddClassifier(new RandomForestTree(0, EjemplosPorHoja, NumAtributos, 
                                                 MultiSplits, AttsToRandomize));
  }

  ProgresoConsola pc("tree");
  b->Build(nd, &pc);

  Clasif = b;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/*void BuildBaggingComMod::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new BuildBaggingComMod());
}
//---------------------------------------------------------------------------
Mandato* BuildBaggingComMod::Ejecutar()
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
void BuildBoosting::DarDeAlta()
{
  ifns().AnadirInterfazFuncion( new BuildBoosting(), "Classifiers");
}
//---------------------------------------------------------------------------
void BuildBoosting::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "DAtaset variable",
                               "Nombre de la base de datos con los ejemplos"));

  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Number of CART trees", "Numero de arboles que compondran el conjunto");
  Parametros.push_back(p);

//  p = new ParametroBooleano(true, "Podar arboles", "Por omision: El arbol se poda");
  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("Stumps");
  CadenasValidas.push_back("Pruned");
  CadenasValidas.push_back("Unpruned");
  ParametroCadena *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(true, 0, "Prune trees?", "By default: Pruned");
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
Mandato* BuildBoosting::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();

  if (Clasif) delete Clasif;

  Ensemble *b;
  MINomData *mind = dynamic_cast<MINomData*>(nd);

  if (mind) {
    b = new MIBoosting();
  }
  else {
    b = new Boosting();
  }

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

  ProgresoConsola pc("tree");
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
void BuildBoostingC45::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new BuildBoostingC45(), "Classifiers");
}
//---------------------------------------------------------------------------
void BuildBoostingC45::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Dataset variable",
                               "Nombre de la base de datos con los ejemplos"));
  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Number of C45 trees", "Numero que compondran el conjunto");
  Parametros.push_back(p);
  p = new ParametroBooleano(true, "Pruned?", "By default: true");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* BuildBoostingC45::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();
  int nArb = Param(1)->ComoEntero();

  if (Clasif) delete Clasif;

  Ensemble *b;
  MINomData *mind = dynamic_cast<MINomData*>(nd);

  if (mind) {
    b = new MIBoosting();
  }
  else {  
    b = new Boosting();
  } 

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

  ProgresoConsola pc("tree");
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
void BuildClassSwitching::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new BuildClassSwitching(), "Classifiers");
}

//---------------------------------------------------------------------------
void BuildClassSwitching::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Dataset variable",
                               "Nombre de la base de datos con los ejemplos"));

  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Number of treess (C45/CART) or neural networks",
                                      "Number of classifiers of the ensemble");
  Parametros.push_back(p);

  p = new ParametroNumero(-1.0, 1.0);
  p->PonPropiedades(false, 0, "Class switching parameter [-1.0;1.0]", 
     "If it is negative the diversity is generated using Kuncheva's (not very useful) method");
  Parametros.push_back(p);

  p = new ParametroEntero(-1, 1);
  p->PonPropiedades(false, 0, "Architecture for the ensemble", 
                                                        "Trees=1 NNetworks=-1");

  Parametros.push_back(p);

  p = new ParametroEntero(-1, 1024);
  p->PonPropiedades(false, 0, "Mixed parameter.", 
    "If trees then (C45=1,CART=-1) else if NN then this is the number of nodes in the hidden layer");
  Parametros.push_back(p);
}

//---------------------------------------------------------------------------
Mandato* BuildClassSwitching::Ejecutar()
{
  
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();

  //se cogen los datos
  Data* nd = (Data*)Param(0)->ComoDatos();

  //coger numero de clasificadores
  int nClasificadores = Param(1)->ComoEntero();
  
  //coger swiching
  bool GenerateOutputsAsKuncheva = false;
  double err = Param(2)->ComoNumero();
  if (err<0.0) {
    err = -err;
    GenerateOutputsAsKuncheva = true;
  }

  //coger si se utilzan arboles (>0) o redes (<0)
  int arboles = 1;
  if (NumParamsAsignados()>3) {
    Param(3)->Ejecutar();
    arboles = Param(3)->ComoEntero();
  }
  
  int TipoDeArbol = 1;
  int nNeuronas = 0;

  //se coge o bien el tipo de arboles o bien el numero de neuronas
  if (NumParamsAsignados()>4) {
    Param(4)->Ejecutar();
    nNeuronas = TipoDeArbol = Param(4)->ComoEntero();
  }

  bool AddC45Trees = arboles > 0 && TipoDeArbol > 0 ? true : false;

  int epocas = 1000;
  if (NumParamsAsignados()>5) {
    Param(5)->Ejecutar();
    epocas = Param(5)->ComoEntero();
  }

  bool UseDummyCls = false;
/*  if (NumParamsAsignados()>5) {
    Param(5)->Ejecutar();
    UseDummyCls = Param(5)->ComoBooleano();
  }*/

  double f_incr = 0.0;
/*  if (NumParamsAsignados()>6) {
    Param(6)->Ejecutar();
    f_incr = Param(6)->ComoNumero();
  }*/

  Data *test = 0;
/*  if (NumParamsAsignados()>7) {
    Param(7)->Ejecutar();
    test = (Data*)Param(7)->ComoDatos();
  }*/

  if (Clasif) delete Clasif;
  IndepClassifiers *b = new IndepClassifiers(err, UseDummyCls, f_incr);

  bool Podar = false;
  int ndatos = nd->GetNTotal();
  nd->SetNTrain(ndatos);
  
  if (AddC45Trees) {
    C45Tree::SetGlobalOpt('m', (char*)"1");
    C45Tree::SetGlobalOpt('d', (char*)"0");
  }

  //se crean los conjuntos de clasificadores
  for(int i=0;i<nClasificadores;i++) {
    Classifier *c;
    c =  arboles < 0 ? (Classifier *) new NNet(nNeuronas, epocas) :
        (AddC45Trees ? (Classifier *) new C45Tree(Podar)  :
                     (Classifier *) new CARTTree(0, 1, Podar));
    b->AddClassifier(c);
  }

  ProgresoEnsemble *pe = test ? new ProgresoEnsemble(b, test) : 0;
  ProgresoConsola pc(pe);

  if (GenerateOutputsAsKuncheva) b->SetGenerateOutputsAsKuncheva(true);

  b->Build(nd, &pc);
  if (AddC45Trees) {
    //Reset Default value without checking what was there before
    C45Tree::SetGlobalOpt('m', (char*)"2");
    C45Tree::SetGlobalOpt('d', (char*)"1");
  }

  Clasif = b;

  return this;
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/*void BuildBoostingWithBaggingInfo::DarDeAlta()
{
  ifns().AnadirInterfazFuncion( new BuildBoostingWithBaggingInfo(),
                                                              "Classifiers");
}
//---------------------------------------------------------------------------
void BuildBoostingWithBaggingInfo::CreateParams()
{
}
//---------------------------------------------------------------------------
Mandato* BuildBoostingWithBaggingInfo::Ejecutar()
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
void BuildCART::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Dataset variable",
                               "Nombre de la base de datos con los ejemplos"));
  Parametro *p = new ParametroBooleano(true, "Prune?", "By default: true");
  Parametros.push_back(p);
  p = new ParametroBooleano(true, "Multisplits?",
                            "By default: false (unidimensional splits)");
  Parametros.push_back(p);
  p = new ParametroBooleano(true, "SE 0 rule?", "By default: false (SE 1.0)");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
void BuildCART::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new BuildCART(), "Classifiers");
}
//---------------------------------------------------------------------------
Mandato* BuildCART::Ejecutar()
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
void BuildC45::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new BuildC45(), "Classifiers");
}
//---------------------------------------------------------------------------
void BuildC45::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Dataset variable",
                               "Nombre de la base de datos con los ejemplos"));
  Parametro *p = new ParametroBooleano(true, "Prune?", "By default: true");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* BuildC45::Ejecutar()
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
void BuildNNet::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new BuildNNet(), "Classifiers");
}

//---------------------------------------------------------------------------
void BuildNNet::CreateParams()
{
  Parametros.push_back(new Parametro(false, TipoData::TData(), "Data set variable",
                               "Nombre de la base de datos con los ejemplos"));

  Parametro *p = new ParametroEntero(1, 1024);
  p->PonPropiedades(false, 0, "Number of nodes in the hidde layer","Numero de nodos en la capa oculta para la red");
  Parametros.push_back(p);

  p = new ParametroEntero(1, 10000);
  p->PonPropiedades(false, 0, "Number of epochs for training","Numero de epocas en el aprendizaje de la red");
  Parametros.push_back(p);

}

//---------------------------------------------------------------------------

Mandato* BuildNNet::Ejecutar()
{


  Param(0)->Ejecutar();

  Data* nd = (Data*)Param(0)->ComoDatos();

  //se coge  el numero de neuronas
  int nNeuronas=0;
  
  if (NumParamsAsignados()>1) {
    Param(1)->Ejecutar();
    nNeuronas = Param(1)->ComoEntero();
  }

  //se coge  el numero de epocas
  int nEpocas=0;
  
  if (NumParamsAsignados()>2) {
    Param(2)->Ejecutar();
    nEpocas = Param(2)->ComoEntero();
  }
 
  if (nnet) delete nnet;

  int ndatos = nd->GetNTotal();
  nd->SetNTrain(ndatos);

  if (nEpocas==0){
    nnet = new NNet(nNeuronas);    
  }else
    nnet = new NNet(nNeuronas,nEpocas);    
 
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
void ClassificationMap2D::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ClassificationMap2D(), "Tools");
}
//---------------------------------------------------------------------------
void ClassificationMap2D::CreateParams()
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
Mandato* ClassificationMap2D::Ejecutar()
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
  DoMapaDeAlturas2D(fm, x1, x2, min1, min2, max1, max2, Ancho, Alto, dato,
                                                 (nom + "margen.bmp").c_str());
  DoMapaDeAlturas2D(fmc, x1, x2, min1, min2, max1, max2, Ancho, Alto, dato,
                                            (nom + "margenclasf.bmp").c_str());

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
      C45Tree::SetGlobalOpt('m', (char*)"1");
      C45Tree::SetGlobalOpt('d', (char*)"0");
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
      C45Tree::SetGlobalOpt('m', (char*)"2");
      C45Tree::SetGlobalOpt('d', (char*)"1");
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
    C45Tree::SetGlobalOpt('m', (char*)"1");
    C45Tree::SetGlobalOpt('d', (char*)"0");
    for(int i=0;i<nArb;i++) 
      b->AddClassifier(new C45Tree(false));

    ProgresoConsola pc;
    b->Build(nd, &pc);
    C45Tree::SetGlobalOpt('m', (char*)"2");
    C45Tree::SetGlobalOpt('d', (char*)"1");

    int err = 0;
    b->SetData(dtr);
    int **Mcd = b->MatClasif0(dtr);
    savematriz(Mcd, nArb, dtr->GetNTotal(), (char*)"mtr.csv");
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
    savematriz(Mcd, nArb, dts->GetNTotal(), (char*)"mts.csv");
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
      C45Tree::SetGlobalOpt('m', (char*)"1");
      C45Tree::SetGlobalOpt('d', (char*)"0");
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
      C45Tree::SetGlobalOpt('m', (char*)"2");
      C45Tree::SetGlobalOpt('d', (char*)"1");
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
  int nClasf = ens->GetClassifiersToUse();

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
    
    //Ensemble *ens = (Ensemble*)Param(1)->ComoDatos();
    //Data* dt = (Data*)Param(2)->ComoDatos();

    //for(int i=0;i<dt->GetNumVar()-1;i++){
    //}
  }
  else if (what.compare("AllTreeTest")==0) { //---  AllTreeTest  ---

    Param(1)->Ejecutar();
    Ensemble *e = (Ensemble*)Param(1)->ComoDatos();
    vector<set<double> > splits(e->GetData()->GetNumVar());

    int NumDivisionesTotal = 0;
    for(int i=0;i<e->GetClassifiersToUse();i++) {
      Tree *t = ((DecisionTree*)e->GetClassifier(i))->GetTree();
      int K = ((DecisionTree*)e->GetClassifier(i))->GetK();
      Node *cur = t->GetRoot();
      int numLeafs = 0;
      while(cur) {
        if (!cur->IsLeaf(K)) {
          NumDivisionesTotal++;
          splits[cur->att].insert(cur->fSplit);
    //      cout << cur->fSplit << "(" << cur->att << ")" << endl;
          cur = cur->child;
        }
        else {
          numLeafs++;
          while(cur->parent && !cur->sib) {
            cur = cur->parent;
          }
          cur = cur->sib;
        }
      }
      cout << numLeafs << "*";
    }
    cout << endl;
    int NumDivisionesDistintas = 0;
    for (unsigned i=0;i<splits.size();i++) {
      NumDivisionesDistintas += splits[i].size();
     // cout << "Atts: " << i << " - Elements: " << splits[i].size() << endl;
     // copy(splits[i].begin(), splits[i].end(), 
     //                                     ostream_iterator<double>(cout, " "));
     // cout << endl << "----------------" << endl;
    }
    cout << "DIVISIONES: Total=" << NumDivisionesTotal << " Distintas=" << 
                                               NumDivisionesDistintas << endl;
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
    SpaceCubes sc(t, mins, maxs);
//M    printf("Acumulad + tree      = resultad -> after merging\n");
//printf("%s", sc.Print().c_str());
//SpaceCubes sc2(sc.ToTree(), mins, maxs);
//sc.ToTree();
//printf("%s", sc2.Print().c_str());
    for(int i=1;i<e->GetClassifiersToUse();i++) {
      int ac, tr, re;
      ac = sc.size();
      t = (CARTTree*)e->GetClassifier(i);
      SpaceCubes spt(t, mins, maxs);
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
    cur = cart->GetTree()->GetRoot();
    while(cur->child) cur = cur->child;
    while (cur) {
      PruneInfo *pi = trees[cur];
      char buf[1000];
      cart->GetData()->WriteNode(cur, buf, -1);
      cout << "-----------------------------------" << "\n";
      cout << buf;
      for(unsigned i=0;i<pi->min_cost_alpha.size();i++) 
        cout << "\t"/* << pi->min_cost_left[i] << ", " << pi->min_cost_right[i] << ")" */<< pi->min_cost_alpha[i];
      cout << "\n";
      cur = cur->nextDown(-1);
    }
  }
  else if (what.compare("ValidacionMCT")==0) { // Validacion cruzada para MCT
    Param(1)->Ejecutar();
    Param(2)->Ejecutar();
    Param(3)->Ejecutar();
	  
    CARTTree *cart = (CARTTree*)Param(1)->ComoDatos();
    Data *data = (Data*)Param(2)->ComoDatos(); // Seguro que se puede sacar de Param(1)...
    bool SE_0 = Param(3)->ComoBooleano(); 
//    map<Node*, PruneInfo*> trees = cart->GetTree()->MinimumCostTrees();	  
	  
//	cart->GetTree()->SelectMinimumCostTrees(trees);
	  	
	  cout << "Llego" << endl;
	 
    int Nmin = 1; 
	cart->GetTree()->PruneMCT(data, Nmin, SE_0);
	  
	// Faltan las funciones de test.
  }
  else if (what.compare("FamilyError")==0) { // ---  Multiclassing  ---
    Param(1)->Ejecutar();
    Param(2)->Ejecutar();
    Param(3)->Ejecutar();
    Param(4)->Ejecutar();

    DecisionTree *tree = (DecisionTree*)Param(1)->ComoDatos();
    Data* dtr = (Data*)Param(2)->ComoDatos();
    Data* dts = (Data*)Param(3)->ComoDatos();
    string nf = Param(4)->ComoCadena();

    FILE *f = fopen(nf.c_str(), "a");
    int KMax = tree->GetTree()->GetRoot()->child->K+1;
    tree->SetData(dtr);
    for(int i=0;i<KMax;i++) {
      DecisionTree::K = i;
      fprintf(f, "%0.3f", 100.0*tree->Error(0, dtr->GetNTotal()-1) );
      fprintf(f, "%s\t", (i==tree->GetTree()->GetK() ? "*" : ""));
    }
    fprintf(f, "\n");
    tree->SetData(dts);
    for(int i=0;i<KMax;i++) {
      DecisionTree::K = i;
      fprintf(f, "%0.3f", 100.0*tree->Error(0, dts->GetNTotal()-1) );
      fprintf(f, "%s\t", (i==tree->GetTree()->GetK() ? "*" : ""));
    }
    fprintf(f, "\n");
    fclose(f);
  }
  else if (what.compare("Twonormcilla")==0) { //---  Twonormcilla  --- 
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
void Clear::CreateParams()
{
  Parametro *p = new Parametro(false, TipoData::TData(), "Data variable to free", 
               "Frees the memory associated with the variable");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
Mandato* Clear::Ejecutar()
{
  Param(0)->Ejecutar();

  Data *d = (Data*)Param(0)->ComoDatos();
  delete d;

  return this;
}
//---------------------------------------------------------------------------
void Clear::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Clear());
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void GValue::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new GValue(), "Classifiers");
}
//---------------------------------------------------------------------------
void GValue::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
Mandato* GValue::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *data = (Data *)Param(1)->ComoDatos();

  gValue = ens->computeGValue(data);

  return this;
}

