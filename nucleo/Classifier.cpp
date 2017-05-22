//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#include <typeinfo>

#include "Classifier.h"
#include "Evaluations.h"

#include <stdexcept>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include "Utils.h"
#include "values.h"
#include "Matriz.h"

using namespace std;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
ClsfRegister *ClsfRegister:: _Register=0;
ClsfRegister *ClsfRegister:: Register()
{
  if (!_Register) _Register = new ClsfRegister();
  return _Register;
}
//---------------------------------------------------------------------------
ClsfReg *ClsfRegister::GetRegByType(const type_info &type)
{
  for(unsigned i=0;i<ClsfRegs.size();i++) {
    if (ClsfRegs[i]->ClassifierTypeIs(type)) {
      return ClsfRegs[i];
    }
  }
  return 0;
}
//---------------------------------------------------------------------------
Classifier *ClsfRegister:: CreateClassifier(string clsfname)
{
  for(unsigned i=0;i<ClsfRegs.size();i++) {
    if (ClsfRegs[i]->ClassifierName()==clsfname) {
      return ClsfRegs[i]->CreateClassifier();
    }
  }
  return 0;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
char Classifier::nf_algo[260] = {"temp0983kj3he.dat"};
FILE* Classifier::f = 0;
bool Classifier::Leyendo = false;
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Classifier::Classifier()
{
  NumVar = 0;
  NumVarOrd = 0;
  NumVarFuz = 0;
  NumVarNom = 0;
}
Classifier::~Classifier()
{
  if (f) fclose(f);
  f = 0;
}

//---------------------------------------------------------------------------

void Classifier::Init(Data * data)
{
  this->data = data;
  int NTrain = data->GetNTrain();

  //Me guardo las caracteristicas del conjunto de datos
  NumVar = data->GetNumVar();
  NumVarOrd = data->GetNumVarOrd();
  NumVarFuz = data->GetNumVarFuz();
  NumVarNom = data->GetNumVarNom();
  DatosUtilizados.clear();
  for(int i=0;i<NTrain;i++) {
    int k = data->GetDatIniPos(i);
    DatosUtilizados.push_back(k);
  }

  //Calculo las probabilidades a priori de las distintas clases
  if (((NomData*)data)->NumClass > 0) {
    APrioriClas.resize(((NomData*)data)->NumClass);
    for(int i=0;i<data->GetNTrain();i++) {
      int cl =  data->GetDatClass(i);
      APrioriClas[cl] += data->GetDatWeight(i);
    }
  }
}
//---------------------------------------------------------------------------

void Classifier::Classify(double **dat, int NumeroDatos)
{
  Data *dt, *holdData;
  int nIndepVar;

  nIndepVar = data->GetNumVar()-1;

  //Generamos un conjunto de datos nuevo
  dt = data->Clone(1);
  dt->redim(NumeroDatos);
  for(int i=0;i<NumeroDatos;i++) {
    for(int j=0;j<nIndepVar;j++)
      dt->SetValueVar(i, j, dat[i][j]);
  }

  //Classify examples
  holdData = data;
  SetData(dt);
  for(int i=0;i<NumeroDatos;i++) {
    dat[i][nIndepVar] = Classify(i) + 0.5;
  }

  //Reset old data
  SetData(holdData);
  delete dt;

/*  int nIndepVar = data->GetNumVar()-1;
  double *m = new double[nIndepVar];
  //Guarda el ejemplo 0
  for(int j=0;j<nIndepVar;j++)
    m[j] = data->GetValueVar(0, j);
  //Classifica los ejemplos
  for(int i=0;i<NumeroDatos;i++) {
    for(int j=0;j<nIndepVar;j++)
      data->SetValueVar(0, j, dat[i][j]);
    dat[i][nIndepVar] = Classify(0) + 0.5;
  }
  //Recoloca el ejemplo 0
  for(int j=0;j<nIndepVar;j++)
    data->SetValueVar(0, j, m[j]);
  delete m;*/
}
//---------------------------------------------------------------------------
int Classifier::Classify(int ElementIndex)
{
  vector<double> distrb = Distribution(ElementIndex);
  return Classify(distrb);
}
//---------------------------------------------------------------------------
int Classifier::Classify(std::vector<double> distrb)
{
  return WhichClass(distrb);
}
//---------------------------------------------------------------------------
vector<double> Classifier::Distribution(int ElementIndex)
{
  vector<double> distrb;
  double TotWeight;

  distrb = UnnormalizedDistribution(ElementIndex);

  //Normalizing
  TotWeight = 0.0;
  for(int i = 0; i < ((NomData*)data)->NumClass; i++) TotWeight += distrb[i];
  for(int i = 0; i < ((NomData*)data)->NumClass; i++) distrb[i] /= TotWeight;

  return distrb;
}
//---------------------------------------------------------------------------
vector<double> Classifier::UnnormalizedDistribution(int ElementIndex)
{
  vector<double> distrb(((NomData*)data)->NumClass, 0.0);

  int Class = Classify(ElementIndex);
  distrb[Class] = 1.0;

  return distrb;
}
Evaluation* Classifier::Evaluate(Data *data, string EvaluationClassName, string params)
{
  if (EvaluationClassName.length() == 0) {
    if (data->GetDepVarType() == NOM) {
      return new ErrorEvaluation(this, data);
    }
    else if ( data->GetDepVarType() == ORD) {
      return new MSEEvaluation(this, data);
    }
  }
  else {
    return Evaluation::EvaluationByClassName(EvaluationClassName, this, data, params);
  }

  return 0;
}
//---------------------------------------------------------------------------
  /**
   *  Selecciona la clase mAs votada (del parAmetro distrb). A igualdad de
   *    votos se queda: o con la primera, o con la de mayor distribucion
   *    a priori (=mayor nUmero de ejemplos de entrenamiento) si se le pasa 
   *    el parAmetro APrioriClas.
   *    AdemAs si se le pasa el parAmetro winners introduce todas las clases
   *    con igual distrib, en la posiciOn 0 de winners se mete el no de clases
   *    ganadoras
   */

int Classifier::WhichClass(vector<double> distrb, int *winners,
                                              std::vector<double> *APrioriClas)
{
  double max = distrb[0];
  int maxj = 0;
  bool mine = false;

  if (!winners) {
    winners = new int[distrb.size()+1];
    mine = true;
  }

  if (winners) {
    winners[0] = 1;
    winners[1] = 0;
  }

  for(unsigned j=1; j<distrb.size(); j++)  {
    if (distrb[j]>max) {
      maxj = j;
      max = distrb[j];
      if (winners) {
        winners[0] = 1;
        winners[1] = j;
      }
    }
    else if (distrb[j]==max) {
      if (APrioriClas) {
//        if ((*APrioriClas)[j]<(*APrioriClas)[maxj]) { //Sele
        if ((*APrioriClas)[j]>(*APrioriClas)[maxj]) {
          maxj = j;
          max = distrb[j];
          if (winners) {
            winners[0] = 1;
            winners[1] = j;
          }
        }
        else if((*APrioriClas)[j]==(*APrioriClas)[maxj] && winners) {
          winners[0]++;
          winners[winners[0]] = j;
        }
      }
      else if(winners) {
        winners[0]++;
        winners[winners[0]] = j;
      }
    }
  }

  // Random class
  bool UseRandomClass = false;
  if (winners && UseRandomClass) {
    maxj = winners[ 1 + (int) (1.0*winners[0]*rand()/(RAND_MAX+1.0)) ];
  }
  if (winners && mine) delete [] winners;

  return maxj;
}
int Classifier::WhichClass(double* dist, int NumClass, int *winners, 
                                              std::vector<double> *APrioriClas)
{
  vector<double> distrb;

  for(int i=0;i<NumClass;i++) distrb.push_back(dist[i]);

  return WhichClass(distrb, winners, APrioriClas);
}
//---------------------------------------------------------------------------
double Classifier::Margin(Instance& dis, int RealClass, bool Normalized)
{
  vector<double> ins = (vector<double>)dis;
  ins.resize(dis.GetNVars());
  return Margin(ins, RealClass, Normalized);
}
double Classifier::Margin(vector<double> dis, int RealClass, bool Normalized)
{ 
  double maxmalo = 0.0;
  double tot = 0.0;

  if ( RealClass < 0 ) {
    RealClass = Classifier::WhichClass(dis);
  }

  for(unsigned id = 0; id < dis.size(); id++) {
    tot += dis[id];
    if ( (int)id == RealClass ) continue;
    if ( dis[id] > maxmalo ) {
      maxmalo = dis[id];
    }
  }
  return Normalized ? (dis[RealClass] - maxmalo)/tot : (dis[RealClass] - maxmalo);
}
double Classifier::Margin(double* dis, int NumClass, int RealClass, bool Normalized)
{
  vector<double> distrb;

  for(int i = 0; i < NumClass; i++) 
    distrb.push_back(dis[i]);

  return Margin(distrb, RealClass, Normalized);
}
//---------------------------------------------------------------------------

double Classifier::ClassificationCertainty(int ElementIndex, int &Class)
{
  vector <double> distrb = Distribution(ElementIndex);

  double TotWeight=0.0;
  for(int i=0;i<((NomData*)data)->NumClass;i++) TotWeight += distrb[i];
  Class = WhichClass(distrb);

  return distrb[Class]/TotWeight;
}
//---------------------------------------------------------------------------

double Classifier::Error(Data *d)
{
  Data *hold = GetData();
  SetData(d);
  double err = Error(0, d->GetNTotal()-1);
  SetData(hold);
  return err;
}
//---------------------------------------------------------------------------
double Classifier::Error(int first, int last)
{
  double err=0.0;
  double TotalWeight=0.0;
  rc.clear();
  pc.clear();
  for (int i=first;i<=last;i++) {
    int RealClass = data->GetDatClass(i);
    int PredictedClass = Classify(i);

    TotalWeight += data->GetDatWeight(i);

    if (PredictedClass!=RealClass){
      err += data->GetDatWeight(i);
    }
    rc.push_back(RealClass);
    pc.push_back(PredictedClass);
  }
  return err/TotalWeight;
}
//---------------------------------------------------------------------------

void Classifier::ReadOrder(Data *data)
{
  if (!f)
    f = fopen(nf_algo,"rb");
  int NTotar = data->GetTotalData();
  int iOrden = data->GetNumVar() + Data::GetIniPosIndex();
  int iAux   = data->GetNumVar();
  int iGrp   = data->GetNumVar() + Data::GetGroupIndex();
  int h;
  data->SortOn(iOrden, 0, NTotar-1);
  fread( &h, sizeof(int), 1, f);
  data->GroupsCreated = h;
  if (data->GroupsCreated) {
    fread( &data->GroupCount[0], sizeof(int), 1, f);
    fread( &data->GroupCount[1], sizeof(int), 1, f);
  }
  for(int i=0;i<NTotar;i++) {
    fread( &h, sizeof(int), 1, f);
    data->SetValueVar(h, iAux, 0.5 + i);
  }
  data->SortOn(iAux, 0, NTotar-1);
  if (data->GroupsCreated) {
    int i;
    for(i=0;i<data->GroupCount[0];i++)
      data->SetValueVar(i, iGrp, 0.5);
    for(;i<data->GroupCount[0]+data->GroupCount[1];i++)
      data->SetValueVar(i, iGrp, 1.5);
  }
}
//---------------------------------------------------------------------------

void Classifier::WriteOrder(Data *data)
{
  if (!f)
    f = fopen(nf_algo,"wb");
  int NTotar = data->GetTotalData();
  int iOrden = data->GetNumVar() + Data::GetIniPosIndex();
  int h = data->GroupsCreated;
  fwrite( &h, sizeof(int), 1, f);
  if (data->GroupsCreated) {
    fwrite( &data->GroupCount[0], sizeof(int), 1, f);
    fwrite( &data->GroupCount[1], sizeof(int), 1, f);
  }
  for(int i=0;i<NTotar;i++) {
    h = (int)data->GetValueVar(i, iOrden);
    fwrite( &h, sizeof(int), 1, f);
  }
}
//---------------------------------------------------------------------------
void Classifier::SetData(Data *dat)
{

  if (NumVar != dat->GetNumVar() || NumVarOrd != dat->GetNumVarOrd() ||
       NumVarFuz != dat->GetNumVarFuz() || NumVarNom != dat->GetNumVarNom()) {
    //int k = 1/0;//throw new invalid_argument("Tronco, esto que me das es distinto");
    throw /*exception(*/"Tronco, esto que me das es distinto"/*)*/;
  }
  data = dat;
}
//---------------------------------------------------------------------------
void Classifier::ResetOriginalTrainingData()
{
  int Total =  data->GetTotalData();
  if (!data || (int)DatosUtilizados.size()>Total) return;
  int iPos = NumVar + Data::IniPosIndex;
  data->SortOn(iPos, 0, Total-1);
  for(int i=0;i<Total;i++) {
    data->SetValueVar(i, NumVar, Total);
  }
  for(unsigned i=0;i<DatosUtilizados.size();i++) {
    data->SetValueVar(DatosUtilizados[i], NumVar, i);
  }
  data->SortOn(NumVar, 0, Total-1);
  data->SetNTrain(DatosUtilizados.size());
  Init(data);//??????????????No sE si es mejor ponerlo o no
}
//---------------------------------------------------------------------------
bool Classifier::UsedInOriginalTrainingData(int value)
{
  for (unsigned i=0;i<DatosUtilizados.size();i++) {
    if (value==DatosUtilizados[i]) return true;
  }
  return false;
}
//---------------------------------------------------------------------------
string Classifier::Info(int value)
{
  string res = string("{Clase: ") + typeid(*this).name() + "}";
  char kk[256];
  sprintf(kk, " -Dataset: # clases: %d NTrain: %d ",
                         ((NomData*)data)->NumClass, DatosUtilizados.size());

  if (value>0) {
    char kk[24];
    sprintf(kk, " -NumVar: %d", NumVar);
    res = res + kk;
    if (data){
      sprintf(kk, " -# clases: %d", ((NomData*)data)->NumClass);
      res = res + kk;
    }
  }
  return res;
}
//---------------------------------------------------------------------------
void Classifier::Guardar(string salida, int version)
{
  ofstream ofs(salida.c_str(), ios_base::out | ios_base::binary);
  Guardar(ofs, version);
  ofs.close();
}
//---------------------------------------------------------------------------
void Classifier::Leer(string in)
{
  ifstream ifs(in.c_str(), ios_base::in | ios_base::binary);
  LeerConCabecera(ifs);
  ifs.close();
}
//---------------------------------------------------------------------------
void Classifier::Guardar(std::ostream &sal, int version)
{
  GuardaCabecera(sal, version);

  //Datos del conjunto de datos
  sal.write((char*)&NumVar, sizeof(NumVar));
  sal.write((char*)&NumVarOrd, sizeof(NumVarOrd));
  sal.write((char*)&NumVarFuz, sizeof(NumVarFuz));
  sal.write((char*)&NumVarNom, sizeof(NumVarNom));

  if (version==0) {
    //Indices de los datos utilizados
    unsigned num = DatosUtilizados.size();
    sal.write((char*)&num, sizeof(unsigned));
    for(unsigned i=0;i<num;i++) {
      int dat = DatosUtilizados[i];
      sal.write((char*)&dat, sizeof(int));
    }
  }
  else if (version==1) {
    //Indices de los datos utilizados (version bit)
    unsigned i, num = DatosUtilizados.size();
    unsigned char *usados = new unsigned char[(7+num)/8];
    memset(usados, 0, (7+num)/8);
    for(i=0;i<num;i++) {
//      int j = ix[Classifiers[k]->OriginalTrainingDataPos(i)];
      int j = DatosUtilizados[i];
      if (j>=(int)DatosUtilizados.size()) break;
      usados[j/8] = usados[j/8] | 1<<(7-j%8);
    }
    if (i==num) {
      sal.write((char*)&num, sizeof(unsigned));
      sal.write((char*)usados, (7+num)/8);
    }
    else {//Como en version 0
      num = -num;
      sal.write((char*)&num, sizeof(unsigned));
      num = -num;
      for(unsigned i=0;i<num;i++) {
        int dat = DatosUtilizados[i];
        sal.write((char*)&dat, sizeof(int));
      }
    }
  }
}
//---------------------------------------------------------------------------
int Classifier::LeerConCabecera(std::istream &in)
{
  int version;

  LeeCabecera(in, version);
  Leer(in, version);

  return version;
}
//---------------------------------------------------------------------------
void Classifier::Leer(std::istream &in, int version)
{
  int dat;
  int tam;

  //Datos del conjunto de datos
  in.read((char*)&NumVar, sizeof(NumVar));
  in.read((char*)&NumVarOrd, sizeof(NumVarOrd));
  in.read((char*)&NumVarFuz, sizeof(NumVarFuz));
  in.read((char*)&NumVarNom, sizeof(NumVarNom));

  if (version==0) {
    in.read((char*)&tam, sizeof(unsigned));
    for(int i = 0; i < tam; i++) {
      in.read((char*)&dat, sizeof(int));
      DatosUtilizados.push_back(dat);
    }
  }
  else if (version==1){
    in.read((char*)&tam, sizeof(unsigned));
    if (tam<0) {
      tam = -tam;
      for(int i = 0; i < tam; i++) {
        in.read((char*)&dat, sizeof(int));
        DatosUtilizados.push_back(dat);
      }
    }
    else {
      //Indices de los datos utilizados (version bit)
      unsigned char *usados = new unsigned char[(7+tam)/8];
      in.read((char*)usados, (7+tam)/8);
      DatosUtilizados.resize(tam, -1);
      int k=0;
      //FILE *f=fopen("oob.csv", "a");
      for(int i=0;i<tam;i++) { 
        if (usados[i/8] & 1<<(7-i%8)) {
          DatosUtilizados[k++] = i;
        //  fprintf(f, "0, ");
        }
        else { 
          //fprintf(f, "1, ");
        }
      }
//      fprintf(f, "\n");fclose(f);
      for(int i=k;i<tam;i++) { 
        DatosUtilizados[i] = DatosUtilizados[0];
      }
      delete []usados;
    }
  }
}
//---------------------------------------------------------------------------
Classifier *Classifier::LeerClasificadorGeneral(std::istream &in)
{
  int version;
  string nombre = LeeCabecera(in, version);
  Classifier *clas;
  /*
  if      (strcmp(nombre.c_str(),typeid(Bagging).name())==0) clas = new Bagging();
  else if (strcmp(nombre.c_str(),typeid(CART).name())==0) clas = new CARTTree();
  else if (strcmp(nombre.c_str(),typeid(Boosting).name())==0) clas = new Boosting();
  else if (strcmp(nombre.c_str(), typeid(BoostingWithBaggingInfo).name())==0)
                                           clas = new BoostingWithBaggingInfo();
  else     clas = 0;
*/
  clas = ClsfRegister::Register()->CreateClassifier(nombre);
  if (clas) clas->Leer(in, version);

  return clas;
}
//---------------------------------------------------------------------------
string Classifier::LeeCabecera(std::istream &in, int &version)
{
  char nombre[65];
  char nom[65];

  in.get(nombre, 65);

  int nl = sscanf(nombre, "%s %d", nom, &version);
  if (nl==1) version = 0;

/*  int i, j;
  for(i=0, j=0;i<64;i++) {
    in.get(a);
    if (a!=' ') {
      nombre[j]=a;
      j++;
    }
  }
  nombre[j] = '\0';
  nom = nombre;*/
  return nom;
}
//---------------------------------------------------------------------------
string Classifier::nombre()
{
  string nom;/* = typeid(*this).name();
#ifdef _WIN32
#else
  int i=0;
  while(isdigit(nom[i])) i++;
  nom.erase(0, i);
#endif*/
  ClsfReg *reg = ClsfRegister::Register()->GetRegByType(typeid(*this));
  if (!reg) throw domain_error("Class not registered");
  return reg->ClassifierName();
}
//---------------------------------------------------------------------------
void Classifier::GuardaCabecera(std::ostream &salida, int version)
{
  if (version==0) 
    salida << setw(64) << nombre();
  else {
    ostringstream buf;
    buf << nombre() << " " << version;
    salida << setw(64) << buf.str();
  }
}
//---------------------------------------------------------------------------
int **MatCorrect(int **MatClasf, int nClsf, Data *data) 
{
  return 0;
}
//---------------------------------------------------------------------------
int **MatClasif(vector<Classifier*> clsfs, Data *data)
{
  int nc = clsfs.size();
  int nd = data->GetNTotal();

  int **Mcd = new int*[2*nc];

  for(int i=0;i<nc;i++) {
    Data *holdData = clsfs[i]->GetData();
    clsfs[i]->SetData(data);
    Mcd[i] = new int[nd];
    Mcd[i+nc] = new int[nd];
    for(int j=0;j<nd;j++) {
      Mcd[i][j] = clsfs[i]->Classify(j);
      Mcd[i+nc][j] = Mcd[i][j]==data->GetDatClass(j) ? 1 : 0;
    }
    clsfs[i]->SetData(holdData);
  }

  return Mcd;
}
//---------------------------------------------------------------------------
double QStatistic(Classifier *c1, Classifier *c2, Data *data) 
{
  int nDatos = data->GetNTotal();
  int *er1 = new int[nDatos];
  int *er2 = new int[nDatos];

  c1->SetData(data);
  c2->SetData(data);
  
  for(int i=0;i<nDatos;i++) {
    er1[i] = c1->Classify(i)==data->GetDatClass(i) ? 1 : 0;
    er2[i] = c2->Classify(i)==data->GetDatClass(i) ? 1 : 0;
  }

  double Q = QStatistic(er1, er2, nDatos);

  delete []er1;
  delete []er2;

  return Q;
}
//---------------------------------------------------------------------------
double QStatistic(int *c1, int *c2, int *clase_ok, int nDatos)
{
  int *er1 = new int[nDatos];
  int *er2 = new int[nDatos];

  for(int i=0;i<nDatos;i++) {
    er1[i] = c1[i]==clase_ok[i] ? 1 : 0;
    er2[i] = c2[i]==clase_ok[i] ? 1 : 0;
  }

  double Q = QStatistic(er1, er2, nDatos);

  delete []er1;
  delete []er2;

  return Q;
}
//---------------------------------------------------------------------------
double QStatistic(int *er1, int *er2, int nDatos) 
{
  int N00, N10, N01, N11;

  N00 = N10 = N01 = N11 = 0;

  for(int i=0;i<nDatos;i++) {
    if (er1[i]==1 && er2[i]==1) N11++;
    else if (er1[i]==1 && er2[i]==0) N10++;
    else if (er1[i]==0 && er2[i]==1) N01++;
    else N00++;
  }
  
  return ((double)N11*N00-N01*N10)/(N11*N00+N01*N10);
}
//---------------------------------------------------------------------------
double Classifier::CommonErrors(int *c1, int *c2, int *trueClass, int size)
{
	
  int commonErrors = 0; 

  // We check wherther or not the two classifiers misclassify a given instance
  
  for (int i=0;i< size;i++) {
    if (c1[ i ] == c2[ i ] && c1[ i ] != trueClass[ i ])
	commonErrors++;
  }
  
  return (double) commonErrors;
}

//---------------------------------------------------------------------------
double Classifier::KappaStatistic(Classifier *c1, Classifier *c2, int ini, 
                                                                        int fin)
{
  int *cc1, *cc2;
  int m;
  NomData *data = (NomData*)c1->GetData();

  ini = ini<0 ? 0 : ini;
  fin = fin<0 ? c1->GetData()->GetNTrain()-1 : fin;
  m = fin-ini+1;
  
  cc1 = new int[m];
  cc2 = new int[m];
  
  c2->SetData(c1->GetData());
  
  for (int i=0;i<m;i++) {
    cc1[i] = c1->Classify(i);
    cc2[i] = c2->Classify(i);
  }
  
  return KappaStatistic(cc1, cc2, m, data->NumClass);
}
//---------------------------------------------------------------------------
double Classifier::KappaStatistic(int *c1, int *c2, int m, int L)
{
  double theta1, theta2;
  int **C;

  C = new int*[L];
  for (int i=0;i<L;i++) {
    C[i] = new int[L];
    for (int j=0;j<L;j++) C[i][j] = 0;
  }

  for (int i=0;i<m;i++)
    C[ c1[i] ][ c2[i] ] ++;
  
  theta1 = 0.0;
  for (int i=0;i<L;i++) 
    theta1 += C[i][i];

  if ( theta1 == m ) return 1.0; //Coinciden en todo
  theta1 /= m;

  theta2 = 0.0;
  for (int i=0;i<L;i++) {
    double f1=0.0, f2=0.0;
    for (int j=0;j<L;j++) { 
      f1 += C[i][j];
      f2 += C[j][i];
    }
    theta2 += f1*f2;
  }
  theta2 /= (m*m);
  return (theta1-theta2)/(1.0 - theta2);
}
//---------------------------------------------------------------------------
Data* ReduceData(Data *data, Classifier *clsf, double Umbral)
{
  int NTrain = data->GetNTrain();
  int iVar = data->GetNumVar()+1;

/*  if (!clsf) {
    Ensemble *ens = new Bagging();
    for(unsigned i=0;i<Classifiers.size();i++) {
      ens->AddClassifier(new C45Tree());
    }
    ens->Build(data);
    clsf = ens;
  }*/

  Data *data0 = clsf->GetData();
  clsf->SetData(data);

  int clase;
  int datosguenos = 0;
  for(int i=0;i<NTrain;i++) {
//    int inipos = data->GetDatIniPos(i);
    double certeza = clsf->ClassificationCertainty(i, clase);
    if (clase!=data->GetDatClass(i) && certeza>Umbral) {
         data->SetValueVar(i, iVar, 1.0);
  //       excluidos.push_back(inipos);
    }
    else {
      data->SetValueVar(i, iVar, 0.0);
      datosguenos++;
    }
  }
  data->SortOn(iVar, 0, NTrain-1);
  clsf->SetData(data0);
 
  Data *dt = data->Clone(0, NTrain-1);
  dt->SetNTrain(datosguenos);

  return dt;
}
//---------------------------------------------------------------------------
// Genera N*N datos adicionales calculando el pto medio de cada dato con
// todos los demas. En el conjunto devuelto se incluyen tambien los
// datos originales.
//---------------------------------------------------------------------------
Data* ExpandData(Data *data, Classifier *clsf)
{
  int N = data->GetNTrain();
  int nVar = data->GetNumVar();

  Data *dt = data->Clone(0, N-1);
  for(int i=0;i<N;i++) 
    dt->SetValueVar(i, nVar+Data::IniPosIndex, i);

  int idt = N;
  int ejem = N-1;//(int)(0.1*N);
  dt->redim(N*ejem+N);
  dt->SetNTrain(dt->GetNTotal());
  for(int i=0;i<N;i++) {
    for(int j=0;j<N;j++) {
      if (i==j) continue;
      for(int k=0;k<nVar+7;k++) {
        double v_ik, v_jk, val;
        v_ik = data->GetValueVar(i, k);
        v_jk = data->GetValueVar(j, k);
        if (v_ik==DBL_MAX)      val = v_jk;
        else if (v_jk==DBL_MAX) val = v_ik;
        else         val = 0.5*(v_ik+v_jk);
        dt->SetValueVar(idt, k, val);
      }
      dt->SetValueVar(idt, nVar+Data::IniPosIndex, idt);
      idt++;
    }
  }

  Data *data0 = clsf->GetData();
  clsf->SetData(dt);
  for(int i=N;i<dt->GetNTotal();i++) 
    dt->SetValueVar(i, nVar-1, clsf->Classify(i));
  clsf->SetData(data0);
 
  return dt;
}
//---------------------------------------------------------------------------
// Igual que arriba pero en vez de usar el conjunto entero (de tamanyo N)
// se usan los datos con mayor incertidumbre, que se supone que
// son los de las fronteras  
//---------------------------------------------------------------------------
Data* ExpandData2(Data *data, Classifier *clsf)
{
  double corte = 0.6;
  int icorte, dummy;
  int N = data->GetNTrain();
  int nVar = data->GetNumVar();

  clsf->SetData(data);

  for(int i=0;i<N;i++) 
    data->SetValueVar(i, nVar, clsf->ClassificationCertainty(i, dummy));

  data->SortOn(nVar, 0, N-1);

  icorte = 0;
  while(icorte<N && data->GetValueVar(icorte, nVar)<corte)
    icorte++;

  if (icorte==0) return data->Clone();

  data->SetNTrain(icorte);

  Data *dt = ExpandData(data, clsf);

  int tot = dt->GetNTotal();
  int idt = 0;
  dt->redim(tot+N-icorte);
  for(int i=icorte;i<N;i++) {
    for(int k=0;k<nVar+7;k++) {
      dt->SetValueVar(tot+idt, k,data->GetValueVar(i, k));
    }
    idt++;
  }

  data->SetNTrain(N);

  //Chapucilla para evitar probremas en C45
  dt->SetNTrain(dt->GetNTotal()); 
  Data *dtt = dt->Clone(0, dt->GetNTotal()-1); 
  delete dt;

  return dtt;
}
//---------------------------------------------------------------------------
//------------------------------------------------------  MIClassifier  -----
//---------------------------------------------------------------------------

int MIClassifier::Classify(int ElementIndex)
{
  
  vector<double> distrb = Distribution(ElementIndex);
  return Classifier::Classify(distrb);
}
//---------------------------------------------------------------------------
vector<double> MIClassifier::Distribution(int ElementIndex)
{
  vector<double> distrb(((NomData*)data)->NumClass, 0.0);

  int Class=Classify(ElementIndex);
  distrb[Class] = 1.0;

  return distrb;
}
//---------------------------------------------------------------------------
double MIClassifier::Error(Data *d)
{
  Data *hold = GetData();
  SetData(d);
  double err = Error(0, d->GetNTotal()-1);
  SetData(hold);
  return err;
}
//---------------------------------------------------------------------------
double MIClassifier::Error(int first, int last)
{
  double err=0.0;
  double TotalWeight=0.0;
  for (int i=first;i<=last;i++) {
    int RealClass = data->GetDatClass(i);
    TotalWeight += data->GetDatWeight(i);
    if (Classify(i)!=RealClass)
      err += data->GetDatWeight(i);
  }
  return err/TotalWeight;
}

