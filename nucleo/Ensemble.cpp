//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#include "Ensemble.h"
#include "tree20.h"
#include "node20.h"
#include "GPTree.h"
#include <math.h>
#include <string.h>
#include <sstream>
#include <iomanip>
#include <map>

using namespace std;
//---------------------------------------------------------------------------
//----------------------------------------------------------  Ensemble  -----
//---------------------------------------------------------------------------

Ensemble::Ensemble():Classifier()
{
  ClassifiersToUse = -1;
  UseWeights = false;//km
}

Ensemble::~Ensemble()
{
  for(unsigned i=0;i<Classifiers.size();i++)
    delete Classifiers.at(i);
}
//---------------------------------------------------------------------------
void Ensemble::ResetOrdenOriginal()
{
  OrdenOriginal.clear();
  for(unsigned i=0;i<Classifiers.size();i++)
    OrdenOriginal.push_back(i);

  PesosOriginales = ClassifierWeights;
}
//---------------------------------------------------------------------------
void Ensemble::Init(Data * data)
{
  Classifier::Init(data);
}

//---------------------------------------------------------------------------

/*int Ensemble::Classify(int ElementIndex)
{
  int nclases = ((NomData*)data)->NumClass;
  int clase;
  double *votes = new double[nclases];
  double *votes2 = new double[nclases];
  for (int j=0;j<nclases;j++) votes2[j]=votes[j]=0;
  int nClassifiersToUse = GetClassifiersToUse();
  for(int i=0;i<nClassifiersToUse;i++) {
    clase = Classifiers.at(i)->Classify(ElementIndex);
    votes[clase]++;
    votes2[clase] += ClassifierWeights.at(i);
  }
  ClassWithoutWeigths = WhichClass(votes, ((NomData*)data)->NumClass);
  ClassWithWeigths    = WhichClass(votes2, ((NomData*)data)->NumClass);

  //Si hay n clases con el mismo numero de votos e igual nUmero de
  //ejemplos se quitan los votos de los Ultimos clasificadores hasta
  //que sOlo haya un ganador
  for(int i=nClassifiersToUse-1;i>0 && ClassWithoutWeigths<0; i--) {
    clase = Classifiers.at(i)->Classify(ElementIndex);
    votes[clase]--;
    ClassWithoutWeigths = WhichClass(votes, ((NomData*)data)->NumClass);
  }
  for(int i=nClassifiersToUse-1;i>0 && ClassWithWeigths<0; i--) {
    clase = Classifiers.at(i)->Classify(ElementIndex);
    votes2[clase] -= ClassifierWeights.at(i);
    ClassWithWeigths = WhichClass(votes2, ((NomData*)data)->NumClass);
  }

  delete []votes;
  delete []votes2;
  clase = (UseWeights) ? ClassWithWeigths : ClassWithoutWeigths;
  return clase;
}*/
//---------------------------------------------------------------------------

/*double Ensemble::ClassificationCertainty(int ElementIndex, int &Class*//*, 
                  bool ExcluirTrain*//*)
{
  int nclases = data->GetNumVar()-1;
  double percent;
  int clase;
  double *votes = new double[nclases];
  double *votesW = new double[nclases];
  double TotPeso=0;
  double TotPesoW=0;
  for (int j=0;j<nclases;j++) votesW[j]=votes[j]=0;
  int nClassifiersToUse = GetClassifiersToUse();
  for(int i=0;i<nClassifiersToUse;i++) {
    Classifier *c = Classifiers.at(i);
//    if (ExcluirTrain && c->UsedInOriginalTrainingData(ElementIndex)) continue;
    clase = c->Classify(ElementIndex);
    votes[clase]++;
    TotPeso++;
    votesW[clase] += ClassifierWeights.at(i);
    TotPesoW += ClassifierWeights.at(i);
  }
  ClassWithoutWeigths = WhichClass(votes, ((NomData*)data)->NumClass);
  ClassWithWeigths    = WhichClass(votesW, ((NomData*)data)->NumClass);

  //Si hay n clases con el mismo numero de votos e igual nUmero de
  //ejemplos se quitan los votos de los Ultimos clasificadores hasta
  //que sOlo haya un ganador
  for(int i=nClassifiersToUse-1;i>0 && ClassWithoutWeigths<0; i--) {
    clase = Classifiers.at(i)->Classify(ElementIndex);
    votes[clase]--;
    TotPeso--;
    ClassWithoutWeigths = WhichClass(votes, ((NomData*)data)->NumClass);
  }
  for(int i=nClassifiersToUse-1;i>0 && ClassWithWeigths<0; i--) {
    clase = Classifiers.at(i)->Classify(ElementIndex);
    votesW[clase] -= ClassifierWeights.at(i);
    TotPesoW -= ClassifierWeights.at(i);;
    ClassWithWeigths = WhichClass(votesW, ((NomData*)data)->NumClass);
  }

  percent = (UseWeights) ? votesW[ClassWithWeigths]/TotPesoW :
                                             votes[ClassWithoutWeigths]/TotPeso;
  Class = (UseWeights) ? ClassWithWeigths : ClassWithoutWeigths;
  //INFO
*//*  FILE *ff;
  ff=fopen("info.txt", "at");
  fprintf(ff, "%d\t%6.3f", Class, percent);
  if (UseWeights) {
    for(int i=0;i<((NomData*)data)->NumClass;i++) {
      fprintf(ff, "\t%4.1lf", votesW[i]);
    }
    fprintf(ff, "\t%6.3f\n", TotPesoW);
  }
  else {
    for(int i=0;i<((NomData*)data)->NumClass;i++) {
      fprintf(ff, "\t%d", (int)votes[i]);
    }
    fprintf(ff, "\t%d", (int)TotPeso);
  }
  fprintf(ff, "\n");
  fclose(ff);*/
  //
  /*delete []votes;
  delete []votesW;
  return percent;
}*/
//---------------------------------------------------------------------------
double Ensemble::Average(int ElementIndex)
{
  vector<double>& weights = GetUsingWeights();
  double suma = 0.0;
  double total = 0.0;

  for(int i=0;i<GetClassifiersToUse();i++) {
    Classifier *c = Classifiers.at(i);
    suma += weights[i]*c->Average(ElementIndex);
    total += weights[i];
  }

  return suma/total;
}
double Ensemble::Average(int ElementIndex, int AttributeIndex)
{
  vector<double>& weights = GetUsingWeights();
  double suma = 0.0;
  double total = 0.0;

  for(int i=0;i<GetClassifiersToUse();i++) {
    Classifier *c = Classifiers.at(i);
    suma += weights[i]*c->Average(ElementIndex, AttributeIndex);
    total += weights[i];
  }

  return suma/total;
}
//---------------------------------------------------------------------------
std::vector<double> Ensemble::MultipleAverage(int ElementIndex)
{
  int numberObjectives = MultiobjectiveInstance::GetNumMultiobjectives();
  vector<double>& weights = GetUsingWeights();
  std::vector<double> sumas (numberObjectives, 0.0);
  double total = 0.0;

  for(int i=0;i<GetClassifiersToUse();i++) {
    Classifier *c = Classifiers.at(i);
    std::vector<double> averages = c->MultipleAverage(ElementIndex);
    for (int nObjective = 0; nObjective < numberObjectives; nObjective++) {
      sumas[nObjective] += weights[i] * averages[nObjective];
    }
    total += weights[i];
  }

  for (int nObjective = 0; nObjective < numberObjectives; nObjective++) {
    sumas[nObjective] /= total;
  }

  return sumas;
}
//---------------------------------------------------------------------------
vector<double> Ensemble::UnnormalizedDistribution(int ElementIndex)
{
  vector<double> distrb;
  vector<double>& weights = GetUsingWeights();
 // double total = 0.0;

  for (int i=0;i<((NomData*)data)->NumClass;i++) distrb.push_back(0.0);

  for(int i=0;i<GetClassifiersToUse();i++) {
    Classifier *c = Classifiers.at(i);
    int clase = c->Classify(ElementIndex);

    distrb[clase] += weights[i];
//    total += weights[i];

/*    if (!UseWeights) {
      distrb[clase] += 1.0;
      total += 1.0;
    }
    else {
      distrb[clase] += ClassifierWeights[i];
      total += ClassifierWeights[i];
    }*/
  }

//  //Normalizamos
//  for (int i=0;i<((NomData*)data)->NumClass;i++) distrb[i] = distrb[i]/total;

  return distrb;
}
//---------------------------------------------------------------------------
vector<double> Ensemble::Margen(int ini, int fin, int div)
{
  if (div<0) div=GetClassifiersToUse();

//  double sum = 0;
  double pesosdatos;
  vector<double> salida;

  for (int i=0;i<div+1;i++) {
    salida.push_back(0.0);
  }

/*  bool UseWeights0 = GetUseWeights();
  SetUseWeights(true);
  vector<double> Pesos0 = ClassifierWeights;

  //Normalizamos los pesos a 1
  for (unsigned i=0;i<ClassifierWeights.size();i++) 
    sum += ClassifierWeights[i];
  for (unsigned i=0;i<ClassifierWeights.size();i++) 
    ClassifierWeights[i] = ClassifierWeights[i]/sum;*/

  pesosdatos = 0;
  for (int i=ini;i<=fin;i++) {
    int clase = data->GetDatClass(i);
    vector<double> dis = Distribution(i);
    double maxmalo = 0.0;
    for(unsigned id=0;id<dis.size();id++) {
      if ((int)id==clase) continue;
      if (dis[id]>maxmalo) 
        maxmalo = dis[id];
    }
    double icer = dis[clase] - maxmalo;
    int pos =(int) ((icer+1.0)*div/2.0 + 0.5);
    salida[pos] = salida[pos] + data->GetDatWeight(i);
    pesosdatos += data->GetDatWeight(i);
  }

  for (int i=0;i<div+1;i++) {
    salida[i] = salida[i]/pesosdatos;
  }

  /*ClassifierWeights = Pesos0;
  SetUseWeights(UseWeights0);*/

  return salida;
}
//---------------------------------------------------------------------------

Classifier* Ensemble::GetClassifier(int Index)
{
  return Classifiers.at(Index);
}
//---------------------------------------------------------------------------

void Ensemble::SetClassifier(int Index, Classifier *c)
{
  Classifiers[Index] = c;
}
//---------------------------------------------------------------------------

void Ensemble::Exchange(int Index1, int Index2)
{
  double T = ClassifierWeights[Index1];
  Classifier *C = Classifiers[Index1];
  ClassifierWeights[Index1] = ClassifierWeights[Index2];
  Classifiers[Index1] = Classifiers[Index2];
  ClassifierWeights[Index2] = T;
  Classifiers[Index2] = C;

  if (OrdenOriginal.size()==0) {
    ResetOrdenOriginal();
  }
  int jal = OrdenOriginal[Index1];
  OrdenOriginal[Index1] = OrdenOriginal[Index2];
  OrdenOriginal[Index2] = jal;
}
//---------------------------------------------------------------------------

void Ensemble::AddClassifier(Classifier *Clasf)
{
  Classifiers.push_back(Clasf);
  ClassifierWeights.push_back(1.0);
}
//---------------------------------------------------------------------------

void Ensemble::RemoveClassifier(int i)
{
  Classifiers.erase(Classifiers.begin() + i);
  ClassifierWeights.erase(ClassifierWeights.begin() + i);
}
//---------------------------------------------------------------------------

double Ensemble::GetClassifierWeight(int Index)
{
  return ClassifierWeights.at(Index);
}
//---------------------------------------------------------------------------

double Ensemble::GetWeightInUse(int Index)
{
  return UseWeights ? ClassifierWeights.at(Index) : 1.0;
}
//---------------------------------------------------------------------------

vector<double> UnitaryWeights(1000, 1.0);
vector<double>& Ensemble::GetUsingWeights()
{
  if (UseWeights) return ClassifierWeights;
  else {
    if (UnitaryWeights.size()<(unsigned)Count()) {
      UnitaryWeights.resize(Count(), 1.0);
    }
    return UnitaryWeights;
  }
}
//---------------------------------------------------------------------------

void Ensemble::SetClassifierWeight(int Index, double Weight)
{
  ClassifierWeights.at(Index) = Weight;
}
//---------------------------------------------------------------------------
    /*
double Ensemble::EstimError()
{
  double esterr=0.0;
  if (UseWeights) {
    double TotWeight=0.0;
    for(unsigned i=0;i<Classifiers.size();i++) {
      double Weight = ClassifierWeights.at(i);
      esterr += Classifiers.at(i)->EstimError()*Weight;
      TotWeight += Weight;
    }
    esterr /= (Classifiers.size()*TotWeight);
  }
  else {
    for(unsigned i=0;i<Classifiers.size();i++) {
      esterr += Classifiers.at(i)->EstimError();
    }
    esterr /= Classifiers.size();
  }
  if (esterr>1 || esterr<0) throw exception();
  return esterr;
}
  */

int Ensemble::Count()
{
  return Classifiers.size();
}
//---------------------------------------------------------------------------

void Ensemble::SecuencialClassify(int ElementIndex, vector<int> *classes,
                                                     std::vector<int> *indclss)
{
  int nclases = ((NomData*)data)->NumClass;
  int clase;
  int *winners;

  vector<double> votes(nclases, 0.0);
  int nClassifiersToUse = GetClassifiersToUse();
  classes->clear();
  if (indclss) indclss->clear();

  winners = new int[nclases+1];

  vector<double> &weights = GetUsingWeights();
  for(int i=0;i<nClassifiersToUse;i++) {
    clase = Classifiers.at(i)->Classify(ElementIndex);
    votes[clase] += weights[i];

    int cl = WhichClass(votes, winners);

    //int cl = WhichClass(votes, ((NomData*)data)->NumClass);
    //Si hay empate el clasificador clasifica como el anterior ensemble
    if (winners[0] > 1) classes->push_back((*classes)[classes->size()-1]);
    else                classes->push_back(cl);
    if (indclss)        indclss->push_back(clase);
  }

  delete []winners;

/*  if (UseWeights) {
    for(int i=0;i<nClassifiersToUse;i++) {
      clase = Classifiers.at(i)->Classify(ElementIndex);
      votes[clase] += ClassifierWeights.at(i);
      int cl = WhichClass(votes, ((NomData*)data)->NumClass);
      //Si hay empate el clasificador clasifica como el anterior ensemble
      if (cl<0) classes->push_back((*classes)[classes->size()-1]);
      else      classes->push_back(cl);
    }
  }
  else {
    for(int i=0;i<nClassifiersToUse;i++) {
      clase = Classifiers.at(i)->Classify(ElementIndex);
      votes[clase] += 1.0;
      int cl = WhichClass(votes, ((NomData*)data)->NumClass);
      //Si hay empate el clasificador clasifica como el anterior ensemble
      if (cl<0) classes->push_back((*classes)[classes->size()-1]);
      else      classes->push_back(cl);
    }
  }*/
//  delete []votes;
}
//---------------------------------------------------------------------------

void Ensemble::SecuencialError(int first, int last, vector<double> *errors,
                            vector<double> *inderrs, vector<int> *final_class)
{
  double TotalWeight=0.0;
  vector<int> *classes = new vector<int>;
  vector<int> *indclss = inderrs || final_class ? new vector<int> : 0;

  
  int nClassifiersToUse = GetClassifiersToUse();

  errors->clear();
  errors->resize(nClassifiersToUse, 0.0);
  if (inderrs) {
    inderrs->clear();
    inderrs->resize(nClassifiersToUse, 0.0);
  }
  if (final_class) {
    final_class->clear();
    final_class->resize(last-first+1, 0);
  }

  unsigned char **clases = new unsigned char*[nClassifiersToUse];
  for(int j = 0; j < nClassifiersToUse; j++) {
    clases[j] = new unsigned char[last-first+1];
  }

  for (int i = first; i <= last; i++) {

    int RealClass = data->GetDatClass(i);
    TotalWeight += data->GetDatWeight(i);
    SecuencialClassify(i, classes, indclss);
    if (final_class) (*final_class)[i-first] = (*classes)[nClassifiersToUse-1];

    for(int j = 0; j < nClassifiersToUse; j++) {
      clases[j][i-first] = (unsigned char) classes->at(j); 
      if (classes->at(j)!=RealClass)           (*errors)[j] += data->GetDatWeight(i);
      if (inderrs && (*indclss)[j]!=RealClass) (*inderrs)[j] += 1.0;
    }

  }

  for(int j=0;j<nClassifiersToUse;j++) { 
    (*errors)[j] /= TotalWeight;
    if (inderrs) (*inderrs)[j] /= (last-first+1);
  }

  int refs[7];
  refs[0] = 100; 
  refs[1] = 500; 
  refs[2] = 1000; 
  refs[3] = 5000; 
  refs[4] = 10000; 
  refs[5] = 50000;
  refs[6] = 100000;
  for(int i = 0; i < 7; i++) {
    if (nClassifiersToUse < refs[i]) continue;
    char nf[256];
    sprintf(nf, "disagreement_%d.txt", refs[i]);
    FILE *f=fopen(nf, "w");
    for(int j = 0; j < nClassifiersToUse; j++) {
      int nd = 0;
      for(int k = 0; k < last-first+1; k++) {
        if (clases[j][k] != clases[refs[i]-1][k]) nd++;
      }
      fprintf(f, "%g ", (double)nd/(last-first+1));
    }
    fprintf(f, "\n");
    fclose(f);
  }
 
  for(int j = 0; j < nClassifiersToUse; j++) {
    delete []clases[j];
  }
  delete []clases;

  delete classes;
  if (indclss) delete indclss;
}
//---------------------------------------------------------------------
void Ensemble::OrdenarClasificadoresPorOrdenOriginal()
{
  if (OrdenOriginal.size()>0) {
    vector<Classifier*> TempClas;

    //Pesos
    ClassifierWeights = PesosOriginales;

    //Clasificadores
    TempClas = Classifiers;
    for(unsigned i=0;i<OrdenOriginal.size();i++) {
      TempClas[OrdenOriginal[i]] = Classifiers[i];
    }
    Classifiers = TempClas;

    ResetOrdenOriginal();
  }
}
//---------------------------------------------------------------------
void Ensemble::OrdenarClasificadoresPorPeso(double *W)
{
  if (OrdenOriginal.size()==0) {
    ResetOrdenOriginal();
  }
  QuickSort(Classifiers.size(), 0, GetClassifiersToUse()-1, W);
}
//---------------------------------------------------------------------
void Ensemble::QuickSort(int const AHigh, int iLo, int iHi, double *W)
{
  bool IsWMine = false;
  if (!W) {
    IsWMine = true;
    W = new double[GetClassifiersToUse()];
    for(int i=0;i<GetClassifiersToUse();i++) {
      W[i] = ClassifierWeights[i];
    }
  }

  DoQuickSort(AHigh, iLo, iHi, W);

  if (IsWMine) delete []W;

}
void Ensemble::DoQuickSort(int const AHigh, int iLo, int iHi, double *A)
{
  int Lo, Hi;
  double /*T,*/ Mid;
//  Classifier *C;
//  vector<double> &A = ClassifierWeights;

  Lo = iLo;
  Hi = iHi;
  Mid = A[(Lo+Hi)/2];

  do
  {
    while (A[Lo] > Mid)
        Lo++;
    while (A[Hi] < Mid)
        Hi--;
    if (Lo <= Hi)
    {
      double aux = A[Lo];
      A[Lo] = A[Hi];
      A[Hi] = aux;
      Exchange(Lo, Hi);
      Lo++;
      Hi--;
    }
  }
  while (Lo <= Hi);

  if (Hi > iLo)
    DoQuickSort(AHigh, iLo, Hi, A);
  if (Lo < iHi)
    DoQuickSort(AHigh, Lo, iHi, A);
}
//---------------------------------------------------------------------
void Ensemble::Guardar(std::ostream &salida, int version)
{
  unsigned int num;
//  unsigned char *usados = 0;

  Classifier::Guardar(salida, version);

  num = APrioriClas.size();
  salida.write((char*)&num, sizeof(unsigned int));
  for(unsigned int i=0;i<num;i++) {
    double Val = APrioriClas[i];
    salida.write((char*)&Val, sizeof(double));
  }

  int ncls = Classifiers.size();
  salida.write((char*)&ncls, sizeof(unsigned int));

/*  if (version==1) {
    //Guardamos quE datos utilizO cada clasf en su entrenamiento	// V1
    unsigned duti = DatosUtilizados.size();				// |
    salida.write((char*)&duti, sizeof(unsigned));			// |
    map<int, int> ix;							// |
    for(unsigned i=0;i<duti;i++) {					// |
      ix[DatosUtilizados[i]] = i;					// |
      int v = DatosUtilizados[i];					// |
      salida.write((char*)&v, sizeof(unsigned));			// |
    }									// V
    usados = new unsigned char[(7+duti)/8];				// V1
  }*/

  for(int k=0;k<ncls;k++) {
/*    if (version==1) {
      //Guardamos quE datos utilizO cada clasf en su entrenamiento	// V1
      memset(usados, 0, (7+duti)/8);					// |
      for(int i=0;i<Classifiers[k]->OriginalTrainingDataCount();i++) {	// |
        int j = ix[Classifiers[k]->OriginalTrainingDataPos(i)];		// |
        usados[j/8] = usados[j/8] | 1<<(7-j%8);				// |
      }									// V
      salida.write((char*)usados, (7+duti)/8);				// V1
    }*/
    Classifiers[k]->Guardar(salida, version);
  }

/*  if (usados)  delete []usados;*/

  num = ClassifierWeights.size();
  salida.write((char*)&num, sizeof(unsigned int));
  for(unsigned int i=0;i<num;i++) {
    double Val = ClassifierWeights[i];
    salida.write((char*)&Val, sizeof(double));
  }

  salida.write((char*)&UseWeights, sizeof(bool));
  if (version==0) {
    salida.write((char*)&ClassWithWeigths, sizeof(int));
    salida.write((char*)&ClassWithoutWeigths, sizeof(int));
  }
  salida.write((char*)&ClassifiersToUse, sizeof(int));
}
//---------------------------------------------------------------------
void Ensemble::Leer(std::istream &in, int version)
{
  unsigned int Tam;
  double Val;
//  unsigned char *usados = 0;

  Classifier::Leer(in, version);

  APrioriClas.clear();
  in.read((char*)&Tam, sizeof(unsigned int));
  for(unsigned int i=0;i<Tam;i++) {
    in.read((char*)&Val, sizeof(double));
    APrioriClas.push_back(Val);
  }

  //Borramos los clasf si los haylos
  for(unsigned i=0;i<Classifiers.size();i++)
    delete Classifiers.at(i);
  Classifiers.clear();

  in.read((char*)&Tam, sizeof(unsigned int));

/*  if (version==1) {
    //Leemos quE datos utilizO cada clasf en su entrenamiento		// V1
    unsigned duti;							// |
    in.read((char*)&duti, sizeof(unsigned));				// |
    for(unsigned i=0;i<duti;i++) {					// |
      int v;								// |
      in.read((char*)&v, sizeof(int));					// |
      DatosUtilizados.push_back(v);					// |
    }									// V
    usados = new unsigned char[(7+duti)/8];				// V1
  }*/

  for(unsigned int i=0;i<Tam;i++) {
/*    if (version==1) {
      //Guardamos quE datos utilizO cada clasf en su entrenamiento	// V1
      memset(usados, 0, (7+duti)/8);					// |
      for(int i=0;i<Classifiers[k]->OriginalTrainingDataCount();i++) {	// |
        int j = ix[Classifiers[k]->OriginalTrainingDataPos(i)];		// |
        usados[j/8] = usados[j/8] | 1<<(7-j%8);				// |
      }									// V
      salida.write((char*)usados, (7+duti)/8);				// V1
    }*/
    Classifier * clas = LeerClasificadorGeneral(in);
    Classifiers.push_back(clas);
    clas->DatosUtilizados.size();
  }

  in.read((char*)&Tam, sizeof(unsigned int));
  for(unsigned int i=0;i<Tam;i++) {
    in.read((char*)&Val, sizeof(double));
    ClassifierWeights.push_back(Val);
  }

  in.read((char*)&UseWeights, sizeof(bool));
  if (version==0) {
    in.read((char*)&ClassWithWeigths, sizeof(int));
    in.read((char*)&ClassWithoutWeigths, sizeof(int));
  }
  in.read((char*)&ClassifiersToUse, sizeof(int));

  ResetOrdenOriginal();
}
//---------------------------------------------------------------------------
void Ensemble::SetData(Data *dat)
{
  for(unsigned i=0;i<Classifiers.size();i++) {
    Classifiers[i]->SetData(dat);
  }
  Classifier::SetData(dat);
}
//---------------------------------------------------------------------------
string Ensemble::Info(int value)
{
  if (value == 735) {//Orden original de los clasificadores
    ostringstream buf;
    vector<int> oo = DameOrdenOriginal();
    for (unsigned i=0;i<oo.size();i++) { 
      buf << setw(6) << oo[i];
    }
    buf << endl;
    for (int i = 0; i < Count(); i++) {
      long int k = (long int)Classifiers[i];
      k = (k<<16)>>16;
      buf << uppercase << hex << setw(6) << k;
    }
    buf << endl;
    return buf.str();
  }
  else if ( value == 33 ) {//Cuenta los nodeos de los arboles
    ostringstream buf;
    for(int i=0;i<GetClassifiersToUse();i++) {
      buf << Classifiers[i]->Info(value);
    }
    return buf.str();
  }
  else if ( value == 843 ) { //Saca estadisticas de atributos elegidos
    int natts = data->GetNumVar() - 1;
    int **atts = new int*[201];
    atts[200] = 0; // Testigo
    for(int i = 0; i < 200; i++ ) {
      atts[i] = new int[natts];
      for(int j = 0; j < natts; j++ ) {
        atts[i][j] = 0;
      }
    }
    int depth = 0;
    for(int i = 0; i < GetClassifiersToUse(); i++) {
      int K = ((DecisionTree*)Classifiers[i])->GetK();
      Node *n = ((DecisionTree*)Classifiers[i])->GetTree()->GetRoot();
      while (n) {
        if(!n->IsLeaf(K)) {
          if ( n->depth > depth ) depth = n->depth;
          atts[ n->depth ][ n->att ]++;
          n = n->child;
        }
        else {
          while ( n && !n->sib) n = n->parent;
          if ( n ) n = n->sib;
        }
      }
    }
    ostringstream buf;
    for(int i = 0; i <= depth; i++) {
      for(int j = 0; j < natts; j++ ) {
        buf << atts[i][j] << " ";
      }
      buf << endl;
    }
    return buf.str();
  }

  string res = Classifier::Info(value);
  char kk[90];

  sprintf(kk, "%d", GetClassifiersToUse());
  res = res + " -UseWeights: " + (UseWeights?"T":"F") +
              " -# of Classifiers: " + kk;
  res = res + " -PriorProbs: ";
  for(unsigned i=0;i< APrioriClas.size();i++){
    sprintf(kk, "%g, ", (float)APrioriClas[i]);
    res = res + kk; 
  }
  if (value>1) {
    for(int i=0;i<GetClassifiersToUse();i++) {
      res = res + "\n" + Classifiers[i]->Info(value);
    }
    string sep = "";
    res = res + "\nClassifiersWeights: ";
    for(int i=0;i<GetClassifiersToUse();i++) {
      sprintf(kk, "%g", ClassifierWeights[i]);
      res = res + sep + kk;
      sep = ", ";
    }
  }

  return res;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int **Ensemble::MatClasif0(Data *data)
{
  SetData(data);
  vector<Classifier*> clsfs = Classifiers;
  clsfs.resize(GetClassifiersToUse());
  return MatClasif(clsfs, data);
}
/*void Ensemble::ResetOriginalClassifierTrainingData(int Indice)
{

} */

/** 
 * Returns the value of the strengh/diversity used in 
 * Zhang 06 Ensemble Prunning via SDP
 */
double Ensemble::computeGValue(Data *data) {
  double **G = computeGMatrix(data);
  int nClasf = this->GetClassifiersToUse();

  // gValue is the sum of all the values in G.
  // As G is symmetric, we can add the diagonal terms
  // a the value of the upper triangular matrix multiplied by 2.
  double gValue = 0;
  for (int i = 0; i < nClasf; ++i) {
    // The diagonal term
    gValue += G[i][i];
    for (int j = i + 1; j < nClasf; ++j) {
       gValue += ( 2 * G[i][j] );
    } 
    delete[] G[i];
  }
  delete[] G;

  return gValue;
}

/**
 * Returns the normalized version of the G matrix defined in 
 * Zhang 06 Ensemble Prunning via SDP
 */
double** Ensemble::computeGMatrix(Data *data) {
  // The G Matrix
  int nClasf = this->GetClassifiersToUse();

  double **G = new double*[nClasf]; //(double**) malloc(sizeof(double*) * nClasf);
  for (int i = 0; i < nClasf; ++i) G[i] = new double[nClasf]; // (double *) malloc(sizeof(double) *nClasf);


  int ini = 0;
  int fin = data->GetNTotal()-1;
  int nDatos = fin-ini+1;
  int *trueClass = new int[ nDatos ];

  for (int i = 0; i < nDatos; i++) {
    trueClass[ i ] = data->GetDatClass(i);
  }

   // We pass the whole data set through each classifier
   int **Ccd = this->MatClasif0(data);

// We compute the diagonal of the G matrix and the off diagonal terms
  for (int i = 0; i < nClasf; ++i) {
        G[i][i] = Classifier::CommonErrors(Ccd[i], Ccd[i], trueClass, nDatos);
  }
  for (int i = 0; i < (nClasf - 1); i++) {
    for (int j = (i + 1); j < nClasf; j++) {
      G[i][j] = Classifier::CommonErrors(Ccd[i], Ccd[j], trueClass, nDatos);
      G[j][i] = G[i][j];
    }
  }

// We compute the normalized ~G matrix
  for (int i = 0; i < (nClasf - 1); ++i) {
    for (int j = (i + 1); j < nClasf; ++j) {
      G[i][j] =
               0.5 * ( (G[j][j] == 0 ? 0 : G[j][i] / G[j][j]) +
                       (G[i][i] == 0 ? 0 : G[i][j] / G[i][i]) );
      G[j][i] = G[i][j];
    }
  }
  for (int i = 0; i < nClasf; ++i) {
    G[i][i] = G[i][i]  / nDatos;
  }

  return G;
}

/**
 * - G: A G matrix obtained with computeGMatrix.
 * - size: The dimension of G
 * - selectedClassifiers: an array of n selected classifiers.
 * - n: The number of classifiers selected.
 */
double Ensemble::evaluateGMatrix(double** G, int size, int *selectedClassifiers, int n) {
//  int size = Count();
  int* isSelected = new int[size];

  for (int i = 0; i < size; ++i) {
    isSelected[i] = 0;
  }

  for (int i = 0; i < n; ++i) {
    isSelected[selectedClassifiers[i]] = 1;
  }

  double g = evaluateGMatrix(G, size, isSelected);

  delete[] isSelected;

  return g;
}

/**
 * - G: A G matrix obtained with computeGMatrix().
 * - size: The dimension of G
 * - isSelected: An array indicating if the i classifier is selected or not (of size 'size')
 */
double Ensemble::evaluateGMatrix(double** G, int size, int *isSelected) {
//  int size = Count();

  double g = 0;
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
       g += isSelected[i] * G[i][j] * isSelected[j];
    }
  }

  return g;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/*void SimpleStatistics::Calculate() 
{
  double sum = 0.0;
  double sum2 = 0.0;

  for(unsigned i=0;i<pop.size();i++) {
    sum += pop[i];
    sum2 += pop[i]*pop[i];
  }

  mn = sum/pop.size();
  stdv = pop.size()==1 ? 0.0 : sqrt(sum2/pop.size()-mn*mn);

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
EnsembleStatistics::EnsembleStatistics(Ensemble *ens, int **Mcd)
{
}
EnsembleStatistics::EnsembleStatistics(Ensemble *ens, Data *data)
{
}
//---------------------------------------------------------------
vector<double> EnsembleStatistics::Distribs()
{
}
//---------------------------------------------------------------
vector<double> EnsembleStatistics::Certezas()
{
}
//---------------------------------------------------------------
vector<double> EnsembleStatistics::Margen()
{
}
//---------------------------------------------------------------
vector<double> EnsembleStatistics::IndivError()
{
  vector<double> err;
  double sum, NDatos;

  for(int i=0;i<ens->Count();i++) {
    for(int j=0;j<NDatos;j++) {
      sum += Mcd[i+ens->Count()][j];
    }
    err.push_back(sum/NDatos);
  }
}
//---------------------------------------------------------------
vector<double> EnsembleStatistics::SecuencialError()
{
}
//---------------------------------------------------------------
double EnsembleStatistics::Error()
{
}
//---------------------------------------------------------------

void EnsembleStatistics::ThStat(double &mean, double &stdev)
{
}
//---------------------------------------------------------------
void EnsembleStatistics::QStat(double &mean, double &stdev)
{
}
//---------------------------------------------------------------
double EnsembleStatistics::PairQStat(int i, int j)
{
}
//---------------------------------------------------------------
void EnsembleStatistics::KStat(double &mean, double &stdev)
{
}
//---------------------------------------------------------------
double EnsembleStatistics::PairKStat(int i, int j)
{
}
//---------------------------------------------------------------
*/
