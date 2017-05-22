//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#include "Arcing.h"
#include "GPTree.h"
#include "tree20.h"
#include "node20.h"
#include "util20.h"
#include "Utils.h"
#include "OutputGenerator.cpp"
#include <math.h>
#include <fstream>
#include <algorithm>


using namespace std;

BaggingReg *BaggingReg::reg = BaggingReg::autoreg();
BoostingWithBaggingInfoReg *BoostingWithBaggingInfoReg::reg = 
                                         BoostingWithBaggingInfoReg::autoreg();
BoostingReg *BoostingReg::reg = BoostingReg::autoreg();

void CalcWeight(Instances& d, int f, int c, int *RealClass, double *weights, int iw);

//---------------------------------------------------------------------------
//----------------------------------------------------  ArcingTemplate  -----
//---------------------------------------------------------------------------

ArcingTemplate::ArcingTemplate():Ensemble()
{
}

ArcingTemplate::~ArcingTemplate()
{
}
//---------------------------------------------------------------------------

void ArcingTemplate::Build(Data *data, FuncionDeProgreso *fp)
{
  bool Continue = true;
  unsigned i;

  Init(data);

  if (fp) fp->Start(Classifiers.size());

  StartArc();
  for(i = 0; Continue && i < Classifiers.size(); i++ ) {
    Data *sample = Resample(i);
    Classifiers.at(i)->Build(sample);
    bool ValidClassifier = true;
    AfterBuilding(i, Continue, ValidClassifier);
    if (!ValidClassifier) i--; //Recalculamos el clasificador de nuevo
    if (fp && !(*fp)(ValidClassifier ? 1 : 0)) break;
    if (sample!=data) {
      if (ValidClassifier) Classifiers.at(i)->SetData(data);
      delete sample;
    }
  }
  //Borra los clasificadores que sobran
  Classifiers.resize(i);

  EndArc();

  if (fp) fp->End();
}
//---------------------------------------------------------------------------
//------------------------------------------------------  GroupBagging  -----
//---------------------------------------------------------------------------

GroupBagging::GroupBagging():ArcingTemplate()
{
}

GroupBagging::~GroupBagging()
{
}
//---------------------------------------------------------------------------

Data* GroupBagging::Resample(int)
{
  data->CreateGroups(2,0,data->GetNTrain()-1);
  return data;
}
//---------------------------------------------------------------------------
//-----------------------------------------------------------  Bagging  -----
//---------------------------------------------------------------------------

Bagging::Bagging(int _ResampleSize, bool _WithReplacement, bool _EqualizeClasses):ArcingTemplate()
{
  ResampleSize = _ResampleSize;
  WithReplacement = _WithReplacement;
  EqualizeClasses = _EqualizeClasses;
  EvidenceWeights = 0;
cout << ResampleSize << " " << WithReplacement << " " << EqualizeClasses << endl;
}

Bagging::~Bagging()
{
}
//---------------------------------------------------------------------------

Data* Bagging::Resample(int)
{
  if (WithReplacement) {
    if (ResampleSize >= 0) { 
      if (EqualizeClasses) {
        if (ResampleSize == 0) { // If zero, set it as the smallest class times # of classes
          ResampleSize = data->GetNTotal();
          for (unsigned int i = 0; i < APrioriClas.size(); i++) {
            if (APrioriClas[i] < ResampleSize) {
              ResampleSize = (int)APrioriClas[i];
            }
          }
cout << "Smallest representant has: " << ResampleSize << " elements" << endl;  
          ResampleSize =  APrioriClas.size() * ResampleSize; 
        }
        return data->GetBootstrapSampleStratified(ResampleSize);
      }
      else {
        return data->GetBootstrapSample(ResampleSize);
      }
    }
    else {
      if (EqualizeClasses) {
        return data->GetBootstrapSampleStratified();
      }
      else {
        return data->GetBootstrapSample(); //100% of data*/
      }
    } 
  }
  else {
    if ( ResampleSize > data->GetNTotal() ) {
      ResampleSize = data->GetNTotal();
    }
    data->Scramble(0, data->GetNTotal() - 1);
    data->SetNTrain(ResampleSize);

    return data;
  }
}
//---------------------------------------------------------------------------
void Bagging::AfterBuilding(int ic, bool &Continue,  bool &ValidClassifier)
{
  double ee = -1;

  Data *bsp = Classifiers[ic]->GetData();

  ee = Classifiers[ic]->Error(bsp->GetNTrain(), bsp->GetNTotal()-1);
  ClassifierWeights[ic] = ee>=0.5 ? 0.0 : ee>0.0 ? log((1.0-ee)/ee) : 10E100;
}
//---------------------------------------------------------------------------
void Bagging::EndArc()
{
  int nc = ((NomData*)data)->NumClass;
  double *w = new double[nc];
  int *RealClass = new int[data->GetNTrain()];

  for (int i = 0; i < nc; i++) w[i] = 1.0;

  NomData *ndout = new NomData(data->GetNTrain(), nc+1);
  Instances& m = ndout->GetM();

  for (int j = 0; j < data->GetNTrain(); j++) {
    m[j].SetWorkingVar(0, 0.0);//(nc, 0.0);
  }

  for (unsigned int i = 0; i < Classifiers.size(); i++) {
    for (int j = 0; j < data->GetNTrain(); j++) {
      int iori = data->GetDatIniPos(j);
      if (Classifiers[i]->UsedInOriginalTrainingData(iori)) continue;

      vector<double> d = Classifiers[i]->UnnormalizedDistribution(j);
      for(int k = 0; k < nc; k++) m[j](k, m[j][k] + d[k]);
      //m[j](nc, m[j][nc] + 1);
      m[j].SetWorkingVar(0, m[j].GetWorkingVar(0) + 1);

    }
  }

  for (int j = 0; j < data->GetNTrain(); j++) {
    for(int k = 0; k < nc; k++) { 
      m[j](k, m[j][k]/m[j].GetWorkingVar(0));//[nc]);
    }
  }

  for (int i = 0; i < data->GetNTrain(); i++) 
    RealClass[i] = data->GetDatClass(i);

  for (int i = 1; i < nc; i++) {
    CalcWeight(ndout->GetM(), ndout->GetNTotal(), i+1, RealClass, w, i);
  }

  for (int k = 0; k < nc; k++) {
    for (int i = 0; i < nc; i++) {
      CalcWeight(ndout->GetM(), ndout->GetNTotal(), nc, RealClass, w, i);
    }
  }

  EvidenceWeights = w;

  delete []RealClass;
  delete ndout;

}
//---------------------------------------------------------------------------
//-----------------------------------------------------------  Wagging  -----
//---------------------------------------------------------------------------

Wagging::Wagging():ArcingTemplate()
{
}

Wagging::~Wagging()
{
}
//---------------------------------------------------------------------------

Data* Wagging::Resample(int)
{
  //New random instance weights conforming to the continuous 
  //Poisson distribution (see webb00multiboosting)
  int NTrain = data->GetNTrain();
  double total = 0.0;
  for(int i=0; i<NTrain;i++) {
    int num = 1+(int) (999.0*rand()/(RAND_MAX+1.0));
    double w = -log((double)num/1000.0);
    data->SetDatWeight(i, w);
    total += w;
  } 
  //Normalizamos a NTrain
  for(int i=0; i<NTrain;i++) {
    double w = data->GetDatWeight(i);
    data->SetDatWeight(i, w*NTrain/total);
  }

  return data; 
}
//---------------------------------------------------------------------------
//----------------------------------------------------------  Boosting  -----
//---------------------------------------------------------------------------

Boosting::Boosting():ArcingTemplate()
{
  SetAsInWebb00multiboosting();
  SetUseWeights(true);
}

Boosting::~Boosting()
{
}
//---------------------------------------------------------------------------
void Boosting::RandomReweight(Data *data, int first, int last)
{
  double TotWeight = 0.0;
  last++;
  double TotReal = last-first;
  double newWeight;
  int iWeight = data->GetNumVar() + Data::WeightIndex;

  //Reweigth as wagging defined in webb_00_multiboosting.pdf
  for (int i=first;i<last;i++) {
    newWeight = -log((1.0 + (double)(rand() % 999))/1000.0);
    data->SetValueVar(i, iWeight, newWeight);
    TotWeight += newWeight;
  }

  //Renormalizar
  for (int i=first;i<last;i++) {
    newWeight = data->GetValueVar(i, iWeight);
    data->SetValueVar(i, iWeight, newWeight*TotReal/TotWeight);
  }
}
//---------------------------------------------------------------------------
string Boosting::Info(int value)
{
  string res = Ensemble::Info(value);
  char kk[1024];

  sprintf(kk, "\n Reweights(%lu): ", (long int)ReweightingPoints.size());
  res = res + kk;
  for (unsigned i=0;i<ReweightingPoints.size();i++) {
    sprintf(kk, "%d, ", ReweightingPoints[i]);
    res = res + kk;
  }

  return res;
}
//---------------------------------------------------------------------------
void Boosting::StartArc()
{
  data->ResetWeights();
  ReweightingPoints.clear();
  NewBootstrap = false;
}
//---------------------------------------------------------------------------
void Boosting::EndArc()
{
  data->ResetWeights();
}
//---------------------------------------------------------------------------

Data* Boosting::Resample(int iClasf)
{
  if (NewBootstrap){
    NewBootstrap = false;
    data->ResetWeights();
    return data->GetBootstrapSample();
  }

  return ArcingTemplate::Resample(iClasf);
}
//---------------------------------------------------------------------------
void Boosting::AfterBuilding(int ic, bool &Continue, bool &ValidClassifer) 
{
  double ee, newWeight;
  int NTrain = data->GetNTrain();
/*
 * A variant of Freund and Schapires (1995, 1996) AdaBoost algorithm is used.
 * The algorithm presented is Bauer and Kohavis (1999) variant of the original
 * AdaBoost.M1. This variant is chosen because:
 *    1.- It uses a one step weight update process that is less subject to
 *        numeric underflow than the original two step process.
 *    2.- It prevents numeric underflow .
 */
    Classifier *Clasf = Classifiers[ic];
    Clasf->SetData(data);
    ee = Clasf->Error(0, NTrain-1);

    /* Classifier weight calculation */
    if      (ee>=0.5)                             ClassifierWeights[ic] = 0.0;
    else if (ee> 0.0)                ClassifierWeights[ic] = log((1.0-ee)/ee);
    else if (ee==0.0 && AsInWebb00multiboosting) 
                                           ClassifierWeights[ic] = log(1E+10);
    else   /*ee==0.0*/                          ClassifierWeights[ic] = 1E100;

    // Check halt condition
    if (ee==0.0) {
      Continue = false;
      if (ForceContinuationByReweighting) {
        Continue       = true;
        NewBootstrap   = true;
        ValidClassifer = true;
        printf("\b0"); 
        ReweightingPoints.push_back(ic+1);
      }
      return;
    }
    else if (ee>=0.5) {
      Continue = false;
      if (ForceContinuationByReweighting) {
        Continue       = true;
        NewBootstrap   = true;
        ValidClassifer = false;
        printf("\b*");
        ReweightingPoints.push_back(ic);
      }
      return;
    }

    // Reweighting of the examples
    double KMal = 1.0/(2.0*ee);         //Mal clasificado
    double KBien = 1.0/(2.0*(1.0-ee));  //Bien clasificado
    double TotWeight = 0.0;
    double TotReal = NTrain;

    //Cambiamos los pesos de los ejemplos
    for (int i=0;i<NTrain;i++) {
      newWeight = data->GetDatWeight(i);
      newWeight *= Clasf->Classify(i) != data->GetDatClass(i) ? KMal : KBien;
//if (Clasf->Classify(i) != data->GetDatClass(i))
//printf("%g ", data->GetDatWeight(i));
      if (newWeight<MinWeight) newWeight = MinWeight;
      data->SetDatWeight(i, newWeight);
      TotWeight += newWeight;
    }
    //Renormalizar for si acaso
//FILE *ff=fopen("weights.csv","a");
//data->OriginalOrder();
    for (int i=0;i<NTrain;i++) {
      newWeight = data->GetDatWeight(i);
      data->SetDatWeight(i, newWeight*TotReal/TotWeight);
//fprintf(ff, "%g\t", data->GetDatWeight(i));
    }
//fprintf(ff, "\n");
//fclose(ff);

    //Options
    if (AsInWebb00multiboosting) {
      //As in webb00multiboosting, to allow finding a new 0 error classifier
      //if (ee==0.0) ClassifierWeights[ic] = log(1E+10);
    }
    if (AsInBauer99empirical) {
    }
}
//---------------------------------------------------------------------------
void Boosting::SetAsInWebb00multiboosting()
{
  AsInWebb00multiboosting = true;
  AsInBauer99empirical = false;
  ForceContinuationByReweighting = true;
  MinWeight = 1e-8;
}
//---------------------------------------------------------------------------
void Boosting::SetAsInBauer99empirical()
{
  AsInWebb00multiboosting = false;
  AsInBauer99empirical = true;
  ForceContinuationByReweighting = false;
  MinWeight = 1e-6;
} 
//---------------------------------------------------------------------------
void Boosting::SetAdaBoost()
{
  AsInWebb00multiboosting = false;
  AsInBauer99empirical = false;
  ForceContinuationByReweighting = false;
  MinWeight = 0.0;
}
//---------------------------------------------------------------------------
void Boosting::SetDefaultClassifiers()
{
  int num = ReweightingPoints.size()>0 ? ReweightingPoints[0] : Count();
  if (ClassifierWeights[num-1]==log(1E+10)) ClassifierWeights[num-1] = 1E100;
  SetClassifiersToUse(num); 
} 
//---------------------------------------------------------------------------
void Boosting::Guardar(std::ostream &salida, int version)
{
  ArcingTemplate::Guardar(salida, version);

  unsigned int num = ReweightingPoints.size();
  salida.write((char*)&num, sizeof(unsigned int));
  for(unsigned int i=0;i<num;i++) {
    int Val = ReweightingPoints[i];
    salida.write((char*)&Val, sizeof(int));
  }
  
}
//---------------------------------------------------------------------------
void Boosting::Leer(std::istream &in, int version)
{
  ArcingTemplate::Leer(in, version);

  unsigned int Tam;
  in.read((char*)&Tam, sizeof(unsigned int));
  for(unsigned int i=0;i<Tam;i++) {
    int Val;
    in.read((char*)&Val, sizeof(int));
    ReweightingPoints.push_back(Val);
  }
}
//---------------------------------------------------------------------------
//--------------------------------------------------------  MIBoosting  -----
//---------------------------------------------------------------------------

MIBoosting::MIBoosting():Boosting()
{
  SetAsInWebb00multiboosting();
  SetUseWeights(true);
  Weights = 0;
  Resampling = true;
}

MIBoosting::~MIBoosting()
{
  if (Weights) delete[] Weights;
}
//---------------------------------------------------------------------------
void MIBoosting::ResetWeights()
{
  if (Weights) delete []Weights;

  Weights = new double[((MINomData*)data)->GetMITrain()];
  for(int i = 0; i < ((MINomData*)data)->GetMITrain(); i++) {
    Weights[i] = 1.0;
  }

  for(int i = 0; i < data->GetNTrain(); i++) {
    data->SetDatWeight(i, Weights[data->GetDatGroup(i)]);
  }
}
//---------------------------------------------------------------------------
void MIBoosting::StartArc()
{
  ResetWeights();
  ReweightingPoints.clear();
  NewBootstrap = false;
}
//---------------------------------------------------------------------------
void MIBoosting::EndArc()
{
  ResetWeights();
}
//---------------------------------------------------------------------------

Data* MIBoosting::Resample(int iClasf)
{
  Data *d = 0;

  if (Resampling) {
    for(int i = 0; i < ((MINomData*)data)->GetMITrain(); i++) {
      ((MINomData*)data)->SetWeight(i, Weights[i]);
    }
    d = data->GetBootstrapSample(0, true, false); 
  }
  else if (NewBootstrap){
    NewBootstrap = false;
    ResetWeights();
    d = data->GetBootstrapSample(); 
  }
  else {
    d = ArcingTemplate::Resample(iClasf);
  }

  return d;
}
//---------------------------------------------------------------------------
void MIBoosting::AfterBuilding(int ic, bool &Continue, bool &ValidClassifer) 
{
    double ee;
    int NTrain = data->GetNTrain();
    int NMITrain = ((MINomData*)data)->GetMITrain();
    int nClasses = ((NomData*)data)->NumClass;
    double TotWeight = 0.0;
 
    vector<double>** dis;

    Classifier *Clasf = Classifiers[ic];
    Clasf->SetData(data);

    int *RealClass = new int[ NMITrain ];
    int *PredictedClass = new int[ NMITrain ];

    dis = new vector<double>*[ NMITrain ];
    for(int i = 0; i < NMITrain; i++) {
      dis[i] = new vector<double>(nClasses, 0.0);
    }

    for(int i = 0; i < NTrain; i++) {
      vector<double> id = Clasf->UnnormalizedDistribution(i);
      vector<double>* di = dis[ data->GetDatGroup(i) ]; 
      for(int j = 0; j < nClasses; j++) {
        (*di)[j] = (*di)[j] + id[j];
      }
      RealClass[data->GetDatGroup(i)] = data->GetDatClass(i);
    }

    ee = 0.0;
int e21=0;
int e12=0;
    for(int i = 0; i < NMITrain; i++) {
      PredictedClass[i] = WhichClass(*(dis[i]));
      delete dis[i];
      TotWeight += Weights[i];
      if ( PredictedClass[i] != RealClass[i] ) {
        ee += Weights[i];
if (RealClass[i] == 0) e12++;
else e21++;
      }
    }
    ee /= TotWeight;
    delete []dis;
cout << "E" << ee << "(e21=" << e21 << ";e12=" << e12 << ")" << endl;
    /* Classifier weight calculation */
    if      (ee>=0.5)                             ClassifierWeights[ic] = 0.0;
    else if (ee> 0.0)                ClassifierWeights[ic] = log((1.0-ee)/ee);
    else if (ee==0.0 && AsInWebb00multiboosting) 
                                           ClassifierWeights[ic] = log(1E+10);
    else   /*ee==0.0*/                          ClassifierWeights[ic] = 1E100;

    // Check halt condition
    if (ee==0.0) {
      Continue = false;
      if (ForceContinuationByReweighting) {
        Continue       = true;
        NewBootstrap   = true;
        ValidClassifer = true;
        printf("\b0"); 
        ReweightingPoints.push_back(ic+1);
      }
    }
    else if (ee>=0.5) {
      Continue = false;
      if (ForceContinuationByReweighting) {
        Continue       = true;
        NewBootstrap   = true;
        ValidClassifer = false;
        printf("\b*");
        ReweightingPoints.push_back(ic);
      }
    }
    else {
      // Reweighting of the examples
      double KMal = 1.0/(2.0*ee);         //Mal clasificado
      double KBien = 1.0/(2.0*(1.0-ee));  //Bien clasificado

      //Cambiamos los pesos de los miejemplos
      TotWeight = 0.0;
      for(int i = 0; i < NMITrain; i++) {
        Weights[i] *=  PredictedClass[i] != RealClass[i] ? KMal : KBien;
        if ( Weights[i] < MinWeight ) Weights[i] = MinWeight;
        TotWeight += Weights[i];
      }
      for (int i = 0; i < NMITrain; i++) {
        Weights[i] *= NMITrain;  
        Weights[i] /= TotWeight;  
      }
      for (int i = 0; i < NTrain; i++) {  
        data->SetDatWeight(i, Weights[data->GetDatGroup(i)]);
        TotWeight += Weights[data->GetDatGroup(i)];
      }

      for (int i = 0; i < NTrain; i++) {  
        double w = data->GetDatWeight(i);
        data->SetDatWeight(i, w * NTrain / TotWeight);
      }
    }

    delete []RealClass;
    delete []PredictedClass;

}
//---------------------------------------------------------------------------
double MIBoosting::Error(int first, int last)
{
  int nclases = ((NomData*)data)->NumClass;

  vector<double> &weights = GetUsingWeights();

  data->MarkOrder();

  double tw = 0.0;
  double err = 0.0;

  for(int i = first; i <= last; i++) {
    vector<double> votes(nclases, 0.0);
    do {
      for(int j = 0; j < GetClassifiersToUse(); j++) {
        vector<double> d = Classifiers.at(j)->UnnormalizedDistribution(i);
        for(int k = 0; k < nclases; k++) {
          votes[k] = votes[k] + d[k] * weights[j];
        } 
      }
      i++;
    } while ( i <= last && data->GetDatGroup(i) == data->GetDatGroup(i - 1) );

    i--;

    int cl = WhichClass(votes);

    if (cl != data->GetDatClass(i-1) ) {
      err += data->GetDatWeight(i-1);
    } 

    tw +=  data->GetDatWeight(i-1);
  }
  data->ResetOrder();

  return err/tw;
}

//---------------------------------------------------------------------------
void MIBoosting::SecuencialClassify(int ElementIndex, vector<int> *classes,
                                                     std::vector<int> *indclss)
{
  int nclases = ((NomData*)data)->NumClass;
  int j;
  int *winners;

  vector<double> votes(nclases, 0.0);
  int nClassifiersToUse = GetClassifiersToUse();
  classes->clear();
  if (indclss) indclss->clear();

  winners = new int[nclases+1];

  vector<double> &weights = GetUsingWeights();
  for(int i = 0; i < nClassifiersToUse; i++) {
    j = ElementIndex;
    do {
      vector<double> d = Classifiers.at(i)->UnnormalizedDistribution(j);
      for(int k = 0; k < nclases; k++) {
        votes[k] = votes[k] + d[k] * weights[i];
      }
      j++;
    }
    while( j < data->GetNTotal() && data->GetDatGroup(j) == data->GetDatGroup(j-1) );

    int cl = WhichClass(votes, winners);

    //Si hay empate el clasificador clasifica como el anterior ensemble
/*    if (winners[0] > 1 && i > 0 ) classes->push_back((*classes)[i-1]);
    else  */                        classes->push_back(cl);

    if (indclss)        indclss->push_back(cl);
  }

  delete[] winners;

}
//---------------------------------------------------------------------------

void MIBoosting::SecuencialError(int first, int last, vector<double> *errors,
                            vector<double> *inderrs, vector<int> *final_class)
{
  double TotalWeight=0.0;
  vector<int> *classes = new vector<int>;
  vector<int> *indclss = inderrs || final_class ? new vector<int> : 0;

  int nClassifiersToUse = GetClassifiersToUse();

  data->SortByGroup(first, last);

  errors->clear();
  errors->resize(nClassifiersToUse, 0.0);
  if (inderrs) {
    inderrs->clear();
    inderrs->resize(nClassifiersToUse, 0.0);
  }
  if (final_class) {
    final_class->clear();
  }

  int last_mi = -1;
  for (int i = first; i <= last; i++) {

    if (last_mi == data->GetDatGroup(i)) continue;
    last_mi = data->GetDatGroup(i);

    int RealClass = data->GetDatClass(i);
    TotalWeight += data->GetDatWeight(i);

    SecuencialClassify(i, classes, indclss);

    if (final_class) {
      (*final_class).push_back((*classes)[nClassifiersToUse-1]);
    }
    for(int j = 0; j < nClassifiersToUse; j++) {
      if (classes->at(j)!=RealClass)
        (*errors)[j] += data->GetDatWeight(i);
      if (inderrs && (*indclss)[j]!=RealClass)
        (*inderrs)[j] += 1.0;
    }
  }

  for(int j=0;j<nClassifiersToUse;j++) { 
    (*errors)[j] /= TotalWeight;
    if (inderrs) (*inderrs)[j] /= (last-first+1);
  }

  delete classes;
  if (indclss) delete indclss;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------  MIBagging  -----
//---------------------------------------------------------------------------

MIBagging::MIBagging(int _ResampleSize):MIBoosting()
{
  ResampleSize = _ResampleSize;
cout <<"RR " << ResampleSize << endl;
  SetUseWeights(false);
  MaxSingleInstanceMargin = false;
  PositiveClass = 1;
  ndout_all = 0;
  ndout = 0;
  OutputLevel = 0x0;
}

MIBagging::~MIBagging()
{
}
//---------------------------------------------------------------------------
int nn(const void* AA, const void *BB)
{
  int ret;
  double A = **((double**) AA);
  double B = **((double**) BB);

  if (A > B) ret = 1;
  else if (A < B) ret = -1;
  else ret = 0;

  return ret;
}
//---------------------------------------------------------------------------
void CalcWeight(Instances& d, int f, int c, int *RealClass, double *weights, int iw)
{
  double **iww;
  int error;
  int *PredictedClass;
  
  PredictedClass = new int[f];
  iww = new double*[f];
  error = 0;
  for(int i = 0; i < f; i++) {
    PredictedClass[i] = -1;
    double mmax = -1;
    for(int j = 0; j < c; j++) {
      if (j == iw) continue;
      double iwww = d[i][j] * weights[j];
      if (iwww > mmax && iwww > 0) {
        mmax = iwww;
        PredictedClass[i] = j;
      }
    }
    iww[i] = new double[2];
    iww[i][0] = mmax / d[i][iw];
    iww[i][1] = i; 
    if (PredictedClass[i] != RealClass[i]) error++;
  }

  qsort(iww, f, sizeof(double*), nn);

  int iiww = -1;
  int iiwwM = -1;
  int minerror = error;
  for(int i = 0; i < f; i++) {
    if (PredictedClass[(int)iww[i][1]] == RealClass[(int)iww[i][1]]) error++;
    else if (iw == RealClass[(int)iww[i][1]]) error--;
    if (error < minerror) {
      minerror = error;
      iiww = i;
      iiwwM = i == f-1 ? i : i + 1;
    } 
    else if (error == minerror) {
      iiwwM = i == f-1 ? i : i + 1;
    }
  }
  if (iiww > -1) {
 //   if (iiww == f-1) {
 //     weights[iw] = iww[iiww][0] + 0.1;
 //   }
 //   else {
      weights[iw] = (iww[iiww][0]+iww[iiwwM][0]) / 2;
 //   }
  }
  delete []PredictedClass;
  for(int i = 0; i < f; i++)
    delete []iww[i];
  delete []iww;

}
//---------------------------------------------------------------------------
void MIBagging::StartArc()
{
  int NTrain = data->GetNTotal();
  int NMITrain = ((MINomData*)data)->GetMITrain();
//ew  int n_atts_out = 1 + ((NomData*)data)->NumClass * (((NomData*)data)->NumClass + 1) / 2;
  int n_atts_out = 1 + ((NomData*)data)->NumClass;

//  ndout_all = new int*[NTrain];
  ndout_all_votes = new int*[NTrain];
  for(int i = 0; i < NTrain; i++) {
//    ndout_all[i] = new int[n_atts_out];
    ndout_all_votes[i] = new int[n_atts_out-1];
    for(int j = 0; j < n_atts_out-1; j++) {
//      ndout_all[i][j] = 0;
      ndout_all_votes[i][j] = 0;
    }
//    ndout_all[i][n_atts_out-1] = 0;
  }
  ndout_all = new NomData(NTrain, n_atts_out+1);//, new FloatInstance(0));

  ndout = new NomData(NMITrain, n_atts_out);
  for(int i = 0; i < NMITrain; i++) {
    ndout->GetM()[i].SetWorkingVar(0, 0.0);
    ndout->GetM()[i].SetWorkingVar(1, 0.0);
    //ndout->SetValueVar(i, n_atts_out+1,  0.0);
  }
}
//---------------------------------------------------------------------------
bool FindIdx(vector<int> &idxs_oob, int to_find)
{
  return binary_search(idxs_oob.begin(), idxs_oob.end(), to_find);
}
//---------------------------------------------------------------------------
void  AllocTreeInfoStructures(DecisionTree* tree, int nqq, int nc)
{
  Node* cur = tree->GetTree()->GetRoot();

  while (cur->child) { 
    cur = cur->child;
  }
  while (cur != tree->GetTree()->GetRoot()) {
    if ( cur->child == 0 ) {
      cur->info = new Matriz(nqq, nc);
    }
    cur = cur->nextDown();
  }
}
//---------------------------------------------------------------------------
void MIBagging::EndArc()
{
  int nc = ((NomData*)data)->NumClass;
//ew  CalcDiffs();
  CalcDiffs(ndout, 1);
//KK2  CalcDiffs(ndout_all, 1);

//  string n = GenerateProcNameFileName("train_",".cre");
//  ndout->SaveToFile((char*)n.c_str());

  //CalcDiffs(ndout, 3); //MEAN
  //n = GenerateProcNameFileName("train_mean_",".cre");
  //ndout->SaveToFile((char*)n.c_str());

#ifdef KK2
  if (MaxSingleInstanceMargin) {
    ndout_all->SortByGroup();
    int igrp = 0;
    for(int i = 0; i < ndout_all->GetNTotal(); i++) {
      double margin = -1000;
      do {
        double mi = Classifier::Margin(ndout_all->GetInstance(i), PositiveClass, false);
        if ( mi > margin ) {
          for(int j = 0; j < nc; j++) {
            ndout->SetValueVar(igrp, j, ndout_all->GetValueVar(i, j));
          }
          ndout->SetValueVar(igrp, nc,  ndout_all->GetDatClass(i));
          margin = mi;
        }
        i++;
      } while(i < ndout_all->GetNTotal() && ndout_all->GetDatGroup(i) == ndout_all->GetDatGroup(i-1));
      igrp++;
      i--;
    }
  }
#endif

  double *w = new double[nc];
  int *RealClass = new int[ndout->GetNTotal()];

  for (int i = 0; i < nc; i++) w[i] = 1.0;
  for (int i = 0; i < nc; i++) { cout << w[i] << " "; } cout << endl;

  if ( ! MaxSingleInstanceMargin ) {
    for (int i = 0; i < ndout->GetNTotal(); i++) RealClass[i] = ndout->GetDatClass(i);
cout << endl << "1.0 ";

    for (int i = 1; i < nc; i++) {
      CalcWeight(ndout->GetM(), ndout->GetNTotal(), i+1, RealClass, w, i);
    }
    for (int i = 0; i < nc; i++) { cout << w[i] << " "; } cout << endl;

    for (int i = 0; i < nc; i++) {
      CalcWeight(ndout->GetM(), ndout->GetNTotal(), nc, RealClass, w, i);
    }
    for (int i = 0; i < nc; i++) { cout << w[i] << " "; } cout << endl;

    for (int i = 0; i < nc; i++) {
      CalcWeight(ndout->GetM(), ndout->GetNTotal(), nc, RealClass, w, i);
    }
    for (int i = 0; i < nc; i++) { cout << w[i] << " "; } cout << endl;

    for (int i = 0; i < nc; i++) {
      CalcWeight(ndout->GetM(), ndout->GetNTotal(), nc, RealClass, w, i);
    }
    for (int i = 0; i < nc; i++) { cout << w[i] << " "; } cout << endl;

  }

  for (int i = 0; i < (int)Classifiers.size(); i++) {
    for (int j = 0; j < nc; j++) {
      DistributionWeights[i][j] = w[j];
    }
  }


  string n = GenerateProcNameFileName("sifts_train_",".txt");
  FILE *f=fopen(n.c_str(), "w");
  n = GenerateProcNameFileName("votes_train_",".txt");
  FILE *f2=fopen(n.c_str(), "w");
  data->OriginalOrder();
//  CalcDiffs(ndout_all, 2);
  for(int i = 0; i < data->GetNTotal(); i++) {
    fprintf(f, "%d %d %d ", data->GetDatIniPos(i), data->GetDatGroup(i), data->GetDatClass(i));
    fprintf(f2, "%d %d %d ", data->GetDatIniPos(i), data->GetDatGroup(i), data->GetDatClass(i));
    for(int k = 0; k < nc; k++) {
      ndout_all->SetValueVar(i, k, ndout_all->GetValueVar(i, k)/ndout_all->GetValueVar(i, nc));
      fprintf(f, "%g ", (double)ndout_all->GetValueVar(i, k));
      fprintf(f2, "%d ", ndout_all_votes[i][k]);
    }
    fprintf(f, "\n ");
    fprintf(f2, "\n ");
  }
  fclose(f);
  fclose(f2);

  if (OutputLevel & 0x8) {
    int desde, hasta, n_ins_clase;
    int nqq = 10;

    q = new Matriz(nc, nqq);

    // Get the quantiles for each class
    desde = 0;
    ndout_all->SortByClass(); 
    for(int i = 0; i < nc; i++) {
      hasta = desde;
      while (hasta < ndout_all->GetNTotal() && 
        ndout_all->GetDatClass(hasta)==ndout_all->GetDatClass(desde)) hasta++;

      ndout_all->SortOn(i, desde, hasta-1);
      n_ins_clase = hasta - desde;
      for(int j = 0; j < nqq; j++) {
        (*q)[i][j] = ndout_all->GetValueVar(desde + round((double)n_ins_clase*j/nqq), i);
cout << (*q)[i][j] << " ";
      }
cout << endl;
      desde = hasta;
    }

    // Compute the node populations
    ndout_all->OriginalOrder();
    data->OriginalOrder();
    for(int j = 0; j < GetClassifiersToUse(); j++) {
      DecisionTree *tree = (DecisionTree*)GetClassifier(j);
      AllocTreeInfoStructures(tree, nqq, nc);
      tree->SetData(data);
      for(int i = 0; i < data->GetNTotal(); i++) {
        if (FindIdx(idxs_oobs[j], data->GetDatGroup(i))) continue;
        do {
          tree->Classify(i);
          Matriz *n = (Matriz*)tree->GetLastClassifyingNode()->info;
          int ic = ndout_all->GetDatClass(i);
          for(int k = 0; k < nqq; k++) {
            if ((*q)[ic][k] <= ndout_all->GetValueVar(i, ic)) {
              (*n)[k][ic]++;
            }
            else {
              break;
            }
          }
          i++;
        } while(i < data->GetNTotal() && data->GetDatGroup(i) == data->GetDatGroup(i-1));
        i--;
      }
    }

    // Compute the oob training set for all quantiles
    FILE *ff[30];
    for (int i = 0; i < nqq; i++) {
      n = GenerateNewFileName("evidence_filtered_train_",".txt");
      ff[i] = fopen(n.c_str(), "w");
    }
    for(int i = 0; i < data->GetNTotal(); i++) {
      int ic = ndout_all->GetDatClass(i);
      Matriz mm(nqq, nc);
      int n_trees = 0;
      for(int j = 0; j < GetClassifiersToUse(); j++) {
        if (!FindIdx(idxs_oobs[j], data->GetDatGroup(i))) continue;
        n_trees++;
        DecisionTree *tree = (DecisionTree*)GetClassifier(j);
        tree->Classify(i);
        Matriz *n = (Matriz*)tree->GetLastClassifyingNode()->info;
        for(int k = 0; k < nqq; k++) {
          if ((*q)[ic][k] <= ndout_all->GetValueVar(i, ic)) {
            for (int l = 0; l < nc; l++) {
              mm[k][l] += (*n)[k][l]; //Para noise=0 no se recupera la estadistica guardada
              //en el fichero de arriba pq aqui no sabemos cuantas veces cada instancia fue 
              //extraida en el bootstrap
            }
          }
          else {
            break;
          }
        }
      }
      for (int k = 0; k < nqq; k++) {
        fprintf(ff[k], "%d %d %d ", data->GetDatIniPos(i), data->GetDatGroup(i), data->GetDatClass(i));
        for(int l = 0; l < nc; l++) {
          fprintf(ff[k], "%g ", (double)mm[k][l]/n_trees);
        }
        fprintf(ff[k], "\n ");
      }
    }

    for (int i = 0; i < nqq; i++) {
      fclose(ff[i]);
    }


  }

  delete []w;
  delete []RealClass;
  delete ndout;
  if (ndout_all) {
    delete ndout_all;
    for (int i = 0; i < data->GetNTotal(); i++) {
//      delete []ndout_all[i];
      delete []ndout_all_votes[i];
    }
//    delete []ndout_all;
    delete []ndout_all_votes;
  }
  ndout_all = 0;
  ndout_all_votes = 0;
  ndout = 0;
}
//---------------------------------------------------------------------------
void MIBagging::CalcDiffs(Data* nd, int locualo)
{
  int n_classes = ((NomData*)data)->NumClass;

  for (int k = 0; k < nd->GetNTotal(); k++) {
    double isk = 0.0;
    if (locualo == 3) {
      for (int i = 0; i < n_classes; i++) {
        isk += nd->GetValueVar(k, i);
      }
    }
    for (int i = 0; i < n_classes; i++) {
      if (nd->GetValueVar(k, n_classes+1) != 0) {
        if (locualo == 1) {//Div by number of classifiers
//cout << nd->GetValueVar(k, i) << " ";
          nd->SetValueVar(k, i, nd->GetValueVar(k, i) / nd->GetValueVar(k, n_classes+1));
        }
        else if (locualo == 2) {//Weighted
          nd->SetValueVar(k, i, nd->GetValueVar(k, i) * DistributionWeights[0][i]);
        }
        else if (locualo == 3) {//Norm
          nd->SetValueVar(k, i, nd->GetValueVar(k, i) / isk);
        }
        else if (locualo == 4) {//Div by the number of detections
          nd->SetValueVar(k, i, nd->GetValueVar(k, i) / nd->GetValueVar(k, n_classes+2));
        }
      }
    }
//if (locualo==3)
//cout << nd->GetValueVar(k, n_classes+2) << " ";
  }
//cout << endl;
/*ew  for (int k = 0; k < nd->GetNTotal(); k++) {
    for (int i = 0; i < n_classes; i++) {
      nd->SetValueVar(k, i, nd->GetValueVar(k,i) / nd->GetValueVar(k, n_classes));
    }
    for (int i = 1; i < n_classes; i++) {
      for (int j = 0; j < i; j++) {
        double d1 = nd->GetValueVar(k, i);
        double d2 = nd->GetValueVar(k, j);
        nd->SetValueVar(k, n_classes + i*(i-1)/2 + j,  d1 - d2);
      }
    }
  }*/
}
//---------------------------------------------------------------------------

Data* MIBagging::Resample(int iClasf)
{
  Data *ret = ResampleSize > 0 ? data->GetBootstrapSample(ResampleSize, false, true, &idxs_oob) :
                            data->GetBootstrapSample(0, false, true, &idxs_oob) ;
  idxs_oobs.push_back(idxs_oob);

  return ret;
}
//---------------------------------------------------------------------------
void MIBagging::AfterBuilding(int ic, bool &Continue, bool &ValidClassifer) 
{
//    vector<double>** dis;

    Classifier *Clasf = Classifiers[ic];

    // int NMITrain = ((MINomData*)bs)->GetMITrain();
   // int NMITotal = ((MINomData*)bs)->GetMITotal();
    int nClasses = ((NomData*)data)->NumClass;
 
//    int *RealClass = new int[ NMITrain ];
  //  int *PredictedClass = new int[ NMITrain ];

  /*  dis = new vector<double>*[ NMITotal - NMITrain ];
    for(int i = 0; i < NMITotal - NMITrain; i++) {
      dis[i] = new vector<double>(nClasses, 0.0);
    }*/

#ifdef KKFUTIS
    Data *bs = Clasf->GetData();
    int NTrain = bs->GetNTrain();
    int NTotal = bs->GetNTotal();
    bs->SortByGroup(NTrain, NTotal - 1);
    int igrp = 0;
    for(int i = NTrain; i < NTotal; i++) {
      do {
        vector<double> id = Clasf->UnnormalizedDistribution(i);
        //vector<double> id = Clasf->Distribution(i);
        if (MaxSingleInstanceMargin || ndout_all) {
          for(int j = 0; j < nClasses; j++) {
            ndout_all->SetValueVar(bs->GetDatIniPos(i), j, ndout_all->GetValueVar(bs->GetDatIniPos(i), j) + id[j]);
          }
          ndout_all->SetValueVar(bs->GetDatIniPos(i), nClasses+1,  ndout_all->GetValueVar(bs->GetDatIniPos(i), nClasses+1) + 1);
          ndout_all->SetValueVar(bs->GetDatIniPos(i), nClasses,  bs->GetDatClass(i));
          ndout_all->SetDatGroup(bs->GetDatIniPos(i), bs->GetDatGroup(i));
        }
        if ( ! MaxSingleInstanceMargin ) {
/*          vector<double>* di = dis[ igrp ]; */
          for(int j = 0; j < nClasses; j++) {
            /*(*di)[j] = (*di)[j] + id[j];*/
            ndout->SetValueVar(bs->GetDatGroup(i), j,  ndout->GetValueVar(bs->GetDatGroup(i), j) + id[j]);
          }
        }
        i++;
      } while(i < NTotal && bs->GetDatGroup(i) == bs->GetDatGroup(i-1));
      i--;
     // RealClass[igrp] = bs->GetDatClass(i);
      if ( ! MaxSingleInstanceMargin ) {
//ew      ndout->SetValueVar(bs->GetDatGroup(i), nClasses*(nClasses+1)/2,  bs->GetDatClass(i));
//ew      ndout->SetValueVar(bs->GetDatGroup(i), nClasses,  ndout->GetValueVar(bs->GetDatGroup(i), nClasses) + 1);
        ndout->SetValueVar(bs->GetDatGroup(i), nClasses,  bs->GetDatClass(i));
        ndout->SetValueVar(bs->GetDatGroup(i), nClasses+1,  ndout->GetValueVar(bs->GetDatGroup(i), nClasses+1) + 1);
      }
      igrp++;
    }
#else
UseOOBForNodePops = false;
    if (UseOOBForNodePops) {
      DecisionTree *tree = (DecisionTree*)Clasf;
      // Recorro el arbol y pongo a 0 NodePops
      Node* cur = tree->GetTree()->GetRoot();

      while (cur->child) { 
        cur = cur->child;
      }
      string n = GenerateProcNameFileName("train_pops_", ".txt");
      FILE *f = ic==0 ? fopen(n.c_str(), "w") : fopen(n.c_str(), "a");
      while (cur != tree->GetTree()->GetRoot()) {
        for (int i = 0; i < nClasses; i++) {
          if ( cur->child == 0 ) {
            fprintf(f, "%g ", cur->NodePop[i]);
          }
          cur->NodePop[i] = 0.0;
        }
        if ( cur->child == 0 ) {
          fprintf(f, "\n");
        }
        cur = cur->nextDown();
      }
      for (int i = 0; i < nClasses; i++) {
        cur->NodePop[i] = 0.0;
      }
      fclose(f);

      // Lanzo las instancias oob por el arbol y actualizo las NodePops
      int igrp = 0;
      Clasf->SetData(data);
      for(int i = 0; i < data->GetNTotal() && igrp < (int)idxs_oob.size(); i++) {
        if (idxs_oob[igrp] != data->GetDatGroup(i)) continue;
        int ngrp = 0;
        do {
          tree->Classify(i);
          tree->GetLastClassifyingNode()->NodePop[data->GetDatClass(i)]++;
          i++;
          ngrp++;
        } while(i < data->GetNTotal() && data->GetDatGroup(i) == data->GetDatGroup(i-1));
        i--;
        igrp++;
      }

      // Recorro el arbol y pongo a 0 NodePops
      cur = tree->GetTree()->GetRoot();

      while (cur->child) { 
        cur = cur->child;
      }

      n = GenerateProcNameFileName("oob_pops_", ".txt");
      f = ic==0 ? fopen(n.c_str(), "w") : fopen(n.c_str(), "a");
      while (cur != tree->GetTree()->GetRoot()) {
        if (cur->child==0) {
          for (int i = 0; i < nClasses; i++) {
            fprintf(f, "%g ", cur->NodePop[i]);
          }
          fprintf(f, "\n");
        }
        cur = cur->nextDown();
      }
      fclose(f);

    }
    int igrp = 0;
    Clasf->SetData(data);
    for(int i = 0; i < data->GetNTotal() && igrp < (int)idxs_oob.size(); i++) {
//      if (idxs_oob[igrp] != data->GetDatGroup(i)) continue;
      if (!UseOOBForNodePops) {
        if (idxs_oob[igrp] != data->GetDatGroup(i)) continue;
        else igrp++;
      }
      else if (idxs_oob[igrp] == data->GetDatGroup(i)) {
        igrp++;
        i--;
        continue; // Usamos train
      }
      int ngrp = 0;
      do {
        vector<double> id = Clasf->UnnormalizedDistribution(i);
        int cls = Clasf->Classify(i);
        if (MaxSingleInstanceMargin || ndout_all) {
          int i0 = data->GetDatIniPos(i);
          ndout_all_votes[i0][cls] += 1;//KK2
          for(int j = 0; j < nClasses; j++) {
            //KK2ndout_all->SetValueVar(i0, j, ndout_all->GetValueVar(data->GetDatIniPos(i), j) + id[j]);
            ndout_all->SetValueVar(i0, j, ndout_all->GetValueVar(i0, j) + (int)id[j]);//KK2
          }
          //KK2 ndout_all->SetValueVar(i0, nClasses+1,  ndout_all->GetValueVar(i0, nClasses+1) + 1);
          ndout_all->SetValueVar(i0, nClasses, ndout_all->GetValueVar(i0, nClasses) + 1);//KK2
          ndout_all->SetDatClass(i0, data->GetDatClass(i));//KK2
          //KK2 ndout_all->SetValueVar(i0, nClasses,  data->GetDatClass(i));
          //KK2 ndout_all->SetDatGroup(i0, data->GetDatGroup(i));
        }
        if ( ! MaxSingleInstanceMargin ) {
          for(int j = 0; j < nClasses; j++) {
            ndout->SetValueVar(data->GetDatGroup(i), j,  ndout->GetValueVar(data->GetDatGroup(i), j) + id[j]);
          }
        }
        i++;
        ngrp++;
      } while(i < data->GetNTotal() && data->GetDatGroup(i) == data->GetDatGroup(i-1));
      i--;
      if ( ! MaxSingleInstanceMargin ) {
        ndout->SetValueVar(data->GetDatGroup(i), nClasses,  data->GetDatClass(i));
        //ndout->SetValueVar(data->GetDatGroup(i), nClasses+1,  ndout->GetValueVar(data->GetDatGroup(i), nClasses+1) + 1);
        ndout->GetM()[data->GetDatGroup(i)].SetWorkingVar(0, ndout->GetM()[data->GetDatGroup(i)].GetWorkingVar(0) + 1);
if (ndout->GetM()[data->GetDatGroup(i)].GetWorkingVar(1) > 0 && ndout->GetM()[data->GetDatGroup(i)].GetWorkingVar(1) != ngrp)
{
cout << "Error diff n detects for " << i << endl;
}
//        ndout->SetValueVar(data->GetDatGroup(i), nClasses+2,  ngrp);//MEAN
        ndout->GetM()[data->GetDatGroup(i)].SetWorkingVar(1,  ngrp);//MEAN
      }
    }
#endif

    vector<double> ww;
/*    if (nClasses == 2) {
      //(solo para dos clases)
      int minerr = NMITotal - NMITrain;
      int w = 0;
      for(int i = 1; i < 100; i++) {
        int err = 0;
        for(int j = 0; j < NMITotal - NMITrain; j++) {
          if ( (*(dis[j]))[0]*i/100.0 > (*(dis[j]))[1]*(100.0-i)/100.0 ) {
            if (RealClass[j] == 1) err++;
          }
          else {
            if (RealClass[j] == 0) err++;
          }
        }
        if (err < minerr) {
          minerr = err;
          w = i;
        }
      }
      ww.push_back(w/100.0);
      ww.push_back((100.0-w)/100.0);
    }
    else {*/
      for(int i = 0; i < nClasses; i++) {
        ww.push_back(1.0 / Clasf->APrioriClas[i]);
      }
  /*  }*/
    DistributionWeights.push_back(ww);

/*    for(int i = 0; i < NMITotal - NMITrain; i++) {
      delete dis[i];
    }
    delete []dis;*/

  //  delete []RealClass;
//    delete []PredictedClass;

}
//---------------------------------------------------------------------------
double MIBagging::Error(int first, int last)
{
  int nclases = ((NomData*)data)->NumClass;

  vector<double> &weights = GetUsingWeights();

  data->MarkOrder();
//  data->SortByGroup(first, last);

  double tw = 0.0;
  double err = 0.0;

//ew  ndout = new NomData(((MINomData*)data)->CountMIInstances(first, last), 1 + nclases * ( nclases + 1) / 2);
  ndout = new NomData(((MINomData*)data)->CountMIInstances(first, last), 1 + nclases );
  int igrp = 0;

  rc.clear();
  pc.clear();

FILE *f, *fv, *fall;

string n = GenerateNewFileName("sifts_test_",".txt");
f=fopen(n.c_str(), "w");
n = GenerateNewFileName("votes_test_",".txt");
fv=fopen(n.c_str(), "w");
if (OutputLevel & 0x4) {
n = GenerateNewFileName("all_reached_nodes_test_",".txt");
fall=fopen(n.c_str(), "w");
}
else {
  fall=0;
}
  FILE *ff[30];
  int nqq = 10;
  if (OutputLevel & 0x8) {
    for (int i = 0; i < nqq; i++) {
      n = GenerateNewFileName("evidence_filtered_test_",".txt");
      ff[i] = fopen(n.c_str(), "w");
      if (ff[i] == 0) {
        cout << "Could not open " << n << endl;

      }
    }
  }


  for(int i = first; i <= last; i++) {
    vector<double> votes(nclases, 0.0);
    int ngrp = 0;
    do {
      Matriz mm(nqq, nclases);
vector<double> votes_2(nclases, 0.0);
vector<double> votes_v(nclases, 0.0);
      for(int j = 0; j < GetClassifiersToUse(); j++) {
        vector<double> d = Classifiers.at(j)->UnnormalizedDistribution(i);
if (fall) {
Node* n = ((DecisionTree*)Classifiers.at(j))->GetLastClassifyingNode();
fprintf(fall, "%d %d %ld ", i, j, (long int)n);
}
        if (OutputLevel & 0x8) {
          Matriz* n = (Matriz*) ((DecisionTree*)Classifiers.at(j))->GetLastClassifyingNode()->info;
          for(int k = 0; k < nclases; k++) {
            for (int iq = 0; iq < nqq; iq++) {
              mm[iq][k] = mm[iq][k] + (*n)[iq][k];
            }
          }
        }
        int cls = Classifiers.at(j)->Classify(i);
        //vector<double> d = Classifiers.at(j)->Distribution(i);
        bool Sustituye = false;
        if (MaxSingleInstanceMargin) {
          double m1 = Classifier::Margin(d, PositiveClass, false);
          if (m1 > Classifier::Margin(ndout->GetInstance(igrp), PositiveClass, false)) {
            Sustituye = true;
          }
        }
votes_v[cls] += 1.0;//weights[j];
        for(int k = 0; k < nclases; k++) {
if (fall) {
fprintf(fall, "%g ", d[k]);
}
votes_2[k] = votes_2[k] + d[k];// * weights[j];
          if (Sustituye) {
            votes[k] = d[k] * weights[j] * DistributionWeights[j][k];
            ndout->SetValueVar(igrp, k, d[k]);
          }
          else {
            votes[k] = votes[k] + d[k] * weights[j] * DistributionWeights[j][k];
            ndout->SetValueVar(igrp, k, ndout->GetValueVar(igrp, k) + d[k]);
          } 
        } 
if (fall) {
fprintf(fall, "\n");
}
      }
      ngrp++;
      if (OutputLevel & 0x8) {
        for (int k = 0; k < nqq; k++) {
          fprintf(ff[k], "%d %d %d ", data->GetDatIniPos(i), data->GetDatGroup(i), data->GetDatClass(i));
  //        if ((*q)[ic][k] <= ndout_all->GetValueVar(i, ic)) {
            for(int l = 0; l < nclases; l++) {
              fprintf(ff[k], "%g ", (double)mm[k][l]/GetClassifiersToUse());
            }
    //      }
      //    else {
        //    for(int l = 0; l < nclases; l++) {
          //    fprintf(ff[k], "%g ", 0.0);
            //}
          //}
          fprintf(ff[k], "\n ");
        }
      }
fprintf(f, "%d %d %d ", data->GetDatIniPos(i), data->GetDatGroup(i), data->GetDatClass(i));
for(int k = 0; k < nclases; k++) {
fprintf(f, "%g ", votes_2[k]/GetClassifiersToUse());
}
fprintf(f, "\n");

fprintf(fv, "%d %d %d ", data->GetDatIniPos(i), data->GetDatGroup(i), data->GetDatClass(i));
for(int k = 0; k < nclases; k++) {
fprintf(fv, "%d ", (int)votes_v[k]);
}
fprintf(fv, "\n");

      i++;
    } while ( i <= last && data->GetDatGroup(i) == data->GetDatGroup(i - 1) );

//ew    ndout->SetValueVar(igrp, nclases, GetClassifiersToUse());
//ew    ndout->SetValueVar(igrp, nclases*(nclases+1)/2, data->GetDatClass(i-1));
    ndout->SetValueVar(igrp, nclases+1, GetClassifiersToUse());//ew
    ndout->SetValueVar(igrp, nclases+2, ngrp);//MEAN
    ndout->SetValueVar(igrp, nclases, data->GetDatClass(i-1));//ew

    int cl = WhichClass(votes);

    if (cl != data->GetDatClass(i-1) ) {
      err += data->GetDatWeight(i-1);
    } 

    rc.push_back(data->GetDatClass(i-1));
    pc.push_back(cl);

    tw +=  data->GetDatWeight(i-1);
    igrp++;
    i--;
  }

fclose(f);
fclose(fv);
if (fall) {
fclose(fall);
}

  if (OutputLevel & 0x8) {
    for (int i = 0; i < nqq; i++) {
      fclose(ff[i]);
    }
  }

//ew  CalcDiffs();
//  n = GenerateProcNameFileName("test_",".cre");
//  CalcDiffs(ndout, 1);
//  ndout->SaveToFile((char*)n.c_str());

//  n = GenerateProcNameFileName("test_mean_",".cre");
//  CalcDiffs(ndout, 3);
//  ndout->SaveToFile((char*)n.c_str());
//  CalcDiffs(ndout, 2);
//  n = GenerateProcNameFileName("test_weighted_",".cre");
//  ndout->SaveToFile((char*)n.c_str());
  delete ndout;
  ndout = 0,

  data->ResetOrder();

  return err/tw;
}

//---------------------------------------------------------------------------
void MIBagging::SecuencialClassify(int ElementIndex, vector<int> *classes,
                                                     std::vector<int> *indclss)
{
  int nclases = ((NomData*)data)->NumClass;
  int j;
  int *winners;

  vector<double> votes(nclases, 0.0);
  int nClassifiersToUse = GetClassifiersToUse();
  classes->clear();
  if (indclss) indclss->clear();

  winners = new int[nclases+1];

  vector<double> &weights = GetUsingWeights();
  for(int i = 0; i < nClassifiersToUse; i++) {
    j = ElementIndex;
    do {
      vector<double> d = Classifiers.at(i)->UnnormalizedDistribution(j);
      //vector<double> d = Classifiers.at(i)->Distribution(j);
      for(int k = 0; k < nclases; k++) {
        votes[k] = votes[k] + d[k] * weights[i] * DistributionWeights[i][k];
      }
      j++;
    }
    while( j < data->GetNTotal() && data->GetDatGroup(j) == data->GetDatGroup(j-1) );

    int cl = WhichClass(votes, winners);

    //Si hay empate el clasificador clasifica como el anterior ensemble
/*    if (winners[0] > 1 && i > 0 ) classes->push_back((*classes)[i-1]);
    else  */                        classes->push_back(cl);

    if (indclss)        indclss->push_back(cl);
  }

  delete[] winners;

}
//---------------------------------------------------------------------------
//-------------------------------------------  BoostingWithBaggingInfo  -----
//---------------------------------------------------------------------------

BoostingWithBaggingInfo::BoostingWithBaggingInfo(double Umbral, Ensemble *ens)
:Boosting()
{
  this->Umbral = Umbral;
  this->ens = ens;
  percents = 0;
}

BoostingWithBaggingInfo::~BoostingWithBaggingInfo()
{
  if (percents) delete []percents;
//  if (ids) delete []ids;
}
//---------------------------------------------------------------------------
/*void BoostingWithBaggingInfo::GetPercent(int pos)
{
} */
//---------------------------------------------------------------------------
void BoostingWithBaggingInfo::StartArc()
{
  int NTrain = data->GetNTrain();
  int NTotal = data->GetNTotal();
//  int NIniPos  = data->GetNumVar()+ Data::IniPosIndex;
//  int iClass = data->GetNumVar()-1;
  int iVar = data->GetNumVar()+1;
  if (percents) delete []percents;
  percents = new double[NTotal];
  if (!ens) {
    ens = new GroupBagging();
    for(unsigned i=0;i<Classifiers.size();i++) {
      ens->AddClassifier(new GPTree(0, 1));
    }
    ens->Build(data);
  }

//  int ini = data->GetNTrain();
//  int fin = data->GetTotalData();
  ens->SetData(data);

  int clase;
  int datosguenos = 0;
//  bool ExcluirTrain = false;
//  FILE *f = fopen("percents", "at");
  for(int i=0;i<NTrain;i++) {
    int inipos = data->GetDatIniPos(i);
    percents[inipos] = ens->ClassificationCertainty(i, clase/*, ExcluirTrain*/);
    //Esto funciona para dos clases para mas hay que pensarlo
    if (clase!=data->GetDatClass(i)) percents[inipos] = -percents[inipos];
  //  fprintf(f, "%d\t%g\n", inipos, percents[inipos]);
    if (clase!=data->GetDatClass(i) && -percents[inipos]>Umbral) {
         data->SetValueVar(i, iVar, 1.0);
         excluidos.push_back(inipos);
    }
    else {
      data->SetValueVar(i, iVar, 0.0);
      datosguenos++;
    }
  }
//  fclose(f);
  data->SortOn(iVar, 0, NTrain-1);//0, 50, 300);//
//  data->SaveToFile("kra.cre",0,datosguenos);
//  delete []data;
//  data = Data::DataFromFile("kra.cre");
  data->SetNTrain(datosguenos);

  //Inicializo
  Boosting::StartArc();
}
//---------------------------------------------------------------------------

void BoostingWithBaggingInfo::AfterBuilding(int ic, bool &Continue, 
                                                           bool &ValidClassifer)
{
  Boosting::AfterBuilding(ic, Continue, ValidClassifer);
}
//---------------------------------------------------------------------------
void BoostingWithBaggingInfo::EndArc()
{
  data->SetNTrain(data->GetNTrain()+excluidos.size());
}
//---------------------------------------------------------------------------
string BoostingWithBaggingInfo::Info(int value)
{
  string res = Boosting::Info(value);
  char kk[90];

  sprintf(kk, "\n Excluidos(%lu): ", (long unsigned int)excluidos.size());
  res = res + kk;
  for (unsigned i=0;i<excluidos.size();i++) {
    sprintf(kk, "%d, ", excluidos[i]);
    res = res + kk;
  }

  return res;
}
//---------------------------------------------------------------------------
//---------------------------------------------------  BoostingGPTrees  -----
//---------------------------------------------------------------------------

BoostingGPTrees::BoostingGPTrees( void (*_salida)(char *),
                                    int _nClassifiers, int _NMin):Boosting()
{
  nMaxClassifiers = _nClassifiers;
  NMin =_NMin;
  salida = _salida;
  ForceContinuationByReweighting = true;
}
//---------------------------------------------------------------------------

void BoostingGPTrees::Init(Data *data)
{
  Boosting::Init(data);
  for(int i=0;i<nMaxClassifiers;i++)
    AddClassifier(new GPTree(salida, NMin));
}
//---------------------------------------------------------------------------
//------------------------------------------------  InvBoostingGPTrees  -----
//---------------------------------------------------------------------------
InvBoostingGPTrees::InvBoostingGPTrees( void (*_salida)(char *),
    int _nClassifiers, int _NMin):BoostingGPTrees(_salida, _nClassifiers, _NMin)
{
}
//---------------------------------------------------------------------------
void InvBoostingGPTrees::AfterBuilding(int ic, bool &Continue,
                                                          bool &ValidClassifer)
{
  double e, newWeight, k;
  int clase;
  int iClass = data->GetNumVar()-1;
  int iWeight = iClass + 5;
  Classifier *Clasf = Classifiers[ic];

  GPTree *tree = (GPTree*)Clasf;
  for (int iGrupo=0;iGrupo<2;iGrupo++) {
    tree->SetGroup(iGrupo);
    e = tree->Error(0, data->GroupCount[iGrupo]-1);
    if (e==0.0 || e>=0.4999) {
      Continue = false;
      return;
    }
    k = (1-e)/e;
    for (int i=0;i<data->GroupCount[iGrupo];i++) {
      clase=tree->Classify(i);
      if (clase!=(int)data->GetValueVar(i, iClass))
        newWeight = data->GetValueVar(i, iWeight)*k;  //Mal clasificado
      else
        newWeight = data->GetValueVar(i, iWeight)/k;  //Bien clasificado
      data->SetValueVar(i, iWeight, newWeight);
    }
  }
}
//---------------------------------------------------------------------------
//--------------------------------------------------  IndepClassifiers  -----
//---------------------------------------------------------------------------

IndepClassifiers::IndepClassifiers(double factor, bool UseDummyClass,
                                         double factor_incr) : ArcingTemplate()
{
  this->factor = factor;
  this->factor_incr = factor_incr;
  this->UseDummyClass = UseDummyClass;
  AsBreiman = false;

  GenerateOutputsAsKuncheva = false;
  outputs = 0;
}

IndepClassifiers::~IndepClassifiers()
{
  if (outputs)
    deleteMatrix(outputs, Classifiers.size());
}
//---------------------------------------------------------------------------
int IndepClassifiers::Classify(std::vector<double> distrb)
{
  vector<double> d=distrb;

  if (UseDummyClass) {
    d[dummy_class] = 0.0;
    //Renormalizamos
    double sum = 0.0;
    for (int i=0;i<dummy_class;i++) sum += d[i];
    for (int i=0;i<dummy_class;i++) d[i] = sum==0.0 ? 1./dummy_class : d[i]/sum;
  }
  
  return ArcingTemplate::Classify(d);
}
//---------------------------------------------------------------------------
void IndepClassifiers::SetData(Data *dat)
{
  if (UseDummyClass) ((NomData*)dat)->AddClass("~~~dummy~~~");
  ArcingTemplate::SetData(dat);
}
//---------------------------------------------------------------------------
void IndepClassifiers::Init(Data *dat)
{
  if (UseDummyClass) ((NomData*)dat)->AddClass("~~~dummy~~~");
  ArcingTemplate::Init(dat);
}
//---------------------------------------------------------------------------
int **dis;
Data *ddd;
void init()
{
  dis = new int*[360000];
  for(int i=0;i<360000;i++) {
    dis[i] = new int[2];
    dis[i][0] = dis[i][1] =0;
  }
  //ddd = Data::DataFromFile((char*)"kk2.cre");
}
void IndepClassifiers::StartArc()
{
init();
  e_ok    = 0.0;
  e_nook  = 0.0;
  e_train = 0.0;
 
  dummy_class = -1; 
}
//---------------------------------------------------------------------------
void IndepClassifiers::EndArc()
{
/*  int x2div=600;
  int x1div=600;
  unsigned char *fila = new unsigned char[600];
  unsigned char v     = 1;
  int AnchoEx = ((3 + (7  + x1div*v)/8)/4)*4;
  FILE *bmp=fopen("sal.bmp", "wb");

  for(int i=0;i<x2div;i++) {
    int bit, byte;
    byte = 0;
    bit = sizeof(*fila)*8;
    fila[byte]=0;
    for(int j=0;j<x1div;j++) {
      unsigned char clase = dis[x2div*i+j][0] > dis[x2div*i+j][1] ? 0 : 1;
      bit-=v;
      fila[byte] = fila[byte] | (clase<<bit);
      if (bit==0) {byte++; fila[byte]=0; bit=sizeof(*fila)*8;}
    }
    fwrite(fila, AnchoEx, 1, bmp);
  }
  fclose(bmp);*/
}
//---------------------------------------------------------------------------
/*int random_class(vector<double> pop)
{
  unsigned i;
  double pos, val;

  pos = (1.0 * rand ()) / (RAND_MAX + 1.0);

  val = 0.0;
  i = 0;
  while (val<pos && i<pop.size()) {
    val += pop[i];
    i++;
  }

  if (i>0) i--;
  return i;
}*/

Data* IndepClassifiers::Resample(int iClasf)
{

  int NClasses = ((NomData*)data)->NumClass;
  int NTrain = data->GetNTrain();
  int NVar = data->GetNumVar();

  Data *sample = data->Clone(0, NTrain-1);
  sample->SetNTrain(NTrain);

  if (AsBreiman) {
    vector< vector<double> > mat(NClasses);
    vector<double> pc(NClasses);

    double w = 0.0;
    for(int i=0;i<NClasses;i++) {
      pc[i] = APrioriClas[i]/NTrain;
      w += pc[i]*pc[i]; 
    }
    w = factor/(1-w);

    for(int i=0;i<NClasses;i++) {
      for(int j=0;j<NClasses;j++) {
        if (i==j) mat[i].push_back(1.0-w*(1.0-pc[i]));
        else mat[i].push_back(w*pc[j]);
      }
    }
    
    for(int i=0; i<NTrain;i++) {
      int ic = sample->GetDatClass(i);
      sample->SetValueVar(i, NVar-1, 0.5 + RandomIntProportional(mat[ic]));
    }
  }
  else if (GenerateOutputsAsKuncheva) {
    if (!outputs) {
      int iN = Classifiers.size()%2 == 0 ? 2 : 3;

      OutputGenerator o(iN, NTrain, 1.0 - factor, -1);
      outputs = getMatrix<int>(Classifiers.size(), NTrain); 
      o.generate(outputs);

      for (unsigned ii=iN;ii<Classifiers.size();ii+=2) {
        OutputGenerator o(2, NTrain, 1.0 - factor, -1);
        o.generate(outputs+ii);
      }
    }

    for(int i=0; i<NTrain;i++) {
      if ( outputs[iClasf][i] == 0 ) {
        int new_class = (int) ( (double) (NClasses-1) * rand() / (RAND_MAX+1.0) );
        if ( new_class >= sample->GetDatClass(i) ) new_class++;
        sample->SetValueVar(i, NVar-1, 0.5 + new_class); // Cambio clase
        sample->SetValueVar(i, NVar+1, 123);  // Indicador
      }
    }
    
  }
  else {
    int NFlips = (int) (factor*NTrain + 0.5);
 
    vector<int> aux_aleat;
    for(int i=0; i<NTrain;i++) {
      aux_aleat.push_back(i);
      sample->SetValueVar(i, NVar+1, 100);
    }

    int new_class = -1;
    if (UseDummyClass) {
      new_class = ((NomData*)sample)->AddClass("~~~dummy~~~");
      dummy_class = new_class;
    }
    for(int i=0; i<NFlips;i++) {
      int elem = (int)((double)(NTrain-i)*rand()/(RAND_MAX+1.0));
      int aux = aux_aleat[elem+i];
      aux_aleat[elem+i] = aux_aleat[i];
      aux_aleat[i] = aux;
      if (!UseDummyClass) {
        new_class = (int) ((double)(NClasses-1)*rand()/(RAND_MAX+1.0));
        if (new_class>=sample->GetDatClass(aux_aleat[i])) new_class++;
      }
      sample->SetValueVar(aux_aleat[i], NVar-1, 0.5 + new_class);
      sample->SetValueVar(aux_aleat[i], NVar+1, 123);
    }
  }

  return sample;  
}
//---------------------------------------------------------------------------
void IndepClassifiers::AfterBuilding(int ic, bool &Continue, 
                                                         bool &ValidClassifier)
{
  double e1, e2, e3;
  int NTrain = data->GetNTrain();
  double NClass = ((NomData*)data)->NumClass;
  int NVar = data->GetNumVar();
  int NFlips = (int) (factor*NTrain + 0.5);

  Data *bsp = Classifiers[ic]->GetData();
  bsp->ResetWeights();

  bsp->SortOn(NVar+1, 0, bsp->GetNTrain()-1);

  Classifiers[ic]->SetData(bsp);
  e1 = Classifiers[ic]->Error(0, bsp->GetNTrain()-1-NFlips);
  e2 = Classifiers[ic]->Error(bsp->GetNTrain()-NFlips, bsp->GetNTrain()-1);
  Classifiers[ic]->SetData(data);
  e3 = Classifiers[ic]->Error(0, NTrain-1);
  if (e1>=(NClass-1)/NClass /*|| e2>=0.5*/ || e3>=(NClass-1)/NClass) { 
//    ValidClassifier = false;
//    printf("\b*");
  }
  else {
 //   printf("%d\terr_ok=%g\terr_kk=%g\terr_dt=%g\n", NFlips, e1, e2, e3); 
  }
  e_ok    += e1;
  e_nook  += e2;
  e_train += e3;
  factor += factor_incr;


//Classifiers[ic]->SetData(ddd);
/*Classifiers[ic]->GetData()->OriginalOrder();
string n = GenerateNewFileName("e",".txt");
FILE *f = fopen(n.c_str(),"w");
//for(int i =0; i<360000;i++) {
for(int i =0; i<729;i++) {
for(int j =0; j<239;j++) {
int ii = i*239+j;
dis[ii][Classifiers[ic]->Classify(ii)]++;
fprintf(f, "%d\t", dis[ii][0]-dis[ii][1]);
}
fprintf(f, "\n");
}
fclose(f);
//Classifiers[ic]->SetData(data); 
if (ic>0) delete Classifiers[ic-1];
*/
}
//---------------------------------------------------------------------------
string IndepClassifiers::Info(int value)
{
  string res = ArcingTemplate::Info(value);
  char kk[1024];
  int NTrain = data->GetNTrain();
  int NFlips = (int) (factor*NTrain + 0.5);

  sprintf(kk, "\n Factor(Flips): %f(%d)", factor, NFlips);
  res = res + kk;
  sprintf(kk, " Dummy class: %s", UseDummyClass ? "YES" : "NO");
  res = res + kk;
  double c= Count();
  sprintf(kk, "\nEr(ok nook train): %f\t%f\t%f", e_ok/c, e_nook/c, e_train/c);
  res = res + kk;

  return res;
}
//---------------------------------------------------------------------------
//-----------------------------------------------------  PastingIVotes  -----
//---------------------------------------------------------------------------

PastingIVotes::PastingIVotes(int ResampleSize) : ArcingTemplate()
{
  this->ResampleSize = ResampleSize;
}

PastingIVotes::~PastingIVotes()
{
}
Data* PastingIVotes::Resample(int iClasf)
{

  /*if (iClasf==0) return data->GetBootstrapSample(ResampleSize);



  int NClasses = ((NomData*)data)->NumClass;
  int NTrain = data->GetNTrain();
  int NVar = data->GetNumVar();

    return data->GetBootstrapSample(ResampleSize);
  Data *sample = data->Clone(0, NTrain-1);
  sample->SetNTrain(NTrain);

    int NFlips = (int) (factor*NTrain + 0.5);
 
    vector<int> aux_aleat;
    for(int i=0; i<NTrain;i++) {
      aux_aleat.push_back(i);
      sample->SetValueVar(i, NVar+1, 100);
    }

    int new_class = -1;
    if (UseDummyClass) {
      new_class = ((NomData*)sample)->AddClass("~~~dummy~~~");
      dummy_class = new_class;
    }
    for(int i=0; i<NFlips;i++) {
      int elem = (int)((double)(NTrain-i)*rand()/(RAND_MAX+1.0));
      int aux = aux_aleat[elem+i];
      aux_aleat[elem+i] = aux_aleat[i];
      aux_aleat[i] = aux;
      if (!UseDummyClass) {
        new_class = (int) ((double)(NClasses-1)*rand()/(RAND_MAX+1.0));
        if (new_class>=sample->GetDatClass(aux_aleat[i])) new_class++;
      }
      sample->SetValueVar(aux_aleat[i], NVar-1, 0.5 + new_class);
      sample->SetValueVar(aux_aleat[i], NVar+1, 123);
    }

  return sample;  */
return 0;
}
//---------------------------------------------------------------------------
void PastingIVotes::AfterBuilding(int ic, bool &Continue, 
                                                         bool &ValidClassifier)
{
/*  double e1, e2, e3;
  int NTrain = data->GetNTrain();
  double NClass = ((NomData*)data)->NumClass;
  int NVar = data->GetNumVar();
  int NFlips = (int) (factor*NTrain + 0.5);

  Data *bsp = Classifiers[ic]->GetData();
  bsp->ResetWeights();

  bsp->SortOn(NVar+1, 0, bsp->GetNTrain()-1);

  Classifiers[ic]->SetData(bsp);
  e1 = Classifiers[ic]->Error(0, bsp->GetNTrain()-1-NFlips);
  e2 = Classifiers[ic]->Error(bsp->GetNTrain()-NFlips, bsp->GetNTrain()-1);
  Classifiers[ic]->SetData(data);
  e3 = Classifiers[ic]->Error(0, NTrain-1);*/
  return;
}
//---------------------------------------------------------------------------
//---------------------------------------------  CompetentModelBagging  -----
//---------------------------------------------------------------------------
double CompetentModelBagging::Error(int first, int last)
{
/*  double err=0.0;
  int iClass     = data->GetNumVar()-1;
  int iWeight    = iClass + 5;
  double TotalWeight=0.0;

  int e=0;
  int eC=0;
  int eW=0;
  int emconf[2][2];
  emconf[0][0] = emconf[1][0] = emconf[0][1] = emconf[1][1] = 0;
  
  for (int i=first;i<=last;i++) {
    int RealClass = (int)data->GetValueVar(i, iClass);
    TotalWeight += data->GetValueVar(i, iWeight);
    if (Classify(i)!=RealClass)
      err += data->GetValueVar(i, iWeight);
    if (clasevotes!=RealClass) e++;
    if (clasevotesC!=RealClass) eC++;
    if (clasevotesW!=RealClass) eW++;
    emconf[0][0] += mconf[0][0];
    emconf[1][0] += mconf[1][0];
    emconf[0][1] += mconf[0][1];
    emconf[1][1] += mconf[1][1];
  }

  FILE *f=fopen("err.txt", "at");
  fprintf(f, "%d\t%d\t%d\t%d\t", last-first+1, e, eW, eC);
  fprintf(f, "%d\t%d\t%d\t%d\n", emconf[0][0], emconf[1][0], emconf[0][1], 
                                                                  emconf[1][1]);
  fclose(f);
      */
  return 0.0;//err/TotalWeight;
}
void CompetentModelBagging::GetVotes(int ElementIndex, vector<int> &WithRefs,
                                                            vector<int> &NoRefs)
{
/*  int nclases = ((NomData*)data)->NumClass;
  int clase;
  WithRefs.clear();
  NoRefs.clear();
  for (int j=0;j<nclases;j++) {
    WithRefs.push_back(0);
    NoRefs.push_back(0);
  }
  int nClassifiersToUse = GetClassifiersToUse();
  for(int i=0;i<nClassifiersToUse;i++) {
    clase = Classifiers.at(i)->Classify(ElementIndex);
    NoRefs[clase]++;
    //Aqui entra el tema del modelo
    Referees.at(i)->SetData(GetData());
    if  (Referees.at(i)->Classify(ElementIndex)==0)  WithRefs[clase]++;
    else {
      int kl = clase==0 ? 1 : 0;//2-class chapuza
      WithRefs[kl]++;
    }
  }*/
}

int CompetentModelBagging::Classify(int ElementIndex)
{
/*  int nclases = ((NomData*)data)->NumClass;
  int clase;
  double *votes = new double[nclases];
  double *votesC = new double[nclases];
  double *votesW = new double[nclases];
  for (int j=0;j<nclases;j++) votesW[j]=votesC[j]=votes[j]=0.0;
  mconf[0][0] = mconf[1][0] = mconf[0][1] = mconf[1][1] = 0;
  int nClassifiersToUse = GetClassifiersToUse();
  for(int i=0;i<nClassifiersToUse;i++) {
    clase = Classifiers.at(i)->Classify(ElementIndex);
    votes[clase]++;
    //Aqui entra el tema del modelo
    Referees.at(i)->SetData(GetData());
    if  (!UseReferees || Referees.at(i)->Classify(ElementIndex)==0) {
      votesC[clase]++;
      if (clase==data->GetDatClass(ElementIndex))
         mconf[0][0]++;
      else mconf[1][0]++;
    }
    else {
//      int kl = clase==0 ? 1 : 0;//2-class chapuza
//      votesC[kl]++;
      if (clase==data->GetDatClass(ElementIndex))
         mconf[0][1]++;
      else mconf[1][1]++;
    }
    //
    votesW[clase] += ClassifierWeights.at(i);
  }
  ClassWithoutWeigths = WhichClass(votes, nclases);
  ClassWithWeigths    = WhichClass(votesW, nclases);

  //Si hay n clases con el mismo numero de votos e igual nUmero de
  //ejemplos se quitan los votos de los Ultimos clasificadores hasta
  //que sOlo haya un ganador
  for(int i=nClassifiersToUse-1;i>0 && ClassWithoutWeigths<0; i--) {
    clase = Classifiers.at(i)->Classify(ElementIndex);
    votes[clase]--;
    ClassWithoutWeigths = WhichClass(votes, nclases);
  }
  for(int i=nClassifiersToUse-1;i>0 && ClassWithWeigths<0; i--) {
    clase = Classifiers.at(i)->Classify(ElementIndex);
    votesW[clase] -= ClassifierWeights.at(i);
    ClassWithWeigths = WhichClass(votesW, nclases);
  }

  clasevotes   = ClassWithoutWeigths;
  clasevotesC  = WhichClass(votesC, nclases);
  clasevotesW  = ClassWithWeigths;
  
  if (clasevotesC<0)  clasevotesC  = clasevotes;
  
  delete []votes;
  delete []votesC;
  delete []votesW;
  
  clase = (UseWeights) ? ClassWithWeigths : ClassWithoutWeigths;
*/  return -1;//clase;
}
//---------------------------------------------------------------------------
#ifdef COMPBAG
void CompetentModelBagging::CalcWeights(int idx)
{
/*  Data *dt = Bagging::Resample(idx);*/

#else
void CompetentModelBagging::CalcWeights(Classifier *Clasf, bool &Continue)
{
/*  Boosting::CalcWeights(Clasf, Continue);
  int idx=Referees.size();
  int NTrain=data->GetNTrain();
  int NCols=data->GetNumVar();  */
#endif  
/*  double ee = -1;
  int nNoUsados = 0;
  int ccls = data->GetNumVar()-1;
  int nBien, nMal;
  
  nBien=nMal=0;

  for(int i=0; i<NTrain;i++) {*/
#ifdef COMPBAG
/*    if (residx[i] == 1) continue;
    for(int j=0; j<NCols;j++) {
      data->SetValueVar(nNoUsados, j, m[i][j]);
    }
    int clscls = Classifiers[idx]->Classify(nNoUsados);*/
#else
/*    int clscls = Clasf->Classify(nNoUsados);*/
#endif
  /*  int clsreal = data->GetDatClass(nNoUsados);
    double err = clscls==clsreal ? 0.5*/ /*Bien*//* : 1.5*/ /*Mal*//*;
    if (clscls==clsreal) nBien++; else nMal++;
    data->SetValueVar(nNoUsados, ccls, err);
    nNoUsados++;
  }
  //
  CARTTree *cart = new CARTTree();
  data->SetNTrain(nNoUsados);
  //Ponderate
  double wBien, wMal;
  double w[2000];
  wBien = (double)nNoUsados/(2.0*nBien);
  wMal = (double)nNoUsados/(2.0*nMal);
  for(int i=0; i<nNoUsados;i++) {
    w[i]=data->GetDatWeight(i);
    if (data->GetValueVar(i, ccls)<1.0) data->SetDatWeight(i, wBien);
    else data->SetDatWeight(i, wMal);
  }
  cart->Build(data);
  for(int i=0; i<nNoUsados;i++) data->SetDatWeight(i, w[i]);
  Referees.push_back(cart);
  data->SetNTrain(NTrain);*/
  
}
//---------------------------------------------------------------------------

//--------------------------------------------------  Main for testing  -----
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#ifdef _DEBUG_ARCING_CPP
#include "C45Tree.h"
int main(int argc, char* argv[])
{
  NomData *dat = new NomData((char*)"train.cre");
  MINomData *data = (MINomData*)NomData::GenerateMIDataSet(dat, 20);
  delete dat;

  NomData *dat2 = new NomData((char*)"test.cre");
  MINomData *data2 = (MINomData*)NomData::GenerateMIDataSet(dat2, 20);
  delete dat2;

  //MIBoosting *b = new MIBoosting();
  MIBagging *b = new MIBagging();
  //Bagging *b = new Bagging();

  cout << data->NumClass << endl;

  int ndatos = data->GetNTotal();
  data->SetNTrain(ndatos);
  for(int i=0;i<3;i++) {
    b->AddClassifier(new RandomForestTree(0, 1, 4, false));
    //  b->AddClassifier(new CARTTree(0, 4, false));
      //b->AddClassifier(new C45Tree(true));
  }

  ProgresoConsola pc("tree");
  data->ResetWeights();
  data->SortOn(data->GetNumVar() + Data::IniPosIndex, 0, data->GetNTrain()-1);
  b->Build(data, &pc);
  data->ResetWeights();

  cout << b->Info(2);

  vector<double> errores;
  b->SecuencialError(0, ndatos-1, &errores);

  for (unsigned i=0; i<errores.size(); i++) {
    cout << errores[i] << " ";
  }
  cout << endl;
  
  cout << b->Error(0, ndatos-1) << endl;

  b->SetData(data2);
  b->SecuencialError(0, data2->GetNTotal()-1, &errores);

  for (unsigned i=0; i<errores.size(); i++) {
    cout << errores[i] << " ";
  }
  cout << endl;
  
  cout << b->Error(0, data2->GetNTotal()-1) << endl;

  delete data;
  delete data2;
  delete b;
}
#endif

