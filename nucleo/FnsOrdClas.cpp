//---------------------------------------------------------------------------


#include "FnsOrdClas.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Matriz.h"
#include "Ensemble.h"
//#include "UCramerIterativo.h"
//#include <Dialogs.hpp>
#include <math.h>
#include <nrr.h>
#include "FnsClasificacion.h"
extern "C" {
  #include <genetic.h>
#ifdef COMP_SDP
  #include <declarations.h> // SDP library
#endif
}
#include "FuncOptim.h"
#include "Utils.h"

#include <limits>
#include <iostream>
#include <sstream>
#include <time.h>

typedef struct {
	double diversity;
	int order;
} StructToSort;

int cmpToSort(const void * p1,  const void *p2) {

	StructToSort *pp1 = (StructToSort *) p1;
	StructToSort *pp2 = (StructToSort *) p2;

	if (pp1->diversity > pp2->diversity)
		return 1;
	else 
		return 0;
}


void DoOrdenaPorRERegularizadoPorDiversidad(Ensemble *ens, NomData *data, double rho);

//---------------------------------------------------------------------------
using namespace std;
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void DarDeAltaFnsOrdClas()
{
  OrdenOriginal::DarDeAlta();
  ClasificadoresNecesarios::DarDeAlta();
  OrdenaPorContribucion::DarDeAlta();
  OrdenaPorReduccionError::DarDeAlta();
  OrdenaPorReduccionDeDistancias::DarDeAlta();
  OrdenaPorAngulos::DarDeAlta();
  OrdenaPorEPIC::DarDeAlta();
  InvertirOrden::DarDeAlta();
  OrdenaPorReduccionErrorPonderado::DarDeAlta();
  PseudoBoosting::DarDeAlta();
  OrdenaPorReduccionErrorDoble::DarDeAlta();
  SeleccionGenetica::DarDeAlta();
#ifdef COMP_SDP
  SeleccionSDP::DarDeAlta();
#endif
  AnalisisDeOrden::DarDeAlta();
  OrdenaEnPruebas::DarDeAlta();
  SeleccionCombinatoria::DarDeAlta();
  Mejores::DarDeAlta();
  OrdenaPorGGreedy::DarDeAlta();
  OrdenaPorRERegularizadoPorDiversidad::DarDeAlta(); 
  OrdenaPorUWA::DarDeAlta(); 
}
//---------------------------------------------------------------------
class BestSolutions
{
    int **Ccd; 
    int nDatos;
    int nClasf;
    int nClases;
    Matriz *res;
    int *ClasReal;
    int start_at;
    vector< vector<int> > all_sols;

  public:
    BestSolutions(Ensemble *ens, Data *dat, int start_at=0) {
      this->start_at = start_at;
      nDatos = dat->GetNTotal();
      nClasf = ens->GetClassifiersToUse();
      nClases = ((NomData*)dat)->NumClass;
      if (nClasf-start_at>32) 
        throw MandatoError("Tu de que vas, tiempo estimado mayor de 12 dias (ademas se me acaba el int)");

      //Obtiene la clsificacion para cada clasificador del conjunto
      //y para cada dato
      Ccd = ens->MatClasif0(dat);

      //Obtiene la clase de cada dato
      ClasReal = new int[nDatos];
      for(int j=0;j<nDatos;j++) {
        ClasReal[j] = (int)dat->GetDatClass(j);
      }

      all_sols.resize(nClasf);

      res = new Matriz(nClasf+1, nClasf+2);
    }

    ~BestSolutions() {
      delete res;
      delete []ClasReal;
      for(int i=0;i<nClasf;i++) {
        delete []Ccd[i];
        delete []Ccd[i+nClasf];
      }
      delete []Ccd;
    }

    Matriz *GetBests() {return res;}
    vector<vector<int> > GetAllBests() {return all_sols;}

    void run() {
      int **acum = new int*[nDatos];
      for(int i=0;i<nDatos;i++) {
        acum[i] = new int[nClases];
        for(int j=0;j<nClases;j++) {
          acum[i][j] = 0;
        }
      }

      for(int j=0;j<start_at;j++) {
        for(int i=0;i<nDatos;i++) {
          acum[i][Ccd[j][i]]++;
        }
      }
      Ccd = Ccd + start_at;

      unsigned comb = 0;//(0xff>>(32-start_at))<<32-start_at;
      //comb <<= nClasf-start_at;
      MinComb(nClasf-1-start_at, acum, comb, start_at);

      Ccd = Ccd - start_at;

      for(int i=0;i<nClasf+1;i++) {
        unsigned comb = (int) (*res)[i][0];
        unsigned ibit = 1;
        for(int j=0;j<nClasf;j++) {
          if (j<start_at) (*res)[i][j] = 1;
          else {
            (*res)[i][j] = ibit & comb ? 1: 0;
            ibit <<= 1;
          }
        }
      }

      for(int i=0;i<nDatos;i++)
        delete []acum[i];
      delete []acum;
    }

    void MinComb(int ibit, int **acum, unsigned comb, int ups) {
      if (ibit==-1) {
        if (ups%2==0) return;
        Error(acum, comb, ups);
        static int k = 0;
        if ((k%10000)==0) cout << hex << comb << endl;
        k++;
        return;
      }
      MinComb(ibit-1, acum, comb, ups);

      int **acum2 = new int*[nDatos];
      for(int i=0;i<nDatos;i++) {
        acum2[i] = new int[nClases];
        for(int j=0;j<nClases;j++) {
          acum2[i][j] = acum[i][j];
        }
        acum2[i][Ccd[ibit][i]]++;
      }

      unsigned bit = 1 << ibit;
      MinComb(ibit-1, acum2, comb | bit, ups+1);

      for(int i=0;i<nDatos;i++)
        delete []acum2[i];
      delete []acum2;
    }

    void Error(int **acum, unsigned comb, int ups) {
      int oks=0;
      double kk[1024];

      for(int i=0;i<nDatos;i++) {
        for(int j=0;j<nClases;j++) {
           kk[j] = acum[i][j];
        }
        int cls = Classifier::WhichClass(kk, nClases);
        if (cls==ClasReal[i]) oks++;
      }
      if (oks>(*res)[ups][nClasf]) {
        (*res)[ups][nClasf+1] = 1;
        (*res)[ups][nClasf] = oks;
        (*res)[ups][0] = comb;
        all_sols[ups-1].clear();
        all_sols[ups-1].push_back(comb);
        //cout << dec << "\t"  << hex << comb << dec << "\t" << oks  << "\t" << ups << endl;
      }
      else if (oks==(*res)[ups][nClasf]) {
        (*res)[ups][nClasf+1] += 1.0;
        all_sols[ups-1].push_back(comb);
        //cout << dec << "\t"  << hex << comb << dec << "\t" << oks  << "\t" << ups << endl;
      }
/*
int pp=0;
for(int jjj=0;jjj<nClasf;jjj++){
int kk = comb & (1 << jjj) ? 1 : 0;
printf("%d", kk);
pp+=kk;
}
printf(" %d %d\n", pp, oks);
*/
    }  

};
//---------------------------------------------------------------------
// Orders adding at each step the most different classifier
void DoDiversityOrdering(Ensemble *ens, NomData *data)
{
  int **Ccd;  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
  int *clase_temp;
  double **distribs_temp;
  int num_clases;
  int nDatos;
  int nClasf;

  int ini = 0;
  int fin = data->GetNTotal()-1;
  num_clases = data->NumClass;
  nDatos = fin-ini+1;
  nClasf = ens->GetClassifiersToUse();

  if (nClasf==1) return;

  //Obtiene la clsificacion para cada clasificador del conjunto
  //y para cada dato
  Ccd = ens->MatClasif0(data);
  //Se busca la pareja de clasificadores mAs distintos
  int i_min=-1, j_min=-1;
  double kappa_min = numeric_limits<double>::infinity();
  for(int i=0;i<nClasf;i++) {
    for(int j=i+1;j<nClasf;j++) {
      double kp =
               Classifier::KappaStatistic(Ccd[i], Ccd[j], nDatos, num_clases);
      if (kp<kappa_min) {
        kappa_min = kp;
        i_min = i;
        j_min = j;
      }
    }
  }

/*  double i_error = nDatos;
  for(int i=0;i<nClasf;i++) {
    int error = 0;
    Classifier *c = ens->GetClassifier(i);
    for(int j=0;j<nDatos;j++) {
      if (c->Classify(ini+j)!=data->GetDatClass(ini+j)) error++;
    }

    if(i_error > error) {
	i_min = i;
	i_error = error;
    }
  }

  clase_temp = new int[nDatos];
  distribs_temp = new int*[nDatos];
  for(int j=0;j<nDatos;j++) {
    distribs_temp[j] = new int[num_clases];
    for(int i=0;i<num_clases;i++) distribs_temp[j][i] = 0;
    distribs_temp[j][Ccd[i_min][j]]++;
    clase_temp[j] = Ccd[i_min][j];
  }

  ens->Exchange(0, i_min);
  int *cc = Ccd[0];
  Ccd[0] = Ccd[i_min];
  Ccd[i_min] = cc; */

  //Distribucion temporal de cada ejemplo
  clase_temp = new int[nDatos];
  distribs_temp = new double*[nDatos];
  for(int j=0;j<nDatos;j++) {
    distribs_temp[j] = new double[num_clases];
    for(int i=0;i<num_clases;i++) distribs_temp[j][i] = 0;
    distribs_temp[j][Ccd[i_min][j]] += ens->GetWeightInUse(i_min);
    distribs_temp[j][Ccd[j_min][j]] += ens->GetWeightInUse(j_min);
    clase_temp[j] = ens->GetWeightInUse(i_min) <= ens->GetWeightInUse(j_min) ? Ccd[i_min][j] : Ccd[j_min][j];
  }

  ens->Exchange(0, i_min);
  ens->Exchange(1, j_min);
  int *cc = Ccd[0];
  Ccd[0] = Ccd[i_min];
  Ccd[i_min] = cc;
  cc = Ccd[1];
  Ccd[1] = Ccd[j_min];
  Ccd[j_min] = cc;

  //Se anyade el clasificador menos semejante a la combinacion
  for(int i=2;i<nClasf;i++) {
    kappa_min = numeric_limits<double>::infinity();
    for(int j=i+1;j<nClasf;j++) {
      double kp = 
            Classifier::KappaStatistic(clase_temp, Ccd[j], nDatos, num_clases);
      if (kp<kappa_min) {
        kappa_min = kp;
        j_min = j;
      }
    }

    for(int j=0;j<nDatos;j++) {
      distribs_temp[j][Ccd[j_min][j]] += ens->GetWeightInUse(j_min);
      if (Ccd[j_min][j]!=clase_temp[j] && 
            distribs_temp[j][Ccd[j_min][j]]>distribs_temp[j][clase_temp[j]]) {
        clase_temp[j] = Ccd[j_min][j];
      }
    }
    ens->Exchange(i, j_min);
    int *cc = Ccd[i];
    Ccd[i] = Ccd[j_min];
    Ccd[j_min] = cc;
  }


  for(int i=0;i<nClasf;i++) {
    delete []Ccd[i];
    delete []Ccd[i+nClasf];
  }
  delete []Ccd;
  for(int i=0;i<nDatos;i++) {
    delete []distribs_temp[i];
  }
  delete []distribs_temp;
  delete []clase_temp;
}
//---------------------------------------------------------------------
int DoPseudoBoosting(Ensemble *ens, int ini, int fin, bool Resampling)
{
  int **Ccd;  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
  double *Wd, *Wd_sample, *Wd_aux;
  int *ClasReal; //Clase de cada dato
  int nClases;
  int nDatos;
  int nClasf;
  int FirstStop=-1;

  //Recupera los datos con que se construyo el clasificador
//  ens->ResetOriginalTrainingData();

  NomData *data = (NomData*)ens->GetData();
  nClases = data->NumClass;
  nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  nClasf = ens->GetClassifiersToUse();

  if (nClasf==1) return 1;

  //Obtiene la clase de cada dato
  ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++)
    ClasReal[j] = data->GetDatClass(ini+j);

  //Obtiene la clsificacion para cada clasificador del conjunto y para cada dato
  //Tambien obtiene los dos clasificadores con menor error
  Ccd = new int*[nClasf];
  for(int i=0;i<nClasf;i++) {
    Classifier *c = ens->GetClassifier(i);
    Ccd[i] = new int[nDatos];
    for(int j=0;j<nDatos;j++) {
      Ccd[i][j] = c->Classify(ini+j)==ClasReal[j] ? 1 : -1;
    }
  }

  //Pesos de datos
  double sumWd = nDatos;
  Wd = new double[nDatos];
  Wd_sample = new double[nDatos];
  for(int i=0;i<nDatos;i++) {
    Wd[i] = 1.0;
  }

  //Se busca el clasificador que tenga menor error usando los pesos Wd para
  //los ejemplos
  double *Dctemp = new double[nClases];
  double MargenTemp, ErrorTemp;
  double MargenMax, error;
  int iMargenMax=-1;
  FILE *f=fopen("pseudoboostinfo.txt", "at");
  for(int iClasf=0;iClasf<nClasf;iClasf++) {
    if (Resampling) {
      for(int i=0;i<nDatos;i++) {
        Wd_sample[i] = 0.0;
      }
      for(int j=0;j<nDatos;j++) {
        int ii = RandomIntProportional(Wd, nDatos, sumWd);
        Wd_sample[ii]++;
      }
    }

    Wd_aux = Resampling ? Wd_sample : Wd;
    
    MargenMax = -sumWd;
    error = 0.5;
    for(int i=iClasf;i<nClasf;i++) {
      MargenTemp = 0.0;
      ErrorTemp = 0.0;
      for(int j=0;j<nDatos;j++) {
        MargenTemp += Wd_aux[j]*Ccd[i][j];
        if (Ccd[i][j]<0) ErrorTemp += Wd_aux[j];
      }
      if (MargenTemp>MargenMax) {
        MargenMax = MargenTemp;
        error = ErrorTemp;
        iMargenMax = i;
      }
    }
    
    //Si el clasificador seleccionado tiene un error >0.5 reinicializamos
    //los pesos y volvemos a seleccionar un clasificador
//    double error = (-1.0/(2.0*nDatos))*MargenMax+0.5;
    error /= sumWd;
    if (error>0.5) { 
      if (FirstStop<0) FirstStop = iClasf;
      for(int j=0;j<nDatos;j++) Wd[j] = 1.0;  
      fprintf(f, "\t%d", iClasf);
      continue;
    } else if (error==0.0) {
      error = 0.5/nDatos;
    }
    
    //Se da por bueno el clasf seleccionado
    ens->SetClassifierWeight(iMargenMax, log((1.0-error)/error)); 
    ens->Exchange(iClasf, iMargenMax);
    int *cc = Ccd[iClasf];
    Ccd[iClasf] = Ccd[iMargenMax];
    Ccd[iMargenMax] = cc;
    double KMal = 1.0/(2.0*error);
    double KBien = 1.0/(2.0*(1.0-error));
    //Se reajustan los pesos
    sumWd = 0.0;
    for(int j=0;j<nDatos;j++) {
      Wd[j] *=  Ccd[iClasf][j] > 0 ? KBien : KMal;
      sumWd += Wd[j];
    }
  }
  fprintf(f, "\n");
  fclose(f);

  delete []Wd;
  for(int i=0;i<nClasf;i++) delete Ccd[i];
  delete []Ccd;
  delete ClasReal;
  delete Dctemp;

  return FirstStop;
}
//---------------------------------------------------------------------
//---------------------------------------------------------------------------
void DoOrdenaPorError(Ensemble *ens, int ini, int fin, bool ExcluirTrain)
{
  NomData *data = (NomData*)ens->GetData();
  int nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  int nClasf = ens->GetClassifiersToUse();

  for(int i=0;i<nClasf;i++) {
    int ndat  = 0;
    int error = 0;
    Classifier *c = ens->GetClassifier(i);
    for(int j=0;j<nDatos;j++) {
      if (ExcluirTrain &&
                     c->UsedInOriginalTrainingData(data->GetDatIniPos(ini+j))) {
        continue;
      }
      ndat++;
      if (c->Classify(ini+j)!=data->GetDatClass(ini+j)) error++;
    }
    ens->SetClassifierWeight(i, (double)ndat/error);
  }

  ens->OrdenarClasificadoresPorPeso();

}
//---------------------------------------------------------------------------
void DoOrdenaPorAlgunAcum(Ensemble *ens, int ini, int fin,
                double mastrain,double menostrain,double masval,double menosval)
{
  NomData *data = (NomData*)ens->GetData();
  int nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  int nClasf = ens->GetClassifiersToUse();

  for(int i=0;i<nClasf;i++) {
    double acum  = 0.0;
    Classifier *c = ens->GetClassifier(i);
    for(int j=0;j<nDatos;j++) {
      bool clasificadoOK = (c->Classify(ini+j)==data->GetDatClass(ini+j));
      if (c->UsedInOriginalTrainingData(data->GetDatIniPos(ini+j))) {
        if (clasificadoOK) acum += mastrain;
        else               acum += menostrain;
      }
      else {
        if (clasificadoOK) acum += masval;
        else               acum += menosval;
      }
    }
    ens->SetClassifierWeight(i, acum);
  }

  ens->OrdenarClasificadoresPorPeso();

}
//---------------------------------------------------------------------------
void DoOrdenaPorContribucion(Ensemble *ens, int ini, int fin)
{
  int **Mcd;  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
  double **Ddc;  //Matriz con el nUmero de votos para cada datos y para cada clase
  int *Nd;  //Numero de clasificadores que tienen el dato d en el conjunto de validacion
  int/*double*/ ErrAux;
  int *ClasReal; //Clase de cada dato
  int nClases;
  int nDatos;
  int nClasf;

  //Recupera los datos con que se construyo el clasificador
//  ens->ResetOriginalTrainingData();

  NomData *data = (NomData*)ens->GetData();
  nClases = data->NumClass;
  nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  nClasf = ens->GetClassifiersToUse();

  //Obtiene la clase de cada dato
  ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++)
    ClasReal[j] = data->GetDatClass(ini+j);

  //Obtiene la clsificacion para cada clasificador del conjunto y para cada dato
  //Tambien obtiene los dos clasificadores con menor error
  bool ExcluirTrain = false;
  Ddc = new double*[nDatos];
  Mcd = new int*[nClasf];
  Nd = new int[nDatos];
  for(int j=0;j<nDatos;j++) {
    Ddc[j] = new double[nClases];
    Nd[j] = 0;
    for(int i=0;i<nClases;i++) {
      Ddc[j][i] = 0.0;
    }
  }
  for(int i=0;i<nClasf;i++) {
    Classifier *c = ens->GetClassifier(i);
    Mcd[i] = new int[nDatos];
    for(int j=0;j<nDatos;j++) {
      if (ExcluirTrain &&
                     c->UsedInOriginalTrainingData(data->GetDatIniPos(ini+j))) {
        Mcd[i][j] = -1;
        continue;
      }
      Mcd[i][j] = c->Classify(ini+j);
      Ddc[j][Mcd[i][j]]++;
      Nd[j]++;
    }
  }

  for(int i=0;i<nClasf;i++) {
    ErrAux = 0;
    for(int j=0;j<nDatos;j++) {
// DONE AS IN:
// http://www.cs.washington.edu/homes/grossman/projects/573projects/learning/
//  
      //if (Mcd[i][j] == -1) continue;
      //if (Mcd[i][j]==ClasReal[j]) 
    //    ErrAux += 1.0 + (((double)Nd[j]-Ddc[j][ClasReal[j]])/Nd[j]);
  //    else
//        ErrAux += 1.0 - ((double)Ddc[j][ClasReal[j]]/Nd[j]);
// FORMA MIA
      if (Mcd[i][j]!=ClasReal[j]) continue;
      ErrAux += (int)((double)Nd[j]/Ddc[j][ClasReal[j]]);
    }
    ens->SetClassifierWeight(i, ErrAux);
  }
  ens->OrdenarClasificadoresPorPeso();

  vector<int> oo = ens->DameOrdenOriginal();
  FILE *f=fopen("kkorden.txt", "at");
  for (unsigned i=0;i<oo.size();i++) {
    fprintf(f, "\t%d", oo[i]);
  }
  fprintf(f, "\n");
  fclose(f);

  for(int i=0;i<nDatos;i++) delete []Ddc[i];
  delete []Ddc;
  for(int i=0;i<nClasf;i++) delete []Mcd[i];
  delete []Mcd;
  delete []Nd;
}
//---------------------------------------------------------------------
int DoOrdenaPorReduccionError(Ensemble *ens, int ini, int fin, int start_at, 
                                                                  bool inverse)
{
  int **Mcd;  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
  double **Ddc;  //Matriz con el nUmero de votos para cada datos y para cada clase
  int ErrMin1, ErrMin2;
  int iErrMin1=-1;//, iErrMin2=-1;
  int ErrAux;
  int ErrAct;
  int *ErrorClas; //Error de cada clasificador
  int *ClasReal; //Clase de cada dato
  int *ClasAct;  //Clase del conjunto para el nUmero de clasificadores actuales
  int nClases;
  int nDatos;
  int nClasf;
  int nmin=0, ErrActMin;

  //Recupera los datos con que se construyo el clasificador
//  ens->ResetOriginalTrainingData();

  NomData *data = (NomData*)ens->GetData();
  nClases = data->NumClass;
  nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  nClasf = ens->GetClassifiersToUse();

  if (nClasf==1) return 1;

  //Obtiene la clase de cada dato
  ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++)
    ClasReal[j] = data->GetDatClass(ini+j);

  //Obtiene la clsificacion para cada clasificador del conjunto y para cada dato
  //Tambien obtiene los dos clasificadores con menor error
  ErrorClas = new int[nClasf];
  Mcd = new int*[nClasf];
  ErrMin1 = ErrMin2 = nDatos;
  int incremento = inverse ? -1 : 1;
  for(int i=0;i<nClasf;i++) {
    Classifier *c = ens->GetClassifier(i);
    Mcd[i] = new int[nDatos];
    ErrAux = 0;
    for(int j=0;j<nDatos;j++) {
      Mcd[i][j] = c->Classify(ini+j);
      if (Mcd[i][j]!=ClasReal[j]) ErrAux+=incremento;
    }
    ErrorClas[i] = nDatos - ErrAux;
//    ens->SetClassifierWeight(i, 
//            log(((double)nDatos-incremento*ErrAux))/(incremento*ErrAux));

/*    if (ErrAux<ErrMin1) {
      ErrMin2  = ErrMin1;
      iErrMin2 = iErrMin1;
      ErrMin1  = ErrAux;
      iErrMin1 = i;
    }
    else if (ErrAux<ErrMin2) {
      ErrMin2  = ErrAux;
      iErrMin2 = i;
    }*/
  }

  int *cc;
/*  if (start_at<=0) {*/
    //Coloca los dos mejores clasificadores en primera posicion
/*    ens->Exchange(0, iErrMin1);
    cc = Mcd[0];
    Mcd[0] = Mcd[iErrMin1];
    Mcd[iErrMin1] = cc;
    ens->Exchange(1, iErrMin2);
    cc = Mcd[1];
    Mcd[1] = Mcd[iErrMin2];
    Mcd[iErrMin2] = cc; *///quita quita que lo descomentes

/*    ErrMin1 = ErrorClas[iErrMin1];
    ErrorClas[iErrMin1] = ErrorClas[0];
    ErrorClas[0] = ErrMin1;
    ErrMin2 = ErrorClas[iErrMin2];
    ErrorClas[iErrMin2] = ErrorClas[1];
    ErrorClas[1] = ErrMin2;

    start_at = 2;
  }*/start_at = 0;

  //Distribuciones de clase
  Ddc = new double*[nDatos];
  for(int i=0;i<nDatos;i++) {
    Ddc[i] = new double[nClases];
    for(int j=0;j<nClases;j++) {
      Ddc[i][j] = 0.0;
    }
  }
  ClasAct = new int[nDatos];
  ErrAct = 0;
  for(int j=0;j<nDatos;j++) {
    for(int k=0;k<start_at;k++) {
      Ddc[j][Mcd[k][j]]++;
    }
    ClasAct[j] = start_at == 0 ? ClasReal[j]+1 : ens->WhichClass(Ddc[j], nClases);
//    if (ClasAct[j]<0) printf("JARLLLLLLLLR\n\n"), ClasAct[j] = Mcd[0][j];
    if (ClasAct[j]!=ClasReal[j]) ErrAct++;
  }

  //Busca el clasificador que reduzca mas el error(o lo aumente menos)
  double *Dctemp = new double[nClases];
  int ErrTempSube, ErrTempBaja, ErrTemp;
  ErrActMin = ErrAct;
  nmin = 1;
  for(int iClasf=start_at;iClasf<nClasf;iClasf++) {
    ErrMin1 = nDatos;
    for(int i=iClasf;i<nClasf;i++) {
      ErrTempSube = 0;
      ErrTempBaja = 0;
      for(int j=0;j<nDatos;j++) {
        if (Mcd[i][j]==ClasAct[j] || ens->GetWeightInUse(i) + Ddc[j][Mcd[i][j]] < Ddc[j][ClasAct[j]])
          continue; //Este ejemplo ni sube ni baja el error (Optimizamos)
        for(int k=0;k<nClases;k++) {
          Dctemp[k] = Ddc[j][k];
        }
        Dctemp[Mcd[i][j]] += ens->GetWeightInUse(i);
        int ClasTemp = ens->WhichClass(Dctemp, nClases);//, 0, &(ens->APrioriClas));
        if (ClasTemp>=0 && ClasTemp!=ClasAct[j]) {
          if (ClasAct[j]==ClasReal[j]) ErrTempSube++;    //Sube el error
          else if (ClasTemp==ClasReal[j]) ErrTempBaja--; //Baja el error
        }
      }
      ErrTemp = incremento*(ErrTempBaja + ErrTempSube);
      if (ErrTemp<ErrMin1 /*|| (ErrTemp==ErrMin1 &&
                                       ErrorClas[i] > ErrorClas[iErrMin1])*/) {
        ErrMin1 = ErrTemp;
    /*    Baja = ErrTempBaja;
        Sube = ErrTempSube;*/
        iErrMin1 = i;
      }
    }
//FILE *fff=fopen("flujo.csv", "a");
//printf("%d\t%d\n", -Baja, Sube);
//fclose(fff);
    ens->Exchange(iClasf, iErrMin1);
    cc = Mcd[iClasf];
    Mcd[iClasf] = Mcd[iErrMin1];
    Mcd[iErrMin1] = cc;
    ErrAct += ErrMin1;
    ErrMin1 = ErrorClas[iErrMin1];
    ErrorClas[iErrMin1] = ErrorClas[iClasf];
    ErrorClas[iClasf] = ErrMin1;
    if (ErrAct<=ErrActMin) {
//cout << nmin << (ErrAct<ErrActMin?"(A)":"(B)") << "..";
      ErrActMin = ErrAct;
      nmin = iClasf+1;
    }
    for(int j=0;j<nDatos;j++) {
      Ddc[j][Mcd[iClasf][j]]++;
      int ct = ens->WhichClass(Ddc[j], nClases);//, 0, &(ens->APrioriClas));
      if (ct>=0) ClasAct[j] = ct;
    }
  }

  //Calculo los diagramas de kappa-error
/*  int e[1000];
  for(int i=0;i<nClasf;i++) {
    e[i] = 0;
    for(int j=0;j<nDatos;j++) {
      if (Mcd[i][j]!=ClasReal[j]) e[i]++;
    }
  }
  double k;
  FILE *fke = fopen("k-error.txt", "wt");
  for(int i=nClasf-1;i>=0;i--) {
    for(int j=i-1;j>=0;j--) {
      k = Classifier::KappaStatistic(Mcd[i], Mcd[j], nDatos, nClases);
      fprintf(fke, "%g\t", k);
      fprintf(fke, "%g\t", (double)e[i]/nDatos);
      fprintf(fke, "%g\t", (double)e[j]/nDatos);
      fprintf(fke, "%g\t", ((double)e[i]+e[j])/(2.0*nDatos));
      fprintf(fke, "%d\t", i);
      fprintf(fke, "%d\n", j); 
    }
  }
  fclose(fke);
*/
/*  vector<int> oo = ens->DameOrdenOriginal();
  FILE *f=fopen("kkorden.txt", "at");
  for (unsigned i=0;i<oo.size();i++) {
    fprintf(f, "\t%d", oo[i]);
  }
  fprintf(f, "\n");
  fclose(f);
*/
  for(int i=0;i<nDatos;i++) delete Ddc[i];
  delete []Ddc;
  for(int i=0;i<nClasf;i++) delete Mcd[i];
  delete []Mcd;
  delete ErrorClas;
  delete ClasReal;
  delete ClasAct;
  delete Dctemp;

//cout << nmin << endl;
  return nmin;
}
//---------------------------------------------------------------------
int DoOrdenaPorUWA(Ensemble *ens, NomData *data)
{
  int **Mcd;  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
  double **Ddc;  //Matriz con el nUmero de votos para cada datos y para cada clase
  int Min1;
  int iMin1=-1;
  int ErrAux;
  int *ErrorClas; //Error de cada clasificador
  int *ClasReal; //Clase de cada dato
  int *ClasAct;  //Clase del conjunto para el nUmero de clasificadores actuales
  int nClases;
  int nDatos;
  int nClasf;

  nClasf = ens->GetClassifiersToUse();

  if (nClasf==1) return 1;

  NomData *dataold = (NomData*)ens->GetData();
  ens->SetData(data);

  nClases = data->NumClass;
  nDatos = data->GetNTotal();

  //Obtiene la clase de cada dato
  ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++)
    ClasReal[j] = data->GetDatClass(j);

  //Obtiene la clsificacion para cada clasificador del conjunto y dato
  //Tambien obtiene el clasificador con menor error
  ErrorClas = new int[nClasf];
  Mcd = new int*[nClasf];
  Min1 = nDatos;
  for(int i = 0; i < nClasf; i++) {
    Classifier *c = ens->GetClassifier(i);
    Mcd[i] = new int[nDatos];
    ErrAux = 0;
    for(int j = 0; j < nDatos; j++) {
      Mcd[i][j] = c->Classify(j);
      if (Mcd[i][j] != ClasReal[j]) ErrAux++;
    }
    ErrorClas[i] = ErrAux;
    if (ErrorClas[i] < Min1) {
      Min1 = ErrorClas[i];
      iMin1 = i;
    }
  }

  int *cc;

  //Ponemos el clasificador ocn menor error en la primera posicion
  ens->Exchange(0, iMin1);
  cc = Mcd[0];
  Mcd[0] = Mcd[iMin1];
  Mcd[iMin1] = cc;
  ErrAux = ErrorClas[iMin1];
  ErrorClas[iMin1] = ErrorClas[0];
  ErrorClas[0] = ErrAux;

  //Distribuciones de clase por ejemplo para el subensemble t
  Ddc = new double*[nDatos];
  for(int i=0;i<nDatos;i++) {
    Ddc[i] = new double[nClases];
    for(int j=0;j<nClases;j++) {
      Ddc[i][j] = 0.0;
    }
  }
  ClasAct = new int[nDatos]; //Clase predicha de cada dato por el subensemble t
  for(int j=0;j<nDatos;j++) {
    Ddc[j][Mcd[0][j]]++;
    ClasAct[j] = Mcd[0][j];
  }

  //Busca el clasificador que con mejor UWA
  double UWA;
  for(int iClasf = 1; iClasf < nClasf; iClasf++) {
    Min1 = -1;
    for(int i = iClasf; i < nClasf; i++) {
      UWA = 0;
      for(int j = 0; j < nDatos; j++) { //Calculamos el UWA para cada clasificador
        if (Mcd[i][j]==ClasReal[j] && ClasAct[j]==ClasReal[j]) {
          UWA += 1.0 - Ddc[j][ClasReal[j]]/iClasf;
        }
        else if (Mcd[i][j]!=ClasReal[j] && ClasAct[j]==ClasReal[j]) {
          UWA -= 1.0 - Ddc[j][ClasReal[j]]/iClasf;
        }
        else if (Mcd[i][j]==ClasReal[j] && ClasAct[j]!=ClasReal[j]) {
          UWA += Ddc[j][ClasReal[j]]/iClasf;
        }
        else if (Mcd[i][j]!=ClasReal[j] && ClasAct[j]!=ClasReal[j]) {
          UWA -= Ddc[j][ClasReal[j]]/iClasf;
        }
        else  {
          exit(1); //We shouldnt be here
        }
      }
      if ( Min1 < 0 || UWA > Min1 ) { 
        Min1 = UWA;            //De momento el mejor
        iMin1 = i;
      }
    }

    //Intercambiamos clasificadores y datos
    ens->Exchange(iClasf, iMin1);
    cc = Mcd[iClasf];
    Mcd[iClasf] = Mcd[iMin1];
    Mcd[iMin1] = cc;
    ErrAux = ErrorClas[iMin1];
    ErrorClas[iMin1] = ErrorClas[iClasf];
    ErrorClas[iClasf] = ErrAux;

    //Actualizamos distribuciones
    for(int j = 0; j < nDatos; j++) {
      Ddc[j][Mcd[iClasf][j]]++;
      int ct = ens->WhichClass(Ddc[j], nClases);
      if (ct>=0) ClasAct[j] = ct;
    }
  }

  for(int i = 0; i < nDatos; i++) delete []Ddc[i];
  delete []Ddc;
  for(int i = 0; i < nClasf; i++) delete []Mcd[i];
  delete []Mcd;
  delete []ErrorClas;
  delete []ClasReal;
  delete []ClasAct;

  ens->SetData(dataold);

  return 0;
}
//---------------------------------------------------------------------
int BuscaClasificadorReduceError(int **Mcd, int nClasf, int nDatos, double **Ddc, int nClases, 
                                 int iClasf, int &ErrMin1, Ensemble *ens, int *ClasReal)
{
    int iErrMin = -1;
    int ErrTemp;
    int ClasAct;
    double *Dctemp = new double[nClases];
 
    ErrMin1 = nDatos;
    for(int i=iClasf;i<nClasf;i++) {
      int ErrTempSube = 0;
      int ErrTempBaja = 0;
      for(int j=0;j<nDatos;j++) {
        ClasAct = ens->WhichClass(Ddc[j], nClases);
        if (Mcd[i][j]==ClasAct || Ddc[j][Mcd[i][j]]+1<Ddc[j][ClasAct])
          continue; //Este ejemplo ni sube ni baja el error
        for(int k=0;k<nClases;k++) {
          Dctemp[k] = Ddc[j][k];
        }
        Dctemp[Mcd[i][j]]++;
        int ClasTemp = ens->WhichClass(Dctemp, nClases);
        if (ClasTemp>=0 && ClasTemp!=ClasAct) {
          if (ClasAct==ClasReal[j]) ErrTempSube++;    //Sube el error
          else if (ClasTemp==ClasReal[j]) ErrTempBaja--; //Baja el error
        }
      }
      ErrTemp = (ErrTempBaja + ErrTempSube);
      if (ErrTemp<ErrMin1 || (ErrTemp==ErrMin1 &&
            ens->GetClassifierWeight(i) > ens->GetClassifierWeight(iErrMin))) {
        ErrMin1 = ErrTemp;
        iErrMin = i;
      }
    }

  
  return iErrMin;
}
Matriz* DoOrdenaPorReduccionErrorBackfitting(Ensemble *ens, Data *dtr, Data *dts)
{
  int **Mcd;  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
  double **Ddc;  //Matriz con el nUmero de votos para cada datos y para cada clase
  int ErrMin1, ErrMin2;
  int iErrMin1=-1, iErrMin2=-1;
  int ErrAux;
  int ErrAct;
  int *ClasReal; //Clase de cada dato
  int *ClasAct;  //Clase del conjunto para el nUmero de clasificadores actuales
  int nClases;
  int nDatos;
  int nClasf;
  int nmin, ErrActMin;

  //Recupera los datos con que se construyo el clasificador
//  ens->ResetOriginalTrainingData();

  ens->SetData(dtr);
  NomData *data = (NomData*)dtr;
  nClases = data->NumClass;
  nDatos = data->GetNTotal();
  nClasf = ens->GetClassifiersToUse();

  if (nClasf==1) return new Matriz(2, 1);

  //Obtiene la clase de cada dato
  ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++)
    ClasReal[j] = data->GetDatClass(j);

  //Obtiene la clsificacion para cada clasificador del conjunto y para cada dato
  //Tambien obtiene los dos clasificadores con menor error
  Mcd = new int*[nClasf];
  ErrMin1 = ErrMin2 = nDatos;
  for(int i=0;i<nClasf;i++) {
    Classifier *c = ens->GetClassifier(i);
    Mcd[i] = new int[nDatos];
    ErrAux = 0;
    for(int j=0;j<nDatos;j++) {
      Mcd[i][j] = c->Classify(j);
      if (Mcd[i][j]!=ClasReal[j]) ErrAux++;
    }
//    ens->SetClassifierWeight(i, 
//            log(((double)nDatos-ErrAux)/ErrAux));
    if (ErrAux<ErrMin1) {
      ErrMin2  = ErrMin1;
      iErrMin2 = iErrMin1;
      ErrMin1  = ErrAux;
      iErrMin1 = i;
    }
    else if (ErrAux<ErrMin2) {
      ErrMin2  = ErrAux;
      iErrMin2 = i;
    }
  }

  int *cc;
  //Coloca los dos mejores clasificadores en primera posicion
  ens->Exchange(0, iErrMin1);
  cc = Mcd[0];
  Mcd[0] = Mcd[iErrMin1];
  Mcd[iErrMin1] = cc;
  ens->Exchange(1, iErrMin2);
  cc = Mcd[1];
  Mcd[1] = Mcd[iErrMin2];
  Mcd[iErrMin2] = cc;//quita quita que lo descomentes

  int start_at = 2;

  //Distribuciones de clase
  Ddc = new double*[nDatos];
  for(int i=0;i<nDatos;i++) {
    Ddc[i] = new double[nClases];
    for(int j=0;j<nClases;j++) {
      Ddc[i][j] = 0.0;
    }
  }
  ClasAct = new int[nDatos];
  ErrAct = 0;
  for(int j=0;j<nDatos;j++) {
    for(int k=0;k<start_at;k++)
      Ddc[j][Mcd[k][j]]++;
//    Ddc[j][Mcd[1][j]]++;
    ClasAct[j] = ens->WhichClass(Ddc[j], nClases);
    if (ClasAct[j]<0) printf("JARLLLLLLLLR\n\n"), ClasAct[j] = Mcd[0][j];
    if (ClasAct[j]!=ClasReal[j]) ErrAct++;
  }

  //Busca el clasificador que reduzca mas el error(o lo aumente menos)
  double *Dctemp = new double[nClases];
  ErrActMin = ErrAct;
  nmin = 1;
  Matriz *errors = new Matriz(2, nClasf);
ens->SetClassifiersToUse(1);
(*errors)[0][0] = ens->Error(0, dtr->GetNTotal()-1);
ens->SetData(dts);
(*errors)[1][0] = ens->Error(0, dts->GetNTotal()-1);
ens->SetData(dtr);
ens->SetClassifiersToUse(nClasf);

  for(int iClasf=start_at;iClasf<nClasf;iClasf++) {
//if (nClasf==200 && iClasf==41) break;
//if (nClasf==100 && iClasf==21) break;

//Calculamos el error en train y test
ens->SetClassifiersToUse(iClasf);
(*errors)[0][iClasf-1] = ens->Error(0, dtr->GetNTotal()-1);
ens->SetData(dts);
(*errors)[1][iClasf-1] = ens->Error(0, dts->GetNTotal()-1);
ens->SetData(dtr);
ens->SetClassifiersToUse(nClasf);

    iErrMin1 = BuscaClasificadorReduceError(Mcd, nClasf, nDatos, Ddc, nClases,
                                 iClasf, ErrMin1, ens, ClasReal);
    ens->Exchange(iClasf, iErrMin1);
    cc = Mcd[iClasf];
    Mcd[iClasf] = Mcd[iErrMin1];
    Mcd[iErrMin1] = cc;
    ErrAct += ErrMin1;
    if (ErrAct<=ErrActMin) {
      ErrActMin = ErrAct;
      nmin = iClasf+1;
    }
    for(int j=0;j<nDatos;j++) {
      Ddc[j][Mcd[iClasf][j]]++;
      int ct = ens->WhichClass(Ddc[j], nClases);
      if (ct>=0) ClasAct[j] = ct;
    }
   //Backfitting 
    bool changed = true;
    for(int k=0;k<100 && changed ;k++) {
      changed = false;
//printf("C");
      for(int ii=0;ii<=iClasf;ii++) {
        //Quito el clasf ii
        int err = 0;
        for(int j=0;j<nDatos;j++) {
          Ddc[j][Mcd[ii][j]]--;
          int ct = ens->WhichClass(Ddc[j], nClases);
          if (ct!=ClasReal[j]) err++;
        }
        //Miro si hay alguno que mejore eñ error
        iErrMin1 = BuscaClasificadorReduceError(Mcd, nClasf, nDatos, Ddc, nClases,
                                 iClasf+1, ErrMin1, ens, ClasReal);
        if (err+ErrMin1 < ErrAct) {
//printf("x");
          changed = true;
          ens->Exchange(ii, iErrMin1);
          cc = Mcd[ii];
          Mcd[ii] = Mcd[iErrMin1];
          Mcd[iErrMin1] = cc;
          ErrAct = err + ErrMin1;
          if (ErrAct<=ErrActMin) {
            ErrActMin = ErrAct;
            nmin = iClasf+1;
          }
          for(int j=0;j<nDatos;j++) {
            Ddc[j][Mcd[ii][j]]++;
            int ct = ens->WhichClass(Ddc[j], nClases);
            if (ct>=0) ClasAct[j] = ct;
          }
        }
        else {
          for(int j=0;j<nDatos;j++) {
            Ddc[j][Mcd[ii][j]]++;
          }
        }
      }
    }
//vector<int> oo = ens->DameOrdenOriginal();
//for (unsigned ii=0;ii<oo.size();ii++) { 
//printf("%3d ", oo[ii]); 
//}
//printf("\n");
  }

(*errors)[0][nClasf-1] = ens->Error(0, dtr->GetNTotal()-1);
ens->SetData(dts);
(*errors)[1][nClasf-1] = ens->Error(0, dts->GetNTotal()-1);
ens->SetData(dtr);

  for(int i=0;i<nDatos;i++) delete Ddc[i];
  delete []Ddc;
  for(int i=0;i<nClasf;i++) delete Mcd[i];
  delete []Mcd;
  delete ClasReal;
  delete ClasAct;
  delete Dctemp;

  cout << nmin << endl;

  return errors;
}
//---------------------------------------------------------------------
/*int DoOrdenaPorReduccionFuncion(Ensemble *ens, Data *data, Funcion f)
{
  
}*/
//---------------------------------------------------------------------
//---------------------------------------------------------------------
double ErrorVal(int *P1, int nDatos)
{
  int pos=0, neg=0;
  for (int i=0;i<nDatos;i++)
    if (P1>0) pos++;
    else if (P1<0) neg++;
  return (double)neg/(neg+pos);
}
void DoOrdenaPorReduccionDeVal(Ensemble *ens, int ini, int fin)
{
  int **Icd;  //Matriz[iclasificador][idato] con 1 si bien clas. y -1 si mal.
  int *ClasReal; //Clase de cada dato
  int *SumAct;
  int *SumTmp;
  double ErrorValMin, ErrorValTemp;
  int nDatos;
  int nClasf;
  vector<int> ErrorValMinIndices;

  NomData *data = (NomData*)ens->GetData();
  nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  nClasf = ens->GetClassifiersToUse();

  //Obtiene la clase de cada dato
  ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++)
    ClasReal[j] = data->GetDatClass(ini+j);

  bool ExcluirTrain = true;
  Icd = new int*[nClasf];
  for(int i=0;i<nClasf;i++) {
    Classifier *c = ens->GetClassifier(i);
    Icd[i] = new int[nDatos];
    for(int j=0;j<nDatos;j++) {
      bool EnTrain = ExcluirTrain &&
               c->UsedInOriginalTrainingData(data->GetDatIniPos(ini+j));
      Icd[i][j] = EnTrain ? 0 : (c->Classify(ini+j)!=ClasReal[j]) ? -1 : 1;
    }
  }

  SumAct = new int[nDatos];
  SumTmp = new int[nDatos];
  for (int i=0;i<nDatos;i++) {
    SumAct[i] = 0;
    SumTmp[i] = 0;
  }

  int deende = 1;
  //Busca el clasificador que se acerque mas al punto PtoRef
  for(int iClasf=0;iClasf<nClasf;iClasf++) {
    ErrorValMin = 1.0;
/*    for(int i=iClasf;i<nClasf;i++) {
      //Calcula distancia a PtoRef
      for(int j=0;j<nDatos;j++) {
        PtoTmp[j] = PtoAct[j] + Icd[i][j];
      }
      DistTemp = Distancia(PtoTmp, PtoRef, nDatos);
      if (DistTemp<DistMin) {
        //De momento es el que mas reduce la distancia
        DistMin = DistTemp;
        iDistMin = i;
      }
    ens->Exchange(iClasf, iDistMin);
    int *cc = Icd[iClasf];
    Icd[iClasf] = Icd[iDistMin];
    Icd[iDistMin] = cc;
    for(int j=0;j<nDatos;j++) {
      PtoAct[j] += Icd[iClasf][j];
   //   PtoRef[j] = dist + sqrt(iClasf+1);
    }
      */
    //Combinaciones de nClasf-iClasf elementos tomados de 1 en 1
    Combinaciones k(nClasf-iClasf, deende<nClasf-iClasf ? deende : nClasf-iClasf);
    while(k.HayMasCombinaciones()) {
      vector<int> comb = k.SiguienteCombinacion();

      for(int j=0;j<nDatos;j++)
        SumTmp[j] = SumAct[j];

      for (unsigned i=0;i<comb.size();i++) {
        //Calcula distancia a PtoRef
        for(int j=0;j<nDatos;j++) {
          SumTmp[j] += Icd[iClasf+comb[i]][j];
        }
      }

      ErrorValTemp = ErrorVal(SumTmp, nDatos);
      if (ErrorValTemp<ErrorValMin) {
        //De momento es el que mas reduce el erro validacion estimado
        ErrorValMin = ErrorValTemp;
        ErrorValMinIndices = comb;
      }
    }
    int iii = iClasf;
    for (unsigned i=0;i<ErrorValMinIndices.size();i++, iii++) {
      ens->Exchange(iii, ErrorValMinIndices[i]+iClasf);
      int *cc = Icd[iii];
      Icd[iii] = Icd[ErrorValMinIndices[i]+iClasf];
      Icd[ErrorValMinIndices[i]+iClasf] = cc;
      for(int j=0;j<nDatos;j++) {
        SumAct[j] += Icd[iii][j];
      }
    }
    iClasf += ErrorValMinIndices.size()-1;
  }

  for(int i=0;i<nClasf;i++) delete Icd[i];
  delete []Icd;
  delete []ClasReal;
  delete []SumTmp;
  delete []SumAct;
}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
double Modulo(double *P, int d)
{
  double mod = 0;

  for (int i=0;i<d;i++)
    mod += P[i]*P[i];

  return sqrt(mod);
}
double ProductoEscalar(double *P1, double *P2, int d)
{
  double pe = 0;

  for (int i=0;i<d;i++)
    pe += P2[i]*P1[i];

  return pe;
}
double Angulo(double *P1, double *P2, int d)
{
  return acos(ProductoEscalar(P1, P2, d)/(Modulo(P1,d)*Modulo(P2,d)));
}
double Distancia(double *P1, double *P2, int d)
{
  double Dist = 0;

  for (int i=0;i<d;i++)
    Dist += pow(P2[i]-P1[i], 2);

  return sqrt(Dist);
}
// Devuelve la cte K tal que los vectores u=(v + K * n) y n sean 
//   perperdiculares. Cte de proyeccion de v sobre el plano
//   definido por n.
double CteDeProyeccion(double *n, double *v, int d)
{
  return -ProductoEscalar(n, v, d)/ProductoEscalar(n, n, d);
}
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
void DoNumClasfCV(Ensemble *ens, int ini, int fin, int ncv, 
		int &iminmax, int &iminmin)
{
  NomData *data = (NomData*)ens->GetData();
  vector<double> err;
  vector<int> maxidxs;
  double minmin, minmax;
  int imin=-1, imax;
  int vuertas = 50;

  ncv=2;
  iminmin = iminmax = 0;
 for(int k=0;k<vuertas;k++) {
  data->ResetCV(ncv);
  data->SetNGrow(data->GetNTotal());
  for(int j=0;j<ncv;j++) {
    data->SetCV(j+1);
    ens->OrdenarClasificadoresPorOrdenOriginal();
//    DoOrdenaPorReduccionError(ens, 0, data->GetNumData()-1);
//    DoOrdenaPorReduccionErrorPonderado(ens, 0, data->GetNumData()-1);
    DoOrdenaPorReduccionDeDistancias(ens, 50, 0, data->GetNumData()-1);
    ens->SecuencialError(data->GetNumData(), fin, &err);
    minmin=1.0;
    minmax=1.0;
    for (unsigned i=0; i<err.size(); i++) {
//      fprintf(fff, "\t%g", err[i]);
      if (minmax>=err[i]) {
        minmax = err[i];
        imax = i;
      }
      if (minmin>err[i]) {
        minmin = err[i];
        imin = i;
      }
    }
//    fprintf(fff, "\nestim %d\t%d\n", imin, imax);
    maxidxs.push_back(imax);
    iminmin += imin;
    iminmax += imax;
  }
 }
  FILE *fff=fopen("errCV", "at");
  for (unsigned i=0; i<maxidxs.size(); i++) 
     fprintf(fff, "\t%d", maxidxs[i]);
  fprintf(fff, "\n");
  fclose(fff);
  iminmin = (int)(0.5 + iminmin/(ncv*vuertas));
  iminmax = (int)(0.5 + iminmax/(ncv*vuertas));
}
void DoOrdenaCV(Ensemble *ens, int ini, int fin, int ncv)
{
/*  NomData *data = (NomData*)ens->GetData();
//  int ncv = data->GetNumCV();
  vector<double> err;
  double min;
  int pos, posmin=0;

  data->ResetCV(ncv);
  data->SetNGrow(data->GetNTrain());
  for(int i=0;i<ncv;i++) {
    data->SetCV(i+1);
    ens->OrdenarClasificadoresPorOrdenOriginal();
    DoOrdenaPorReduccionError(ens, 0, data->GetNumData()-1);
//    DoOrdenaPorReduccionDeDistancias(ens, 10, 0, data->GetNumData()-1);
    ens->SecuencialError(data->GetNumData(), fin, &err);
    min=1.0;
    for (unsigned i=0; i<err.size(); i++) {
      if (min>=err[i]) {
        min = err[i];
        pos = i;
      }
    }
    posmin += pos;
  }*/
  int iminmax, iminmin;
  DoNumClasfCV(ens, ini, fin, ncv, iminmax, iminmin);
  ens->OrdenarClasificadoresPorOrdenOriginal();
  DoOrdenaPorReduccionError(ens, ini, fin);
//  DoOrdenaPorReduccionDeDistancias(ens, 10, 0, fin);
  ens->SetClassifiersToUse(iminmin);

}
void DoOrdenaPorReduccionDeDistanciasCV(Ensemble *ens, int dist,
                                                      int ini, int fin, int nCV)
{

}
int TestDistancias(Ensemble *ens, Data *data, bool SetNegativeToZero, bool UseOOB)
{
  int nClasf            = ens->GetClassifiersToUse();
  int nDatos            = data->GetNTotal();
  double **Icd          = CreaMatrizDistancias(ens, data, UseOOB);
  double *W;
  vector<double> pesos(nClasf);
  vector<Classifier*> clsfs(nClasf);
  
  W = new double[nClasf];

  double *PtoAcum = new double[nDatos];
//  double *PtoMrg  = new double[nDatos];
  double *PtoRef  = new double[nDatos];
  double *PtoDiag  = new double[nDatos];
  double *errs  = new double[nClasf];
  for(int j=0;j<nClasf;j++) {
    pesos[j] = ens->GetClassifierWeight(j);
    clsfs[j] = ens->GetClassifier(j);
    errs[j] = 0.0; 
  }
  for (int i=0;i<nDatos;i++) {
    PtoAcum[i] = 0.0;
    PtoRef[i] = 1.0;
    PtoDiag[i] = 1.0;
//    PtoMrg[i] = nClasf*Classifier::Margin(ens->Distribution(i), data->GetDatClass(i)); 
    for(int j=0;j<nClasf;j++) {
      PtoAcum[i] += ens->GetWeightInUse(j) * Icd[j][i];
      errs[j] += Icd[j][i];
    }
  }
  
/*double *kk = PtoMrg;//This changes the PtoAcum to the margin point of the ensemble
PtoMrg = PtoAcum;
PtoAcum = kk;*/

  double lambda = 0.0;//, denom = 0.0;
/*  for (int i=0;i<nDatos;i++) {
    lambda += PtoAcum[i]*PtoRef[i];
    denom  += PtoAcum[i]*PtoAcum[i];
  }
  lambda /= denom;*/
  lambda = CteDeProyeccion(PtoAcum, PtoRef, nDatos);

/*  for (int i=0;i<20;i++) {
    cout << PtoAcum[i] << ", ";
  }
  cout << "..." << endl;
  for (int i=0;i<20;i++) {
    cout << PtoMrg[i] << ", ";
  }
  cout << "..." << endl;*/

  int kllk = 0;
  for (int i=0;i<nDatos;i++) {
    PtoRef[i] = PtoRef[i] + lambda*PtoAcum[i];
    if (PtoRef[i]<0.0) {
      kllk++;
      if (SetNegativeToZero)
        PtoRef[i] = 0.0; 
    }
  //  if (i<20)  cout << PtoRef[i] << ", ";
  }
/*  cout << endl << "neg =                             ";
  cout << (SetNegativeToZero ? "(zeroed)" : "        ");
  cout << "               NEGATIVES: " << kllk << endl;
  cout << "PE = " <<  ProductoEscalar(PtoAcum, PtoRef, nDatos) << endl;*/
 
  int ok = 0; 
  int ok_45 = 0; 
  double mb = Modulo(PtoRef, nDatos);
  //double md = Modulo(PtoDiag, nDatos);
  for(int j=0;j<nClasf;j++) {
    double pe = ProductoEscalar(Icd[j], PtoRef, nDatos);
    double ma = Modulo(Icd[j], nDatos);
//    double a = pe/(ma*mb);
    double a = pe/(ma);
    W[j] = a;
//////////////////    ens->SetClassifierWeight(j, a);
  //  cout.precision(3);
//    cout << "{"<<(180.0*acos(a/mb)/3.14159265358979323844);
    //if (a>0.0) {
    if (pe>0.0) {
      ok++;
//      cout << "*";
    }
    if (a/mb>.08715574274765817357) {
      ok_45++;
//      cout << "*";
    }
    //cout << "}";
//Proyection VVVVVVVVVVVVVVVVVVVV
/*    double lambda_1 = CteDeProyeccion(PtoAcum, Icd[j], nDatos);
    double z = -Modulo(PtoAcum, nDatos)*lambda_1;
    for (int i=0;i<nDatos;i++)
      Icd[j][i] = Icd[j][i] + lambda_1*PtoAcum[i];
    double lambda_2 = CteDeProyeccion(PtoRef, Icd[j], nDatos);
    double x = -mb*lambda_2;
    for (int i=0;i<nDatos;i++)
      Icd[j][i] = Icd[j][i] + lambda_2*PtoRef[i];

    double lambda_3 = CteDeProyeccion(Icd[0], Icd[j], nDatos);
//    cout << Angulo(Icd[0], Icd[j], nDatos) << endl;
    double y = -Modulo(Icd[0], nDatos)*lambda_3;
//    double y = Modulo(Icd[j], nDatos) * (Icd[j][1]>=0.0 ? 1.0 : -1.0);
    cout << x << "\t" << y << "\t" << z << endl;
*/
  }
  //cout << endl << "NO. 85º: " << ok_45 << endl;
  char nf[256];
  sprintf(nf, "stop85_%d.txt", nClasf);
  FILE *f = fopen(nf, "a");
  fprintf(f, "%d\t%d\n", ok_45, kllk);
  fclose(f);

  ens->OrdenarClasificadoresPorPeso(W);

  if (ens->GetUseWeights()) {// Parte cuadratica (que se puede quitar o hacer mejor)
    for(int j=0;j<nClasf;j++) {//que recupera los pesos originales de
      for(int k=0;k<nClasf;k++) {//los clasificadores
        if (clsfs[k]==ens->GetClassifier(j))
          ens->SetClassifierWeight(j, pesos[k]);
      }
    }
  }
  else {
    sprintf(nf, "angles_diag_%d.txt", nClasf);
    f = fopen(nf, "a");
    for(int j=0;j<nClasf;j++) {
      fprintf(f, "%g\t", acos(ens->GetClassifierWeight(j)/mb));
      ens->SetClassifierWeight(j, 1.0);
    }
    fprintf(f, "\n");
    fclose(f);
  }

  delete []PtoRef;
  delete []PtoDiag;
  delete []PtoAcum;
//  delete []PtoMrg;
  delete []errs;

  for(int j=0;j<nClasf;j++) delete []Icd[j];
  delete []Icd;
 
  return ok;
}
double **CreaMatrizDistancias(Ensemble *ens, Data *data, bool ExcluirTrain)
{
  int nClasf            = ens->GetClassifiersToUse();
  int nDatos            = data->GetNTotal();
  vector<double> &pesos = ens-> GetUsingWeights();
  double **Icd             = new double*[nClasf];

  ens->SetData(data);
  for(int i=0;i<nClasf;i++) {
    Classifier *c = ens->GetClassifier(i);
    Icd[i] = new double[nDatos];
    double w = pesos[i];
    for(int j=0;j<nDatos;j++) {
      bool EnTrain = ExcluirTrain &&
                         c->UsedInOriginalTrainingData(data->GetDatIniPos(j));
      Icd[i][j] = EnTrain ? 0 : (c->Classify(j)!=data->GetDatClass(j) ? -w : w);
    }
  }

  return Icd;
}
void DoOrdenaPorReduccionDeDistancias(Ensemble *ens, int dist, int ini, int fin)
{
  double **Icd;  //Matriz[iclasificador][idato] con 1 si bien clas. y -1 si mal.
  int *ClasReal; //Clase de cada dato
  double *margen;
  double *PtoAct;
  double *PtoTmp;
  double *PtoRef;
  double DistMin, DistTemp;
  int nDatos;
  int nClasf;
  vector<int> DistMinIndices;

  NomData *data = (NomData*)ens->GetData();
  nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  nClasf = ens->GetClassifiersToUse();

  if (nClasf==1) return;

  //Obtiene la clase de cada dato
  ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++) {
    ClasReal[j] = data->GetDatClass(ini+j);
  }

  //Obtiene la clsificacion para cada clasificador del conjunto y para cada dato
  //Tambien obtiene los dos clasificadores con menor error
  bool ExcluirTrain = false;
  Icd = CreaMatrizDistancias(ens, data, ExcluirTrain);

  margen = new double[nDatos];
  for(int i=0;i<nClasf;i++) {
    for(int j=0;j<nDatos;j++) {
      margen[j]+=Icd[i][j];
    }
  }
if (dist==0) cout << "------------->>>>>>>" << dist << "<<<<<<<" << endl;
  //Inicializa puntos en el espacio de datos.
  PtoAct = new double[nDatos];
  PtoTmp = new double[nDatos];
  PtoRef = new double[nDatos];
  bool metodo2 = false;
  for (int i=0;i<nDatos;i++) {
    PtoAct[i] = 0.0;
    PtoTmp[i] = 0.0;
    if (metodo2) PtoRef[i] = margen[i] > dist ? 1 : dist;
    else PtoRef[i] = dist==0 ? 1.0 : dist;
  }

  int nestimmin=-1;
  double distestimmin=1000000.0;
  int deende = 1;
  //Busca el clasificador que se acerque mas al punto PtoRef
  for(int iClasf=0;iClasf<nClasf;iClasf++) {
    DistMin = 1000000.0;
/*    for(int i=iClasf;i<nClasf;i++) {
      //Calcula distancia a PtoRef
      for(int j=0;j<nDatos;j++) {
        PtoTmp[j] = PtoAct[j] + Icd[i][j];
      }
      DistTemp = Distancia(PtoTmp, PtoRef, nDatos);
      if (DistTemp<DistMin) {
        //De momento es el que mas reduce la distancia
        DistMin = DistTemp;
        iDistMin = i;
      }
    ens->Exchange(iClasf, iDistMin);
    int *cc = Icd[iClasf];
    Icd[iClasf] = Icd[iDistMin];
    Icd[iDistMin] = cc;
    for(int j=0;j<nDatos;j++) {
      PtoAct[j] += Icd[iClasf][j];
   //   PtoRef[j] = dist + sqrt(iClasf+1);
    }
      */
    //Combinaciones de nClasf-iClasf elementos tomados de 1 en 1
    Combinaciones k(nClasf-iClasf, deende<nClasf-iClasf ? deende : nClasf-iClasf);
    while(k.HayMasCombinaciones()) {
      vector<int> comb = k.SiguienteCombinacion();

      for(int j=0;j<nDatos;j++)
        PtoTmp[j] = PtoAct[j];

      for (unsigned i=0;i<comb.size();i++) {
        //Calcula distancia a PtoRef
        for(int j=0;j<nDatos;j++) {
          PtoTmp[j] += ens->GetWeightInUse(iClasf) * Icd[iClasf+comb[i]][j];
        }
      }

      DistTemp = Distancia(PtoTmp, PtoRef, nDatos);
      if (DistTemp<=DistMin) {
        //De momento es el que mas reduce la distancia
        DistMin = DistTemp;
        DistMinIndices = comb;
      }
    }
    int iii = iClasf;
    for (unsigned i=0;i<DistMinIndices.size();i++, iii++) {
      ens->Exchange(iii, DistMinIndices[i]+iClasf);
      double *cc = Icd[iii];
      Icd[iii] = Icd[DistMinIndices[i]+iClasf];
      Icd[DistMinIndices[i]+iClasf] = cc;
      for(int j=0;j<nDatos;j++) {
        PtoAct[j] += ens->GetWeightInUse(iClasf) * Icd[iii][j];
        if (metodo2) PtoRef[j] = margen[i] > dist ? PtoAct[j]+1 : dist;
        if (dist==0) PtoRef[j] =  2.0 * sqrt(2.0*(iClasf+1));
      }
    }
    iClasf += DistMinIndices.size()-1;
    //Se guarda la distancia minima alcanzada para ver si es un buen
    //estimador del nUmero necesario de clasificadores
    if (distestimmin>DistMin) {
        distestimmin = DistMin;
        nestimmin = iClasf+1;
    }
  }

  //guardamos a fichero la distancia minima alcanzada i el num
  //de clasificadores utilizados
//  FILE *fff = fopen("distmin.txt", "at");
 // fprintf(fff, "%g\t%d\n", distestimmin, nestimmin);
  //fclose(fff);
  //char filename[256];
  //sprintf(filename, "dist%04dinfo.txt", 1000*dist/nClasf);
  //fff = fopen(filename, "at");
  //fprintf(fff, "%d\n", nestimmin);
  //fclose(fff);

  for(int i=0;i<nClasf;i++) delete Icd[i];
  delete []margen;
  delete []Icd;
  delete []ClasReal;
  delete []PtoTmp;
  delete []PtoRef;
  delete []PtoAct;
}
//---------------------------------------------------------------------
void DoOrdenaPorReduccionErrorPonderado(Ensemble *ens, int ini, int fin)
{
  int **Mcd;  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
//  double **Mpd;  //Matriz con todas los pesos de los clasificadores del ens para cada dato de ens->data
  double **Ddc;  //Matriz con el numero de votos para cada dato y para cada clase
  int ErrMin1, ErrMin2;
  int iErrMin1=-1/*, iErrMin2*/;
  int ErrAux;
  int ErrAct;
  int *ClasReal; //Clase de cada dato
  int *ClasAct;  //Clase del conjunto para elumero de clasificadores actuales
  int nClases;
  int nDatos;
  int nClasf;

  //Recupera los datos con que se construyo el clasificador
//  ens->ResetOriginalTrainingData();

  NomData *data = (NomData*)ens->GetData();
  nClases = data->NumClass;
  nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  nClasf = ens->GetClassifiersToUse();

  if (nClasf==1) return;

  //Obtiene la clase de cada dato
  ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++)
    ClasReal[j] = data->GetDatClass(ini+j);

  //Obtiene la clsificacion para cada clasificador del conjunto y para cada dato
  //Tambien obtiene los dos clasificadores con menor error
  Mcd = new int*[nClasf];
  ErrMin1 = ErrMin2 = nDatos;
  bool ExcluirTrain = false;
  for(int i=0;i<nClasf;i++) {
    Classifier *c = ens->GetClassifier(i);
    Mcd[i] = new int[nDatos];
    ErrAux = 0;
    for(int j=0;j<nDatos;j++) {
      bool EnTrain = ExcluirTrain &&
               c->UsedInOriginalTrainingData(data->GetDatIniPos(ini+j));
      Mcd[i][j] = EnTrain ? -1/*Asi lo ignora en los calculos*/ : c->Classify(ini+j);
      if (Mcd[i][j]!=-1 && Mcd[i][j]!=ClasReal[j]) ErrAux++;
//      Mcd[i][j] = c->Classify(ini+j);
//      if (Mcd[i][j]!=ClasReal[j]) ErrAux++;
    }
//    ens->SetClassifierWeight(i, log(((double)nDatos-(double)ErrAux)/(double)ErrAux));
    if (ErrAux<ErrMin1) {
      ErrMin2  = ErrMin1;
  //    iErrMin2 = iErrMin1;
      ErrMin1  = ErrAux;
      iErrMin1 = i;
    }
    else if (ErrAux<ErrMin2) {
      ErrMin2  = ErrAux;
    //  iErrMin2 = i;
    }
  }

  //Coloca los dos mejores clasificadores en primera posicion
  ens->Exchange(0, iErrMin1);
  int *cc = Mcd[0];
  Mcd[0] = Mcd[iErrMin1];
  Mcd[iErrMin1] = cc;
/*  ens->Exchange(1, iErrMin2);
  cc = Mcd[1];
  Mcd[1] = Mcd[iErrMin2];
  Mcd[iErrMin2] = cc;*/

  //Distribuciones de clase
  Ddc = new double*[nDatos];
  for(int i=0;i<nDatos;i++) {
    Ddc[i] = new double[nClases];
    for(int j=0;j<nClases;j++) {
      Ddc[i][j] = 0.0;
    }
  }
  ClasAct = new int[nDatos];
  ErrAct = 0;
  for(int j=0;j<nDatos;j++) {
    if(Mcd[0][j]>=0) Ddc[j][Mcd[0][j]]++;
//    if(Mcd[1][j]>=0) Ddc[j][Mcd[1][j]]++;
    ClasAct[j] = ens->WhichClass(Ddc[j], nClases);
    if (ClasAct[j]<0) ClasAct[j] = Mcd[0][j];
    if (ClasAct[j]!=ClasReal[j]) ErrAct++;
  }

  //Busca el clasificador que contribuya mas a mejorar el error del conjunto.
  //Para ello se suma uno por cada dato bien clasificado por un clasificador y
  //que mal clasificando por el conjunto (clasificadores anyadidos hasta el momento)
  double *Dctemp = new double[nClases];
//  int *w = new int[nClases+1];
  double ErrTemp, ErrMinimo;
  vector<int> ComErrMin;
  int deende = 1;
  for(int iClasf=1;iClasf<nClasf;iClasf++) {
    ErrMinimo = -1.0;
/*    for(int i=iClasf;i<nClasf;i++) {
      ErrTemp = 0.0;
      //Calcula lo que el clasf. i variaria el error de conjunto
      for(int j=0;j<nDatos;j++) {
        if (ClasReal[j]==ClasAct[j] || ClasReal[j]!=Mcd[i][j]) continue;
//        double tot=0;
//        for(int k=0;k<nClases;k++) {
//          tot += Ddc[j][k];
//        }
        ErrTemp += 1;//(Ddc[j][ClasReal[j]]==0.0) ? 5 : tot/Ddc[j][ClasReal[j]];
      }
      if (ErrTemp>ErrMinimo) {
        //De momento es el que mas baja el error
        ErrMinimo = ErrTemp;
        iErrMin1 = i;
      }
    }
    ens->Exchange(iClasf, iErrMin1);
    cc = Mcd[iClasf];
    Mcd[iClasf] = Mcd[iErrMin1];
    Mcd[iErrMin1] = cc;
    for(int j=0;j<nDatos;j++) {
      Ddc[j][Mcd[iClasf][j]]++;
      int ct = ens->WhichClass(Ddc[j], nClases);
      if (ct>=0) ClasAct[j] = ct;
    }*/
    Combinaciones k(nClasf-iClasf, deende<nClasf-iClasf ? deende : nClasf-iClasf);
    while(k.HayMasCombinaciones()) {
      vector<int> comb = k.SiguienteCombinacion();
      ErrTemp = 0.0;
      //Calcula lo que el clasf. i variaria el error de conjunto
      for (unsigned i=0;i<comb.size();i++) {
        for(int j=0;j<nDatos;j++) {
          if (ClasReal[j]==ClasAct[j] || ClasReal[j]!=Mcd[iClasf+comb[i]][j])
            continue;
          ErrTemp += ens->GetWeightInUse(iClasf+comb[i]);
        }
      }
      if (ErrTemp>ErrMinimo) {
        //De momento es el que mas baja el error
        ErrMinimo = ErrTemp;
        ComErrMin = comb;
      }
    }
    int iii = iClasf;
    for (unsigned i=0;i<ComErrMin.size();i++, iii++) {
      ens->Exchange(iii, ComErrMin[i]+iClasf);
      cc = Mcd[iii];
      Mcd[iii] = Mcd[ComErrMin[i]+iClasf];
      Mcd[ComErrMin[i]+iClasf] = cc;
      for(int j=0;j<nDatos;j++) {
        if (Mcd[iii][j]>=0) Ddc[j][Mcd[iii][j]]++;
        int ct = ens->WhichClass(Ddc[j], nClases);//, w, &(ens->APrioriClas));
        ClasAct[j] =  ct;//w[0] == 1 ? ct : -1;
      }
    }
    iClasf += ComErrMin.size()-1;
  }

//  delete []w;
  for(int i=0;i<nDatos;i++) delete Ddc[i];
  delete []Ddc;
  for(int i=0;i<nClasf;i++) delete Mcd[i];
  delete []Mcd;
  delete ClasReal;
  delete ClasAct;
  delete Dctemp;
}
//---------------------------------------------------------------------
void DoOrdenaPorEPIC(Ensemble *ens, NomData *data)
{
  int **Mcd;  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
//  double **Mpd;  //Matriz con todas los pesos de los clasificadores del ens para cada dato de ens->data
  double **Ddc;  //Matriz con el numero de votos para cada dato y para cada clase
  int *ClasReal; //Clase de cada dato
  int nClases;
  int nDatos;
  int nClasf;
  double *W;

  //Recupera los datos con que se construyo el clasificador
//  ens->ResetOriginalTrainingData();

  NomData *dataHold = (NomData*)ens->GetData();
  ens->SetData(data);
  nClases = data->NumClass;
  nDatos = data->GetNTotal();
  nClasf = ens->GetClassifiersToUse();
  W = new double[nClasf];

  if (nClasf==1) return;

  //Obtiene la clase de cada dato
  ClasReal = new int[nDatos];
  for(int j = 0; j < nDatos; j++) 
    ClasReal[j] = data->GetDatClass(j);

  //Obtiene la clsificacion para cada clasificador del conjunto y para cada dato
  //Tambien obtiene los dos clasificadores con menor error
  Mcd = new int*[nClasf];
  for(int i = 0; i < nClasf; i++) {
    Classifier *c = ens->GetClassifier(i);
    Mcd[i] = new int[nDatos];
    for(int j = 0; j < nDatos; j++) {
      Mcd[i][j] = c->Classify(j);
    }
  }

  //Distribuciones de clase
  Ddc = new double*[nDatos];
  for(int i = 0; i < nDatos; i++) {
    Ddc[i] = new double[nClases+2];
    for(int j = 0; j < nClases; j++) {
      Ddc[i][j] = 0.0;
    }
  }
  for(int j = 0; j < nDatos; j++) {
    for(int i = 0; i < nClasf; i++) {
      Ddc[j][Mcd[i][j]]++;
    }
    Ddc[j][nClases] = ens->WhichClass(Ddc[j], nClases); //Majority class index;
    Ddc[j][nClases+1] = 0;
    for(int k = 0; k < nClases; k++) {//Second common class votes
      if (k==Ddc[j][nClases]) continue;
      if (Ddc[j][k] >= Ddc[j][nClases+1]) Ddc[j][nClases+1] = Ddc[j][k];
    }
  }

  for(int iClasf = 0; iClasf < nClasf; iClasf++) {
    W[iClasf] = 0;
    for(int j = 0; j < nDatos; j++) {
      if ( ClasReal[j] == Mcd[iClasf][j] && Mcd[iClasf][j] != Ddc[j][nClases]) { //Classifier ok && in minority group
        int kk=Ddc[j][nClases];
        W[iClasf] += 2*Ddc[j][kk]; 
        W[iClasf] -= Ddc[j][ Mcd[iClasf][j] ] ;
      }
      else if ( ClasReal[j] == Mcd[iClasf][j] && Mcd[iClasf][j] == Ddc[j][nClases]) //Classifier ok && in mayority group
        W[iClasf] += Ddc[j][ nClases+1 ];
      else if ( ClasReal[j] != Mcd[iClasf][j]) { //Classifier wrong
        W[iClasf] += Ddc[j][ ClasReal[j] ];
        W[iClasf] += -Ddc[j][ Mcd[iClasf][j] ];
        int kk=Ddc[j][nClases];
        W[iClasf] += -Ddc[j][kk];
      }
      else
        exit(1); //We could not get here
    }
  }

  ens->OrdenarClasificadoresPorPeso(W);

  ens->SetData(dataHold);

  delete []W;
  for(int i=0;i<nDatos;i++) delete []Ddc[i];
  delete []Ddc;
  for(int i=0;i<nClasf;i++) delete []Mcd[i];
  delete []Mcd;
  delete []ClasReal;
}
//---------------------------------------------------------------------
void DoOrdenaPorReduccionErrorDoble(Ensemble *ens, int ini, int fin)
{
  int **Mcd;  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
  double **Ddc;  //Matriz con el numero de votos para cada datos y para cada clase
  int ErrMin;
  int iErrMin=-1, iErrMin2=-1;
  int ErrAux;
  int ErrAct;
  int *ClasReal; //Clase de cada dato
  int *ClasAct;  //Clase del conjunto para elumero de clasificadores actuales
  int nClases;
  int nDatos;
  int nClasf;

  //Recupera los datos con que se construyo el clasificador
//  ens->ResetOriginalTrainingData();

  NomData *data = (NomData*)ens->GetData();
  nClases = data->NumClass;
  nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  nClasf = ens->GetClassifiersToUse();

  //Obtiene la clase de cada dato
  ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++)
    ClasReal[j] = data->GetDatClass(ini+j);

  //Obtiene la clsificacion para cada clasificador del conjunto y para cada dato
  //Tambien obtiene los dos clasificadores con menor error
  Mcd = new int*[nClasf];
  ErrMin = nDatos;
  for(int i=0;i<nClasf;i++) {
    Classifier *c = ens->GetClassifier(i);
    Mcd[i] = new int[nDatos];
    ErrAux = 0;
    for(int j=0;j<nDatos;j++) {
      Mcd[i][j] = c->Classify(ini+j);
      if (Mcd[i][j]!=ClasReal[j]) ErrAux++;
    }
//    ens->SetClassifierWeight(i, log(((double)nDatos-(double)ErrAux)/(double)ErrAux));
    if (ErrAux<ErrMin) {
      ErrMin  = ErrAux;
      iErrMin = i;
    }
  }

  //Coloca el mejor clasificador en primera posicion
  ens->Exchange(0, iErrMin);
  int *cc = Mcd[0];
  Mcd[0] = Mcd[iErrMin];
  Mcd[iErrMin] = cc;

  //Distribuciones de clase
  Ddc = new double*[nDatos];
  for(int i=0;i<nDatos;i++) {
    Ddc[i] = new double[nClases];
    for(int j=0;j<nClases;j++) {
      Ddc[i][j] = 0.0;
    }
  }
  ClasAct = new int[nDatos];
  ErrAct = 0;
  for(int j=0;j<nDatos;j++) {
    Ddc[j][Mcd[0][j]]++;
    ClasAct[j] = Mcd[0][j];
    if (ClasAct[j]!=ClasReal[j]) ErrAct++;
  }

  //Busca el clasificador que reduzca mas el error(o lo aumente menos)
  double *Dctemp = new double[nClases];
  int ErrTemp;
  for(int iClasf=1;iClasf<nClasf;iClasf+=2) {
    if (iClasf==nClasf-1) break;
    ErrMin = nDatos;
    for(int i=iClasf;i<nClasf-1;i++) {
     for(int ii=i+1;ii<nClasf;ii++) {
      ErrTemp = 0;
      for(int j=0;j<nDatos;j++) {
        if (Mcd[i][j]==ClasAct[j] || Mcd[ii][j]==ClasAct[j]) continue;
    //    else if (Ddc[j][Mcd[i][j]]+1<Ddc[j][ClasAct[j]])
      //    continue; //Este ejemplo ni sube ni baja el error
        for(int k=0;k<nClases;k++) {
          Dctemp[k] = Ddc[j][k];
        }
        Dctemp[Mcd[i][j]]++;
        Dctemp[Mcd[ii][j]]++;
        int ClasTemp = ens->WhichClass(Dctemp, nClases);
        if (ClasTemp>=0 && ClasTemp!=ClasAct[j]) {
          if (ClasAct[j]==ClasReal[j]) ErrTemp++;    //Sube el error
          else if (ClasTemp==ClasReal[j]) ErrTemp--; //Baja el error
        }
      }
      if (ErrTemp<ErrMin || (ErrTemp==ErrMin &&
            ens->GetClassifierWeight(i) > ens->GetClassifierWeight(iErrMin))) {
        ErrMin = ErrTemp;
        iErrMin = i;
        iErrMin2 = ii;
      }
     }
    }
    ens->Exchange(iClasf, iErrMin);
    cc = Mcd[iClasf];
    Mcd[iClasf] = Mcd[iErrMin];
    Mcd[iErrMin] = cc;
    ens->Exchange(iClasf+1, iErrMin2);
    cc = Mcd[iClasf+1];
    Mcd[iClasf+1] = Mcd[iErrMin2];
    Mcd[iErrMin2] = cc;
    ErrAct += ErrMin;
    for(int j=0;j<nDatos;j++) {
      Ddc[j][Mcd[iClasf][j]]++;
      Ddc[j][Mcd[iClasf+1][j]]++;
      int ct = ens->WhichClass(Ddc[j], nClases);
      if (ct>=0) ClasAct[j] = ct;
    }
  }

  for(int i=0;i<nDatos;i++) delete Ddc[i];
  delete []Ddc;
  for(int i=0;i<nClasf;i++) delete Mcd[i];
  delete []Mcd;
  delete ClasReal;
  delete ClasAct;
  delete Dctemp;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Intento 2 de maximizacion de una ec. lineal con restricciones lineales
// La idea es maximizar la clasificacion del dato mas dificil de clasificar
// teniendo en cuenta que el resto de datos han de estar bien clasificados
//---------------------------------------------------------------------------
void DoOrdenaSimplex2(Ensemble *ens, int ini, int fin)
{
  NomData *data = (NomData*)ens->GetData();
  //int nClases = data->NumClass;
  int nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  int nClasf = ens->GetClassifiersToUse();
  int M = nDatos-1+1; //Numero de restricciones
  int N = nClasf; //Numero de variables
  int clase;

  //Se ordenan los datos por orden de peor clasificacion
  int iColAux = data->GetNumVar();
  for (int i=0;i<nDatos;i++)
    data->SetValueVar(i, iColAux, ens->ClassificationCertainty(ini+i, clase));
  data->SortOn(iColAux, ini, ini+nDatos-1);

  //Obtiene la matriz de clasificaciones
  Matriz *Cmn = CreaMatrizCmn(ens, ini, fin);

//  double **a = new double*[M+3];
  double **a = nr::matrix(1, M+2, 1, N+1);//new double*[M+3];
  //for(int i=0;i<M+3;i++)
  //  a[i] = new double[N+2];
   int icase;
  int *izrov = new int[N+1];
  int *iposv = new int[M+1];
 do {
  for(int i=1;i<M+3;i++)
    for(int j=1;j<N+2;j++)
      a[i][j] = 0;

  for(int i=0;i<nClasf;i++)
    a[1][i+2] = (*Cmn)[i][0];  //Funcion objetivo
  a[1][1] = 0;

//Suma de las pesos de los clasificadores ha de ser <= a 1
  for(int i=1;i<=N;i++) a[2][i+1] = -1; //Restricciones m1

  for(int j=1;j<M/*nDatos*/;j++) {
    a[j+2][1] = 0.5;
    for(int i=0;i<N/*nClasf*/;i++)
      a[j+2][i+2] = -(*Cmn)[i][j];  //Restricciones tipo m2
  }


  FILE *f=fopen("c:\\kk.txt", "wt");
  for(int i=1;i<M+3;i++) {
    for(int j=1;j<N+2;j++) {
      fprintf(f, "\t%f", a[i][j]);
    }
    fprintf(f, "\n");
  }

  nr::simplx(a,M,N,1,M-1,0,&icase,izrov,iposv);

  for(int i=1;i<M+3;i++) {
    for(int j=1;j<N+2;j++) {
      fprintf(f, "\t%f", a[i][j]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n\n");
  for(int j=1;j<M+1;j++) {float kl=iposv[j]   ;
    fprintf(f, "\tk%f", kl);
  }
  fprintf(f, "\n\n");
  for(int j=1;j<N+1;j++) {    float kj= izrov[j];
    fprintf(f, "\tk%f", kj);
  }
  fclose(f);
  M--;
 } while (icase!=0);

 M++;
  if (icase==0) {
    int j;
    for(j=1;iposv[j]<=N && j<=M;j++){
      ens->SetClassifierWeight(iposv[j]-1, a[j+1][1]);
    }
    for(j=1;izrov[j]<=N && j<=N;j++){
      ens->SetClassifierWeight(izrov[j]-1, 0.0);
    }
    ens->OrdenarClasificadoresPorPeso();
  }

 // for(int i=0;i<M+3;i++)
   // delete []a[i];
  //delete []a;
  nr::free_matrix(a, 1, M+2, 1, N+1);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void DoOrdenaSimplex(Ensemble *ens, int ini, int fin)
{
//1ro
//NO FUNCIONA. DA COMO SOLUCION TODOS LOS ARBOLES CON PESO 1
  NomData *data = (NomData*)ens->GetData();
  //int nClases = data->NumClass;
  int nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  int nClasf = ens->GetClassifiersToUse();
  int M = nClasf; //Numero de restricciones
  int N = nClasf; //Numero de variables

  int *ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++)
    ClasReal[j] = data->GetDatClass(ini+j);

  double **a = nr::matrix(1, M+2, 1, N+1);//new double*[M+3];
  for(int i=1;i<=M+2;i++) {
//    a[i] = new double[N+2];
//    a[i][0] = 0;
    a[i][1] = 1;
    for(int j=2;j<N+2;j++) {
      a[i][j] = 0;
    }
  }

//Los clasificadores tienen un peso <=1
  for(int i=1;i<=M;i++) a[i+1][i+1] = -1; //Restricciones m1
//Suma de las pesos de los clasificadores ha de ser <= a 1
//  for(int i=1;i<=N;i++) a[2][i+1] = -1; //Restricciones m1

  int **Mcd;  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
  Mcd = new int*[nClasf];
  for(int i=0;i<nClasf;i++) {
    double PesoClasif = 0;
    Classifier *c = ens->GetClassifier(i);
      ens->SetClassifierWeight(i, 0);///////////
    Mcd[i] = new int[nDatos];
    for(int j=0;j<nDatos;j++) {
      int clase = c->Classify(ini+j);
      Mcd[i][j] = (clase==ClasReal[j]) ? 1 : -1;//Falta comprobar si el dato esta en validacion o en train
      PesoClasif += Mcd[i][j];
    }
    a[1][i+2] = PesoClasif;  //Funcion a maximizar
  }
  a[1][1] = 0;

  int *izrov = new int[N+1];
  int *iposv = new int[M+1];
  int icase;
  FILE *f=fopen("c:\\kk.txt", "wt");
  for(int i=1;i<=M+2;i++) {
    for(int j=1;j<=N+1;j++) {
      fprintf(f, "\t%f", a[i][j]);
    }
    fprintf(f, "\n");
  }

  nr::simplx(a,M,N,M,0,0,&icase,izrov,iposv);

  for(int i=1;i<=M+2;i++) {
    for(int j=1;j<=N+1;j++) {
      fprintf(f, "\t%f", a[i][j]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n\n");
  for(int j=1;j<=M+1;j++) {float kl=iposv[j]   ;
    fprintf(f, "\tk%f", kl);
  }
  fprintf(f, "\n\n");
  for(int j=1;j<=N;j++) {    float kj= izrov[j];
    fprintf(f, "\tk%f", kj);
  }
  fclose(f);

  if (icase==0) {
    int j;
    for(j=1;iposv[j]<=N && j<=M;j++){
      ens->SetClassifierWeight(iposv[j]-1, a[j+1][1]);
    }
    for(j=1;izrov[j]<=N && j<=N;j++){
      ens->SetClassifierWeight(izrov[j]-1, 0.0);
    }
    ens->OrdenarClasificadoresPorPeso();
/*    for(j=1;iposv[j]<=N;j++){
      ens->SetClassifierWeight(iposv[j]-1, a[iposv[j]+1][1]);
    }*/

  }

/*  for(int i=0;i<M+3;i++)
    delete []a[i];
  delete []a;*/
  nr::free_matrix(a, 1, M+2, 1, N+1);
  delete []ClasReal;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class F : public MatrizVariable
{
  Matriz *Cmn;
  Matriz *CmnAbs;
  public:
    F(Matriz *_Cmn) : MatrizVariable(_Cmn->numeroFilasMatriz()-1)  {
      Cmn = new Matriz(*_Cmn);
      CmnAbs = new Matriz(*_Cmn);
      for (int i=0;i<Cmn->numeroFilasMatriz();i++)
        for (int j=0;j<Cmn->numeroColumnasMatriz();j++)
          if (Cmn->valorElemento(i, j)<0.0)
            CmnAbs->valorElemento(i, j, -Cmn->valorElemento(i, j));
    }
    virtual ~F(){}

    virtual void asignarValores(Matriz *O) {
      double UnoMenosSum = 1;
      for(int m=0;m<numeroColumnas;m++) {
        UnoMenosSum = UnoMenosSum - O->m(m, 0);
      }
      int N = Cmn->numeroColumnasMatriz();
      int K = numeroFilas;
      int M = numeroColumnas;
//      FILE *fff=fopen("klon123.txt", "wt");
      for(int k=0;k<K;k++) {
        for(int m=0;m<M;m++) {
          double val = 0;
          for(int n=0;n<N;n++) {
            double div=0;
            double nval = CmnAbs->m(m,n)*(Cmn->m(k,n) - Cmn->m(M,n));
            nval += Cmn->m(m,n)*(CmnAbs->m(M,n) - CmnAbs->m(k,n));
            nval += Cmn->m(M,n) * (CmnAbs->m(k,n));
            nval -= Cmn->m(k,n) * (CmnAbs->m(M,n));
  //          fprintf(fff, "%lf\t", nval);
            for(int l=0;l<Cmn->numeroFilasMatriz()-1;l++) {
              div += CmnAbs->m(l,n) * (O->m(l, 0));
            }
            div = (div + CmnAbs->m(M,n)*UnoMenosSum) *
                                             (div + CmnAbs->m(M,n)*UnoMenosSum);
            if (div==0.0) continue;
            nval = nval / div;
            val += nval;
          }
          //valorElemento(k, m, val);//
          mat[k][m] = val;
    //      fprintf(fff, "\n");
        }
      }
      //fclose(fff);
    }
};
class C : public MatrizVariable
{
  Matriz *Cmn;
  Matriz *CmnAbs;
  public:
    C(Matriz *_Cmn) : MatrizVariable(_Cmn->numeroFilasMatriz()-1, 1) {
      Cmn = _Cmn;
      CmnAbs = new Matriz(*_Cmn);
      for (int i=0;i<Cmn->numeroFilasMatriz();i++)
        for (int j=0;j<Cmn->numeroColumnasMatriz();j++)
          if (Cmn->valorElemento(i, j)<0.0)
            CmnAbs->valorElemento(i, j, -Cmn->valorElemento(i, j));
    }
    virtual ~C(){}

    virtual void asignarValores(Matriz *O) {
      double UnoMenosSum = 1;
      for(int m=0;m<numeroColumnas;m++) {
        UnoMenosSum = UnoMenosSum - O->m(m, 0);
      }
      int N = Cmn->numeroColumnasMatriz();
      int K = numeroFilas;
      int M = numeroColumnas;
      for(int k=0;k<K;k++) {
        double val = 0;
        for(int n=0;n<N;n++) {
          double div=0;
          double nval = Cmn->m(k,n) * (CmnAbs->m(M,n)) -
                        Cmn->m(M,n) * (CmnAbs->m(k,n));
          for(int l=0;l<Cmn->numeroFilasMatriz()-1;l++) {
            div += CmnAbs->m(l,n) * (O->m(l, 0));
          }
          div = (div + CmnAbs->m(M,n)*UnoMenosSum) *
                                             (div + CmnAbs->m(M,n)*UnoMenosSum);
          if (div==0.0) continue;
          nval = nval / div;
          val += nval;
        }
        mat[k][0] = val;
      }
    }
};
Matriz* CreaMatrizCmn(Ensemble *ens, int ini, int fin, bool ExcluirTrain)
{
  NomData *data = (NomData*)ens->GetData();
  int nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  int nClasf = ens->GetClassifiersToUse();
  int M = nClasf;
  int N = nDatos;

  //Obtiene la clase real de cada elemento
  int *ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++)
    ClasReal[j] = data->GetDatClass(ini+j);

  //Calculala matriz Cmn. Matriz con todas las clasificaciones de los
  // clasificadores del ens para cada dato de ens->data
//  Mandato *vVal = Variable::VariablePorNombre("vVal", true, Mandato);
  double vValNeg = -1;
  double vValPos = 1;

  Matriz *Cmn = new Matriz(M, N);
//  int iColPos  = data->GetNumVar()+ Data::IniPosIndex;
  for(int m=0;m<M;m++) {
    Classifier *c = ens->GetClassifier(m);
    for(int n=0;n<N;n++) {
      int clase = c->Classify(ini+n);
      bool Val = ExcluirTrain &&
               c->UsedInOriginalTrainingData(data->GetDatIniPos(ini+n));
      (*Cmn)[m][n] = Val ? 0 : (clase==ClasReal[n] ? vValPos : vValNeg);
    }
  }

  delete[] ClasReal;

  return Cmn;
}
void DoMinimoIterativo(Ensemble *ens, int ini, int fin)
{
  //NomData *data = (NomData*)ens->GetData();
  //int nClases = data->NumClass;
  //int nDatos = (fin<0) ? data->GetNTrain() : fin-ini+1;
/*UNIX  AnsiString Value = IntToStr(ens->GetClassifiersToUse());
  if (InputQuery("ACaption", "APrompt", Value)) {
    int val;
    try {
      val = Value.ToInt();
      ens->SetClassifiersToUse(val);
    }
    catch(...){
    }
  }
*/
/*  int nClasf = ens->GetClassifiersToUse();
  int M = nClasf;
  int N = nDatos;

  int *ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++)
    ClasReal[j] = data->GetDatClass(ini+j);
*/
  //Calculala matriz Cmn
  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
//UNIX  Matriz *Cmn = CreaMatrizCmn(ens, ini, fin);
	/*new Matriz(M, N);
//  Matriz *Om = new Matriz(M-1, 1);
  int iColPos  = data->GetNumVar()+ Data::IniPosIndex;
  for(int m=0;m<M;m++) {
    Classifier *c = ens->GetClassifier(m);
  //  Om->valorElemento(m, 0, 1.0/(M-1));
    ens->SetClassifierWeight(m, 0);///////////
    for(int n=0;n<N;n++) {
      int clase = c->Classify(ini+n);
      bool Val=c->UsedInOriginalTrainingData(data->GetValueVar(ini+n, iColPos));
      Cmn->valorElemento(m, n, Val ? 0 : (clase==ClasReal[n] ? 1 : -1));
    }
  }
*/
 /*UNIX F *f = new F(Cmn);
  C *c = new C(Cmn);
  int M = f->numeroFilasMatriz() + 1;
  Matriz *Vars = new Matriz(M-1, 1);
  for(int i=0;i<Vars->numeroFilasMatriz();i++) (*Vars)[i][0] = 1.0/M;
*/
  //    f->asignarValores(Vars);
 //   c->asignarValores(Vars);
 /* for(int ii=0;ii<10;ii++) {
    f->asignarValores(Vars);
    c->asignarValores(Vars);
    Matriz::ResolverEcLineal(f, Vars, c);
//    Synchronize(Pinta);
  }*/
/*UNIX
  Form3->Execute(Vars, f, c);

  for(int j=0;j<M-1;j++){
    ens->SetClassifierWeight(j, (*Vars)[j][0]);
  }
  ens->OrdenarClasificadoresPorPeso();


  delete Cmn;*/
//  delete []ClasReal;
}
//---------------------------------------------------------------------------
class AuxData {
 public:
  static Matriz *Cmn;
  static Matriz *CmnNorm;
  static int N;
  static int M;
  static double gamma;
  static double theta;
  static double k;
  static VariableNumeroHistorico *vhfhist;
  static VariableNumeroHistorico *fvaltr;
  static VariableNumeroHistorico *fvalts;
  static VariableNumeroHistorico *pesos;
  static VariableNumeroHistorico *orden;
  static void init(Matriz *_Cmn) {
    if (!vhfhist) vhfhist = new VariableNumeroHistorico("fhist", 0.0);
    if (!fvaltr) fvaltr = new VariableNumeroHistorico("fvaltr", 0.0);
    if (!fvalts) fvalts = new VariableNumeroHistorico("fvalts", 0.0);
    if (!pesos) pesos = new VariableNumeroHistorico("pesos", 0.0);
    if (!orden) orden = new VariableNumeroHistorico("orden", 0.0);
    vhfhist->Inicializar();
    fvaltr->Inicializar();
    fvalts->Inicializar();
    pesos->Inicializar();
    orden->Inicializar();

    Cmn = _Cmn;
    M = Cmn->numeroFilasMatriz();
    N = Cmn->numeroColumnasMatriz();
    if (CmnNorm) delete CmnNorm;
    CmnNorm = new Matriz(*Cmn);
    for (int i=0;i<Cmn->numeroFilasMatriz();i++)
      for (int j=0;j<Cmn->numeroColumnasMatriz();j++)
        if (Cmn->valorElemento(i, j)<0.0)
          CmnNorm->valorElemento(i, j, -Cmn->valorElemento(i, j));
  }
};
Matriz *AuxData::Cmn = 0;
Matriz *AuxData::CmnNorm = 0;
int AuxData::N = 0;
int AuxData::M = 0;
double AuxData::gamma = 0.1;
double AuxData::theta = 0.183939720585721;// = 0.5/exp(1)
double AuxData::k = 50;
VariableNumeroHistorico *AuxData::vhfhist = 0;
VariableNumeroHistorico *AuxData::fvaltr = 0;
VariableNumeroHistorico *AuxData::fvalts = 0;
VariableNumeroHistorico *AuxData::pesos = 0;
VariableNumeroHistorico *AuxData::orden = 0;

//FUNCIONES DE COSTE

double f(double p[])
{
  double res = 0.0;

  for (int n=0;n<AuxData::N;n++) {
    double num = 0.0;
    double den = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      num += (*AuxData::Cmn)[m][n] * (p[m+1]);
      den += (*AuxData::CmnNorm)[m][n] * p[m+1];
    }
    if (den!=0.0) {
      double val = num/den;
      if (val<=1.001 && val>=-1.001) res += num/den; // Para evitar desbordamientos
      //con p pequenyos.
    }
  }

  AuxData::vhfhist->SetNumero(-res);

  return -res;
}
double f(double p[], int ini, int fin)
{
  double res = 0.0;

  for (int n=ini;n<=fin;n++) {
    double num = 0.0;
    double den = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      num += (*AuxData::Cmn)[m][n] * (p[m+1]);
      den += (*AuxData::CmnNorm)[m][n] * p[m+1];
    }
    if (den!=0.0) res += num/den;
  }

  AuxData::vhfhist->SetNumero(-res);

  return -res;
}
double recta(double x, double a=-1.0, double b=0.0)
{
  return a*x+b;
}
double frecta(double p[], int ini, int fin, double _a=-1.0, double _b=0.0
                                                          , bool Variable=false)
{
  double res = 0.0;

  for (int n=ini;n<=fin;n++) {
    double num = 0.0;
    double den = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      num += (*AuxData::Cmn)[m][n] * (p[m+1]);
      den += (*AuxData::CmnNorm)[m][n] * p[m+1];
    }
    if (!Variable && den==0.0) continue;
    res += Variable ? recta(num, 1.0/den, 0.0) : recta(num, _a, _b) ;
  }

  return res;
}
double exp02(double p[], int ini, int fin)
{
  double res = 0.0;

  for (int n=ini;n<=fin;n++) {
    double num = 0.0;
    double den = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      num += (*AuxData::Cmn)[m][n] * (p[m+1]);
      den += (*AuxData::CmnNorm)[m][n] * p[m+1];
    }
    if (den!=0.0) res += exp((double)-num/den);
  }

  return -res;
}
double exp01(double p[], int ini, int fin, double _min=-1.0, double _max=1.0
                                                          , bool Variable=false)
{
  double res = 0.0;

  for (int n=ini;n<=fin;n++) {
    double num = 0.0;
    double max = 0.0;
    double min = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      num += (*AuxData::Cmn)[m][n] * (p[m+1]);
      if ((*AuxData::Cmn)[m][n]>0)  max+=(*AuxData::Cmn)[m][n];
      else min+=(*AuxData::Cmn)[m][n];
    }
    min/=AuxData::M;
    max/=AuxData::M;
    double x = Variable ? num*2.0/(min-max)+1.0-2.0*min/(min-max) :
                          num*2.0/(_min-_max)+1.0-2.0*_min/(_min-_max);
    res += exp(x);
  }

  return res;
}
double log01(double p[], int ini, int fin, double _min=-1.0, double _max=1.0
                                                          , bool Variable=false)
{
  double res = 0.0;

  for (int n=ini;n<=fin;n++) {
    double num = 0.0;
    double max = 0.0;
    double min = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      num += (*AuxData::Cmn)[m][n] * (p[m+1]);
      if ((*AuxData::Cmn)[m][n]>0)  max+=(*AuxData::Cmn)[m][n];
      else min+=(*AuxData::Cmn)[m][n];
    }
    min/=AuxData::M;
    max/=AuxData::M;
    double x = Variable ? num*2.0/(min-max)+1.0-2.0*min/(min-max) :
                          num*2.0/(_min-_max)+1.0-2.0*_min/(_min-_max);
    if (x>=1.0) x=0.99999;
    res += log(1/(1-x));
  }

  return res;
}
double pwsderivable(double p[])
{
  double res = 0.0;
  double th = AuxData::theta;
  double gm = AuxData::gamma;
  double a1 = -gm;
  double b1 = 1.2-gm;
  double a2 = -(1.2-2.0*gm)/th;
  double b2 = (1.2-gm);
  double a3 = -gm/(1.0-th);
  double b3 = gm/(1.0-th);
  double k = AuxData::k;

  for (int n=0;n<AuxData::N;n++) {
    double num = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      num += (*AuxData::Cmn)[m][n] * (p[m+1]);
    }
    double y = 0.0;
    double x = num;
    y  = y + (a1*x+b1)*1.0/(1+exp((double)k*x));
    y  = y + (1.0/(1+exp((double)-k*x)) - 1.0/(1+exp((double)k*(th-x))))*(a2*x+b2);
    y  = y + (a3*x+b3)*1.0/(1+exp((double)k*(th-x)));
    res += y;
  }

  return res;
}
void dpwsderivable(double p[], double g[])
{
  double th = AuxData::theta;
  double gm = AuxData::gamma;
  double a1 = -gm;
  double b1 = 1.2-gm;
  double a2 = -(1.2-2.0*gm)/th;
  double b2 = (1.2-gm);
  double a3 = -gm/(1.0-th);
  double b3 = gm/(1.0-th);

  vector<double> Mn; //Margen del dato n
  for (int n=0;n<AuxData::N;n++) {
    double num = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      num += (*AuxData::Cmn)[m][n] * (p[m+1]);
    }
    Mn.push_back(num);
  }

  for (int k=0;k<AuxData::M;k++) {
    g[k+1] = 0.0;
    for (int n=0;n<AuxData::N;n++) {
      if ((*AuxData::Cmn)[k][n]==0.0) continue;
      double x = Mn[n];
      double e1 = 1.0/(1+exp((double)k*x));
      double e21 = 1.0/(1+exp((double)-k*x));
      double e22 = 1.0/(1+exp((double)k*(th-x)));
      double e3 = 1.0/(1+exp((double)k*(th-x)));
      double y2 = 0;
      y2 = y2 + a1*e1+(a1*x+b1)*(e1*e1)*(-k*exp((double)x*k));
      y2 = y2 + a2*(e21-e22)+(a2*x+b2)*((e21*e21)*(k*exp((double)-x*k)) - (e22*e22)*(k*exp((double)(th-x)*k)));
      y2 = y2 + a3*e3+(a3*x+b3)*(e3*e3)*(k*exp(((double)th-x)*k));
      g[k+1] += y2*(*AuxData::Cmn)[k][n];
    }
  }
}
double pwsderivable(double p[], int dato, double k)
{
  double res = 0.0;
  double th = AuxData::theta;
  double gm = AuxData::gamma;
  double a1 = -gm;
  double b1 = 1.2-gm;
  double a2 = -(1.2-2.0*gm)/th;
  double b2 = (1.2-gm);
  double a3 = -gm/(1.0-th);
  double b3 = gm/(1.0-th);

  for (int n=dato;n<=dato;n++) {
    double num = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      num += (*AuxData::Cmn)[m][n] * (p[m+1]);
    }
    double y = 0.0;
    double x = num;
    y  = y + (a1*x+b1)*1.0/(1+exp((double)k*x));
    y  = y + (1.0/(1+exp((double)-k*x)) - 1.0/(1+exp((double)k*(th-x))))*(a2*x+b2);
    y  = y + (a3*x+b3)*1.0/(1+exp((double)k*(th-x)));
    res += y;
  }

  return res;
}
double pws01(double p[], int ini, int fin, double gamma=0.1, double theta=0.5)
{
  double res = 0.0;

  for (int n=ini;n<=fin;n++) {
    double num = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      num += (*AuxData::Cmn)[m][n] * (p[m+1]);
    }
    num = (num<0.0)   ? (1.2-gamma) - gamma*num :
          (num<theta) ? (1.2-gamma) - (1.2-2.0*gamma)*num/theta :
                        gamma/(1.0-theta) - gamma*num/(1.0-theta);
    res += num;
  }

  return res;
}
double x(double p[], int pos)
{

    double num = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      num += (*AuxData::Cmn)[m][pos] * (p[m+1]);
    }

  return num;
}
void df(double p[], double g[])
{
  for (int k=0;k<AuxData::M;k++) {
    g[k+1] = 0.0;
    for (int n=0;n<AuxData::N;n++) {
      double num1 = 0;
      double num2 = 0;
      for (int m=0;m<AuxData::M;m++) {
        num1 += (*AuxData::Cmn)[m][n] * p[m+1];
        num2 += (*AuxData::CmnNorm)[m][n] * p[m+1];
      }
      g[k+1] -= ((*AuxData::Cmn)[k][n]*num2 - (*AuxData::CmnNorm)[k][n]*num1)
                                                                  / (num2*num2);
    }
  }
}
double f_double(double p[])
{
  return f(p-1);
/*  double res = 0.0;

  for (int n=0;n<AuxData::N;n++) {
    double num = 0.0;
    double den = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      num += (*AuxData::Cmn)[m][n] * (p[m]);
      den += (*AuxData::CmnNorm)[m][n] * p[m];
    }
    if (den!=0.0) res += num/den;
  }

  AuxData::vhfhist->ComoNumero(-res);

  return -res;*/
}
void df_double(double p[], double g[])
{
  df(p-1, g-1);
/*  for (int k=0;k<AuxData::M;k++) {
    g[k] = 0.0;
    for (int n=0;n<AuxData::N;n++) {
      double num1 = 0;
      double num2 = 0;
      for (int m=0;m<AuxData::M;m++) {
        num1 += (*AuxData::Cmn)[m][n] * p[m];
        num2 += (*AuxData::CmnNorm)[m][n] * p[m];
      }
      if (num2!=0.0) g[k] -= ((*AuxData::Cmn)[k][n]*num2 -
                                  (*AuxData::CmnNorm)[k][n]*num1) / (num2*num2);
    }
  }*/
}

double fexp(double p[])
{
  double res = 0.0;

  for (int n=0;n<AuxData::N;n++) {
    double num = 0.0;
    double den = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      double expc = exp((double)-p[m+1]*p[m+1]);
      num += (*AuxData::Cmn)[m][n] * expc;
      den += (*AuxData::CmnNorm)[m][n] * expc;
    }
    if (den!=0.0) res += num/den;
  }

  AuxData::vhfhist->SetNumero(-res);

  return -res;
}
void dfexp(double p[], double g[])
{
  for (int k=0;k<AuxData::M;k++) {
    g[k+1] = 0.0;
    for (int n=0;n<AuxData::N;n++) {
      double expc;
      double num1 = 0;
      double num2 = 0;
      for (int m=0;m<AuxData::M;m++) {
        expc = exp((double)-p[m+1]*p[m+1]);
        num1 += (*AuxData::Cmn)[m][n] * expc;
        num2 += (*AuxData::CmnNorm)[m][n] * expc;
      }
      expc = exp((double)-p[k+1]*p[k+1]);
      g[k+1] -= ((*AuxData::Cmn)[k][n]*(-2.0)*p[k+1]*expc*num2 -
                 (*AuxData::CmnNorm)[k][n]*(-2.0)*p[k+1]*expc*num1) / (num2*num2);
    }
  }
}
double fexp_d(double p[])
{
  double res = 0.0;

  for (int n=0;n<AuxData::N;n++) {
    double num = 0.0;
    double den = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      double expc = exp((double)-p[m]*p[m]);
      num += (*AuxData::Cmn)[m][n] * expc;
      den += (*AuxData::CmnNorm)[m][n] * expc;
    }
    if (den!=0.0) res += num/den;
  }

  AuxData::vhfhist->SetNumero(-res);

  return -res;
}
void dfexp_d(double p[], double g[])
{
  for (int k=0;k<AuxData::M;k++) {
    g[k] = 0.0;
    for (int n=0;n<AuxData::N;n++) {
      double expc;
      double num1 = 0;
      double num2 = 0;
      for (int m=0;m<AuxData::M;m++) {
        expc = exp((double)-p[m]*p[m]);
        num1 += (*AuxData::Cmn)[m][n] * expc;
        num2 += (*AuxData::CmnNorm)[m][n] * expc;
      }
      expc = exp((double)-p[k]*p[k]);
      g[k] -= ((*AuxData::Cmn)[k][n]*(-2.0)*p[k]*expc*num2 -
                 (*AuxData::CmnNorm)[k][n]*(-2.0)*p[k]*expc*num1) / (num2*num2);
    }
  }
}
double fexp2(double p[])
{
  double res = 0.0;

  for (int n=0;n<AuxData::N;n++) {
    double num = 0.0;
    double den = 0.0;
    for (int m=0;m<AuxData::M;m++) {
      double expc = 2.0*exp((double)-p[m+1]*p[m+1])-1.0;
      num += (*AuxData::Cmn)[m][n] * expc;
      den += (*AuxData::CmnNorm)[m][n] * expc;
    }
    if (den!=0.0) res += num/den;
  }

  AuxData::vhfhist->SetNumero(-res);
  
  return -res;
}
void dfexp2(double p[], double g[])
{
  for (int k=0;k<AuxData::M;k++) {
    g[k+1] = 0.0;
    for (int n=0;n<AuxData::N;n++) {
      double expc;
      double num1 = 0;
      double num2 = 0;
      for (int m=0;m<AuxData::M;m++) {
        expc = 2.0*exp((double)-p[m+1]*p[m+1])-1.0;
        num1 += (*AuxData::Cmn)[m][n] * expc;
        num2 += (*AuxData::CmnNorm)[m][n] * expc;
      }
      expc = exp((double)-p[k+1]*p[k+1]);
      g[k+1] -= ((*AuxData::Cmn)[k][n]*(-4.0)*p[k+1]*expc*num2 -
                 (*AuxData::CmnNorm)[k][n]*(-4.0)*p[k+1]*expc*num1) / (num2*num2);
    }
  }
}
//---------------------------------------------------------------------------
class Kk
{
  public:
  double f(double p[])
  {
    return 3+p[1];
  }
};
void DoOrdenaEnPruebas(Ensemble *ens, int ini, int fin, Mandato*m)
{

  int iter;
  double *p/*, *ph, *pH, *g*/, *u, *l/*, r1[101], r3[101], r4[101]*/;
//  double dp[101], dph[101], dpH[101], r2[101];
//  double e[102], eh[102], eH[102];
  int *nbd;
  double fret;
  double gtol;
  int M;

  Mandato *miter, *mfret, *mgtol;

  if (0==(mgtol=Variable::VariablePorNombre("gtol"))) {
    Mandato::valueOf("gtol=0.1");
    mgtol=Variable::VariablePorNombre("gtol");
  }

  gtol=mgtol->ComoNumero();


  M = ens->GetClassifiersToUse();

  //Calculala matriz Cmn
  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
  Matriz *Cmn = CreaMatrizCmn(ens, ini, fin);
  FILE *kr=fopen("c:\\ron.txt", "wt");
  for (int j=0;j<Cmn->columnas();j++) {//Sobre los datos
    for (int i=0;i<Cmn->filas();i++) {//Sobre los clasificadores
      fprintf(kr, "%f ", (*Cmn)[i][j]);
    }
    fprintf(kr, "%d\n", (int)ens->GetData()->GetDatClass(j));
  }
  fclose(kr);


  AuxData::init(Cmn);

  p = nr::vector(1, M);
//  ph = nr::vector(1, M);
//  pH = nr::vector(1, M);
//  g = nr::vector(1, M);
  u = nr::vector(1, M);
  l = nr::vector(1, M);
  nbd = new int[M];
  /*Mandato *var =*/ Variable::VariablePorNombre("v");
//  double v = var ? var->ComoNumero() : 2;
  srand(time(0));
  for(int i=1;i<=M;i++) {
//    p[i] = dp[i-1] = v;//(float)rand()/RAND_MAX;//v;//1.0/M;//exp(-(1.0/M)*(1.0/M));
  //  e[i]=dp[i-1] = sqrt(-log(p[i]));
    u[i] = 1.0;
    l[i] = 0.0;
    nbd[i-1] = 2;
  }

/*  int kk, iar=0;
  double psitrtr, psitstr, psitrts, psitsts;
  FILE*ff=fopen("c:\\ReErr_OrtrainMin.txt", "rt");
  var = Variable::VariablePorNombre("V");
  float V = var ? var->ComoNumero() : 0.0;
  while(1==fscanf(ff,"%d",&kk)) {
    p[kk+1] = dp[kk] = V;//(float)rand()/RAND_MAX;//V;
    e[kk+1] = dp[kk] = sqrt(-log(p[kk+1]));
    iar++;
    AuxData::fvaltr->ComoNumero(fexp(p));
  }
  fclose(ff);
  */

  //Medicion de la posicion de todos los datos de ejemplo para distintas
  //funciones de coste
  /*FILE*ffff=fopen("c:\\sal5.txt", "wt");
  fprintf(ffff, "ERR_TS\tERR_TR\tCOS\tCOS_EXP\tCOS_EXP2\tCOS_LOG\tCOS_PWS\n");
  double mod = 0.0;
  for(int j=0;j<M;j++) {
    p[j+1] = (float)1.0/M;
    ens->SetClassifierWeight(j, p[j+1]);
  }
  ens->SetUseWeights(true);
  double eeee = exp((double)1);
  for (int i=ini;i<=fin;i++) {
    double xx = x(p, i);
    double coste1 = pwsderivable(p, i, 10);//recta(xx, -eeee, 0.0);
    double coste2 = pwsderivable(p, i, 25);//dpwsderivable(p, i);//exp01(p, i, i, -1.0/eeee, 1.0/eeee);
    double coste3 = pwsderivable(p, i, 74);//exp02(p, i, i);
    double coste4 = 0;//log01(p, i, i, -1.0/eeee, 1.0/eeee);
    double coste5 = 0;//pws01(p, i, i, 0.1, 0.5/eeee);
    fprintf(ffff, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", xx, coste1, coste2, coste3,
                                                                coste4, coste5);
  }
  fclose(ffff);
  delete Cmn;
  return;
*/

/*  FILE*fff=fopen("c:\\sal7.txt", "wt");
  fprintf(fff, "ERR_TS\tERR_TR\tCOS\tCOS_EXP\tCOS_EXP2\tCOS_LOG\tCOS_PWS\n");


  //Medicion de distintas funciones de coste en N puntos (pesos) aleatorios
  double eeee = exp(1);
  for (int i=0;i<1000;i++) {
    double mod = 0.0;
    for(int j=0;j<M;j++) {
      p[j+1] = (float)rand()/RAND_MAX;
      mod += p[j+1];
    }
    for(int j=0;j<M;j++) {
      p[j+1] /= mod;
      ens->SetClassifierWeight(j, p[j+1]);
    }
    ens->SetUseWeights(true);
    double err_ts = 0;//ens->Error(fin+1, Total-1);
    double err_tr = ens->Error(ini, fin);
    double coste1 = frecta(p, ini, fin, -eeee, 0.0);
    double coste2 = exp01(p, ini, fin, -1.0/eeee, 1.0/eeee);
    double coste4 = log01(p, ini, fin, -1.0/eeee, 1.0/eeee);
    double coste5 = pws01(p, ini, fin, 0.1, 0.5/eeee);
    fprintf(fff, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", err_ts, err_tr, coste1,
                                                        coste2, coste4, coste5);
  }
  fclose(fff);

  delete Cmn;
  return;            */
  Cmn = CreaMatrizCmn(ens, ini, fin);

  AuxData::init(Cmn);
  FILE*fff=fopen("sal2.txt", "at");
  fprintf(fff, "ERROR_TEST\tCOSTE_1\terrtrainERROR_TEST\tCOSTE_1\terrtrain\n");
  gtol=0.0;
  fclose(fff);
  for (int i=0;i<10000;i++) {
    double mod = 0.0;
    for(int j=0;j<M;j++) {
      p[j+1] = ((int)((101.0*rand()/(RAND_MAX+1.0)))/100.0);
      mod += p[j+1];
    }
    for(int j=0;j<M;j++) {
      p[j+1] /= mod;
      ens->SetClassifierWeight(j, p[j+1]);
    }
    double coste1 =f(p);
    ens->SetUseWeights(true);
    double errortrain = ens->Error(0, 299);
    double error = ens->Error(300, 1299);
    fff=fopen("sal2.txt", "at");
    fprintf(fff, "%g\t%g\t%g\t", error, errortrain, coste1);
    fclose(fff);
    try {
    //Optimizacion::LBFGSB(p, M, gtol, &iter, &fret, f_double, df_double, l, u, nbd);
    }
    catch(...) {
    fff=fopen("sal2.txt", "at");
    fprintf(fff, "ERROR pesos: ");
    for(int j=0;j<M;j++) fprintf(fff, "%g", p[j+1]);
    fprintf(fff, "\n");
    fclose(fff);
    }
    for(int j=0;j<M;j++) {
      ens->SetClassifierWeight(j, p[j+1]);
    }
    coste1 = f(p);
    errortrain = ens->Error(0, 299);
    error = ens->Error(300, 1299);
    fff=fopen("sal2.txt", "at");
    fprintf(fff, "%d\t%g\t%g\t%g\t%g\n", iter, fret, error, errortrain, coste1);
    fclose(fff);
  }

return;

/*  iar=0;
  for(int i=1;i<=M;i++) p[i] = v;//1.0/M;//exp(-(1.0/M)*(1.0/M));
  ff=fopen("c:\\ReErr_Ortest.txt", "rt");
  while(1==fscanf(ff,"%d",&kk)) {
    p[kk+1] = V;
    iar++;
    AuxData::fvalts->ComoNumero(fexp(p));
  }
  fclose(ff);
  */
//  nr::dfpmin(p, M, gtol, &iter, &fret, fexp2, dfexp2);
  Optimizacion::LBFGSB(p, M, gtol, &iter, &fret, f_double, df_double, l, u, nbd);

  Mandato::valueOf("iter=1");
  Mandato::valueOf("fret=1");
  miter=Variable::VariablePorNombre("iter");
  mfret=Variable::VariablePorNombre("fret");
  ((Variable*)miter)->SetNumero(iter);
  ((Variable*)mfret)->SetNumero(fret);

  for(int j=0;j<M;j++){
    ens->SetClassifierWeight(j, exp((double)-p[j+1]*p[j+1]));
  }
  ens->OrdenarClasificadoresPorPeso();

  vector<int> oo = ens->DameOrdenOriginal();
  for(unsigned j=0;j<oo.size();j++){
    AuxData::orden->SetNumero(oo[j]);
    AuxData::pesos->SetNumero( ens->GetClassifierWeight(j) );
  }

  nr::free_vector(p,1,M);
  delete Cmn;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void OrdenOriginal::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());

  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("Recuperar orden inicial");
  CadenasValidas.push_back("Poner este orden como inicial");
  ParametroCadena *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(true, 0, "Orden inicial", "");
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
void OrdenOriginal::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new OrdenOriginal(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* OrdenOriginal::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();

  if (ParamInfo(1)->Asignado() && Param(1)->ComoNumero()==2)
    ens->PonerEsteOrdenComoOrdenOriginal();
  else
    ens->OrdenarClasificadoresPorOrdenOriginal();
  

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void ClasificadoresNecesarios::CreateParams()
{
  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("OOB");
  CadenasValidas.push_back("Ceros");
  Parametro *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(true, 0, "Orden inicial", "");
  Parametros.push_back(p2);

  p2 = new ParametroEnsemble();
  Parametros.push_back(p2);

  p2 = new ParametroData();
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
void ClasificadoresNecesarios::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new ClasificadoresNecesarios(), "Tools");
}
//---------------------------------------------------------------------------
int DeCeros(Ensemble *ens, Data *dat)
{
  ens->SetData(dat);
  vector<double> errors;

  ens->SecuencialError(0, dat->GetNTotal()-1, &errors);
  double min = errors[0];
  int imin = 0;
  for(int i=1;i<(int)errors.size();i++) {
    if (errors[i]<min) {
      min = errors[i];
      imin = i;
    }
  }

  int nClasf = ens->GetClassifiersToUse();
  int nDatos = dat->GetNTotal();
  int nClases = ((NomData*)dat)->NumClass;
  Matriz *mrg_sum = new Matriz(2, nClasf);
  Matriz cTrain(nDatos, nClases + 1, 0.0);
  vector<double> w = ens->GetWeights();
  int imax=-1;
  double max=-1;
  for(int i=0;i<nClasf;i++) {
    Classifier *c = ens->GetClassifier(i);
    (*mrg_sum)[0][i] = 0.0;
    for(int j=0;j<nDatos;j++) {
      int clase_ij = c->Classify(j);
      cTrain[j][clase_ij] += w[i];
      cTrain[j][nClases] += w[i];
    }
    if (i%2==1) {//no. de clasf. es par
      for(int j=0;j<nDatos;j++) {//Contamos los empates
        int clase_real = dat->GetDatClass(j);
        if (cTrain[j][clase_real]>(i+1)/2) continue;//no hay empates
        double maxmalo = 0.0;
        for(int k=0;k<nClases;k++) {
          if (k==clase_real) continue;
          if (cTrain[j][k]>maxmalo) 
            maxmalo = cTrain[j][k];
        }
        if (cTrain[j][clase_real] == maxmalo)
          (*mrg_sum)[0][i] += 1.0;
      }
    }
    else (*mrg_sum)[0][i] = i==0 ? 0.0 : (*mrg_sum)[0][i-1];
    //
    if (i==imin || (i>imin && max<=(*mrg_sum)[0][i])) {
      imax = i;
      max  = (*mrg_sum)[0][i];
    }
  }

  return imax+1;
}
Mandato* ClasificadoresNecesarios::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();

  string what   =              Param(0)->ComoCadena();
  Ensemble *ens = (Ensemble *) Param(1)->ComoDatos();
  Data *data    = (Data *)     Param(2)->ComoDatos();

  ens->SetData(data);

  if (what=="OOB") 
    DoClasificadoresNecesarios(ens, 0, data->GetNTotal()-1);
  else if (what=="Ceros") {
    posiciones.resize(1);
    posiciones[0] = DeCeros(ens, data);
  }

  return this;
}
//---------------------------------------------------------------------------
string ClasificadoresNecesarios::ComoCadena()
{
  string cad="";
/*  for(unsigned i=0;i<posiciones.size();i++)
    cad = cad + "\t" + Mandato::ACad(posiciones[i]);
  cad = cad + "\n\t";
  for(unsigned i=0;i<porcensref.size();i++)
    cad = cad + "\t" + Mandato::ACad(porcensref[i]);
    */
  for(unsigned i=0;i<estimerror.size();i++)
    cad = cad + "\t" + Mandato::ACad(estimerror[i]);
  return cad;
}
//---------------------------------------------------------------------------
void ClasificadoresNecesarios::DoClasificadoresNecesarios(Ensemble *ens, int ini, int fin)
{
  //Obtiene la matriz de clasificaciones excluyendo train (solo validacion)
  Matriz *Cmn = CreaMatrizCmn(ens, ini, fin, true);

  NomData *data = (NomData*)ens->GetData();
  int nDatos = (fin<0) ? data->GetNTotal() : fin-ini+1;
  int nClasf = ens->GetClassifiersToUse();

  int *acum = new int[nDatos];
  int *n_acum = new int[nDatos];

  int *mases = new int[nClasf];
  int *menoses = new int[nClasf];
  int *ceros = new int[nClasf];
  int *nulos = new int[nClasf];
  double *porcent = new double[nClasf];
  int maxmas = 0;
  posiciones.clear();
  for (unsigned i=0;i<porcensref.size()+1;i++)
      posiciones.push_back(-1);

  for(int i=0;i<nDatos;i++)
    acum[i] = n_acum[i] = 0;


  for(int j=0;j<nClasf;j++) {
    mases[j] = menoses[j] = ceros[j] = nulos[j] = 0;
    for(int i=0;i<nDatos;i++) {
      acum[i] += (int)(*Cmn)[j][i];
      n_acum[i] += (int)fabs((*Cmn)[j][i]);
      //
      if (n_acum[i]==0)   nulos[j]++;
      else if (acum[i]>0) mases[j]++;
      else if (acum[i]<0) menoses[j]++;
      else                ceros[j]++;
    }
    if (mases[j]>maxmas) {
      posiciones[0] = j;
      maxmas = mases[j];
    }
  }

  estimerror.clear();
  for(int j=0;j<nClasf;j++) {
    estimerror.push_back((double)mases[j]/(mases[j]+menoses[j]+ceros[j]));
//    if (estimerror<minerror) {minerror=estimerror; iminerror=j;}
    porcent[j] = (double)mases[j]/mases[nClasf-1];
    for (unsigned i=0;i<porcensref.size();i++) {
      if (porcent[j]+porcensref[i]>=1.0 && posiciones[i+1]==-1)
        posiciones[i+1] = j;
    }

  }

  delete Cmn;
  delete []porcent;
  delete []mases;
  delete []menoses;
  delete []ceros;
  delete []nulos;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void OrdenaPorContribucion::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
void OrdenaPorContribucion::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new OrdenaPorContribucion(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* OrdenaPorContribucion::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  ens->SetData(dat);

  int ini       = 0;
  int fin       = dat->GetNTotal() - 1;

  ens->SetUseWeights(false);

  DoOrdenaPorContribucion(ens, ini, fin);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void OrdenaEnPruebas::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
void OrdenaEnPruebas::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new OrdenaEnPruebas(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* OrdenaEnPruebas::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
//  Param(2)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();
//  int ini = 0;
  //int fin = dat->GetNTotal() - 1;
//  Data *dat2 = (Data *)Param(2)->ComoDatos();

/* ORDENACION ALEATORIA
  FILE *f = fopen("random_ord.csv", "w");
  for (int i=0;i<20;i++) {
    for(int j=0;j<ens->GetClassifiersToUse();j++) {
      ens->Exchange(j, j+RandomInteger(ens->GetClassifiersToUse()-1-j));
    }
    vector<double> errores;
    ens->SetData(dat);
    ens->SecuencialError(0, dat->GetNTotal()-1, &errores);
    for (unsigned k=0; k<errores.size(); k++) {
      fprintf(f, "%g\t", errores[k]);
    }
    fprintf(f, "\n");
    ens->SetData(dat2);
    ens->SecuencialError(0, dat2->GetNTotal()-1, &errores);
    for (unsigned k=0; k<errores.size(); k++) {
      fprintf(f, "%g\t", errores[k]);
    }
    fprintf(f, "\n");
  }
  fclose(f);*/
  
//  int ini       = Param(1)->ComoEntero();
 // int fin       = Param(2)->ComoEntero();

//  ens->SetUseWeights(false);
//if (ParamInfo(2)->Asignado())
//  iminmin = TestDistancias(ens, dat, true);
/*else*/ 
/*if (true) {
  Data *dts = (Data *)Param(2)->ComoDatos();
  if (res) delete res;
  res = new Matriz(2, ens->GetClassifiersToUse());
  ens->SetData(dat);
  Matriz *d = DoOrdenaPorReduccionErrorBackfitting(ens, dat, dts);
  (*res) = (*d);
  delete d;
}
else
{
  int nc = ens->GetClassifiersToUse();
  double t0 = clock();
  for (int i=0;i<10;i++) {
    for (int j=0;j<nc;j++)//Permutacion
      ens->Exchange(j, j + (int)(((double)(nc-j) * rand ()) / (RAND_MAX + 1.0)));
    iminmin = TestDistancias(ens, dat);
  }
  double t1 = clock();
  for (int i=0;i<10;i++) {
    for (int j=0;j<nc;j++)//Permutacion
      ens->Exchange(j, j + (int)(((double)(nc-j) * rand ()) / (RAND_MAX + 1.0)));
    //DoOrdenaPorReduccionDeDistancias(ens, (int)(0.5+0.075*nc), 0, dat->GetNTotal() - 1);
    DoOrdenaPorReduccionError(ens, 0, dat->GetNTotal() - 1);
  }
  double t2 = clock();
  printf("%d\t%g\t%g\n", nc, (t1-t0)/(10.0*CLOCKS_PER_SEC), (t2-t1)/(10.0*CLOCKS_PER_SEC));
}*/
//  DoOrdenaPorReduccionDeVal(ens, ini, fin);
//                double mastrain,double menostrain,double masval,double menosval)
//  DoOrdenaPorAlgunAcum(ens, ini, fin, 0, -1, 1, -1);
//  DoOrdenaPorError(ens, ini, fin, false);
  DoOrdenaPorEPIC(ens, (NomData*)dat);
//  DoOrdenaEnPruebas(ens, ini, fin, this);
//  DoOrdenaSimplex(ens, ini, fin);
//  DoOrdenaCV(ens, ini, fin);
//  DoNumClasfCV(ens, ini, fin, 3, iminmax, iminmin);
//  DoPseudoBoosting(ens, ini, fin);
//  DoDiversityOrdering(ens, (NomData *)dat);
//    DoOrdenaPorContribucion(ens, ini, fin);
  return this;
}
//---------------------------------------------------------------------------
string OrdenaEnPruebas::ComoCadena()
{
  string cad = "";
  return cad;
/*  std::ostringstream buf;
  res->saveToStream(buf);
  return buf.str();*/
  /*return Mandato::ACad(iminmin);

  string cad = "";
  cad = cad + "\t" + Mandato::ACad(iminmin);
  cad = cad + "\t" + Mandato::ACad(iminmax);
  return cad;*/
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void OrdenaPorReduccionError::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
void OrdenaPorReduccionError::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new OrdenaPorReduccionError(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* OrdenaPorReduccionError::Ejecutar()
{
  int ini, fin;

  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data*)Param(1)->ComoDatos();
  ens->SetData(dat);
  ini = 0;
  fin = dat->GetNTotal()-1;

  int start_at = 0;
  bool inverse = false;
  if (NumParamsAsignados()>=3) {
    Param(2)->Ejecutar();
    inverse = Param(2)->ComoBooleano();
  }

  num = DoOrdenaPorReduccionError(ens, ini, fin, start_at, inverse);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void OrdenaPorReduccionDeDistancias::CreateParams()
{

  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());

  Parametro *p = new ParametroNumero(0.0, 1.0);
  p->PonPropiedades(true, 0, "Distancia de referencia", "");
  p->SetPorOmision("0.075");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
void OrdenaPorReduccionDeDistancias::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new OrdenaPorReduccionDeDistancias(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* OrdenaPorReduccionDeDistancias::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();
  int dist      = (int)(0.5+Param(2)->ComoNumero()*ens->Count());

  ens->SetData(dat);

  int ini       = 0;
  int fin       = dat->GetNTotal() - 1;

//  ens->SetUseWeights(false);
//
  DoOrdenaPorReduccionDeDistancias(ens, dist, ini, fin);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void OrdenaPorAngulos::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
  Parametros.push_back(new ParametroBooleano());
}
//---------------------------------------------------------------------------
void OrdenaPorAngulos::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new OrdenaPorAngulos(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* OrdenaPorAngulos::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  bool UseOOB = false;
  if (NumParamsAsignados()>=3) {
    Param(2)->Ejecutar();
    UseOOB = Param(2)->ComoBooleano();
  }

  num = TestDistancias(ens, dat, false, UseOOB);
cout << num << endl;
  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void OrdenaPorEPIC::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
void OrdenaPorEPIC::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new OrdenaPorEPIC(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* OrdenaPorEPIC::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  NomData *dat = (NomData *)Param(1)->ComoDatos();

  DoOrdenaPorEPIC(ens, dat);

  num = 0;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void OrdenaPorUWA::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
void OrdenaPorUWA::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new OrdenaPorUWA(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* OrdenaPorUWA::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  DoOrdenaPorUWA(ens, (NomData*)dat);

  num = 0;

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void InvertirOrden::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
}
//---------------------------------------------------------------------------
void InvertirOrden::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new InvertirOrden(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* InvertirOrden::Ejecutar()
{
  Param(0)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();

  int nc = ens->GetClassifiersToUse();
  for(int j=0;j<nc/2;j++) {
    ens->Exchange(j, nc-j-1);
  }

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void OrdenaPorReduccionErrorPonderado::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
void OrdenaPorReduccionErrorPonderado::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new OrdenaPorReduccionErrorPonderado(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* OrdenaPorReduccionErrorPonderado::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
//  Param(2)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  ens->SetData(dat);

  int ini       = 0;
  int fin       = dat->GetNTotal() - 1;

//  ens->SetUseWeights(false);
//
  DoOrdenaPorReduccionErrorPonderado(ens, ini, fin);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void PseudoBoosting::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());

  Parametro *p = new ParametroBooleano();
  p->PonPropiedades(true, 0, "Remuestreo", "true/false");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
void PseudoBoosting::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new PseudoBoosting(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* PseudoBoosting::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
//  Param(2)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();
  bool Resample = ParamInfo(2)->Asignado() ? 
                                 Param(2)->Ejecutar()->ComoBooleano() : false;
  ens->SetData(dat);

  int ini       = 0;
  int fin       = dat->GetNTotal() - 1;

  ens->SetUseWeights(false);

  res = DoPseudoBoosting(ens, ini, fin, Resample);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void OrdenaPorReduccionErrorDoble::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());
  Parametros.push_back(new ParametroData());
}
//---------------------------------------------------------------------------
void OrdenaPorReduccionErrorDoble::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new OrdenaPorReduccionErrorDoble(), "Tools");
}
//---------------------------------------------------------------------------
Mandato* OrdenaPorReduccionErrorDoble::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
//  Param(2)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dat = (Data *)Param(1)->ComoDatos();

  ens->SetData(dat);

  int ini       = 0;
  int fin       = dat->GetNTotal() - 1;

  ens->SetUseWeights(false);

  DoOrdenaPorReduccionErrorDoble(ens, ini, fin);

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void SeleccionGenetica::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());

  Parametro *p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(false, 0, "Epocas", "");
  p->SetPorOmision("20");
  Parametros.push_back(p);

  p = new ParametroNumero(0.0, 1.0);
  p->PonPropiedades(false, 0, "Probabilidad de mutacion", "");
  p->SetPorOmision("0.01");
  Parametros.push_back(p);

  p = new ParametroData();
  p->PonPropiedades(false, 0, "Datos de ejemplo", "Ejemplos con los que se realizara la seleccion");
  Parametros.push_back(p);

  p = new ParametroEntero(1, 2100000000);
  p->PonPropiedades(true, 0, "No. clasf de la solucion inicial", "");
  Parametros.push_back(p);

  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("Error");
  CadenasValidas.push_back("ErrorSinBias");
  ParametroCadena *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(true, 0, "Funcion de coste", "");
  Parametros.push_back(p2);
}
//---------------------------------------------------------------------------
void SeleccionGenetica::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new SeleccionGenetica(), "Tools");
}
//---------------------------------------------------------------------------
class GAObjective
{
  public:
    virtual ~GAObjective(){}
    
  public:
    virtual int Calls()=0;
    virtual double operator()(GA *ag, Chromosome x, long length)=0;
};
class GAOError : public GAObjective
{
  protected:
    int **genMcd;
    int nc, nd;
    int nclasses;
    vector<double> weights;
    Data *dt;
    int nDraws;
    int nCalls;
    bool bigger_size_bias;
    vector<double> APrioriClas;

    
  public:
    GAOError(Ensemble *ens, Data *data, bool bigger_size_bias=true){
      data->SetNTrain(data->GetNTotal());
      genMcd = ens->MatClasif0(data);
      nd = data->GetNTotal();
      nc = ens->GetClassifiersToUse();
      nclasses = ((NomData*)data)->NumClass;
      dt = data;
      weights = ens->GetUsingWeights();
      nDraws = nCalls = 0;
      this->bigger_size_bias = bigger_size_bias;
      APrioriClas = ens->APrioriClas;
    }
    virtual ~GAOError(){
      for(int i=0;i<nc;i++) {
        delete []genMcd[i];
        delete []genMcd[i+nc];
      }
      delete []genMcd;
    }

    
  public:
    virtual int Calls(){return nCalls;}
    int Draws(){return nDraws;}
    
  public:
    virtual double operator()(GA *ag, Chromosome x, long length) {
      double error;
      int winners[1024];
 //     double drawpts = 1.0;
      double winpts = 1.0;
      bool morethanone = false;
      double no_clsf_factor = 1./(2.*nc);

      nCalls++;

      int no_clasf_excluidos = 0;
      for(int i=0;i<nc;i++) if (x[i]==0) no_clasf_excluidos++;

      //Esto hace que en caso de empate se seleccione el conjunto con
      // mas clasificadores
      error = bigger_size_bias ? no_clsf_factor*no_clasf_excluidos : 0.0;
     
      for(int j=0;j<nd;j++) {
        vector<double> dis(nclasses, 0.0);
        for(int i=0;i<nc;i++) {
          if (x[i]==1) {
            dis[genMcd[i][j]] += weights[i];
            morethanone = true;
          }
        }
        int cls = Classifier::WhichClass(dis, winners);//, &APrioriClas);
    //    if (winners[0]>1) { //MAs de una clase ganadora
      /*    int k;
          for(k=0;k<winners[0];k++) {*/
            //Buscamos si la buena estA entre las seleccionadas
/*            if (winners[k+1]==dt->GetDatClass(j)) {
              error += drawpts;
              nDraws++;
              break;
            }
          }
          if (k==winners[0]) error += winpts;
        }
        else*/ if (cls!=dt->GetDatClass(j)) error += winpts;
      }

      return morethanone ? error : nd;
    }
};
class GAODistancia : public GAObjective
{
  protected:
    int **genMcd;
    int nc, nd;
    int nclasses;
    vector<double> weights;
    Data *dt;
    int nDraws;
    int nCalls;

    
  public:
    GAODistancia(Ensemble *ens, Data *data){
      data->SetNTrain(data->GetNTotal());
      genMcd = ens->MatClasif0(data);
      nd = data->GetNTotal();
      nc = ens->Count();
      nclasses = ((NomData*)data)->NumClass;
      dt = data;
      weights = ens->GetUsingWeights();
      nDraws = nCalls = 0;
    }
    virtual ~GAODistancia(){
      for(int i=0;i<nc;i++) {
        delete []genMcd[i];
        delete []genMcd[i+nc];
      }
      delete []genMcd;
    }

    
  public:
    virtual int Calls(){return nCalls;}
    int Draws(){return nDraws;}

  public:
    virtual double operator()(GA *ag, Chromosome x, long length) {
      double error;
      int winners[1024];
 //     double drawpts = 1.0;
      double winpts = 1.0;
      bool morethanone = false;
      double no_clsf_factor = 1./(2.*nc);

      nCalls++;

      int no_clasf_excluidos = 0;
      for(int i=0;i<nc;i++) if (x[i]==0) no_clasf_excluidos++;

      //Esto hace que en caso de empate se seleccione el conjunto con 
      // mas clasificadores
      error = no_clsf_factor*no_clasf_excluidos;
     
      for(int j=0;j<nd;j++) {
        vector<double> dis(nclasses, 0.0);
        for(int i=0;i<nc;i++) {
          if (x[i]==1) {
            dis[genMcd[i][j]] += weights[i];
            morethanone = true;
          }
        }
        int cls = Classifier::WhichClass(dis, winners);
//        if (winners[0]>1) { //MAs de una clase ganadora
  //        int k;
    //      for(k=0;k<winners[0];k++) {
            //Buscamos si la buena estA entre las seleccionadas
      //      if (winners[k+1]==dt->GetDatClass(j)) {
        //      error += drawpts;
          //    nDraws++;
            //  break;
//            }
  //        }
    //      if (k==winners[0]) error += winpts;
      //  }
/*        else*/ if (cls!=dt->GetDatClass(j)) error += winpts;
      }

      return morethanone ? error : nd;
    }
};
/*double GA_objective                    (GA *ag, Chromosome x, long length)
{
  double error;
  int winners[1024];
  double drawpts = .5;
  double winpts = 1.0;
  bool morethanone = false;
  double no_clsf_factor = 1./(2.*nc);

  nCalls++;

  int no_clasf_excluidos = 0;
  for(int i=0;i<nc;i++) if (x[i]==0) no_clasf_excluidos++;

  //Esto hace que en caso de empate se seleccione el conjunto con 
  // mas clasificadores
  error = no_clsf_factor*no_clasf_excluidos;
 
  for(int j=0;j<nd;j++) {
    vector<double> dis(nclasses, 0.0);
    for(int i=0;i<nc;i++) {
      if (x[i]==1) {
        dis[genMcd[i][j]] += weights[i];
        morethanone = true;
      }
    }
    int cls = Classifier::WhichClass(dis, winners);
    if (winners[0]>1) { //MAs de una clase ganadora
      int k;
      for(k=0;k<winners[0];k++) {
        //Buscamos si la buena estA entre las seleccionadas
        if (winners[k+1]==dt->GetDatClass(j)) {
          error += drawpts;
          ndraws++;
          break;
        }
      }
      if (k==winners[0]) error += winpts;
    }
    else if (cls!=dt->GetDatClass(j)) error += winpts;
  }

  return morethanone ? error : nd;
}*/
GAObjective *gao = 0;
double GA_objective                    (GA *ag, Chromosome x, long length)
{
  return (*gao)(ag, x, length);
}
/*double GA_objective                    (GA *ag, Chromosome x, long length)
{
  double val = 0.0;

  nCalls++;

  int no_clasf_incluidos = 0;
  for(int i=0;i<nc;i++) if (x[i]==1) no_clasf_incluidos++;

  for(int j=0;j<nd;j++) {
    vector<double> dis(nclasses, 0.0);
    for(int i=0;i<nc;i++) {
      if (x[i]==1) {
        dis[genMcd[i][j]] += weights[i];
      }
    }
    val -= dis[dt->GetDatClass(j)];
    dis[dt->GetDatClass(j)] = 0.0;
    val += dis[Classifier::WhichClass(dis)];
  }
  return (no_clasf_incluidos==0 ? nd : val/(nd*no_clasf_incluidos)) +
                   (no_clasf_incluidos>10 ? 1.0 : 15.0-no_clasf_incluidos);
}*/
Mandato* SeleccionGenetica::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();
  Param(3)->Ejecutar();

  Ensemble *ens   = (Ensemble *)Param(0)->ComoDatos();
  int epochs      = Param(1)->ComoEntero();
  double prob_mut = Param(2)->ComoNumero();
  Data *data      = (Data *)Param(3)->ComoDatos();

  int nclasf = 1;
  if (NumParamsAsignados()>4) {
    // Nos viene dada una solucion inicial ma o meno guena
    Param(4)->Ejecutar();
    nclasf = Param(4)->ComoEntero();
  }

  int cost_fun    = 1;
  if (NumParamsAsignados()>5) {
    Param(5)->Ejecutar();
    string cf = Param(5)->ComoCadena();
    if (cf=="Error") cost_fun = 1;
    else if (cf=="ErrorSinBias") cost_fun = 2;
  }

  GAObjective *valGao = 0;
  if (NumParamsAsignados()>6) {
    Data *val;
    Param(6)->Ejecutar();
    val = (Data*)Param(6)->ComoDatos();
    val->SetNTrain(val->GetNTotal());
    valGao = new GAOError(ens, val);
  }

//  // Buscamos una solucion inicial ma o meno guena
//  int nclasf = DoOrdenaPorReduccionError(ens, 0, data->GetNTrain()-1);

  if (gao) delete gao;

  if (cost_fun==1)
    gao = new GAOError(ens, data);
  else if (cost_fun==2)
    gao = new GAOError(ens, data, false);

  int nc = ens->GetClassifiersToUse();
printf("\n--------------------->%d\n", nclasf);
  GA *g = GA_initialization(0, (nc<=20 ? 30 : nc) /*+ nclasf - 1*/, nc, 0.65, 
                                                                     prob_mut);
  g->elitism = 1;
  g->crossover_type = 0;
 
  if (nclasf>1) {//Adds solution given twice 
    //for (int i=1; i<nc; i++) {
    for (int j=0; j<g->length; j++) {
      g->old->ind [nclasf-1].chrom [j] = nclasf>j ? 1 : 0;
      g->old->ind [nclasf].chrom [j] = nclasf>j ? 1 : 0;
    }
//    }
  }

  /* Calculation of fitness */
  GA_calculate_population_fitness (g, g->old);
  /* Calculate best individual  */
  GA_calculate_best_individual (g);
  /* Variance  */
  GA_fitness_variance (g);
  /* Calculation of selection probabilities */
  GA_calculate_probabilities (g->old);
  /* Report results 
  GA_report (g);*/

  int i=0;
  ProgresoConsola pc;
  pc.Start(epochs);
  double fitnessbest = g->best_fitness;
//      printf("\n%12.10g", (1.0-g->best_fitness)/g->best_fitness);
  //    for (int j=0; j<g->length; j++)
     //   printf("%ld", g->old->ind [g->best].chrom [j]);
//  int last_size_increase = 0; 
//  int last_error_reduction = 0;
  do {
    i++;
    GA_generational_step(g);
//    FILE *f=fopen("sal.out", "at");
 //   double f1 = g->old->best_fitness;
    //double f2 =0.0;
    if (valGao)
      /*f2 = */(*valGao)(g, g->old->ind[g->best].chrom, nc);
//    fprintf(f, "%d\t%g\t%g\n", i, (1.0-f1)/f1, f2);
//    fclose(f);
    if (g->best_fitness>fitnessbest) {
/*      printf("\n%12.10g: ", (1.0-g->best_fitness)/g->best_fitness);
      for (int j=0; j<g->length; j++)
        printf("%ld", g->old->ind [0].chrom [j]);
      printf("\n");*/
      if (1.0/fitnessbest - 1.0/g->best_fitness < .9) {
//        printf("\b>");
  //      last_size_increase = i + 1;
      }
      else {
  //      printf("\b*");
//        last_error_reduction = i + 1;
      }
      fitnessbest = g->best_fitness;
    }
    pc();
  } while (i<epochs);
  pc.End();
  printf("ObjCalls: %d, Draws: %d \n", gao->Calls(), 
                                       ((GAOError*)gao)->Draws());

//  FILE *f=fopen("info.ga","a");
//  fprintf(f, "mut%g ini%d %d %d\n", prob_mut, nclasf, 
//                                     last_error_reduction, last_size_increase);
//  fclose(f);

  Chromosome best = g->old->ind[g->best].chrom;
  int nsel = 0;
  for (int i=0;i<nc;i++)
    if (best[i]==1)
      ens->Exchange(nsel++, i);

  ens->SetClassifiersToUse(nsel);
  res = nsel;

  if (valGao) delete valGao;

  GA_free(g);

  return this;
}
#ifdef COMP_SDP
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void SeleccionSDP::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());

  Parametro *p = new ParametroNumero(0.0, 1.0);
  p->PonPropiedades(false, 0, "Procentaje de poda.", "");
  p->SetPorOmision("0.20");
  Parametros.push_back(p);

  p = new ParametroData();
  p->PonPropiedades(false, 0, "Datos de ejemplo", "Ejemplos con los que se realizara la seleccion");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
void SeleccionSDP::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new SeleccionSDP(), "Tools");
}
Mandato* SeleccionSDP::Ejecutar()
{
	// We read the parameters given to the function SeleccionSDP
	
	Param(0)->Ejecutar();
	Param(1)->Ejecutar();
	Param(2)->Ejecutar();

	Ensemble *ens   = (Ensemble *)Param(0)->ComoDatos();
	double pruningRate = Param(1)->ComoNumero();
	Data *data      = (Data *)Param(2)->ComoDatos();
	int **Ccd; 

	// We compute the H matrix and the D matrix

	int nClasf = ens->GetClassifiersToUse();
	double **X, **H = (double**) malloc(sizeof(double*) * (nClasf + 1 ));
	for (int i = 0; i < (nClasf + 1); i++) H[ i ] = (double *) malloc(sizeof(double) * (nClasf + 1));
	double **D = (double**) malloc(sizeof(double*) * (nClasf + 1 ));
	for (int i = 0; i < (nClasf + 1); i++) D[ i ] = (double *) malloc(sizeof(double) * (nClasf + 1));

	int ini = 0;
	int fin = data->GetNTotal()-1;
	int nDatos = fin-ini+1;
	int *trueClass = new int[ nDatos ];

	for (int i = 0; i < nDatos; i++)
		trueClass[ i ] = data->GetDatClass(i);
	
	// We pass the whole data set through each classifier

	Ccd = ens->MatClasif0(data);

	// We compute the diagonal of the G matrix and the off diagonal terms

	for (int i = 0; i < nClasf; i++)
		H[ i + 1 ][ i + 1 ] = Classifier::CommonErrors(Ccd[i], Ccd[i], trueClass, nDatos);
	for (int i = 0; i < (nClasf - 1); i++)
		for (int j = (i + 1); j < nClasf; j++) {
			H[ i + 1 ][ j + 1 ] = Classifier::CommonErrors(Ccd[i], Ccd[j], trueClass, nDatos);
			H[ j + 1 ][ i + 1 ] = H[ i + 1 ][ j + 1 ];
		}

	// We compute the normalized ~G matrix

	for (int i = 0; i < (nClasf - 1); i++)
		for (int j = (i + 1); j < nClasf; j++) {
			H[ i + 1 ][ j + 1 ] = 
				0.5 * (( H[ j + 1 ][ j + 1 ] == 0 ? 0 : H[ j + 1 ][ i + 1 ] / H[ j + 1 ][ j + 1 ]) + 
					(H[ i + 1 ][ i + 1 ] == 0 ? 0 :  H[ i + 1 ][ j + 1 ] / H[ i + 1 ][ i + 1 ]));
			H[ j + 1 ][ i + 1 ] = H[ i + 1 ][ j + 1 ];
		}
	for (int i = 0; i < nClasf; i++)
		H[ i + 1 ][ i + 1 ] = H[ i + 1 ][ i + 1 ]  / nDatos;

	// We compute e~G terms in the H matrix

	for (int i = 0; i < nClasf; i++) {
		double sum = 0;
		for (int j = 0; j < nClasf; j++) sum += H[ i + 1 ][ j + 1 ];
		H[ 0 ][ i + 1 ] = sum;
		H[ i + 1 ][ 0 ] = H[ 0 ][ i + 1 ];
	}

	// We compute e^t~Ge terms in the H matrix

	H[ 0 ][ 0 ] = 0;
	for (int i = 0; i < nClasf; i++) H[ 0 ][ 0 ] += H[ 0 ][ i + 1 ];

	// We compute the D matrix

	for (int i = 0; i < nClasf + 1; i++) 
		for (int j = 0; j < nClasf + 1; j++) 
			D[ i ][ j ] = 0;			

	for (int i = 0; i < nClasf; i++) {
		D[ i + 1 ][ i + 1 ] = 1;
		D[ 0 ][ i + 1 ] = 1;
		D[ i + 1 ][ 0 ] = 1;
	}

	D[ 0 ][ 0 ] = nClasf;

	// We obtain the solution of the SDP opt problem called X

	if (solve_SDP(H,D, nClasf + 1, (int) (0.5 + pruningRate * nClasf)) != 1) {
		printf("Error while solving SDP problem.\n");
		exit(0);
	}

	X = D;

	// Now we perform the rounding process

	int *S, *Stmp;
	double vS, vStmp;
	S = approx_SDP(H, X, nClasf + 1, (int) (0.5 + pruningRate * nClasf), &vS);
	for (int i = 0 ; i < 10 ; i++) {
		
		Stmp = approx_SDP(H, X, nClasf + 1, (int) (0.5 + pruningRate * nClasf), &vStmp);

		if (vS > vStmp) {
			vS = vStmp;
			free(S);
			S = Stmp;
		} else {
			free(Stmp);
		}	
	}

	int nsel;
	nsel  = 0;
	for (int i = 0 ; i < nClasf ; i++)
		if (S[ i ] != 0)
			ens->Exchange(nsel++, i);

	free(S);

	printf("Numero de Clasificadores en conjunto final: %d.\n", nsel);

	ens->SetClassifiersToUse(nsel);
	res = nsel;

	// We free all the memory 

	for (int i = 0; i < (nClasf + 1); i++) {
		free(H[ i ]); free(D[ i ]);
	}

	free(D); free(H);

	for(int i=0;i<nClasf;i++) {
		delete []Ccd[i];
		delete []Ccd[i+nClasf];
	}

	delete []Ccd;
	delete []trueClass;

	return this;
}

// Function that solves the SDP optimization problem by using the CSPD library

int *approx_SDP(double **H, double **X, int dim, int k, double *value) {

double unif1, unif2, minValue;
int *S = (int*) malloc(sizeof(int) * (dim - 1));
double *U = new double[ dim ];
int *x = new int[ dim ];
int min;

	// We generate a gaussian random variable N(0, X)
	
	for (int i = 0 ; i < dim; i++) {
		unif1 = 1.0 * rand () / (RAND_MAX + 1.0);
		unif2 = 1.0 * rand () / (RAND_MAX + 1.0);
		U[ i ] = cos(2 * 3.1415 * unif1) * sqrt(2 * log(1 / unif2));
	}


	for (int i = 0 ; i < dim; i++) {
		double sum = 0;
		for (int j = 0 ; j < dim ; j++) 
            sum += X[ j ][ i ] * U[ j ];
		x[ i ] = sum >= 0 ? 1 : -1;
	}

	//  We set to one those values equal to x[ 0 ]

	for (int i = 1 ; i < dim; i++)
		if (x[ i ] == x[ 0 ])
			S[ i - 1 ] = 1;
		else 
			S[ i - 1 ] = 0;

	int k_tmp = 0;
	for (int i = 0 ; i < dim - 1 ; i++) 
		if (S[ i ] != 0) k_tmp ++;

//	printf("Ensemble size before geedy algorithm: %d. Target size: %d.\n", k_tmp, k);

	// If k_tmp == k we finish

	if (k_tmp == k) {
		delete []U;
		delete []x;

		double sum = 0;
		for (int i = 0 ; i < dim - 1 ; i++)
			if (S[ i ] != 0)
				for (int j = 0 ; j < dim -1 ; j++)
					sum += S[ j ] * H[ j + 1 ][ i + 1 ];
		*value = sum;
		return S;
	}

	double sum;
		
	// Now we complete the k terms by a greedy search procedure

	do {
		k_tmp = 0;
		for (int i = 0 ; i < dim - 1 ; i++) 
			if (S[ i ] != 0) k_tmp ++;

		if (k_tmp == k) break;

		min = 0; 
		minValue = 1e300;

		// We check wherther we have to remove or add a classifier

		for (int l = 0 ; l < dim -1 ; l++) {

			// We set a new S vector
			
			if (k_tmp < k) {
				if (S[ l ] == 1)
					continue;
				else
					S[ l ] = 1;
			} else
				if (S[ l ] == 0)
					continue;
				else
					S[ l ] = 0;

			// We check the viability of the new S

			sum = 0;
			for (int i = 0 ; i < dim - 1 ; i++)
				if (S[ i ] != 0)
					for (int j = 0 ; j < dim -1 ; j++)
						sum += S[ j ] * H[ j + 1 ][ i + 1 ];

			// We restore the S vector and update the best candidate

			if (k_tmp < k)
				S[ l ] = 0;
			 else 
				S[ l ] = 1;

			if (minValue > sum)  {
				min = l;
				minValue = sum;
				*value = sum;
			}
		}

		if (k_tmp < k) 
			S[ min ] = 1;
		else 
			S[ min ] = 0;

	} while (k_tmp != k);

	k_tmp = 0;
	for (int i = 0 ; i < dim - 1 ; i++) 
		if (S[ i ] != 0) k_tmp ++;

//	printf("Reached size: %d.\n", k_tmp);

	delete []U;
	delete []x;
	return S;
}

// Function that solves the SDP optimization problem by using the CSPD library

int solve_SDP(double **H, double **D, int dim, int k) {

// The problem and solution data

struct blockmatrix C;
double *b;
struct constraintmatrix *constraints;

// Storage for the initial and final solutions.

struct blockmatrix X,Z;
double *y;
double pobj,dobj;

// blockptr will be used to point to blocks in constraint matrices.

struct sparseblock *blockptr;

// A return code for the call to easy_sdp().

int ret;

  /*
   * The first major task is to setup the C matrix and right hand side b.
   */

 /*
   * First, allocate storage for the C matrix.  We have three blocks, but
   * because C starts arrays with index 0, we have to allocate space for
   * four blocks- we'll waste the 0th block.  Notice that we check to
   * make sure that the malloc succeeded.
   */

	C.nblocks=1;
	C.blocks=(struct blockrec *)malloc(2 * sizeof(struct blockrec));
	if (C.blocks == NULL) {
		printf("Couldn't allocate storage for C!\n");
		return -1;
	}

  /*
   * Setup the first block.
   */

	C.blocks[1].blockcategory=MATRIX;
	C.blocks[1].blocksize=dim;
	C.blocks[1].data.mat=(double *)malloc(dim * dim * sizeof(double));
	if (C.blocks[1].data.mat == NULL) {
		printf("Couldn't allocate storage for C!\n");
		return -1;
	}


  /*
   * Put the entries into the first block. We exchange signs because we are minimizing
   */
	
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			C.blocks[1].data.mat[ ijtok(i + 1 ,j + 1, dim) ]= -H[ i ][ j ];

  /*
   * Allocate storage for the right hand side, b.
   */

	b = (double *) malloc((dim + 2) * sizeof(double));
	if (b==NULL) {
		printf("Failed to allocate storage for a!\n");
		return -1;
	}

	/*
	* Fill in the entries in b
	*/

	for (int j = 2; j < dim + 2; j++) b[ j ] = 1;
	b[ 1 ] = 4 * k;

  /*
   * The next major step is to setup the two constraint matrices A1 ... Adim.
   * Again, because C indexing starts with 0, we have to allocate space for
   * one more constraint.  constraints[0] is not used.
   */

	constraints=(struct constraintmatrix *)malloc((dim+1+1)*sizeof(struct constraintmatrix));
	if (constraints==NULL) {
		printf("Failed to allocate storage for constraints!\n");
		return -1;
	}

  /*
   * Setup the A1 matrix.  Note that we start with block 1 of A1 and then
   * do block 1 of A1.  We do this in this order because the blocks will
   * be inserted into the linked list of A1 blocks in reverse order.
   */

  /*
   * Terminate the linked list with a NULL pointer.
   */

	constraints[1].blocks=NULL;

   // We process each constrain. The first one is tr(HX) = b[ 1 ]

	/*
	 * Allocate space for block 1.
	 */

	blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
	if (blockptr==NULL) {
		printf("Allocation of constraint block failed!\n");
		return -1;
	}

	/*
	 * Initialize block 1.
	 */

	blockptr->blocknum=1;
	blockptr->blocksize=dim;
	blockptr->constraintnum=1;
	blockptr->next=NULL;
	blockptr->nextbyblock=NULL;
	blockptr->entries=(double *) malloc((dim * (dim - 1) / 2 + dim + 1)*sizeof(double));

	if (blockptr->entries==NULL) {
		printf("Allocation of constraint block failed!\n");
		return -1;
	}

	blockptr->iindices=(int *) malloc((dim * (dim - 1) / 2 + dim + 1)*sizeof(int));
	if (blockptr->iindices==NULL) {
		printf("Allocation of constraint block failed!\n");
		return -1;
	}

	blockptr->jindices=(int *) malloc((dim * (dim - 1) / 2 + dim + 1)*sizeof(int));
	if (blockptr->jindices==NULL) {
		printf("Allocation of constraint block failed!\n");
		return -1;
	}
	
	/*
	 * We have  dim * (dim - 1) / 2 nonzero entries in the upper triangle of block 1 of A0.
	 */

	blockptr->numentries= dim * (dim - 1) / 2 + dim;

	/*
	 * We set the entries of the block that are given by D
	 */

	int n = 1;
	for (int i = 0 ; i < dim ; i++)
		for (int j = i ; j < dim ; j++) {
			blockptr->iindices[ n ]= i + 1;
			blockptr->jindices[ n ]= j + 1;
			blockptr->entries[ n ]= D[ i ][ j ];
			n++;
		}

	/*
	 * Insert block 1 into the linked list of A1 blocks.
	 */

	blockptr->next=constraints[ 1 ].blocks;
	constraints[ 1 ].blocks=blockptr;

	// We set the constrains diag(V) = e

	for (int i = 1 ; i < dim  + 1; i++) {

		constraints[ i + 1 ].blocks=NULL;

		/*
		 * Allocate space for block 1
		 */
		
		blockptr=(struct sparseblock *)malloc(sizeof(struct sparseblock));
		if (blockptr==NULL) {
			printf("Allocation of constraint block failed!\n");
			return -1;
		}

		/*
		 * Initialize block 1.
		 */

		blockptr->blocknum=1;
		blockptr->blocksize=dim;
		blockptr->constraintnum=i + 1;
		blockptr->next=NULL;
		blockptr->nextbyblock=NULL;
		blockptr->entries=(double *) malloc((1+1)*sizeof(double));

		if (blockptr->entries==NULL) {
			printf("Allocation of constraint block failed!\n");
			return -1;
		}

		blockptr->iindices=(int *) malloc((1+1)*sizeof(int));
		if (blockptr->iindices==NULL) {
			printf("Allocation of constraint block failed!\n");
			return -1;
		}

		blockptr->jindices=(int *) malloc((1+1)*sizeof(int));
		if (blockptr->jindices==NULL) {
			printf("Allocation of constraint block failed!\n");
			return -1;
		}
		
		/*
		 * We have 1 nonzero entry in the upper triangle of block 3 of A1.
		 */

		blockptr->numentries=1;

		/*
		 * The entry in the 1,1 position of block 1 of A1 to Adim is 1.0
		 */

		blockptr->iindices[1]=i;
		blockptr->jindices[1]=i;
		blockptr->entries[1]=1.0;

		/*
		 * Note that the entry in the 2,2 position of block 3 of A1 is 0,
		 * So we don't store anything for it.
		 */

		/*
		 * Insert block 1 into the linked list of A1 blocks.
		 */

		blockptr->next=constraints[ i + 1 ].blocks;
		constraints[ i + 1 ].blocks=blockptr;
	}

	/*
	 * Write the problem out in SDPA sparse format.
	 */

	write_prob((char*)"prob.dat-s", dim, dim + 1, C, b, constraints);

	/*
	 * Create an initial solution.  This allocates space for X, y, and Z,
	 * and sets initial values.
	 */

	initsoln(dim, dim + 1,C,b,constraints,&X,&y,&Z);

	/*
	 * Solve the problem.
	 */
 
	ret = easy_sdp(dim, dim + 1, C, b, constraints, 0.0, &X, &y, &Z, &pobj, &dobj);

    chol(X);

	if (ret == 0)
		printf("The objective value is %.7e \n",(dobj + pobj) / 2);
	else
		printf("SDP failed.\n");

	/*
	 * Write out the problem solution.
	 */

	write_sol((char*)"prob.sol",dim , dim + 1, X, y, Z);

	// We set H to the solution X

	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			D[ i ][ j ] = X.blocks[1].data.mat[ ijtok(i + 1 ,j + 1, dim) ];

	/*
	 * Free storage allocated for the problem and return.
	 */

	free_prob(dim, dim + 1, C, b, constraints, X, y, Z);

	return 1;
}
#endif
//---------------------------------------------------------------------------
void AnalisisDeOrden::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());

  ParametroData *pd = new ParametroData();
  pd->PonPropiedades(false, 0, "Datos de ordenacion", "Ejemplos con los que se realizara la ordenacion");
  Parametros.push_back(pd);

  pd = new ParametroData();
  pd->PonPropiedades(false, 0, "Datos de test", "Datos para la comprobacion de resultados");
  Parametros.push_back(pd);

  std::vector<std::string> CadenasValidas;
  CadenasValidas.push_back("OrdenaPorReduccionError");
  CadenasValidas.push_back("OrdenaPorReduccionDeDistancias");
  CadenasValidas.push_back("OrdenaPorReduccionErrorPonderado");
  CadenasValidas.push_back("PseudoBoosting");
  ParametroCadena *p2 = new ParametroCadena(&CadenasValidas);
  p2->PonPropiedades(false, 0, "Ordenar por metodo", "Metodo de ordenacion a utilizar");
  Parametros.push_back(p2);

  Parametro *pfs = new ParametroFicheroSalida();
  pfs->PonPropiedades(false, 0, "Fichero de salida",
                "Nombre del fichero para guardar los resultados");
  Parametros.push_back(pfs);

//  Parametro *p = new ParametroNumero(0.0, 1.0);
//  p->PonPropiedades(true, 0, "Distancia de referencia", "solo necesario para el metodo OrdenaPorReduccionDeDistancias");
//  p->SetPorOmision("0.075");
//  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
void AnalisisDeOrden::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new AnalisisDeOrden(), "Tools");
}
//---------------------------------------------------------------------------
void AFichero(string out, Ensemble *ens, Data *dtr, Data *dts, string id)
{
  FILE *f = fopen(out.c_str(), "at");
  if (!f) return;
  vector<double> errs;

  ens->SetData(dtr);
  ens->SecuencialError(0, dtr->GetNTotal()-1, &errs);
  fprintf(f, "%05lu_tr_%s:", (long)errs.size(), id.c_str());
  for(unsigned i=0;i<errs.size();i++) {
    fprintf(f, "\t%g", errs[i]);
  }
  fprintf(f, "\n");

  errs.clear();
  ens->SetData(dts);
  ens->SecuencialError(0, dts->GetNTotal()-1, &errs);
  fprintf(f, "%05lu_ts_%s:", (long)errs.size(), id.c_str());
  for(unsigned i=0;i<errs.size();i++) {
    fprintf(f, "\t%g", errs[i]);
  }
  fprintf(f, "\n");

  fclose(f);
}
Mandato* AnalisisDeOrden::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();
  Param(3)->Ejecutar();
  Param(4)->Ejecutar();

  Ensemble *ens = (Ensemble *) Param(0)->ComoDatos();
  Data *dor     = (Data*) Param(1)->ComoDatos();
  Data *dts     = (Data*) Param(2)->ComoDatos();
  string alg    = Param(3)->ComoCadena();
  string out    = Param(4)->ComoCadena();

  int inipar = 5;
  double dist = 0.0;
  string id;
  if (alg == "OrdenaPorReduccionDeDistancias") {
    Param(5)->Ejecutar();
    dist = Param(5)->ComoNumero();
    inipar++;
    id = "d";
  }
  else if (alg == "OrdenaPorReduccionError")
    id = "re";
  else if (alg == "OrdenaPorReduccionErrorPonderado")
    id = "rp";
  else if (alg == "PseudoBoosting")
    id = "pb";
  else if (alg == "OrdenaPorAngulos")
    id = "dd";

  vector<int> tam;
  tam.push_back(1);
  for(int i=inipar;i<NumParamsAsignados();i++)
    tam.push_back(Param(i)->Ejecutar()->ComoEntero());
  tam.push_back(ens->Count());

//  ens->SetUseWeights(false);
//
  int ini = 0;
  int fin = dor->GetNTotal()-1;
  char sep = ':';
  printf("Estamos en");
  for(unsigned i=0;i<tam.size();i++) {
    if (tam[i]>=ens->Count() && (i+1)!=tam.size()) continue;
    ens->SetData(dor);
    ens->OrdenarClasificadoresPorOrdenOriginal();
    ens->SetClassifiersToUse(tam[i]);
    printf("%c %d", sep, tam[i]);
    sep = i==tam.size()-2 ? 'y' : ',';
    fflush(stdout);
    if (alg == "OrdenaPorReduccionDeDistancias")
      DoOrdenaPorReduccionDeDistancias(ens, (int)(0.5 + dist*tam[i]), ini, fin);
    else if (alg == "OrdenaPorReduccionError")
      DoOrdenaPorReduccionError(ens, ini, fin);
    else if (alg == "OrdenaPorReduccionErrorPonderado")
      DoOrdenaPorReduccionErrorPonderado(ens, ini, fin);
    else if (alg == "PseudoBoosting")
      DoPseudoBoosting(ens, ini, fin);
    else if (alg == "OrdenaPorAngulos")
      TestDistancias(ens, dor);

    //Guardamos en fichero
    AFichero(out, ens, dor, dts, id);
  }
  printf("\n");

  return this;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
SeleccionCombinatoria::~SeleccionCombinatoria()
{
  if(res) delete res;
}
void SeleccionCombinatoria::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new SeleccionCombinatoria(), "Tools");
}
//---------------------------------------------------------------------------
void SeleccionCombinatoria::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());

  Parametro *p = new ParametroData();
  p->PonPropiedades(false, 0, "Datos de train", 
              "A partir de estos datos se obtienen las mejores combinaciones");
  Parametros.push_back(p);

  p = new ParametroData();
  p->PonPropiedades(true, 0, "Datos de test",
              "Datos para comprobar el resultado (NO se usan en la búsqueda)");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
int DoSelectBest(Ensemble *ens, Data *dtr, int rango, int start_at)
{
  ens->OrdenarClasificadoresPorOrdenOriginal();

  int ini = 0;
  int nt = dtr->GetNTotal();
  int fin = nt-1;
  vector<double> errores;
  ens->SecuencialError(ini, fin, &errores);

  ens->SetClassifiersToUse(rango+start_at);

  BestSolutions kk(ens, dtr, start_at);
  kk.run();
  Matriz *d = kk.GetBests();

//  if (res) delete res;
//  res = new Matriz(d->filas(), d->columnas() + (dts ? 2 : 1) );
//  (*res) = (*d);

  int nClasfs = ens->GetClassifiersToUse();
  double max_dif = 0.0;
  int i_max_dif = -1;
  for(int i=0;i<nClasfs;i++) {
    if ((i+1)%2==0) continue;
    ens->OrdenarClasificadoresPorOrdenOriginal();
    int k = 0;
    for(int j=0;j<ens->GetClassifiersToUse();j++) {
      if ((*d)[i+1][j]==1) ens->Exchange(k++, j);
    }
    ens->SetClassifiersToUse(i+1);
    double ee = ens->Error(dtr);//(*res)[i+1][d->columnas()] = ens->Error(dtr);
    //if (errores[i]-(*res)[i+1][d->columnas()]>=max_dif) {
    if (errores[i]-ee>=max_dif) {
      max_dif = errores[i]-ee;
      i_max_dif = i;
      cout << max_dif << " " << i << "  ";
    }
//    if (dts) (*res)[i+1][d->columnas()+1] = ens->Error(dts);
    ens->SetClassifiersToUse(nClasfs);
  }

  //Deja el orden con mas dif.
  ens->OrdenarClasificadoresPorOrdenOriginal();
  int k = 0;
  for(int j=0;j<ens->GetClassifiersToUse();j++) {
    if ((*d)[i_max_dif+1][j]==1) ens->Exchange(k++, j);
  }

  ens->SetClassifiersToUse(ens->Count());
  return i_max_dif+1;
}
Mandato* SeleccionCombinatoria::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dtr = (Data *)Param(1)->ComoDatos();
//  Data *dts = ParamInfo(2)->Asignado() ? (Data *)Param(2)->ComoDatos() : 0;

  ens->SetUseWeights(false);
  ens->OrdenarClasificadoresPorOrdenOriginal();

  ens->SetData(dtr);

  int ini = 0;
  int nt = dtr->GetNTotal();
  int fin = nt-1;
  int start_at = 0;
  int rango = 15;
  do {
    DoOrdenaPorReduccionError(ens, ini, fin, start_at);
    ens->PonerEsteOrdenComoOrdenOriginal();
    start_at = DoSelectBest(ens, dtr, rango, start_at);
    ens->PonerEsteOrdenComoOrdenOriginal();
  } while(start_at<ens->Count()-1);


  return this;
}
//---------------------------------------------------------------------------
string SeleccionCombinatoria::ComoCadena()
{
  string cad="";
  for(int i=0;i<res->columnas();i++)
    cad = cad + "\t" + Mandato::ACad((*res)[0][i]);
  return cad;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Mejores::~Mejores()
{
  if(res) delete res;
}
void Mejores::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new Mejores(), "Tools");
}
//---------------------------------------------------------------------------
void Mejores::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());

  Parametro *p = new ParametroData();
  p->PonPropiedades(false, 0, "Datos de train", 
              "A partir de estos datos se obtienen las mejores combinaciones");
  Parametros.push_back(p);

  p = new ParametroData();
  p->PonPropiedades(true, 0, "Datos de test",
              "Datos para comprobar el resultado (NO se usan en la búsqueda)");
  Parametros.push_back(p);
}
Mandato* Mejores::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *dtr = (Data *)Param(1)->ComoDatos();
  Data *dts = ParamInfo(2)->Asignado() ? (Data*)Param(2)->ComoDatos() : 0;

  ens->SetUseWeights(false);
  ens->SetData(dtr);
  int nClasfs = ens->GetClassifiersToUse();

  ens->OrdenarClasificadoresPorOrdenOriginal();

  BestSolutions kk(ens, dtr, 0);
  kk.run();

  Matriz *d = kk.GetBests();

  if (res) delete res;
  res = new Matriz(d->filas(), d->columnas() + (dts ? 2 : 1));
  (*res) = (*d);

  for(int i=0;i<nClasfs;i++) {
    if ((i+1)%2==0) continue;
    ens->OrdenarClasificadoresPorOrdenOriginal();
    int k = 0;
    for(int j=0;j<ens->GetClassifiersToUse();j++) {
      if ((*d)[i+1][j]==1) ens->Exchange(k++, j);
    }
    ens->SetClassifiersToUse(i+1);
    (*res)[i+1][d->columnas()] = ens->Error(dtr);
    if (dts) (*res)[i+1][d->columnas()+1] = ens->Error(dts);
    ens->SetClassifiersToUse(nClasfs);
  }

  FILE *ftr = fopen("train_bests.txt", "a");
  FILE *fts = fopen("test_bests.txt", "a");
  for(int i=0;i<nClasfs;i++) {
    double error_tr, error_ts;
    if ((i+1)%2==0) {
      error_tr = 0.5*((*res)[i][d->columnas()]+(*res)[i+2][d->columnas()]);
      error_ts = 0.5*((*res)[i][d->columnas()+1]+(*res)[i+2][d->columnas()+1]);
    }
    else {
      error_tr = (*res)[i+1][d->columnas()];
      error_ts = (*res)[i+1][d->columnas()+1];
    }
    fprintf(ftr, "%g\t", error_tr);
    fprintf(fts, "%g\t", error_ts);
  }
  fprintf(ftr, "\n");
  fprintf(fts, "\n");
  fclose(ftr);
  fclose(fts);

  vector<vector<int> > bb = kk.GetAllBests();

  for(int i=0;i<nClasfs;i+=2) {
    FILE *fal = fopen("all_bests.txt", "a");
    for(unsigned j=0;j<bb[i].size();j++) {
      int sol =  bb[i][j];
fprintf(fal, "S%d ", (i+1));
for(int jjj=0;jjj<nClasfs;jjj++){
int kk = sol & (1 << jjj) ? 1 : 0;
fprintf(fal, "%d ", kk);
}
fprintf(fal, "\n");
    }
    fclose(fal);
  }


  return this;
}
//---------------------------------------------------------------------------
string Mejores::ComoCadena()
{
  std::ostringstream buf;
  res->saveToStream(buf);
  return buf.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void OrdenaPorGGreedy::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());

  Parametro *p = new ParametroData();
  p->PonPropiedades(false, 0, "Datos de ejemplo", "Ejemplos con los que se realizara la seleccion");
  Parametros.push_back(p);

  p = new ParametroNumero(0.0, 500.0);
  p->PonPropiedades(true, 0, "Parar en (Tamano del subconjunto).", "Por omision: se ordena hasta el final");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
void OrdenaPorGGreedy::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new OrdenaPorGGreedy(), "Tools");
}
//---------------------------------------------------------------------------
// We can do an optimization of the following method.
// if we compute the G value based on the previous evaluations
// we only need to sum the new column and row 
/*
Mandato* OrdenaPorGGreedy::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *data = (Data *)Param(1)->ComoDatos();

double start = clock();

  int nClasif = ens->GetClassifiersToUse();

  // If the last parameter is not specified,
  // then the sort is performed until the ensemble is completed
  int stop = nClasif;
  if (NumParamsAsignados() > 2) {
     Param(2)->Ejecutar();
     stop = Param(2)->ComoEntero();
  }

  double** gMatrix = ens->computeGMatrix(data);
  
  // Inits the buffer where the selected clasifiers will
  // be stored at each iteration.
  int* selectedClassifiers = new int[nClasif];
  for (int i = 0; i < nClasif; ++i) {
    selectedClassifiers[i] = i;
  }

  // The array will contain the selected
  // classifiers until the i-th position,
  // and will put the new trials in the i-th position.
  for (int i = 0; i < stop; ++i) {
    double minG = ens->evaluateGMatrix(gMatrix, nClasif, selectedClassifiers, i+1);
    int bestClasif = i; 

    cout << "Trying " << selectedClassifiers[i] << "(" << minG << ")" <<endl;   

    // This loop searches for the next best clasifier
    // Doesn't need to test the i classifier again:
    for (int j = i+1; j < nClasif; ++j) {
      // Try the next not-selected classifier.
      FnsOrdClasUtils::swap(selectedClassifiers, i, j);

      // Evaluate a G natrix of size 'nClasif' and 'i+1' selectedClassifiers
      double G = ens->evaluateGMatrix(gMatrix, nClasif, selectedClassifiers, i+1);
      cout << "Trying " << selectedClassifiers[i] << "(" << G << ")";
      if (G < minG) {
         minG = G;
         bestClasif = j;
         cout << " selected!";
      }
      cout << endl;

      // Undo the trial
      FnsOrdClasUtils::swap(selectedClassifiers, i, j);
    }

    // Selects the best classifier, that is:
    // swaps the bestClasif with the clasifiers which is in the i-th position:
    FnsOrdClasUtils::swap(selectedClassifiers[bestClasif], selectedClassifiers[i]);

    ens->Exchange(i, bestClasif);
    cout << "Current subset: "<< endl;
    FnsOrdClasUtils::print(selectedClassifiers, i+1, cout);
  }

  ens->SetClassifiersToUse(stop);

double end = clock();

  // This block to open the file to store the number of evaluations
  string filename2("time_greedy.txt");
  ofstream os2(filename2.c_str(), ios::out | ios::app);
  os2 << (end - start) / CLOCKS_PER_SEC << endl;
  os2.close();


  delete[] selectedClassifiers;
  return this;
}
*/


Mandato* OrdenaPorGGreedy::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  Data *data = (Data *)Param(1)->ComoDatos();

double start = clock();

  int nClasif = ens->GetClassifiersToUse();

  // If the last parameter is not specified,
  // then the sort is performed until the ensemble is completed
  int stop = nClasif;
  if (NumParamsAsignados() > 2) {
     Param(2)->Ejecutar();
     stop = Param(2)->ComoEntero();
  }

  double** gMatrix = ens->computeGMatrix(data);
  
  // Inits the buffer where the selected clasifiers will
  // be stored at each iteration.
  int* selectedClassifiers = new int[nClasif];
  for (int i = 0; i < nClasif; ++i) {
    selectedClassifiers[i] = i;
  }

  // The array will contain the selected
  // classifiers until the i-th position,
  // and will put the new trials in the i-th position.
  double previousG = 0;
  for (int i = 0; i < stop; ++i) {
    double minOverhead = computeRow(gMatrix, nClasif, selectedClassifiers, i, selectedClassifiers[i]);
    double minG = previousG + minOverhead;
    int bestClasif = i; 

//    cout << "Trying " << selectedClassifiers[i] << "(" << minG << ")" <<endl;   

    // This loop searches for the next best clasifier
    // Doesn't need to test the i classifier again:
    for (int j = i+1; j < nClasif; ++j) {
      // Try the next not-selected classifier.
     double overhead = computeRow(gMatrix, nClasif, selectedClassifiers, i, selectedClassifiers[j]);

      // Evaluate a G natrix of size 'nClasif' and 'i+1' selectedClassifiers
      double G = previousG + overhead;
  //    cout << "Trying " << selectedClassifiers[i] << "(" << G << ")";
      if (G < minG) {
         minG = G;
         bestClasif = j;
         cout << " selected!";
      }
      cout << endl;
    }

    previousG = minG;   

    // Selects the best classifier, that is:
    // swaps the bestClasif with the clasifiers which is in the i-th position:
    FnsOrdClasUtils::swap(selectedClassifiers[bestClasif], selectedClassifiers[i]);

    ens->Exchange(i, bestClasif);
    cout << "Current subset: "<< endl;
    FnsOrdClasUtils::print(selectedClassifiers, i+1, cout);
  }

  ens->SetClassifiersToUse(stop);

double end = clock();

  // This block to open the file to store the number of evaluations
  string filename2("time_greedy.txt");
  ofstream os2(filename2.c_str(), ios::out | ios::app);
  os2 << (end - start) / CLOCKS_PER_SEC << endl;
  os2.close();


  delete[] selectedClassifiers;
  return this;
}


double OrdenaPorGGreedy::computeRow(double** G, int size, int* previousSelected, int previousSize, int trial) {
    double overhead = 0;
    int i = trial;
    for (int j = 0 ; j < previousSize ; ++j) {
      overhead +=  2 * G[ previousSelected[ j ] ][ i ];
    }

    overhead += G[i][i];
    return overhead;
}

// ------------------------------ DHL

void OrdenaPorRERegularizadoPorDiversidad::CreateParams()
{
  Parametros.push_back(new ParametroEnsemble());

  Parametro *p = new ParametroData();
  p->PonPropiedades(false, 0, "Datos de ejemplo", "Ejemplos con los que se realizara la seleccion");
  Parametros.push_back(p);

  p = new ParametroNumero(0.0, 1.0);
  p->PonPropiedades(false, 0, "Parametro rho del método.", "");
  p->SetPorOmision("0.20");
  Parametros.push_back(p);

  p = new ParametroNumero(0.0, 500.0);
  p->PonPropiedades(true, 0, "Parar en (Tamano del subconjunto).", "Por omision: se ordena hasta el final");
  Parametros.push_back(p);
}
//---------------------------------------------------------------------------
void OrdenaPorRERegularizadoPorDiversidad::DarDeAlta()
{
  ifns().AnadirInterfazFuncion(new OrdenaPorRERegularizadoPorDiversidad(), "Tools");
}

Mandato* OrdenaPorRERegularizadoPorDiversidad::Ejecutar()
{
  Param(0)->Ejecutar();
  Param(1)->Ejecutar();
  Param(2)->Ejecutar();

  Ensemble *ens = (Ensemble *)Param(0)->ComoDatos();
  NomData *data = (NomData *)Param(1)->ComoDatos();
  double rho = (double )Param(2)->ComoNumero();

  DoOrdenaPorRERegularizadoPorDiversidad(ens, data, rho);

  return this;
}

void DoOrdenaPorRERegularizadoPorDiversidad(Ensemble *ens, NomData *data, double rho)
{
  int **Mcd;  //Matriz con todas las clasificaciones de los clasificadores del ens para cada dato de ens->data
  double **Ddc;  //Matriz con el nUmero de votos para cada datos y para cada clase
  int Min1;
  int iMin1=-1;
  int ErrAux;
  int *ErrorClas; //Error de cada clasificador
  int *ClasReal; //Clase de cada dato
  int *ClasAct;  //Clase del conjunto para el nUmero de clasificadores actuales
  int nClases;
  int nDatos;
  int nClasf;
  StructToSort * DivPerClassifier;

  nClasf = ens->GetClassifiersToUse();

  if (nClasf==1) return;

  NomData *dataold = (NomData*)ens->GetData();
  ens->SetData(data);

  nClases = data->NumClass;
  nDatos = data->GetNTotal();

  //Obtiene la clase de cada dato
  ClasReal = new int[nDatos];
  for(int j=0;j<nDatos;j++)
    ClasReal[j] = data->GetDatClass(j);

  //Obtiene la clsificacion para cada clasificador del conjunto y dato
  //Tambien obtiene el clasificador con menor error
  ErrorClas = new int[nClasf];

  DivPerClassifier = new StructToSort[nClasf];

  Mcd = new int*[nClasf];
  Min1 = nDatos;
  for(int i = 0; i < nClasf; i++) {
    Classifier *c = ens->GetClassifier(i);
    Mcd[i] = new int[nDatos];
    ErrAux = 0;
    for(int j = 0; j < nDatos; j++) {
      Mcd[i][j] = c->Classify(j);
      if (Mcd[i][j] != ClasReal[j]) ErrAux++;
    }
    ErrorClas[i] = ErrAux;
    if (ErrorClas[i] < Min1) {
      Min1 = ErrorClas[i];
      iMin1 = i;
    }
  }

  int *cc, toEvaluate;

  //Ponemos el clasificador ocn menor error en la primera posicion
  ens->Exchange(0, iMin1);
  cc = Mcd[0];
  Mcd[0] = Mcd[iMin1];
  Mcd[iMin1] = cc;
  ErrAux = ErrorClas[iMin1];
  ErrorClas[iMin1] = ErrorClas[0];
  ErrorClas[0] = ErrAux;

  //Distribuciones de clase por ejemplo para el subensemble t
  Ddc = new double*[nDatos];
  for(int i=0;i<nDatos;i++) {
    Ddc[i] = new double[nClases];
    for(int j=0;j<nClases;j++) {
      Ddc[i][j] = 0.0;
    }
  }
  ClasAct = new int[nDatos]; //Clase predicha de cada dato por el subensemble t
  for(int j=0;j<nDatos;j++) {
    Ddc[j][Mcd[0][j]]++;
    ClasAct[j] = Mcd[0][j];
  }

  // Recorremos los clasificadores para ordenarlos

  for(int iClasf = 1; iClasf < nClasf; iClasf++) {

	// Calculamos la diversidad  de cada clasificador restante y el sub-conjunto actual

	for (int j = iClasf, u = 0; j < nClasf ; j++, u++) {
		DivPerClassifier[ u ].diversity = 0;
		DivPerClassifier[ u ].order = u;

		for (int k = 0 ; k < nDatos;k++)
			if (ClasAct[ k ] == Mcd[ j ][ k ])
				DivPerClassifier[ u ].diversity++;
			else
				DivPerClassifier[ u ].diversity--;

		DivPerClassifier[ u ].diversity /= nDatos;
	}

	// Ordenamos el array para dar prioridad a los clasificadores que incrementarian la diversidad

	qsort(DivPerClassifier, nClasf - iClasf, sizeof(StructToSort), cmpToSort);

	// We check how many classifiers we need to check

	toEvaluate = (int) ceil(rho * (nClasf - iClasf));
	
	int iErrMin = -1;
	int ErrTemp;
	double *Dctemp = new double[nClases];
 
	int ErrMin1 = nDatos;

	for(int i = 0 ;i < toEvaluate;i++) {

		int ErrTempSube = 0;
		int ErrTempBaja = 0;
		int iClasfToTest = iClasf + DivPerClassifier[ i ].order;

		for(int j = 0; j < nDatos; j++) {

			if (Mcd[ iClasfToTest ][ j ] == ClasAct[ j ] || Ddc[ j ][ Mcd[ iClasfToTest ][ j ] ] + 1 < Ddc[ j ][ ClasAct[ j ] ])
				continue; //Este ejemplo ni sube ni baja el error

			for(int k=0;k<nClases;k++) {
				Dctemp[k] = Ddc[j][k];
			}

			Dctemp[ Mcd[ iClasfToTest ][ j ] ]++;

			int ClasTemp = ens->WhichClass(Dctemp, nClases);

			if (ClasTemp >= 0 && ClasTemp != ClasAct[ j ]) {
				if (ClasAct[ j ] == ClasReal[ j ]) 
					ErrTempSube++;    //Sube el error
				else if (ClasTemp == ClasReal[ j ]) 
					ErrTempBaja--; //Baja el error
			}
		}

		ErrTemp = (ErrTempBaja + ErrTempSube);

		if (ErrTemp < ErrMin1) {
			ErrMin1 = ErrTemp;
			iErrMin = iClasfToTest;
		}
	}

	delete []Dctemp;

	// iErrMin contains the classifier to be included in the ensemble
	// Now we include that element in the ensemble
 
	//Intercambiamos clasificadores y datos

	ens->Exchange(iClasf, iErrMin);
	cc = Mcd[iClasf];
	Mcd[iClasf] = Mcd[iErrMin];
	Mcd[iErrMin] = cc;
	ErrAux = ErrorClas[iErrMin];
	ErrorClas[iErrMin] = ErrorClas[iClasf];
	ErrorClas[iClasf] = ErrAux;

	//Actualizamos distribuciones
	for(int j = 0; j < nDatos; j++) {
		Ddc[j][Mcd[iClasf][j]]++;
		int ct = ens->WhichClass(Ddc[j], nClases);
		if (ct>=0) ClasAct[j] = ct;
	}

  } // Fin del bucle que ordena

  for(int i = 0; i < nDatos; i++) delete []Ddc[i];
  delete []Ddc;
  for(int i = 0; i < nClasf; i++) delete []Mcd[i];
  delete []Mcd;
  delete []ErrorClas;
  delete []ClasReal;
  delete []ClasAct;
  delete []DivPerClassifier;

  ens->SetData(dataold);

  return;
}


//---------------------------------------------------------------------------

void FnsOrdClasUtils::swap(int &a, int &b) {
  int temp = a;
  a = b;
  b = temp;
}

void FnsOrdClasUtils::swap(int *array, int i, int j) {
  if (i == j)
    return;
  swap(array[i], array[j]);
}

void FnsOrdClasUtils::print(int *data, int n, ostream &os) {
  for (int i = 0; i < n; ++i) {
    os << data[i] << " ";
  }
  os << endl;
}

void FnsOrdClasUtils::print(double **data, int n, ostream &os) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
       os << data[i][j] << " ";
    }
    os << endl;
  }
  os << endl;
}

