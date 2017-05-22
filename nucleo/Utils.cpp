//---------------------------------------------------------------------------

#include "Utils.h"
#include "Ensemble.h"
#include "data20.h"
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <cstdlib>

//---------------------------------------------------------------------------
//-----------------------------------------------  Funciones generales  -----
//---------------------------------------------------------------------------
std::string proc_name;

std::string GenerateNewFileName(std::string prev, std::string ext, 
                                                              std::string dir)
{
  int k;
  FILE *f;
  char pid[240];
  char num[24];

  if (dir[dir.size()-1]!=SEP) dir = dir + SEP;

  if ( proc_name.length() > 0 ) {
    sprintf(pid, "%s", proc_name.c_str());
  }
  else {    
    sprintf(pid, "%d", getpid());
  }

  std::string nom = dir + prev + pid + ext;
  k = 0;
  while (0!=(f=fopen(nom.c_str(), "r")) && k<10000) {
    fclose(f);
    sprintf(num, "_%d", k);
    nom = dir + prev + pid + num + ext;
    k++;
  }
  return nom;
}

std::string GeneratePIDFileName(std::string prev, std::string ext, 
                                                              std::string dir)
{
  char pid[24];

  if (dir[dir.size()-1]!=SEP) dir = dir + SEP;

  sprintf(pid, "%d", getpid());
  std::string nom = dir + prev + pid + ext;

  return nom;
}

std::string GenerateProcNameFileName(std::string prev, std::string ext, 
                                                              std::string dir)
{
  char prn[1024];

  if (dir[dir.size()-1]!=SEP) dir = dir + SEP;

  sprintf(prn, "%s", proc_name.c_str());
  std::string nom = dir + prev + prn + ext;

  return nom;
}
void SetProcName(std::string prn)
{
  proc_name = prn;
}
std::string GetProcName()
{
  return proc_name;
}
//---------------------------------------------------------------------------
/* Returns a random integer in [0,n]  */
int RandomInteger (int n)
{
  return (int) floor ( (n+1) * RandomDouble() );
}
//---------------------------------------------------------------------------
/*
 * returns a double in the range [0.0;1.0)
 */
double RandomDouble()
{
  return (1.0 * rand () / (RAND_MAX + 1.0));
}
//---------------------------------------------------------------------------
/*
 */
int RandomIntProportional(double *pop, int npop, double sumPop)
{
  std::vector<double> p(npop);
  for(int i=0;i<npop;i++) p[i] = pop[i];
  return RandomIntProportional(p, sumPop);
}

int RandomIntProportional(std::vector<double> pop, double sumPop)
{
  unsigned i;
  double pos, val;

  pos = (sumPop * rand ()) / (RAND_MAX + 1.0);

  val = 0.0;
  i = 0;
  while (val<pos && i<pop.size()) {
    val += pop[i];
    i++;
  }

  if (i>0) i--;
  return i;
}
//---------------------------------------------------------------------------
/* Returns true with probability prob  */
bool Flip (double prob)
{
  return RandomDouble() < prob ? true : false;
}
//---------------------------------------------------------------------------
/*
 * funcion wue permuta aleatoriamente "npermuts" veces un array de tipo T de 
 *  tamanyo "size".
 */
template<typename T> void RandomPermutation(T *elems, int size, int npermuts)
{
  for(int i=0;i<npermuts;i++) {
    int pos = RandomInteger(size-1);
    T hold = elems[i];
    elems[i] = elems[pos];
    elems[pos] = hold;
  }

}

//---------------------------------------------------------------------------
//-----------------------------------------------------  Combinaciones  -----
//---------------------------------------------------------------------------
Combinaciones::Combinaciones(int nelementos, int tomados_de)
{
  this->nelementos = nelementos;
  for (int i=0;i<tomados_de;i++) {
    indices.push_back(i);
  }
  icomb = 0;
  maxcomb = 1;
  int n, m;
  n = nelementos-tomados_de;
  m = tomados_de;
  if (n<m) {int aux = n; n = m; m = aux;}
  for (__int64 i=nelementos;i>n;i--)
    maxcomb *= i;
  for (__int64 i=1;i<=m;i++)
    maxcomb /= i;
}
//---------------------------------------------------------------------------
__int64 Combinaciones::GetNComb()
{
  return maxcomb;
}
//---------------------------------------------------------------------------
bool Combinaciones::HayMasCombinaciones()
{
  return icomb<maxcomb;
}
//---------------------------------------------------------------------------
std::vector<int> Combinaciones::SiguienteCombinacion()
{
  int iact = 0;
  int tomados_de = indices.size();
  std::vector<int> comb = indices;

  if (HayMasCombinaciones()) {
    for (int i=tomados_de-1;i>=0;i--) {
      if (indices[i]<nelementos-tomados_de+i) {
        indices[i]++;
        iact = i;
        break;
      }
    }
    for (int i=iact+1;i<tomados_de;i++)
      indices[i] = indices[iact] + i - iact;
  }

  icomb++;

  return comb;
}

//---------------------------------------------------------------------------
//-------------------------------------------------  FuncionDeProgreso  -----
//---------------------------------------------------------------------------
void FuncionDeProgreso::Start(int numeroDeLlamadas)
{
  this->numeroDeLlamadas = numeroDeLlamadas;
  if (f) f->Start(numeroDeLlamadas);
}
//---------------------------------------------------------------------------

void FuncionDeProgreso::End()
{
  if (f) f->End();
}

//---------------------------------------------------------------------------
//---------------------------------------------------  ProgresoConsola  -----
//---------------------------------------------------------------------------
void ProgresoConsola::Start(int numeroDeLlamadas)
{
  FuncionDeProgreso::Start(numeroDeLlamadas);
  llamadasPorTic = (double)numeroDeLlamadas/longDeBarra;
  char *kk = new char[longDeBarra];
  if ( llamadasPorTic >= 1.0 ) {
    sprintf(kk, "%g %ss per tic -------->", llamadasPorTic, texto.c_str());
  }
  else {
    sprintf(kk, "%g tics per %s -------->", 1.0/llamadasPorTic, texto.c_str());
  }
  printf("|%48s|\n", kk);
  iLlamadas = iTics = 0;
  t0 = clock();
  int t1 = time(0);
  printf("Starting time: %s", ctime((const time_t*)&t1));
  delete []kk;
}
//---------------------------------------------------------------------------
void ProgresoConsola::End()
{
  FuncionDeProgreso::End();
  time_t t;
  time(&t);
  printf("| Finish time: %s", ctime(&t));
}
//---------------------------------------------------------------------------
bool ProgresoConsola::operator()(int incr)
{
   bool algo = next(incr);
   iLlamadas += incr;

   while (iLlamadas>llamadasPorTic*(iTics+1)) {
       if ( t0 >= 0 ) {
         int t1 = clock();
         double t = ((double)t1-t0)/CLOCKS_PER_SEC;
         int s = (int)(0.5 + (double)(numeroDeLlamadas-iLlamadas)/((double)iLlamadas/t));
 /*       int h = s/3600;
         s = s - s/3600;
         int m = s/60;
         s = s - s/60;*/
         t1 = time(0) + s;
         printf("Estimated finish time: %s", ctime((const time_t*)&t1));
         t0 = -1;
       }

       printf("=");
       fflush(stdout);
       iTics++;
   }

   return true && algo;
}

//---------------------------------------------------------------------------
//--------------------------------------------------  ProgresoEnsemble  -----
//---------------------------------------------------------------------------
void ProgresoEnsemble::Start(int numeroDeLlamadas)
{
  FuncionDeProgreso::Start(numeroDeLlamadas);
  votos = new double*[data->GetNTrain()];
  for(int i=0;i<data->GetNTrain();i++) {
    votos[i] = new double[((NomData*)data)->NumClass];
    for(int j=0;j<((NomData*)data)->NumClass;j++) {
      votos[i][j] = 0.0;
    }
  }
  nClasf = 0;
  sal = GenerateNewFileName("ens_test_", ".txt");
  FILE *f=fopen(sal.c_str(), "w");
  fclose(f);
}
//---------------------------------------------------------------------------
void ProgresoEnsemble::End()
{
  FuncionDeProgreso::End();
}
//---------------------------------------------------------------------------
bool ProgresoEnsemble::operator()(int incr)
{
  bool algo = next(incr);
  nClasf += incr;
  if (incr>0) {
    Classifier *cls = ens->GetClassifier(nClasf-1);
    Data *holdData = cls->GetData();
    cls->SetData(data);
    int errs = 0;
    for(int i=0;i<data->GetNTrain();i++) {
      votos[i][cls->Classify(i)] += ens->GetClassifierWeight(nClasf-1); 
      if (data->GetDatClass(i)!=
                  Classifier::WhichClass(votos[i], ((NomData*)data)->NumClass)) 
        errs++;
    }
    cls->SetData(holdData);
    FILE *f=fopen(sal.c_str(), "at");
    fprintf(f, "%g\n", (double)errs/data->GetNTrain());
    fclose(f);
  }
  return true && algo;
}
//---------------------------------------------------------------------------

