
#include "Matriz.h"
#include <sys/dir.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <limits>

using namespace std;

#define OPT_ALL 0xFFFFFFFF
#define OPT_MEAN 0x00000001
#define OPT_MIN 0x00000002
#define OPT_MAX 0x00000004
#define OPT_STD 0x00000008
#define OPT_HORZ 0x00000010
#define OPT_VAL 0x00000020
#define OPT_WDL 0x00000040
#define OPT_SLD 0x00000080
#define OPT_RNK 0x00000100

#define OPTIONS "wfvsiImMkhpXTr"

//void process(FILE *f, int opt, bool first);
void windrawloss(Matriz *m, int total, Matriz *mtt, double pvalor, char *sep=0, char *fin=0);
Matriz *slide_mean(Matriz *mt, int kn);
Matriz *ranks(Matriz *mt);

void display_help();

int main(int argc,char**argv)
{
  FILE *f;
//  char dir[300];
  extern char *optarg;
//  extern int optind;
  char filename[260], fn_ttest[260], o;
  int kn = 5;
  int opt = OPT_MEAN | OPT_HORZ;
  bool first = false;
  int lopt, total;
  bool positions = false;
  double pvalor=0.01;
  char *sep, *fin;

  static struct option long_options[] = {
      {"slide_mean", 2, 0, 0},
      {0, 0, 0, 0}
  };

  total = 0;
  sep = "\t";
  fin = "";

  filename[0] = '\0';
  fn_ttest[0] = '\0';

  /*  Process options  */
//  while ( (o = getopt_long(argc, argv, "wf:vsim::M::", 
//                                               long_options, &lopt)) != -1 ) {
    while ( (o = getopt(argc, argv, "w::f:vsiIm::M::k::hp:XT::r")) != -1 ) {
      switch (o) {
        case 0:   
        case 'h':
          display_help();
          exit(0); 
        case 'k':   
            opt |= OPT_SLD;
            opt &= ~OPT_MEAN;
            kn = optarg ? atoi(optarg) : 5;
          break;
        case 'f':   
          strcpy(filename, optarg);
          break;
        case 'v':   
          opt &= ~OPT_HORZ;
          break;
        case 'M':   
          if (optarg && optarg[0]=='1') first = true;
          opt |= OPT_MAX;
          opt &= ~OPT_MEAN;
          break;
        case 'm':
          if (optarg && optarg[0]=='1') first = true;
          opt |= OPT_MIN;
          opt &= ~OPT_MEAN;
          break;
        case 's':   
          opt |= OPT_STD;
          break;
        case 'I':
          opt |= OPT_VAL;
          opt &= ~OPT_MEAN;
          positions = true;
          break;
        case 'i':
          opt |= OPT_VAL;
          opt &= ~OPT_MEAN;
          break;
        case 'w':   
          if (optarg) strcpy(fn_ttest, optarg);
          opt |= OPT_WDL;
          opt &= ~OPT_MEAN;
          break;
        case 'p':   
            pvalor = atof(optarg);
          break;
        case 'T':   
          total = optarg ? 2 : 1;
          break;
        case 'X':   
            sep = " & ";
            fin = "\\\\";
          break;
        case 'r':   
          opt |= OPT_RNK;
          opt &= ~OPT_MEAN;
          break;
        case '?':   printf("unrecognised option\n");
          exit(1);
      }
    }

    f = strlen(filename) > 0 ? fopen(filename, "r") : stdin;
    if (!f) {
      printf("Error openening file: %s\n", filename);
      exit(1);
    }

     Matriz *mt, *mres=0;

     mt = Matriz::leerDeFichero(f);

     fclose(f);

     if (!mt) exit(1);
   
     if (!(opt & OPT_HORZ)) {
         Matriz *mm = mt->transpuesta();
         delete mt;
         mt = mm;
     }

    if (opt & OPT_MEAN) {   //media
        mres = mt->mean(true/*opt & OPT_HORZ*/, opt & OPT_STD);
    }
    else if (opt & OPT_MAX) {
        mres = mt->maximos(!first);
    }
    else if (opt & OPT_MIN) {
        mres = mt->minimos(!first);
    }
    else if (opt & OPT_VAL) {
        Matriz *indexs = Matriz::leerDeFichero(stdin);
        if (positions) {
            for(int i=0;i<indexs->filas();i++) {
                (*indexs)[i][0] = (*indexs)[i][0] - 1;
            }
        }
        mres = mt->vals(indexs);
        delete indexs;
    }
    else if (opt & OPT_WDL) {
        Matriz *mtt=0;

        if (strlen(fn_ttest) > 0 ) {
            f = fopen(fn_ttest, "r");
            mtt = Matriz::leerDeFichero(f);
            fclose(f);
            if (mtt->columnas() < mt->columnas()*(mt->columnas()-1)/2) {
                delete mtt;
                mtt = 0;
            }
        }

        windrawloss(mt, total, mtt, pvalor, sep, fin);

        if (mtt) delete mtt;

    }
    else if (opt & OPT_SLD) {
        mres = slide_mean(mt, kn);
    }
    else if (opt & OPT_RNK) {
        mres = ranks(mt);
    }

    if (mres) {
         if (!(opt & OPT_HORZ)) {
             Matriz *mm = mres->transpuesta();
             delete mres;
             mres = mm;
         }
         mres->guardarAFichero(stdout, sep, fin); 
     }

     //min

     delete mt;
     if (mres) delete mres;

     return 0;
}

void windrawloss(Matriz *m, int total, Matriz *mtt, double pvalor, char *sep, char *fin)
{

  if (!sep) sep = "\t";
  if (!fin) fin = "";

  for(int i=0;i<m->columnas();i++) {
    int twins, tloss, tdraws;
    twins = tloss = tdraws = 0;
    for(int j=0;j<m->columnas();j++) {
      if (i==j) printf("X%s", sep);//j < m->columnas()-1 ? sep : fin);
      else {
        int wins, draws, loss;
        wins = draws = loss = 0;
        for(int k=0;k<m->filas();k++) {
          if (mtt) {
            int it = i < j ? i * m->columnas() - i * (i+1) / 2 + j - i - 1 :
                             j * m->columnas() - j * (j+1) / 2 + i - j - 1 ;
            if ((*mtt)[k][it] > pvalor || (*mtt)[k][it] < -pvalor ) draws++;
            else if ((*m)[k][i] < (*m)[k][j]) wins++;
            else loss++;
          }
          else {
            if ((*m)[k][i]<(*m)[k][j] ) wins++;
            else if ((*m)[k][i]==(*m)[k][j]) draws++;
            else if ((*m)[k][i]>(*m)[k][j]) loss++;
          }
        }
        twins += wins;
        tloss += loss;
        tdraws += draws;
        printf("%d/%d/%d%s", wins, draws, loss, total==0 && j == m->columnas()-1 ? fin : sep);
      }
    }
    int t = twins + tloss + tdraws;
    if (total==1)      printf("%d/%d/%d%s\n", twins, tdraws, tloss, fin);
    else if (total==2) printf("%.1f/%.1f/%.1f%s\n", 100.0*twins/t, 100.0*tdraws/t, 100.0*tloss/t, fin);
    else               printf("%s\n", fin);
  }
}
Matriz *slide_mean(Matriz *mt, int kn)
{
  Matriz *m = new Matriz(*mt);
  int i_av = kn/2;
  int i_re = -((1-kn)/2);

  for(int i=0;i<m->filas();i++) {
    double acum = 0.0;
    for(int j=0;j<i_av;j++) acum += (*mt)[i][j];
    for(int j=0;j<i_re;j++) acum += (*mt)[i][0];
    for(int j=0;j<m->columnas();j++) {
      acum += (*mt)[i][j + i_av >= m->columnas() ? m->columnas()-1 : j + i_av];
      (*m)[i][j] = acum/kn;
      acum -= (*mt)[i][j - i_re < 0 ? 0 : j - i_re];
    }
  }

  return m;
}
Matriz *ranks(Matriz *mt)
{
  Matriz *m = new Matriz(*mt);

  // Muy poco eficiente pero que mas da
  for (int i = 0; i < m->filas(); i++) {
    for (int j = 0; j < m->columnas(); j++ ) {
      int menores = 0;
      int iguales = 0;
      for (int k = 0; k < m->columnas(); k++ ) {
        if ( k == j ) continue;
        if ( (*mt)[i][j] == (*mt)[i][k] ) iguales++;
        else if ( (*mt)[i][j] > (*mt)[i][k] ) menores++;
      }
      (*m)[i][j] = ( 2.0 * menores + iguales + 2 ) / 2.0;
    }
  }

  return m;
}
void display_help()
{
  char o;
  unsigned int i;

  for(i=0;i<strlen(OPTIONS);i++) {
    o = OPTIONS[i];
    printf("\t-%c: ", o);
    switch (o) {
      case 0:   
      case 'h':
        printf("%s", "Display this help");
        break;
      case 'k':   
        printf("%s", "Sliding mean [VALUE](Default 5)");
        break;
      case 'f':   
        printf("%s", "file name VALUE");
        break;
      case 'v':   
        printf("%s", "Traspose");
        break;
      case 'M':   
        printf("%s", "Maximum [VALUE](if VALUE==1 =>first else =>last(default))");
        break;
      case 'm':
        printf("%s", "Minimum [VALUE](if VALUE==1 =>first else =>last(default))");
        break;
      case 's':   
        printf("%s", "Standard deviation");
        break;
      case 'I':
        printf("%s", "Values by position");
        break;
      case 'i':
        printf("%s", "Values by index");
        break;
      case 'p':
        printf("%s", "P valor para opcion -w");
        break;
      case 'X':   
        printf("%s", "Salida de tabla en formato LaTeX");
        break;
      case 'T':   
        printf("%s", "Pone totales en tablas de win/draw/loss [Con valor op 2: los pone en %]");
        break;
      case 'r':   
        printf("%s", "Calcula el orden de cada celda dentro de cada fila (ranks)");
        break;
      case 'w':   
        printf("%s", "Win-draw-loss table [Nombre fichero con pvalores]");
        break;
    }
    printf("\n");
  }
}
