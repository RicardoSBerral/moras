
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

//void process(FILE *f, int opt, bool first);
Matriz *mean(Matriz *v, bool horizontal, bool desv);
Matriz *minimos(Matriz *m, bool last);
Matriz *maximos(Matriz *m, bool last);
Matriz *vals(Matriz *m, Matriz *indexs);

int main(int argc,char**argv)
{
  FILE *f;
//  char dir[300];
  extern char *optarg;
//  extern int optind;
  char *filename=0, o;
//  bool vertical = false;
//  bool desv_std = false;
  int opt = OPT_MEAN | OPT_HORZ;
  bool first = false;

  /*  Process options  */
  while ( (o = getopt(argc, argv, "f:vsim::M::")) != EOF ) {
    switch (o) {
      case 'f':   
        filename = optarg;
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
      case 'i':   
        opt |= OPT_VAL;
        opt &= ~OPT_MEAN;
        break;
      case '?':   printf("unrecognised option\n");
        exit(1);
    }
  }

  f= filename ? fopen(filename, "r") : stdin;
  if (!f) {
    printf("Error openening file: %s\n", filename);
    exit(1);
  }

/*  process(f, opt, first);

  return 0; 
}

void process(FILE *f, int opt, bool first)
{*/
     Matriz *mt, *mres;
//     int *mins;
//     double minmean, mintrain;

     mt = Matriz::leerDeFichero(f);
     if (!mt) exit(1);
   
     if (!(opt & OPT_HORZ)) {
         Matriz *mm = mt->transpuesta();
         delete mt;
         mt = mm;
     }

    if (opt & OPT_MEAN) {   //media
        mres = mean(mt, true/*opt & OPT_HORZ*/, opt & OPT_STD);
    }
    else if (opt & OPT_MAX) {
        mres = maximos(mt, !first);
    }
    else if (opt & OPT_MIN) {
        mres = minimos(mt, !first);
    }
    else if (opt & OPT_VAL) {
        Matriz *indexs = Matriz::leerDeFichero(stdin);
        mres = vals(mt, indexs);
        delete indexs;
    }

    if (mres) {
         if (!(opt & OPT_HORZ)) {
             Matriz *mm = mres->transpuesta();
             delete mres;
             mres = mm;
         }
         mres->guardarAFichero(stdout); 
     }

     //min

     delete mt;
     if (mres) delete mres;

     return 0;
}

Matriz *maximos(Matriz *m, bool last)
{
  int i, j, imax;
  double max;
  Matriz *maxs = new Matriz(m->filas(), 2);

  for(i=0;i<m->filas();i++) {
    max = (*m)[i][0];
    imax = 0;
    for(j=1;j<m->columnas();j++) {
      if ((max<=(*m)[i][j] && last) || max<(*m)[i][j]) {
        max = (*m)[i][j];
        imax = j;
      }
    }
    (*maxs)[i][0] = max;
    (*maxs)[i][1] = imax;
  }

  return maxs;
}
Matriz *minimos(Matriz *m, bool last)
{
  int i, j, imin;
  double min;
  Matriz *mins = new Matriz(m->filas(), 2);

  for(i=0;i<m->filas();i++) {
    min = (*m)[i][0];
    imin = 0;
    for(j=1;j<m->columnas();j++) {
      if ((min>=(*m)[i][j] && last) || min>(*m)[i][j]) {
        min = (*m)[i][j];
        imin = j;
      }
    }
    (*mins)[i][0] = min;
    (*mins)[i][1] = imin;
  }

  return mins;
}
Matriz *mean(Matriz *v, bool horizontal, bool desv)
{
  int f, c;

  if (horizontal) {
    f = desv ? 2 : 1;
    c = v->columnas();
  }
  else {
    f = v->filas();
    c = desv ? 2 : 1;
  }
  Matriz *val = new Matriz(f, c);
  
  int i, j;
  if (!horizontal) {
    for(i=0;i<v->filas();i++){
      for(j=0;j<v->columnas();j++){
        (*val)[i][0] += (*v)[i][j];
        if (desv) (*val)[i][1] += ((*v)[i][j]*(*v)[i][j]);
      }
      (*val)[i][0] /= v->columnas();
      if (desv) (*val)[i][1] = sqrt((*val)[i][1]/v->columnas() -
                                                    (*val)[i][0]*(*val)[i][0]);
    }
  }
  else {
    for(j=0;j<v->columnas();j++){
      for(i=0;i<v->filas();i++){
        (*val)[0][j] += (*v)[i][j];
        if (desv) (*val)[1][j] += ((*v)[i][j]*(*v)[i][j]);
      }
      (*val)[0][j] /= v->filas();
      if (desv) (*val)[1][j] = sqrt((*val)[1][j]/v->filas() - 
                                                    (*val)[0][j]*(*val)[0][j]);
    }
  }
  return val;
}

Matriz *vals(Matriz *m, Matriz *ixs)
{
  Matriz *val = new Matriz(m->filas(), 1);
  int cols = m->columnas();

  for(int i=0;i<m->filas();i++) 
    (*val)[i][0] = (*ixs)[i][0]>=cols || (*ixs)[i][0]<0 ? 
        numeric_limits<double>::quiet_NaN() : (*m)[i][(int)(*ixs)[i][0]];

  return val;
}
