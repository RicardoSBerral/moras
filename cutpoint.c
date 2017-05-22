
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

int main(int argc,char**argv)
{
  FILE *f;
  extern char *optarg;
  char *filename=0, o;
  double value = 0.0;
  bool otra_curva = false;
  double ref_arriba;

  while ( (o = getopt(argc, argv, "f:v:i")) != -1 ) {
    switch (o) {
      case 'f':   
        filename = optarg;
        break;
      case 'v':
        value = atof(optarg);
        break;
      case 'i':
        otra_curva = true;
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

   Matriz *mt, *mref;

   mt = Matriz::leerDeFichero(f);
   if (!mt) exit(1);
   if (otra_curva) {
     mref = Matriz::leerDeFichero(stdin);
     if (mref->filas()<mt->filas()) exit(1);
     if (mref->columnas()<mt->columnas()) {
       int old_col = mref->columnas();
       mref->redim(mref->filas(), mt->columnas());
       for(int i=0;i<mref->filas();i++) {
         for(int j=old_col;j<mref->columnas();j++) {
           (*mref)[i][j] = (*mref)[i][old_col-1];
         }
       }
     }
   }
   else {
     mref = new Matriz(mt->filas(), mt->columnas(), value);
   }
   if (!mref) exit(1);
 
   for(int i=0;i<mt->filas();i++) {
     ref_arriba = (*mref)[i][0] - (*mt)[i][0];
     for(int j=1;j<mt->columnas();j++) {
       if (((*mref)[i][j]-(*mt)[i][j])*ref_arriba<0.0)
         printf("%d\t", j);
       if ( (*mref)[i][j]-(*mt)[i][j]!=0.0) ref_arriba = (*mref)[i][j]-(*mt)[i][j];
     }
     printf("\n");
   }

   delete mt;
   delete mref;

   return 0;
}

