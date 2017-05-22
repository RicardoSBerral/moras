
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

#define PI 3.1415926535897932384626433832795

int main(int argc,char**argv)
{
  FILE *f;
  extern char *optarg;
  char *filename=0, o;
  double value = 2.0/3.0;
  double corte;
  double stdde;
  int min_values = 1;
  int start_at = 0;
  bool average = false;
  bool fixed = false;
  bool desv = false;

  while ( (o = getopt(argc, argv, "f:v:m:p:ksx:")) != -1 ) {
    switch (o) {
      case 'f':   
        filename = optarg;
        break;
      case 'v':
        value = atof(optarg);
        break;
      case 'm':
        min_values = atoi(optarg);
        break;
      case 'p':
        start_at = atoi(optarg);
        break;
      case 's':
        desv = true;
        break;
      case 'k':
        average = true;
        break;
      case 'x':
        fixed = true;
        corte = atof(optarg);
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

   Matriz *mt, *mres=0;

   mt = Matriz::leerDeFichero(f);
   if (!mt) exit(1);
 
   for(int i=0;i<mt->filas();i++) {
     int j;
     if (average) {
       stdde = 0.0;
       corte = 0.0;
       for(j=0;j<mt->columnas();j++) { 
         if ((*mt)[i][j]<PI/2.0) {
           corte += (*mt)[i][j];
           if (desv) stdde += (*mt)[i][j]*(*mt)[i][j];
         }
         else break;
       }
       corte /= j;
       if (desv) stdde = sqrt(stdde/j - corte*corte);
     }
     else if (fixed) {
     }
     else {
       stdde = 0.0;
       corte = 0.0;
       for(j=0;j<min_values;j++) corte += (*mt)[i][j+start_at];
       corte /= min_values;
       corte = corte + value*(PI/2.0-corte);
     }
     for(j=0;j<mt->columnas();j++) {
       if ((*mt)[i][j]>corte+stdde) {
         printf("%d\n", j);
         break;
       }
     }
     if (j==mt->columnas()) printf("%d\n", j);
   }
   delete mt;
   if (mres) delete mres;

   return 0;
}

