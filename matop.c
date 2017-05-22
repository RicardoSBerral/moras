
#include "Matriz.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#define OPTIONS "fsrcCxph"

void display_help();

int main(int argc,char**argv)
{
  FILE *f;
  extern char *optarg;
  char *filename=0, o;
  char operations[200];
  Matriz *mats[200];
  double nums[200];
  int iop=0;
  bool corteMin = false; //Iguala las matrices al tamaño minimo
  bool corteMax = false; //Iguala las matrices al tamaño maximo
  int mf, mc, Mf, Mc;
  char *fin;

  Mf = Mc = 0;
  mf = mc = 2100000000;
  /*  Process options  */
  while ( (o = getopt(argc, argv, "f:s:r:cCx:p:h")) != EOF ) {
    mats[iop] = 0;
    switch (o) {
      case 'h':
        display_help();
        exit(0); 
      case 'f':
        filename = optarg;
        break;
      case 's':
        mats[iop] = Matriz::leerDeFichero(optarg);
        if (!mats[iop]) exit(1);
        break;
      case 'r':
        mats[iop] = Matriz::leerDeFichero(optarg);
        if (!mats[iop]) exit(1);
        break;
      case 'x':
        nums[iop] = strtod(optarg, &fin);
        if (fin==optarg) exit(1);
        break;
      case 'p':
        nums[iop] = strtod(optarg, &fin);
        if (fin==optarg) exit(1);
        break;
      case 'C':
        corteMin = false;
        corteMax = true;
        break;
      case 'c':
        corteMin = true;
        corteMax = false;
        break;
      case 't':
        break;
      case '?':
        printf("unrecognised option\n");
        exit(1);
    }
    if (mats[iop]) {//Maximos y mins de cols para hacer el corte de opts c y C
      if (mats[iop]->filas()>Mf) Mf = mats[iop]->filas();
      if (mats[iop]->filas()<mf) mf = mats[iop]->filas();
      if (mats[iop]->columnas()>Mc) Mc = mats[iop]->columnas();
      if (mats[iop]->columnas()<mc) mc = mats[iop]->columnas();
    }
    operations[iop++] = o;   
  }

  f= filename ? fopen(filename, "r") : stdin;
  if (!f) {
    printf("Error openening file: %s\n", filename);
    exit(1);
  }


  Matriz *base, *aux;

  base = Matriz::leerDeFichero(f);
  if (!base) exit(1);

  if (base->filas()>Mf) Mf = base->filas();
  if (base->filas()<mf) mf = base->filas();
  if (base->columnas()>Mc) Mc = base->columnas();
  if (base->columnas()<mc) mc = base->columnas();

  for(int i=0;i<iop;i++) {
    aux = 0;
    if (corteMin) {
      base->redim(mf, mc);
      if (mats[i]) mats[i]->redim(mf, mc);
    }
    else if (corteMax) {
      base->redim(Mf, Mc);
      if (mats[i]) mats[i]->redim(Mf, Mc);
    }
    switch (operations[i]) {   
      case 'a':
        for (int ii=0;ii<base->filas();ii++)
          for (int jj=0;jj<base->columnas();jj++)
            (*base)[ii][jj] = fabs((*base)[ii][jj]);
        break;
      case 's':
        aux = base->sumar(*mats[i]);
        if (!aux) exit(1);
        break;
      case 'r':
        aux = base->restar(*mats[i]);
        if (!aux) exit(1);
        break;
      case 't':
        aux = base->transpuesta(); 
        if (!aux) exit(1);
        break;
      case 'x':
        for (int ii=0;ii<base->filas();ii++)
          for (int jj=0;jj<base->columnas();jj++)
            (*base)[ii][jj] = (*base)[ii][jj]*nums[i];
        break;
      case 'p':
        for (int ii=0;ii<base->filas();ii++)
          for (int jj=0;jj<base->columnas();jj++)
            (*base)[ii][jj] = pow((*base)[ii][jj], nums[i]);
        break;
    }
    if (mats[iop]) delete mats[iop];
    if (aux) {
      delete base;
      base = aux;
    }
  }

  base->guardarAFichero(stdout); 
  delete base;

//    operations[iop] ='\0';
//    printf("%s", operations);

  return 0;
}

void display_help()
{
  char o;
  int i;

  for(i=0;i<strlen(OPTIONS);i++) {
    o = OPTIONS[i];
    printf("\t-%c: ", o);
    switch (o) {
      case 0:   
      case 'h':
        printf("%s\n", "Display this help");
        break;
      case 'f':
        printf("%s\n", "file name VALUE");
        break;
      case 's':
        printf("%s\n", "Suma matriz pasada como argumento (nombre de fichero)");
        break;
      case 'r':
        printf("%s\n", "Resta matriz pasada como argumento (nombre de fichero)");
        break;
      case 'x':
        printf("%s\n", "Multiplica por un numero cada elemento de la matriz");
        break;
      case 'p':
        printf("%s\n", "Eleva a un numero cada elemento de la matriz");
        break;
      case 'C':
        printf("%s\n", "Iguala las matrices al tamaño maximo");
        break;
      case 'c':
        printf("%s\n", "Iguala las matrices al tamaño minimo");
        break;
      case 't':
        break;
    }
  }
}
