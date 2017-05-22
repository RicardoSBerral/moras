
#include "Matriz.h"
#include <sys/dir.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <stdlib.h>

using namespace std;

void print(Matriz *vals, vector<string> titles, string output, string title);
void oddify(Matriz *m);
Matriz *processdir(string dir, vector<string> &names);
Matriz *processfiles(char *test, char *train, char *aux, int naux);
Matriz *minimos(Matriz *m, bool last);
//void freematrix(Matriz *m); 
//Matriz *transpose(Matriz *m1);
Matriz *mean(Matriz *v, bool vertical);
double meanvals(Matriz *v, int *pos);
//void writemat(char *filename, Matriz *m) ;
//Matriz *readmat(char *filename);
void writelatexheader(FILE *f);
void tolatex(Matriz *m, vector<string> &names, FILE *f);

int main(int argv,char**argc)
{
  char name[300];
  char cmd[300];
  char dir[300];
  char *res;
  
  //Se crean nombres de ficheros temporales
  strcpy(name, "./XXXXXX");
  mkstemp(name);
  
  //Busco dirs
  sprintf(cmd, "ls -F -1 -d %s | grep \"/\" > %s", argc[1], name);
  system(cmd);
  FILE *f=fopen(name, "rt");
  res=fgets(dir, 300, f);
  FILE *ftex=fopen("salida.tex", "wt");
  writelatexheader(ftex);

  while(res) {
    vector<string> names;
    printf("entering dir %s", dir);
    dir[strlen(dir)-1] = '\0';
    Matriz *m=processdir(dir, names);
    tolatex(m, names, ftex);
    res=fgets(dir, 300, f);
  }
  
  fprintf(ftex, "\\end{document}\n");
  fclose(ftex);
  fclose(f);
  remove(name);
  system("latex salida.tex");
  system("dvips -o salida.ps salida.dvi");
  system("ggv salida.ps &");

//  processfiles(argc[1], argc[2], argv>3 ? argc[3] : 0, argv>4 ? atoi(argc[4]) : 1);
}
void tolatex(Matriz *m, vector<string> &names, FILE *f)
{
  fprintf(f, "\\noindent\n");
  fprintf(f, "\\begin{tabular}[]{|cc|}\n");
  fprintf(f, "\\hline\n");
  printf("%d\n", names.size());
  printf("%s\n", names[names.size()-3].c_str());
  fprintf(f, "%s & \\\\\n", names[names.size()-3].c_str());
  fprintf(f, "train & test\\\\\n");
  fprintf(f, "\\mbox{\\epsfxsize=7.5cm\\epsfbox{%s}}&\n", names[names.size()-2].c_str());
  fprintf(f, "\\mbox{\\epsfxsize=7.5cm\\epsfbox{%s}}\\\\\n", names[names.size()-1].c_str());
  fprintf(f, "\\hline\n");
  fprintf(f, "\\end{tabular}\n");
  fprintf(f, "\n");
  fprintf(f, "\n");
  fprintf(f, "\\noindent\n");
  fprintf(f, "\\begin{tabular}[]{|lrrrrrr|}\n");
  fprintf(f, "&min ave.&ave.min.&train est.&aux1&aux2&aux3\\\\\n");
  fprintf(f, "\\hline\n");
  Matriz &n = *m;
  for (int i = 0; i < names.size()-3; i++) {
    fprintf(f, "%s & %0.1f[%0.1f] & %0.1f[%0.1f] & %0.1f[%0.1f]", names[i].c_str(),
                    n[i][0]*100, n[i][1], n[i][2]*100, n[i][3],  n[i][4]*100, n[i][5]);
    for (int j =0; j < 3; j++) {
      if (n[i][6+2*j]>0.0) fprintf(f, "& %0.1f[%0.1f]", n[i][6+j*2]*100, n[i][7+j*2]);
      else fprintf(f, "& -[-]");
    }
    fprintf(f, "\\\\\n");
  }
  fprintf(f, "\\hline\n");
  fprintf(f, "\\end{tabular}\n");
  fprintf(f, "\n");
}
void writelatexheader(FILE *f)
{
  fprintf(f, "\\documentclass[a4paper,12pt]{article}\n");
  fprintf(f, "\\usepackage{epsfig}\n");
  fprintf(f, "\\topmargin      =0.mm\n");
  fprintf(f, "\\oddsidemargin  =0.mm\n");
  fprintf(f, "\\evensidemargin =0.mm\n");
  fprintf(f, "\\headheight     =0.mm\n");
  fprintf(f, "\\headsep        =0.mm\n");
  fprintf(f, "\\textheight     =9.7in\n");
  fprintf(f, "\\textwidth      =6.3in\n");
  fprintf(f, "\n");
  fprintf(f, "\\pagestyle{empty}\n");
  fprintf(f, "\\input{psfig.sty}\n");
  fprintf(f, "\n");
  fprintf(f, "\\begin{document}\n");
}

Matriz *processdir(string dir, vector<string> &names)
{
  char test[300], train[300], info[300], cmd[300], name2[300], base[300];
  char *res;
  Matriz *meanstrain = 0;
  Matriz *meanstest = 0;
  Matriz *mret = new Matriz(10, 16);
  
  strcpy(name2, "./XXXXXX");
  mkstemp(name2);

  sprintf(cmd, "ls -1 -rt %s*test.txt > %s", dir.c_str(), name2);
  system(cmd);
  FILE *f2=fopen(name2, "rt");
  res=fgets(test, 300, f2);
  while(res) {
    printf("processing file %s", test);
    test[strlen(test)-1] = '\0';
    strncpy(base, test, strlen(test)-strlen("test.txt"));
    base[strlen(test)-strlen("test.txt")] = '\0';
    strcpy(train, base); strcat(train, "train.txt");
    strcpy(info, base); strcat(info, "info.txt");
    Matriz *mmean = processfiles(test, train, info, 5);
//    if (!mmeans) continue
    if (!meanstrain) {
      meanstrain = new Matriz(2, mmean->columnas());
      meanstest = new Matriz(2, mmean->columnas());
      for (int i = 0; i < mmean->columnas(); i++) {
        (*meanstrain)[0][i] = (*meanstest)[0][i] = i+1;
      }
    }
    else if(names.size()==meanstrain->filas()-1) {
      meanstrain->redim(meanstrain->filas()*2, mmean->columnas());
      meanstest->redim(meanstest->filas()*2, mmean->columnas());
    }
    meanstrain->sustituirFil(names.size()+1, mmean, 0);
    meanstest->sustituirFil(names.size()+1, mmean, 1);
    res=fgets(test, 300, f2);
    names.push_back(base+strlen(dir.c_str()));
    mret->sustituirFil(names.size()-1, mmean, 2);
  }
  fclose(f2);
  remove(name2);
  meanstrain->redim(names.size()+1, meanstrain->columnas());
  meanstest->redim(names.size()+1, meanstest->columnas());
  Matriz *m = meanstrain->transpuesta();
  oddify(m);
  if (dir[dir.size()-1]=='/') 
    dir.resize(dir.size()-1);
  print(m, names, dir+"tr", dir+"train");
  delete m;
  m = meanstest->transpuesta();
  oddify(m);
  print(m, names, dir+"ts", dir+"test");
  delete m;
  delete meanstrain, meanstest;  
  
  names.push_back(dir);
  names.push_back(dir+"tr.eps");
  names.push_back(dir+"ts.eps");
  
  return mret;
}
void oddify(Matriz *m)
{
  int i, j;
  for (i = 1, j=2; j < m->filas(); i++,j+=2) {
    m->sustituirFil(i, m, j);
  }
  if (m->filas()==j) m->sustituirFil(i, m, j-1);
  m->redim(i, m->columnas());
}
void print(Matriz *vals, vector<string> titles, string output, string title)
{
  char name[300], name2[300], cmd[100];
  int i;
  
  strcpy(name, "./XXXXXX");
  mkstemp(name);
  strcpy(name2, "./XXXXXX");
  mkstemp(name2);

  FILE *f=fopen(name, "wt");
  fprintf(f, "set xlabel \"Number of classifiers\"\n");
  fprintf(f, "set ylabel \"Error\"\n");
  fprintf(f, "set title \"%s\"\n", title.c_str());
  fprintf(f, "set term postscript eps\n");
  fprintf(f, "set output \"%s.eps\"\n", output.c_str());
  fprintf(f, "set data style lp\n");
  fprintf(f, "plot [1:%g] ", (*vals)[vals->filas()-1][0]);
  for (i = 0; i < titles.size()-1; i++) {
    fprintf(f, "\"%s\" using 1:%d t \"%s\" with lines, \\\n", name2, i+2, titles[i].c_str());
  }
  fprintf(f, "\"%s\" using 1:%d t \"%s\" with lines", name2, i+2, titles[i].c_str());
  fclose(f);
  
  vals->guardarAFichero(name2);
  sprintf(cmd, "gnuplot %s", name);
  system(cmd);
  remove(name);
  remove(name2);
}

Matriz *processfiles(char *test, char *train, char *aux, int naux)
{
     Matriz *mtr, *mts, *mmins, *maux=0;
     Matriz *mret;
     int *mins;
     double minmean, mintrain;

     mts = Matriz::leerDeFichero(test);
     mtr = Matriz::leerDeFichero(train);
     if (aux) maux = Matriz::leerDeFichero(aux);
   
     //Min con train
     mmins = minimos(mtr, true);
     Matriz *mmintrain = mean(mmins, false);
     mins = (int*)malloc(sizeof(int)*mmins->filas());
     for (int i = 0; i < mmins->filas(); i++) mins[i] = (int)(*mmins)[i][1];
     (*mmintrain)[0][0] = meanvals(mts, mins);
     free(mins);
     delete mmins;

     //media de los min
     mmins = minimos(mts, false);
     Matriz *mminmean = mean(mmins, false);
     delete mmins;

     //min de la media
     Matriz *mmean = mean(mts, false);
     Matriz *mmeanmin = minimos(mmean, false);

   //Meadis de train y test
   mret = new Matriz(3, mtr->columnas());
     mret->sustituirFil(1, mmean, 0);
     delete mmean;
     mmean = mean(mtr, false);
     mret->sustituirFil(0, mmean, 0);
     delete mmean;

     printf("Mean Min: %g[%0.1f] - Min. mean: %g[%0.1f] "
             "- Min train: %g[%0.1f]\n", 
             (*mmeanmin)[0][0], (*mmeanmin)[0][1]+1.0,
             (*mminmean)[0][0], (*mminmean)[0][1]+1.0,
             (*mmintrain)[0][0], (*mmintrain)[0][1]+1.0);

  (*mret)[2][0] = (*mmeanmin)[0][0];
  (*mret)[2][1] = (*mmeanmin)[0][1]+1.0;
  (*mret)[2][2] = (*mminmean)[0][0];
  (*mret)[2][3] = (*mminmean)[0][1]+1.0;
  (*mret)[2][4] = (*mmintrain)[0][0];
  (*mret)[2][5] = (*mmintrain)[0][1]+1.0;
    
      //Mean con aux
      if (maux) {
        double minaux;
        Matriz *mmean = mean(maux, false);
        mins = (int*)malloc(sizeof(int)*maux->filas());
        for (int i = 0; i < naux && i < maux->columnas(); i++) {
          for (int j = 0; j < maux->filas(); j++) mins[j] = (int)(*maux)[j][i]-1;
          minaux = meanvals(mts, mins);
          printf("Mean Aux%d: %g[%0.1f] - ", i+1, minaux, (*mmean)[0][i]);
      (*mret)[2][6+i*2] =minaux;
      (*mret)[2][7+i*2] = (*mmean)[0][i];
        }
        printf("\n");
        free(mins);
        delete mmean;
      }
      

       delete mts;
       delete mtr;
       delete mmintrain;
       delete mminmean;
       delete mmeanmin;

       return mret;
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
Matriz *mean(Matriz *v, bool vertical)
{
  Matriz *val = new Matriz(vertical ? v->filas() : 1, 
               vertical ? 1 : v->columnas());
  int i, j;
  if (vertical) {
    for(i=0;i<v->filas();i++){
      for(j=0;j<v->columnas();j++){
        (*val)[i][0] += (*v)[i][j];
      }
      (*val)[i][0] /= v->columnas();
    }
  }
  else {
    for(j=0;j<v->columnas();j++){
      for(i=0;i<v->filas();i++){
        (*val)[0][j] += (*v)[i][j];
      }
      (*val)[0][j] /= v->filas();
    }
  }
  return val;
}
/*void freematrix(Matriz *m) 
{
  int i;
  
  for(i=0;i<m->h;i++) free(m->m[i]);
  free(m);
}*/
/*Matriz *transpose(Matriz *m1) {
  int i, j;
  Matriz *m2;
  
  m2 = (Matriz*)malloc(sizeof(Matriz));
  m2->h = m1->w;
  m2->w = m1->h;
  m2->m = (double**)malloc(sizeof(double*)*m2->h);
  for(i=0;i<m2->h;i++) {
    m2->m[i] = (double*)malloc(sizeof(double)*m2->w);
    for(j=0;j<m2->w;j++) {
      m2->m[i][j] = m1->m[j][i];
    }
  }
  return m2;
}*/
double meanvals(Matriz *v, int *pos){
  double val=0.0;
  int i;
  for(i=0;i<v->filas();i++){
    val += (*v)[i][pos[i]];
  }
  return val/v->filas();
}

/*void writemat(char *filename, Matriz *m) 
{
  int j, i;
  
  FILE *f=fopen(filename, "wt");
  if (!f) return;

  for(i=0;i<m->h;i++) {
  for(j=0;j<m->w;j++) {
    fprintf(f, "\t%g", m->m[i][j]);
  }
  fprintf(f, "\n");
  }
  fclose(f);

  return;
}
Matriz *readmat(char *filename) {
  char leido[16384];
  char *start, *end;
  int cols, fils, i;
  Matriz *v = (Matriz*)malloc(sizeof(Matriz));
  
  FILE *f=fopen(filename, "rt");
  if (!f) return 0;

  v->h = v->w = 1;
  v->m = (double**)malloc(sizeof(double*));
  v->m[0] = (double*)malloc(sizeof(double));
  
  fils = 0;
  while (!feof(f)) {
    if (fils==v->h) {
      v->h += 10; 
      v->m = (double**)realloc(v->m, sizeof(double*)*v->h);
      for(i=fils;i<v->h;i++) 
        v->m[i] = (double*)malloc(sizeof(double)*v->w);
    }
    fscanf(f, "%s\t", leido);
    fgets(leido, 16384, f);
    if (feof(f)) break;
    start = leido;
    cols=0;
    while (strlen(start)!=1) {
      double val = strtod(start, &end);
      start = end;
      if (v->w==cols) {
        v->w += 10; 
        v->m[fils] = (double*)realloc(v->m[fils], sizeof(double)*v->w);
      }
      v->m[fils][cols] = val;
      cols++;
    };
    fils++;
  };
  
  fclose(f);

  v->h = fils;
  v->w = cols;

  return v;
}*/
