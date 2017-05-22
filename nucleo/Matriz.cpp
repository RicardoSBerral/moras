//---------------------------------------------------------------------------

#include "Matriz.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <nrr.h>
#include <math.h>
#include <limits>


//---------------------------------------------------------------------------

using namespace std;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Matriz::Matriz(int dimension, double val)
{
  mat=0;
  construir(dimension, dimension, val);
}
Matriz::Matriz(int numeroFilas, int numeroColumnas, double val)
{
  mat=0;
  construir(numeroFilas, numeroColumnas, val);
}
Matriz::Matriz(const Matriz &A)
{
  mat=0;
  construir(A.numeroFilas, A.numeroColumnas);
  (*this) = A;
}
//Constructor tipo realloc
Matriz::Matriz(Matriz *A, int numeroFilas, int numeroColumnas)
{
  mat=0;
  construir(numeroFilas, numeroColumnas);
  if (A) (*this) = (*A);
}

Matriz::~Matriz()
{
  for(int i=0;i<numeroFilas;i++)
    free(mat[i]);//delete []mat[i];
  free(mat);//delete []mat;
}
void Matriz::construir(int fil, int col, double val)
{
  if (mat) this->~Matriz();
  mat = (double**)malloc (sizeof(double*)*fil);//new double*[fil];
  for(int i=0;i<fil;i++) {
    mat[i] = (double*)malloc(sizeof(double)*col);//new double[col];
    for(int j=0;j<col;j++){
      mat[i][j] = val;
    }
  }
  this->numeroFilas = fil;
  this->numeroColumnas = col;
}
void Matriz::Asignar(Matriz *A)
{
  if (A) (*this) = (*A);
}
void Matriz::operator =(const Matriz &A)
{
  int filas = A.numeroFilas    < numeroFilas    ? A.numeroFilas : numeroFilas;
  int colms = A.numeroColumnas < numeroColumnas ? A.numeroColumnas :
                                                                 numeroColumnas;
  for(int i=0;i<filas;i++) {
    for(int j=0;j<colms;j++) {
      mat[i][j] = A.mat[i][j];
    }
  }
}
/**
  *  Indica si la matriz es cuadrada.
  */
bool Matriz::esCuadrada()
{
  return (numeroFilas == numeroColumnas);
}

/**
  *  Calcula la traza de la matriz. Si la matriz no es cuadrada
  *  devuelve NaN.
  */
double Matriz::traza()
{
  if (!esCuadrada()) throw exception();//return Double.NaN; //
  double resultadoTraza=0;
  for(int i=0;i<numeroFilas;i++)
    resultadoTraza += mat[i][i];
  return resultadoTraza;
}
/**
  *  Devuelve la matriz transpoesta
  */
Matriz *Matriz::transpuesta()
{
        int i, j;
        Matriz *m2;

        m2 = new Matriz(numeroColumnas, numeroFilas);
        for(i=0;i<m2->filas();i++) {
                for(j=0;j<m2->columnas();j++) {
                        m2->mat[i][j] = this->mat[j][i];
                }
        }
        return m2;
}

/**
  *  Devuelve el adjunto de la matriz de la fila fil y columan col
  */
Matriz *Matriz::adj(int fil, int col)
{
  Matriz *m = new Matriz(numeroFilas-1, numeroColumnas-1);
  for(int i=0, iAdj=0;i<numeroFilas;i++) {
    if (i==fil) continue;
    for(int j=0, jAdj=0;j<numeroColumnas;j++) {
      if (j==col) continue;
      m->valorElemento(iAdj, jAdj, this->valorElemento(i, j));
      jAdj++;
    }
    iAdj++;
  }
  return m;
}
/**
  *  Calcula el determinente de la matriz. Si la matriz no es cuadrada
  *  devuelve NaN.
  */
double Matriz::det()
{
  if (!esCuadrada()) throw exception();//return Double.NaN; //
  double res = 0.0;
  if (numeroFilas==1)
    res = mat[0][0];
  else if (numeroFilas==2)
    res = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
  else {
    res = 0;
    for(int i=0;i<numeroFilas;i++){
      Matriz *a = adj(i, 0);
      res += ((i+1)%2-0.5)*2*mat[i][0]*a->det();
      delete a;
    }
  }
  return res;
}

/**
 *  M�todo que te devuelve el valor de un elemento de la matriz.
 *
 *  @param fila no. de fila de la matriz de la que quieres obtener el valor
 *  @param columna no. de columna de la matriz de la que quieres obtener el valor
 *  @return el valor del elemento (fila,columna)
 */
double Matriz::valorElemento(int fila, int columna){
  return mat[fila][columna];
}
double Matriz::m(int i, int j){
  return mat[i][j];
}
/*double * Matriz::operator[](int fila)
{
  return mat[fila];
}   */
//    double operator ()(int fil, int col){ return valorElemento(fil, col);}
/**
 *  M�todo que le asigna un valor a un elemento de la matriz.
 *
 *  @param fila no. de fila de la matriz.
 *  @param columna no. de columna de la matriz.
 *  @param nuevoValor Valor que se le asigna al elemento (fila,columna).
 */
void Matriz::valorElemento(int fila, int columna, double nuevoValor){
  mat[fila][columna] = nuevoValor;
}

/**
 *  M�todo que devuelve el n�mero de filas de la matriz.
 *
 *  @return el n�mero de filas de la matriz.
 */
int Matriz::numeroFilasMatriz(){
  return numeroFilas;
}
/**
 *  M�todo que devuelve el n�mero de columnas de la matriz.
 *
 *  @return el n�mero de columnas de la matriz.
 */
int Matriz::numeroColumnasMatriz(){
  return numeroColumnas;
}
/**
 *  Da una representaci�n en formadecadena de la matriz
 */
/*String toString(){
  String matriz = "";
  for(int i=0;i<numeroFilas;i++) {
    matriz = matriz + "( ";
    for(int j=0;j<numeroColumnas;j++) {
      matriz = matriz + this->valorElemento(i, j) + " ";
    }
    matriz = matriz + ")\n";
  }
  return matriz;
}*/
/**
 *  M�todo que multiplica la matriz con otra y devuelve una nueva matriz
 *
 *  @return la matriz resultado de la multiplicaci�n.
 */
Matriz *Matriz::multiplicar(Matriz &b){
  if (numeroColumnasMatriz()!=b.numeroFilasMatriz()) return 0;
  int filasC = numeroFilasMatriz();
  int columnasC = b.numeroColumnasMatriz();
  Matriz *c = new Matriz(filasC, columnasC);
  for(int i=0;i<filasC;i++) {
    for(int j=0;j<columnasC;j++) {
      double valor = 0.0;
      for(int k=0;k<numeroColumnasMatriz();k++) {
        valor += valorElemento(i, k)*b.valorElemento(k, j);
      }
      c->valorElemento(i, j, valor);
    }
  }
  return c;
}
/**
 *  M�todo que suma la matriz con otra y devuelve una nueva matriz
 *
 *  @return la matriz resultado de la suma.
 */
Matriz *Matriz::sumar(Matriz &b){
  if (numeroColumnasMatriz()!=b.numeroColumnasMatriz()) return 0;
  if (numeroFilasMatriz()!=b.numeroFilasMatriz()) return 0;
  int filasC = numeroFilasMatriz();
  int columnasC = b.numeroColumnasMatriz();
  Matriz *c = new Matriz(filasC, columnasC);
  for(int i=0;i<filasC;i++) {
    for(int j=0;j<columnasC;j++) {
      double valor = valorElemento(i, j) + b.valorElemento(i, j);
      c->valorElemento(i, j, valor);
    }
  }
  return c;
}
Matriz *Matriz::restar(Matriz &b){
  if (numeroColumnasMatriz()!=b.numeroColumnasMatriz()) return 0;
  if (numeroFilasMatriz()!=b.numeroFilasMatriz()) return 0;
  int filasC = numeroFilasMatriz();
  int columnasC = b.numeroColumnasMatriz();
  Matriz *c = new Matriz(filasC, columnasC);
  for(int i=0;i<filasC;i++) {
    for(int j=0;j<columnasC;j++) {
      double valor = valorElemento(i, j) - b.valorElemento(i, j);
      c->valorElemento(i, j, valor);
    }
  }
  return c;
}
/**
  *  M�todo que lee los datos de una matriz de un fichero
  */
Matriz* Matriz::leerDeFichero(char* filename)
{
  FILE *f=fopen(filename, "rt");
  if (f) {
    Matriz *m = Matriz::leerDeFichero(f);
    fclose(f);
    return m;
  }
  return 0;
}
Matriz* Matriz::leerDeFichero(FILE *f)
{
  char leido[1638400];
  char aux[1638400];
  char *start, *end;
  int cols, fils, i;
  Matriz *v;

  if (!f) return 0;

  v = new Matriz(1, 1);

/*  v->numeroFilas = v->numeroColumnas = 1;
  v->mat = (double**)malloc(sizeof(double*));
  v->mat[0] = (double*)malloc(sizeof(double));
*/
  cols = fils = 0;
  while (!feof(f)) {
    if ( fils == v->numeroFilas ) {
      v->numeroFilas += 10;
      v->mat = (double**)realloc(v->mat, sizeof(double*)*v->numeroFilas);
      for(i=fils;i<v->numeroFilas;i++)
        v->mat[i] = (double*)malloc(sizeof(double)*v->numeroColumnas);
    }
//    fscanf(f, "%s\t", leido);
    fgets(leido, 1638400, f);
    if (feof(f)) break;
    start = leido;
    cols=0;
    while (strlen(start)>1) {
      double val = strtod(start, &end);
      if (start==end) {
        aux[0] = '\0';
        if (0==sscanf(start, "%s \t", aux) || strlen(aux)==0) start[0]='\0';
        else start = start + strlen(aux);
        continue;
      }
      start = end;
      if (v->numeroColumnas==cols) {
        if (fils>0) break;

        v->numeroColumnas += 10;
//        double *fila = new double[v->numeroColumnas];
//        for (int ii=0;ii<cols;ii++) fila[ii] = v->mat[fils][ii];
//        delete [] v->mat[fils];
//        v->mat[fils] = fila;
        v->mat[fils] = (double*)realloc(v->mat[fils], 
                                            sizeof(double)*v->numeroColumnas);
      }
      v->mat[fils][cols] = val;
      cols++;
    };
    fils++;
  };

//  fclose(f);
//
  v->numeroFilas = fils;
  v->numeroColumnas = cols;

  return v;

}

/**
  *  Método que guarda los datos de la matriz a un fichero
  */

void Matriz::guardarAFichero(char* nf, char *sep, char *fin)
{
  FILE *f=fopen(nf, "wt");
  guardarAFichero(f, sep, fin);
  fclose(f);
}
void Matriz::guardarAFichero(FILE* f, char *sep, char *fin)
{
  if (!f) return;
  if (!sep) sep = (char*)"\t";
  if (!fin) fin = "";
  for(int i=0;i<numeroFilas;i++) {
    for(int j=0;j<numeroColumnas;j++) {
      fprintf(f, "%.10g%s", mat[i][j], j < numeroColumnas-1 ? sep : fin);
    }
    fprintf(f, "\n");
  }
}
void Matriz::saveToStream(std::ostream &out)
{
  int oldp = out.precision(15);
  for(int i=0;i<numeroFilas;i++) {
    for(int j=0;j<numeroColumnas;j++) {
      out << mat[i][j] << "\t";
    }
    out << endl;
  }
  out.precision(oldp);
}

/**
 *  M�todo resulave una ec lineal
 *
 */
void Matriz::ResolverEcLineal(Matriz *A, Matriz *Vars, Matriz *Coef){
  if (A->numeroColumnas!=A->numeroFilas) return;
  double **a = new double*[A->numeroFilas+1];
  int *indx = new int[A->numeroFilas+1];
  double *d  = new double[A->numeroFilas+1];
  double *b  = new double[A->numeroFilas+1];
  for(int i=0;i<A->numeroFilas;i++) {
    a[i+1] = new double[A->numeroColumnas+1];
    b[i+1] = Coef->mat[i][0];
    for(int j=0;j<A->numeroColumnas;j++) {
      a[i+1][j+1] = A->mat[i][j];
    }
  }
  nr::ludcmp(a, A->numeroFilas, indx, d);
  nr::lubksb(a, A->numeroFilas, indx, b);
  for(int i=0;i<A->numeroFilas;i++) {
    delete [] a[i+1];
    Vars->mat[i][0] = b[i+1];
  }
  delete []b;
  delete []d;
  delete []indx;
  delete []a;
}

void Matriz::sustituirCol(int col, Matriz *B, int colB)
{
  for(int i=0;i<numeroFilas;i++)
  mat[i][col] = B->mat[i][colB];
}

void Matriz::sustituirFil(int fil, Matriz *B, int filB)
{
  for(int i=0;i<numeroColumnas;i++)
  mat[fil][i] = B->mat[filB][i];
}
void Matriz::redim(int numeroFilas, int numeroColumnas)
{
  Matriz *aux = new Matriz(*this);
    construir(numeroFilas, numeroColumnas);
  (*this) = (*aux);
  delete aux;
}
/**
 *  Mas cosas
 */

Matriz *Matriz::maximos(bool last) 
{
  int i, j, imax;
  double max;
  Matriz *maxs = new Matriz(this->filas(), 2);

  for(i=0;i<this->filas();i++) {
    max = (*this)[i][0];
    imax = 0;
    for(j=1;j<this->columnas();j++) {
      if ((max<=(*this)[i][j] && last) || max<(*this)[i][j]) {
        max = (*this)[i][j];
        imax = j;
      }
    }
    (*maxs)[i][0] = max;
    (*maxs)[i][1] = imax;
  }

  return maxs;
}
Matriz *Matriz::minimos(bool last)
{
  int i, j, imin;
  double min;
  Matriz *mins = new Matriz(this->filas(), 2);

  for(i=0;i<this->filas();i++) {
    min = (*this)[i][0];
    imin = 0;
    for(j=1;j<this->columnas();j++) {
      if ((min>=(*this)[i][j] && last) || min>(*this)[i][j]) {
        min = (*this)[i][j];
        imin = j;
      }
    }
    (*mins)[i][0] = min;
    (*mins)[i][1] = imin;
  }

  return mins;
}
Matriz *Matriz::mean(bool horizontal, bool desv) 
{
  int f, c, N;

  if (horizontal) {
    f = desv ? 2 : 1;
    c = this->columnas();
    N = this->filas();
  }
  else {
    f = this->filas();
    c = desv ? 2 : 1;
    N = this->columnas();
  }
  Matriz *val = new Matriz(f, c);
  
  int i, j;
  if (!horizontal) {
    for(i=0;i<this->filas();i++){
      for(j=0;j<this->columnas();j++){
        (*val)[i][0] += (*this)[i][j];
        if (desv) (*val)[i][1] += ((*this)[i][j]*(*this)[i][j]);
      }
      if (desv) {
        double aux = (*val)[i][1]*N - (*val)[i][0]*(*val)[i][0];
        (*val)[i][1] = aux < 0.0 ? 0.0 : sqrt(aux/((double)N*(N-1)));
      }
      (*val)[i][0] /= N;
    }
  }
  else {
    for(j=0;j<this->columnas();j++){
      for(i=0;i<this->filas();i++){
        (*val)[0][j] += (*this)[i][j];
        if (desv) (*val)[1][j] += ((*this)[i][j]*(*this)[i][j]);
      }
      if (desv) {
        double aux = (*val)[1][j]*N - (*val)[0][j]*(*val)[0][j];
        (*val)[1][j] = aux < 0.0 ? 0.0 : sqrt(aux/((double)N*(N-1)));
      }
      (*val)[0][j] /= N;
    }
  }
  return val;
}

Matriz *Matriz::vals(Matriz *ixs) 
{
  Matriz *val = new Matriz(this->filas(), 1);
  int cols = this->columnas();

  for(int i=0;i<this->filas();i++) 
    (*val)[i][0] = (*ixs)[i][0]>=cols || (*ixs)[i][0]<0 ? 
        numeric_limits<double>::quiet_NaN() : (*this)[i][(int)(*ixs)[i][0]];

  return val;
}
