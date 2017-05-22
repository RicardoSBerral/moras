//---------------------------------------------------------------------------

#ifndef MatrizH
#define MatrizH

#include <stdio.h>
#include <exception>
#include<ostream>
//---------------------------------------------------------------------------
class Matriz
{
  protected:
    /**
     *  N�mero de filas de la matriz.
     */
    int numeroFilas;

    /**
     *  N�mero de columnas de la matriz.
     */
    int numeroColumnas;

    /**
     *  matriz de datos.
     */
    double** mat;


  protected:
    void construir(int fil, int col, double val = 0.0);

  public:

    /**
     *  Contruye un objeto Matriz de una dimensi�n determinada.
     *
     *  @param  numeroFilas N�mero de filas de la matriz a crear.
     *  @param  numeroColumnas N�mero de columnas de la matriz a crear.
     */
    Matriz(int dimension, double val = 0.0);
    Matriz(int numeroFilas, int numeroColumnas, double val = 0.0);
    Matriz(const Matriz &A);
    Matriz(Matriz *A, int numeroFilas, int numeroColumnas);

    ~Matriz();

    void Asignar(Matriz *A);
    void operator =(const Matriz &A);

    /**
     *  Indica si la matriz es cuadrada.
     */
    bool esCuadrada();

    /**
     *  Calcula la traza de la matriz. Si la matriz no es cuadrada
     *  devuelve NaN.
     */
    double traza();

    /**
     *  Devuelve la matriz transpoesta
     */
    Matriz *transpuesta();

    /**
     *  Devuelve el adjunto de la matriz de la fila fil y columan col
     */
    Matriz *adj(int fil, int col);

    /**
     *  Calcula el determinente de la matriz. Si la matriz no es cuadrada
     *  devuelve NaN.
     */
    double det();

    /**
     *  M�todo que te devuelve el valor de un elemento de la matriz.
     *
     *  @param fila no. de fila de la matriz de la que quieres obtener el valor
     *  @param columna no. de columna de la matriz de la que quieres obtener el valor
     *  @return el valor del elemento (fila,columna)
     */
    double valorElemento(int fila, int columna);
    double m(int i, int j);
    inline double * operator[](int fila){ return mat[fila]; }
    inline double** datos(){return mat;}
//    double operator ()(int fil, int col){ return valorElemento(fil, col);}
    /**
     *  M�todo que le asigna un valor a un elemento de la matriz.
     *
     *  @param fila no. de fila de la matriz.
     *  @param columna no. de columna de la matriz.
     *  @param nuevoValor Valor que se le asigna al elemento (fila,columna).
     */
    void valorElemento(int fila, int columna, double nuevoValor);

    /**
     *  M�todo que devuelve el n�mero de filas de la matriz.
     *
     *  @return el n�mero de filas de la matriz.
     */
    int numeroFilasMatriz();
    int filas(){return numeroFilasMatriz();}
    /**
     *  M�todo que devuelve el n�mero de columnas de la matriz.
     *
     *  @return el n�mero de columnas de la matriz.
     */
    int numeroColumnasMatriz();
    int columnas(){ return numeroColumnasMatriz();}

    /**
     *  M�todo que multiplica la matriz con otra y devuelve una nueva matriz
     *
     *  @return la matriz resultado de la multiplicaci�n.
     */
    Matriz &operator * (Matriz &m1) {return *(this->multiplicar(m1));}
    Matriz *multiplicar(Matriz &b);

    /**
     *  M�todo que suma la matriz con otra y devuelve una nueva matriz
     *
     *  @return la matriz resultado de la suma.
     */
    Matriz &operator + (Matriz &m1) {return *(this->sumar(m1));}
    Matriz *sumar(Matriz &b);
    Matriz &operator - (Matriz &m1) {return *(this->restar(m1));}
    Matriz *restar(Matriz &b);

    /**
     *  M�todo que lee los datos de una matriz de un fichero
     */
    static Matriz *leerDeFichero(char* nf);
    static Matriz *leerDeFichero(FILE *f);

    /**
     *  M�todo que guarda los datos de la matriz a un fichero
     */
    void guardarAFichero(char* nf, char *sep=0, char *fin=0);
    void guardarAFichero(FILE* f, char *sep=0, char *fin=0);
    void saveToStream(std::ostream &out);

    /**
     *  M�todo resulave una ec lineal
     *
     */
    static void ResolverEcLineal(Matriz *A, Matriz *Vars, Matriz *Coef);

    void sustituirCol(int col, Matriz *B, int colB);
    void sustituirFil(int fil, Matriz *B, int filB);
    void redim(int numeroFilas, int numeroColumnas);

   /**
    *  Mas cosas
    */
   Matriz *minimos(bool last);
   Matriz *maximos(bool last);
   Matriz *mean(bool horizontal, bool desv);
   Matriz *vals(Matriz *ixs);
};
//---------------------------------------------------------------------------
class MatrizVariable : public Matriz
{
  public:
    MatrizVariable(int dimension) : Matriz(dimension) {}
    MatrizVariable(int fil, int col) : Matriz(fil, col) {}
    virtual ~MatrizVariable(){}

    virtual void asignarValores(Matriz *O)=0;
};
//---------------------------------------------------------------------------
#endif
