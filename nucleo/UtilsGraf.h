//---------------------------------------------------------------------------

#ifndef UtilsGrafH
#define UtilsGrafH

//---------------------------------------------------------------------------
#include <vector>
#include <string>
#include <limits>
//---------------------------------------------------------------------------
#include "Classifier.h" 
//---------------------------------------------------------------------------
float f(float p[]);
void df(float p[], float g[]);
//---------------------------------------------------------------------------
#ifndef WIN32
//#pragma push option -a1
#define DWORD int
#define LONG long
#define WORD short
#define BYTE char
#define BI_RGB 0
typedef struct tagRGBQUAD {
  BYTE    rgbBlue; 
  BYTE    rgbGreen; 
  BYTE    rgbRed; 
  BYTE    rgbReserved; 
} RGBQUAD; 
typedef struct tagBITMAPINFOHEADER{
  DWORD  biSize; 
  LONG   biWidth; 
  LONG   biHeight; 
  WORD   biPlanes;
  WORD   biBitCount; 
  DWORD  biCompression; 
  DWORD  biSizeImage; 
  LONG   biXPelsPerMeter; 
  LONG   biYPelsPerMeter; 
  DWORD  biClrUsed; 
  DWORD  biClrImportant; 
} BITMAPINFOHEADER, *PBITMAPINFOHEADER;
typedef struct tagBITMAPINFO { 
  BITMAPINFOHEADER bmiHeader; 
  RGBQUAD          bmiColors[1]; 
} BITMAPINFO, *PBITMAPINFO; 
typedef struct tagBITMAPFILEHEADER { 
  WORD    bfType; 
  DWORD   bfSize; 
  WORD    bfReserved1; 
  WORD    bfReserved2; 
  DWORD   bfOffBits; 
} BITMAPFILEHEADER, *PBITMAPFILEHEADER; 
#endif
//---------------------------------------------------------------------------
class Classifier;
class Data;
//---------------------------------------------------------------------------

class Function
{
  protected:
    int dim;
  
  public:
//    Function() {dim = 2;}
    Function(int dim=2) {this->dim = dim;}
    virtual ~Function() {}

    virtual double f(std::vector<double> xs) {
      double *x = new double[dim];
      for (int i=0;i<dim;i++) x[i] = xs[i];
      double ret = f(x);
      delete []x;
      return ret;
    }
    
    //Aquellas funciones que puedan ser optimizadas que
    //redefinan esta funcion
    virtual void f(double **xs, int n, double *res){
      for(int i=0;i<n;i++) {
        res[i] = f(xs[i]);
      }
    }
    virtual double f(double *xs)=0;

    inline int GetDim() { return dim; }
    //Just try to indicate limits of the fuction, not necessary
    virtual double max(){return std::numeric_limits<double>::infinity();}
    virtual double min(){return -std::numeric_limits<double>::infinity();}
};
class Recta : public Function
{
  double a, b;
  public:

    Recta(double _a, double _b) : Function(1) {a=_a; b=_b;}
    
    virtual double f(double *xs) {
      return a*(xs[0]) + b;
    }
};
class ComposedFunction : public Function
{
  protected:
    std::vector<Function*> funs;
    std::vector<double> des;
    std::vector<double> as;
  
  public:
    ComposedFunction() : Function(1) {}
    virtual ~ComposedFunction() {}

  public:
    virtual double f(double *xs) {
      for (unsigned i=0;i<funs.size();i++) {
        if(*xs>=des[i] && *xs<as[i]) 
          return funs[i]->f(xs);
      }
      return 0.0;
    }
    
    void AddFunction(Function *f, double de, double a){
   //   if (f->dim!=this->dim) return;
      funs.push_back(f);
      des.push_back(de);
      as.push_back(a<de ? std::numeric_limits<double>::infinity() : a);
    }
    void AddFunction(Function *f, double a){
      AddFunction(f, -std::numeric_limits<double>::infinity(), a);
    }
};
class ProbFunction : public Function
{
  public:
    ProbFunction(int dim=2) : Function(dim) {}

    virtual double max(){return 1.0;}
    virtual double min(){return 0.0;}
};
class Gaussiana : public ProbFunction
{
  int dim;
  double **Sigma;
  double Inv[2][2];
  double *despl;
  double constante;

  public:
    Gaussiana(int _dim, double **mcovarianza=0, double *_despl=0);
    virtual ~Gaussiana();

    virtual double f(double *xs);
};
class GaussianaMultiple : public ProbFunction
{
  std::vector<Gaussiana *> gaussianas;
  std::vector<double> pesos;
  double sumapesos;

  public:
    GaussianaMultiple(){sumapesos = 0.0;}

    void AnadirGaussiana(Gaussiana *g, double peso);
    virtual double f(double *xs); 
};
class Gaussiana4 : public ProbFunction
{
  double despl;
  Gaussiana *g1;
  Gaussiana *g2;
  Gaussiana *g3;
  Gaussiana *g4;

  public:
    Gaussiana4(double _despl);

    virtual double f(double *xs);
};
class Escalon : public ProbFunction
{
  public:
    Escalon() {
    }
    virtual double f(double *xs);
};
class Banana : public ProbFunction
{
  double despl;
  Gaussiana *g;
  Escalon *e;

  public:
    Banana(double _despl); 

    virtual double f(double *xs);
};
class BayesClassifier : public Classifier
{
  std::vector<double> APrioriProb;
  std::vector<ProbFunction*> ProbFuns;

public:
  BayesClassifier(std::vector<double> _APrioriProb, 
                                        std::vector<ProbFunction*> _ProbFuns);
  virtual ~BayesClassifier(){;}

  ProbFunction* GetProbFun(int i);
  double GetAPrioriProb(int i);
  int Classes();

  virtual int Classify(int ElementIndex);
  int Classify(double *ejemplo); 
  virtual void Classify(double **dat, int NumeroDatos);

  virtual void Build(Data *data, FuncionDeProgreso *fp=0);
};
class BayesErrorFunction : public ProbFunction
{
  BayesClassifier *bc;
  Classifier *c;
  public:
  BayesErrorFunction(Classifier *_c, BayesClassifier *_bc);

  virtual double f(double *xs);

};
class MarginDistribution : public ComposedFunction
{
  protected:
    void Build(Classifier *_e, Data *dat, bool interpol);

  public:
    MarginDistribution(Classifier *_e, Data *dat, bool inter) : 
                      ComposedFunction()
    {
      Build(_e, dat, inter);
    }

  public:

    double fm(double x) {return f(&x);}

    virtual double max(){return 1.0;}
    virtual double min(){return 0.0;}
};
class MargenClasfFunction : public Function
{
  Classifier *e;
  public:
  MargenClasfFunction(Classifier *_e, int dim=2);

  virtual void f(double **xs, int n, double *res);
  virtual double f(double *xs);
  virtual double max(){return 1.0;}
  virtual double min(){return -1.0;}
};
class MargenFunction : public Function
{
  Classifier *e;
  public:
  MargenFunction(Classifier *_e, int dim=2);

  virtual void f(double **xs, int n, double *res);
  virtual double f(double *xs);
  virtual double max(){return 1.0;}
  virtual double min(){return .5;}
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void DoMapaDeClasif2D(Classifier *c, int x1, int x2,
                                                double x1min, double x2min,
                                                double x1max, double x2max,
                                                int x1div, int x2div,
                                               std::vector<double> vector_const,
                                                std::string nomfichsal);
void DoMapaDeAlturas2D(Function *f, int x1, int x2, double x1min, double x2min,
                                                double x1max, double x2max,
                                                int x1div, int x2div,
                                               std::vector<double> vector_const,
                                                std::string nomfichsal,
                                              double zmin=2.0, double zmax=1.0);
//---------------------------------------------------------------------------
#endif
