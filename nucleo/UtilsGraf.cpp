//---------------------------------------------------------------------------


#include "UtilsGraf.h"

#include <math.h>
#include <fstream>
#include <algorithm>
#include "Classifier.h"
#include "FnsClasificacion.h"
#include "Ensemble.h"
#include "FnsOrdClas.h"
using namespace std;

Gaussiana::Gaussiana(int _dim, double **mcovarianza, double *_despl) 
                             : ProbFunction(_dim) 
{
  dim = _dim;

  //Reserva memoria
  despl = new double[dim];
  Sigma = new double*[dim];
  for (int i=0;i<dim;i++) Sigma[i] = new double[dim];

  //Inicialización
  for (int i=0;i<dim;i++) { 
    despl[i] = _despl ? _despl[i] : 0.0;
  }

  if (mcovarianza) {
    for (int i=0;i<dim;i++) {
      for (int j=0;j<dim;j++) {
        Sigma[i][j] = mcovarianza[i][j];
      }
    }
  }
  else {
    for (int i=0;i<dim;i++) {
      for (int j=0;j<dim;j++) {
        Sigma[i][j] = i==j ? 1.0: 0.0;
      }
    }
  }

  constante = pow(2.0*3.1415926535, (double)dim/2.0);
  if (Sigma) {
    double det = Sigma[0][0]*Sigma[1][1] - Sigma[1][0]*Sigma[0][1];
    constante *= pow(det, 0.5);
    Inv[0][0] = Sigma[1][1]/det;
    Inv[0][1] = -Sigma[0][1]/det;
    Inv[1][0] = -Sigma[1][0]/det;
    Inv[1][1] = Sigma[0][0]/det;
  }
  else {
    Inv[0][0] = 1.0;
    Inv[0][1] = 0.0;
    Inv[1][0] = 0.0;
    Inv[1][1] = 1.0;
  }
}
Gaussiana::~Gaussiana() 
{
  delete []despl;
  for (int i=0;i<dim;i++) delete []Sigma[i];
  delete []Sigma;
}
double Gaussiana::f(double *xs) 
{
  //de momento solo 2D y sin sigma
  double val;
  double dif[2];
  dif[0] = (xs[0]-despl[0]);
  dif[1] = (xs[1]-despl[1]);
  val = exp(-0.5*(dif[0]*Inv[0][0]*dif[0]+dif[1]*Inv[1][0]*dif[0]+
                  dif[0]*Inv[0][1]*dif[1]+dif[1]*Inv[1][1]*dif[1]));
  return val/constante;
}
//----------------------------------------------------------------------------
void GaussianaMultiple::AnadirGaussiana(Gaussiana *g, double peso) 
{
  gaussianas.push_back(g);
  pesos.push_back(peso);
  sumapesos += peso;
}
double GaussianaMultiple::f(double *xs) {
  double res = 0.0;
  for(unsigned i=0;i<gaussianas.size();i++) {
    res += pesos[i]*gaussianas[i]->f(xs);
  }
  return res/sumapesos;
}
//----------------------------------------------------------------------------
Gaussiana4::Gaussiana4(double _despl) {
  despl = _despl;
  double *dg1 = new double[2];
  double *dg2 = new double[2];
  double *dg3 = new double[2];
  double *dg4 = new double[2];
  dg1[0] = dg1[1] = dg2[0] = dg3[1] = despl;
  dg4[0] = dg4[1] = dg2[1] = dg3[0] = -despl;
  g1 = new Gaussiana(2, 0, dg1);
  g2 = new Gaussiana(2, 0, dg2);
  g3 = new Gaussiana(2, 0, dg3);
  g4 = new Gaussiana(2, 0, dg4);
}
double Gaussiana4::f(double *xs) 
{
  //de momento solo 2D y sin sigma
  return (0.25*g1->f(xs) + 0.25*g2->f(xs) + 0.25*g3->f(xs) + 0.25*g4->f(xs));
}
//----------------------------------------------------------------------------
double Escalon::f(double *xs) {
  for (int i=0;i<dim;i++)//??????????
    if (xs[0]<=0 || xs[0]>=1) return 0.0;
  return 1.0;
}
//----------------------------------------------------------------------------
Banana::Banana(double _despl) 
{
  despl = _despl;
  double *dg1 = new double[2];
  dg1[0] = dg1[1] = 0;
  g = new Gaussiana(2, 0, dg1);
  delete []dg1;
  e = new Escalon();
}
double Banana::f(double *xs) 
{
  //de momento solo 2D y sin sigma
/*    xn=sqrt(0.3)*(0.1*g->f(xs)+5*e->f(xs)-2.5)/3+[1.45;-0.6] ;
  xn(2)=xn(2)-abs(xn(1)-1.45)^2*4 ;*/
  return 0.0;
}
//----------------------------------------------------------------------------

BayesClassifier::BayesClassifier(std::vector<double> _APrioriProb, 
                                  std::vector<ProbFunction*> _ProbFuns) 
{
  APrioriProb = _APrioriProb;
  ProbFuns = _ProbFuns;
}

ProbFunction* BayesClassifier::GetProbFun(int i) 
{
  return ProbFuns[i];
}
double BayesClassifier::GetAPrioriProb(int i) 
{
  return APrioriProb[i];
}
int BayesClassifier::Classes() 
{
  return (int)ProbFuns.size();
}

void BayesClassifier::Classify(double **dat, int NumeroDatos)
{
  double d[1024];
  int dim = APrioriProb.size();

  for(int i=0;i<NumeroDatos;i++) {
    for(int j=0; j<dim;j++) {
      d[j] = dat[i][j];
    }
    dat[i][dim] = Classify(d);
  }
}
int BayesClassifier::Classify(double *ejemplo) 
{
  std::vector<double> dat;
  int dim = APrioriProb.size();

  for(int j=0; j<dim;j++)
    dat.push_back(ejemplo[j]);

  double val = -1;
  int i_probmax = 0;
  for (unsigned i=0;i<APrioriProb.size();i++) {
    double val_i = APrioriProb[i] * ProbFuns[i]->f(dat);
    if (val_i>val) {
      val = val_i;
      i_probmax = i;
    }

  }

  return i_probmax;
}
int BayesClassifier::Classify(int ElementIndex) 
{
  double dat[1024];
  int NCols  = data->GetNumVar() - 1;
  for(int j=0; j<NCols;j++) {
    dat[j] = data->GetValueVar(ElementIndex, j);
  }
  return Classify(dat);
/*  double val = -1;
  int i_probmax = 0;
  for (unsigned i=0;i<APrioriProb.size();i++) {
    double val_i = APrioriProb[i] * ProbFuns[i]->f(dat);
    if (val_i>val) {
      val = val_i;
      i_probmax = i;
    }

  }

  return i_probmax;*/
}
void BayesClassifier::Build(Data *data, FuncionDeProgreso *fp) 
{
  Init(data);
}
//----------------------------------------------------------------------------
BayesErrorFunction::BayesErrorFunction(Classifier *_c, BayesClassifier *_bc) 
{
  bc = _bc;
  c = _c;
}
double BayesErrorFunction::f(double *xs) 
{
  double *x = new double[3];
  double val=0.0;
  x[0] = xs[0]; x[1] = xs[1]; x[2] = 0;
    ((Classifier*)bc)->Classify(&x, 1);
  c->Classify(&x, 1);
  int clase = (int)x[2];
  for(int i=0;i<bc->Classes();i++) {
    if (clase==i) continue;
    val += bc->GetAPrioriProb(i)*bc->GetProbFun(i)->f(xs);
  }
  delete []x;
  return val;
}
//---------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
void MarginDistribution::Build(Classifier *clf, Data *dat, bool interpol)
{
  int N = dat->GetNTotal();
  vector<double> mrgs(N, 0.0);

  //Calculamos el margen de cada dato
  clf->SetData(dat);
  for(int i=0;i<N;i++) {
    vector<double> dis = clf->Distribution(i);
    int clase_real = dat->GetDatClass(i);
    double maxmalo = 0.0;
    for(int j=0;j<(int)dis.size();j++) {
      if (j==clase_real) continue;
      if (dis[j]>maxmalo) maxmalo = dis[j];
    }
    mrgs[i] = dis[clase_real] - maxmalo;
printf("%g\n", mrgs[i]);
  }

  //ordenamos y generamos las rectas
  sort(mrgs.begin(), mrgs.end());
  double de = mrgs[0]==-1.0 ? -100.0 : -1.0;
  double y = 0.0;
  for(int i=0;i<N;i++) {
    int ini = i;
    while(i<N-1 && mrgs[i]==mrgs[i+1]) i++;
    if (interpol) {
      double a = (1.0+i-ini)/((mrgs[ini]-de)*N);
      double b = y - de*a;
      AddFunction(new Recta(a, b), de, mrgs[ini]);
//printf("(%g, %g, %g, %g)\t", a, b, de, y);
    }
    else {
      AddFunction(new Recta(0.0, y), de, mrgs[ini]);
    }
    de = mrgs[ini];
    y += (1.0+i-ini)/N;
  }
  AddFunction(new Recta(0.0, y), de, numeric_limits<double>::infinity());
//printf("(%g %g)\t", de, y);
//printf("\n");

}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
MargenClasfFunction::MargenClasfFunction(Classifier *_e, int dim):Function(dim)
{
  e = _e;
}
void MargenClasfFunction::f(double **xs, int n, double *res)
{
  double val=0.0;
  int Clase;

  Data *holdData = e->GetData();
  Data *dat = holdData->Clone(1);
  dat->redim(n);
  for(int i=0;i<n;i++) {
    dat->SetValueVar(i,0,xs[i][0]);
    dat->SetValueVar(i,1,xs[i][1]);
  }
  e->SetData(dat);
  for(int i=0;i<n;i++) {
    val = e->ClassificationCertainty(i, Clase);
    res[i] = Clase==0 ? -val : val;
  }

  e->SetData(holdData);
  delete dat;
}
double MargenClasfFunction::f(double *xs)
{
  double val=0.0;
  int Clase;
  Data *holdData = e->GetData();
  Data *dat = holdData->Clone(1);
  dat->SetValueVar(0,0,xs[0]);dat->SetValueVar(0,1,xs[1]);
  e->SetData(dat);
  val = e->ClassificationCertainty(0, Clase);
  if (Clase==0) val = -val;

  e->SetData(holdData);
  delete dat;

  return val; 
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
MargenFunction::MargenFunction(Classifier *_e, int dim):Function(dim)
{
  e = _e;
}
void MargenFunction::f(double **xs, int n, double *res)
{
  double val=0.0;
  int Clase;

  Data *holdData = e->GetData();
  Data *dat = holdData->Clone(1);
  dat->redim(n);
  for(int i=0;i<n;i++) {
    dat->SetValueVar(i,0,xs[i][0]);
    dat->SetValueVar(i,1,xs[i][1]);
  }
  e->SetData(dat);
  for(int i=0;i<n;i++) {
    val = e->ClassificationCertainty(i, Clase);
    res[i] = val;
  }

  e->SetData(holdData);
  delete dat;
}
double MargenFunction::f(double *xs)
{
//  double *x = new double[3];
  double val=0.0;
  int Clase;
  Data *holdData = e->GetData();
  Data *dat = holdData->Clone(1);
//  x[0] = dat->GetValueVar(0,0); x[1] = dat->GetValueVar(0,1);
  dat->SetValueVar(0,0,xs[0]);dat->SetValueVar(0,1,xs[1]);
  e->SetData(dat);
  val = e->ClassificationCertainty(0, Clase);
//  if (Clase==0) val = -val;
  //val = ((Ensemble*)e)->Margen(0, 0)[0];
//  dat->SetValueVar(0,0,x[0]);dat->SetValueVar(0,1,x[1]);

  e->SetData(holdData);
  delete dat;

//  delete []x;
  return val; 
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
class InversaFunction : public ProbFunction
{
  public:
  InversaFunction() {
  }
  virtual double f(double *xs) {
    return 1.0/xs[0];
  }
};
class XorClassifier : public Classifier
{
public:
  XorClassifier() :Classifier(){ }
  virtual ~XorClassifier(){;}
  virtual int Classify(int ElementIndex) {
      bool x1 = data->GetValueVar(ElementIndex, 0) > 0.0;
      bool x2 = data->GetValueVar(ElementIndex, 1) > 0.0;
      return x1 && x2 ? 1 : 0;
  }
  virtual void Build(Data *data){ Init(data); }
};
double Integrate(ProbFunction *f, vector<double> &constvector, double min,
                                                   double max, int div, int var)
{
  double paso = (max-min)/div;
  double val=0.0, ant;
  double valvar = constvector[var];

  constvector[var] = min;
  ant = f->f(constvector);
  for (int i=1;i<=div;i++) {
    val += ant;
    constvector[var] = min + paso*i;
    ant = f->f(constvector);
    val += ant;
  }
  constvector[var] = valvar;
  return paso*val/2.0;
}
double HiperIntegrate(ProbFunction *f, vector<double> &cte, vector<double> &min,
                                          vector<double> &max, vector<int> &div)
{
  if (min.size()==1) {
     return Integrate(f, cte, min[0], max[0], div[0], 0);
  }

  double valz = 0.0, antz;
  int dim = (int)min.size();
  double pasoz = (max[dim-1]-min[dim-1])/div[dim-1];
  int divz = div[dim-1]; div.pop_back();
  int minz = (int)min[dim-1]; min.pop_back();
  int maxz = (int)max[dim-1]; max.pop_back();

  //cte[dim-1] = minz;
  antz += HiperIntegrate(f, cte, min, max, div);
  for (int k=1;k<=divz;k++) {
    valz += antz;
    cte[dim-1] = minz + pasoz*k;
    antz += HiperIntegrate(f, cte, min, max, div);
    valz += antz;
  }
  div.push_back(divz);
  min.push_back(minz);
  max.push_back(maxz);

  return pasoz*valz/2.0;
}
double HiperIntegrateEsp(ProbFunction *f, vector<double> &min, vector<double> &max,
                                                               vector<int> &div)
{
  vector<double> cte = min;
  return HiperIntegrate(f, cte, min, max, div);
}


//---------------------------------------------------------------------------
float f(float p[])
{//y=10*sin(x)+x.*x;

  return 10*sin(p[1]) + p[1]*p[1];
}
void df(float p[], float g[])
{//y'=(10*sin(x)+x.*x)' = 10*cos(x) + 2*x

  g[1] = 10*cos(p[1]) + 2*p[1];
}
//---------------------------------------------------------------------------
void DoMapaDeClasif2D(Classifier *c, int x1, int x2,
                                                double x1min, double x2min,
                                                double x1max, double x2max,
                                                int x1div, int x2div,
                                                vector<double> vector_const,
                                                string nomfichsal)
{
  double div1 = (x1max-x1min)/x1div;
  double div2 = (x2max-x2min)/x2div;
  int nc = c->GetData() ? ((NomData*)c->GetData())->NumClass : 2;
  unsigned char v     = (nc<3 ? 1 : (nc<17 ? 4 : 8));
  int ncolores = (int)pow(2.0, v);
  int AnchoEx = ((3 + (7  + x1div*v)/8)/4)*4;
  unsigned char *fila = new unsigned char[AnchoEx];

  BITMAPFILEHEADER bfh;
  BITMAPINFOHEADER bih;
  RGBQUAD color;
  int sizeofbfh = 14;  //sizeof(bfh)
  
  bfh.bfType      = 0x4d42;
  bfh.bfSize      = sizeof(bih) + sizeofbfh + AnchoEx*x2div + ncolores*4;
  bfh.bfReserved1 = 0;
  bfh.bfReserved2 = 0;
  bfh.bfOffBits   = sizeof(bih) + sizeofbfh + ncolores*4;

  bih.biSize          = sizeof(bih);
  bih.biWidth         = x1div;
  bih.biHeight        = x2div;
  bih.biPlanes        = 1;
  bih.biBitCount      = v;
  bih.biCompression   = BI_RGB;
  bih.biSizeImage     = bfh.bfSize - bfh.bfOffBits;
  bih.biXPelsPerMeter = 100;
  bih.biYPelsPerMeter = 100;
  bih.biClrUsed       = 0;
  bih.biClrImportant  = 0;

  FILE *bmp=fopen(nomfichsal.c_str(), "wb");
  if (!bmp) return;

//  fwrite(&bfh, sizeof(bfh), 1, bmp);
  fwrite(&bfh.bfType, 2, 1, bmp);
  fwrite(&bfh.bfSize, 4, 1, bmp);
  fwrite(&bfh.bfReserved1, 2, 1, bmp);
  fwrite(&bfh.bfReserved1, 2, 1, bmp);
  fwrite(&bfh.bfOffBits, 4, 1, bmp);
  fwrite(&bih, sizeof(bih), 1, bmp);

  //Paleta
  double a = 255.0/(nc-1);
  for (int i=0;i<ncolores;i++) {
    color.rgbBlue     = (char) (a*i);
    color.rgbGreen    = (char) (a*i);
    color.rgbRed      = (char) (a*i);
    color.rgbReserved = 0x00;
    fwrite(&color, sizeof(color), 1, bmp);
  }

  double **dato = new double*[x1div];
  for (int j=0;j<x1div;j++) {
    dato[j] = new double[vector_const.size() + 1];
    for (unsigned i=0;i<vector_const.size();i++) {
      dato[j][i] = vector_const[i];
    }
  }

//Data *nd = Data::DataFromFile("kk.cre");
//nd->redim(x2div*x1div);


  for(int i=0;i<x2div;i++) {
    int bit, byte;
    byte = 0;
    bit = sizeof(*fila)*8;
    fila[byte]=0;
    for(int j=0;j<x1div;j++) {
      dato[j][x2] = x2min + div2*i + div2*0.5; //X2
      dato[j][x1] = x1min + div1*j + div1*0.5; //X1
//nd->SetValueVar(x2div*i+j, 1, dato[j][x2]);    
//nd->SetValueVar(x2div*i+j, 0, dato[j][x1]);    
//nd->SetValueVar(x2div*i+j, 2, 0);    
    }
    c->Classify(dato, x1div);
    for(int j=0;j<x1div;j++) {
      unsigned char clase = (char)dato[j][2];
      bit-=v;
      fila[byte] = fila[byte] | (clase<<bit);
      if (bit==0) {byte++; fila[byte]=0; bit=sizeof(*fila)*8;}
    }
    fwrite(fila, AnchoEx, 1, bmp);
  }
//nd->SaveToFile("kk2.cre");
/*  for(int i=0;i<x2div;i++) {
    dato[x2] = x2min + div2*i + div2*0.5; //X2
    int bit, byte;
    bit = byte = 0;
    bit = sizeof(*fila)*8;
    fila[byte]=0;
    for(int j=0;j<x1div;j++) {
      dato[x1] = x1min + div1*j + div1*0.5; //X1
      c->Classify(&dato, 1);
      char clase = (char)dato[2];
      bit--;
      fila[byte] = fila[byte] | (clase<<bit);
      if (bit==0) {byte++; fila[byte]=0; bit=sizeof(*fila)*8;}
    }
    fwrite(fila, AnchoEx, 1, bmp);
  }*/

  fclose(bmp);

  for(int j=0;j<x1div;j++)
    delete []dato[j];
  delete []dato;
  delete []fila;
  return;
}
//---------------------------------------------------------------------------
void DoMapaDeAlturas2D(Function *f, int x1, int x2, double x1min, double x2min,
                                                double x1max, double x2max,
                                                int x1div, int x2div,
                                                vector<double> vector_const,
                                                string nomfichsal,
                                                double zmin, double zmax)
{
/*
int x1, x2;
 double x1min, x2min, x1max, x2max;
 int x1div, x2div;   */
 ;
 ;
 ;

  double div1 = (x1max-x1min)/x1div;
  double div2 = (x2max-x2min)/x2div;
  int AnchoEx = (x1div/4 + (x1div%4==0 ? 0 : 1)) * 4;
  unsigned char *fila = new unsigned char[AnchoEx];

  BITMAPFILEHEADER bfh;
  BITMAPINFOHEADER bih;
  RGBQUAD color;
  int ncolores = 256;
  int sizeofbfh = 14;  //sizeof(bfh)

  bfh.bfType      = 0x4d42;
  bfh.bfSize      = sizeof(bih) + sizeofbfh + AnchoEx*x2div + ncolores*4;
  bfh.bfReserved1 = 0;
  bfh.bfReserved2 = 0;
  bfh.bfOffBits   = sizeof(bih) + sizeofbfh + ncolores*4;

  bih.biSize          = sizeof(bih);
  bih.biWidth         = x1div;
  bih.biHeight        = x2div;
  bih.biPlanes        = 1;
  bih.biBitCount      = 8; 
  bih.biCompression   = BI_RGB;
  bih.biSizeImage     = bfh.bfSize - bfh.bfOffBits;
  bih.biXPelsPerMeter = 100;
  bih.biYPelsPerMeter = 100;
  bih.biClrUsed       = 0;
  bih.biClrImportant  = 0;

  FILE *bmp=fopen(nomfichsal.c_str(), "wb");
  if (!bmp) return;

//  fwrite(&bfh, sizeof(bfh), 1, bmp);
  fwrite(&bfh.bfType, 2, 1, bmp);
  fwrite(&bfh.bfSize, 4, 1, bmp);
  fwrite(&bfh.bfReserved1, 2, 1, bmp);
  fwrite(&bfh.bfReserved1, 2, 1, bmp);
  fwrite(&bfh.bfOffBits, 4, 1, bmp);
  fwrite(&bih, sizeof(bih), 1, bmp);
  for (int i=0;i<ncolores;i++) {
    color.rgbBlue     = i;
    color.rgbGreen    = i;
    color.rgbRed      = i;
    color.rgbReserved = 0x00;
    fwrite(&color, sizeof(color), 1, bmp);
  }


  double **dato = new double*[x1div];
  double *res = new double[x1div];
  for(int j=0;j<x1div;j++) {
    dato[j] = new double[vector_const.size() + 1];
    for (unsigned i=0;i<vector_const.size();i++) {
      dato[j][i] = vector_const[i];
    }
  }

  if (zmin>=zmax) {
    zmin = f->min();
    zmax = f->max();
  }
  
  double a, b;
  a = 256.0/(zmax-zmin);
  b = -a*zmin;
  /*for(int i=0;i<x2div;i++) {
    dato[x2] = x2min + div2*i + div2*0.5; //X2
    for(int j=0;j<x1div;j++) {
      dato[x1] = x1min + div1*j + div1*0.5; //X1
      double val = a*f->f(dato) + b;
      if (val>=256.0) val = 255.5;
      else if (val<0.0) val = 0.5;
      fila[j] = (unsigned char)val;
    }
    fwrite(fila, AnchoEx, 1, bmp);
  }*/
  for(int i=0;i<x2div;i++) {
    for(int j=0;j<x1div;j++) {
      dato[j][x2] = x2min + div2*i + div2*0.5; //X2
      dato[j][x1] = x1min + div1*j + div1*0.5; //X1
    }
    f->f(dato, x1div, res);
    for(int j=0;j<x1div;j++) {
      double val = a*res[j] + b;
      if (val>=256.0) val = 255.5;
      else if (val<0.0) val = 0.5;
      fila[j] = (unsigned char)val;
    }
    fwrite(fila, AnchoEx, 1, bmp);
  }

  fclose(bmp);

  for(int j=0;j<x1div;j++)
    delete []dato[j];
  delete []dato;
  delete []fila;
  delete []res;
  return;
}
//---------------------------------------------------------------------------

