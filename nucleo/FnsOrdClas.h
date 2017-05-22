//---------------------------------------------------------------------------

#ifndef FnsOrdClasH
#define FnsOrdClasH
#include "Lector.h"

//---------------------------------------------------------------------------
class Matriz;
class Ensemble;
class Data;
class NomData;
//---------------------------------------------------------------------------
void DarDeAltaFnsOrdClas();
//---------------------------------------------------------------------------
double Distancia(int *P1, int *P2, int nDatos);
int DoPseudoBoosting(Ensemble *ens, int ini, int fin, bool Resampling=false);
void DoDiversityOrdering(Ensemble *ens, NomData *data);
int solve_SDP(double **H, double **D, int dim, int k);
int *approx_SDP(double **H, double **X, int dim, int k, double *value);
void DoOrdenaPorContribucion(Ensemble *ens, int ini=0, int fin=-1);
// *** DHL
int DoOrdenaPorReduccionError(Ensemble *ens, int ini=0, int fin=-1,
                                           int start_at=0, bool inverse=false);
int DoOrdenaPorRERegularizadoPorDiversidad(Ensemble *ens, int ini=0, int fin=-1, 
                                           int start_at=0, bool inverse=false);

void DoOrdenaPorReduccionDeDistanciasCV(Ensemble *ens, int dist=10,
                                            int ini=0, int fin=-1, int nCV=10);
void DoOrdenaPorReduccionDeDistancias(Ensemble *ens, int dist=10,
                                                        int ini=0, int fin=-1);

void DoOrdenaPorReduccionErrorPonderado(Ensemble *ens, int ini=0, int fin=-1);
void DoOrdenaPorReduccionErrorDoble(Ensemble *ens, int ini=0, int fin=-1);
void DoOrdenaSimplex(Ensemble *ens, int ini=0, int fin=-1);
void DoOrdenaSimplex2(Ensemble *ens, int ini=0, int fin=-1);
void DoOrdenaEnPruebas(Ensemble *ens, int ini=0, int fin=-1, Mandato*m=0);
void DoOrdenaCV(Ensemble *ens, int ini=0, int fin=-1, int ncv=-1);
void DoMinimoIterativo(Ensemble *ens, int ini=0, int fin=-1);
double **CreaMatrizDistancias(Ensemble *ens, Data *data, bool ExcluirTrain);
Matriz* CreaMatrizCmn(Ensemble *ens, int ini=0, int fin=-1, 
                                                       bool ExcluirTrain=true);
void DoOrdenaPorError(Ensemble *ens, int ini, int fin, bool ExcluirTrain=true);
int TestDistancias(Ensemble *ens, Data *data, bool SetNegativeToZero=false, bool UseOOB=false);
//---------------------------------------------------------------------------
class OrdenOriginal : public IntefazFuncion, Funcion
{
  public:
    OrdenOriginal(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    OrdenOriginal() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new OrdenOriginal();}
    virtual std::string NombreFuncion() {return "OrdenOriginal";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class ClasificadoresNecesarios : public IntefazFuncion, Funcion
{
  std::vector<int> posiciones;
  std::vector<double> porcensref;
  std::vector<double> estimerror;
  void DoClasificadoresNecesarios(Ensemble *ens, int ini=0, int fin=-1);
  void init() {
    porcensref.push_back(0.00);
    porcensref.push_back(0.01);
    porcensref.push_back(0.02);
    porcensref.push_back(0.03);
    posiciones.push_back(0);
    CreateParams();
  }

  public:
    ClasificadoresNecesarios(std::vector<Mandato *> mands) : Funcion(mands) 
       {init();}
    ClasificadoresNecesarios() : Funcion() {init();}

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new ClasificadoresNecesarios();}
    virtual std::string NombreFuncion() {return "ClasificadoresNecesarios";}

    virtual bool EsEntero(){return true;}
    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena();
    virtual double ComoNumero(){return (double)posiciones[0];}
    virtual int ComoEntero(){return posiciones[0];}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class OrdenaPorContribucion : public IntefazFuncion, Funcion
{
  public:
    OrdenaPorContribucion(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    OrdenaPorContribucion() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new OrdenaPorContribucion();}
    virtual std::string NombreFuncion() {return "OrdenaPorContribucion";}

    virtual Mandato* Ejecutar();

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class OrdenaEnPruebas : public IntefazFuncion, Funcion
{
  protected:
    Matriz *res;
    int iminmin, iminmax;

  public:
    OrdenaEnPruebas(std::vector<Mandato *> params) : Funcion(params) { CreateParams(); res=0;}
    OrdenaEnPruebas() : Funcion() {CreateParams();res=0;  }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new OrdenaEnPruebas();}
    virtual std::string NombreFuncion() {return "OrdenaEnPruebas";}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual Mandato* Ejecutar();
    virtual bool EsEntero(){return true;}
    virtual std::string ComoCadena();// {return "";}
    virtual double ComoNumero(){return iminmin;}
    virtual int ComoEntero(){return iminmin;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return res;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class OrdenaPorReduccionError : public IntefazFuncion, Funcion
{
  int num;

  public:
    OrdenaPorReduccionError(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    OrdenaPorReduccionError() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new OrdenaPorReduccionError();}
    virtual std::string NombreFuncion() {return "OrdenaPorReduccionError";}

    virtual Mandato* Ejecutar();
    virtual bool EsEntero(){return true;}

    virtual std::string ComoCadena() {return Mandato::ACad(num);}
    virtual double ComoNumero(){return num;}
    virtual int ComoEntero(){return num;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class OrdenaPorReduccionDeDistancias : public IntefazFuncion, Funcion
{
  public:
    OrdenaPorReduccionDeDistancias(std::vector<Mandato *> mands) : Funcion(mands) {  CreateParams();}
    OrdenaPorReduccionDeDistancias() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new OrdenaPorReduccionDeDistancias();}
    virtual std::string NombreFuncion() {return "OrdenaPorReduccionDeDistancias";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class OrdenaPorAngulos : public IntefazFuncion, Funcion
{
  protected:
    int num;

  public:
    OrdenaPorAngulos(std::vector<Mandato *> mands) : Funcion(mands) 
                     {  CreateParams();}
    OrdenaPorAngulos() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new OrdenaPorAngulos();}
    virtual std::string NombreFuncion() {return "OrdenaPorAngulos";}

    virtual Mandato* Ejecutar();
    virtual bool EsEntero(){return true;}

    virtual std::string ComoCadena() {return Mandato::ACad(num);}
    virtual double ComoNumero(){return num;}
    virtual int ComoEntero(){return num;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class OrdenaPorEPIC : public IntefazFuncion, Funcion
{
  protected:
    int num;

  public:
    OrdenaPorEPIC(std::vector<Mandato *> mands) : Funcion(mands) 
                     {  CreateParams();}
    OrdenaPorEPIC() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new OrdenaPorEPIC();}
    virtual std::string NombreFuncion() {return "OrdenaPorEPIC";}

    virtual Mandato* Ejecutar();
    virtual bool EsEntero(){return true;}

    virtual std::string ComoCadena() {return Mandato::ACad(num);}
    virtual double ComoNumero(){return num;}
    virtual int ComoEntero(){return num;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class OrdenaPorUWA : public IntefazFuncion, Funcion
{
  protected:
    int num;

  public:
    OrdenaPorUWA(std::vector<Mandato *> mands) : Funcion(mands) 
                     {  CreateParams();}
    OrdenaPorUWA() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new OrdenaPorUWA();}
    virtual std::string NombreFuncion() {return "OrdenaPorUWA";}

    virtual Mandato* Ejecutar();
    virtual bool EsEntero(){return true;}

    virtual std::string ComoCadena() {return Mandato::ACad(num);}
    virtual double ComoNumero(){return num;}
    virtual int ComoEntero(){return num;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class InvertirOrden : public IntefazFuncion, Funcion
{
  public:
    InvertirOrden(std::vector<Mandato *> mands) : Funcion(mands) 
                     {  CreateParams();}
    InvertirOrden() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 1;}
    virtual Funcion* CrearFuncion() {return new InvertirOrden();}
    virtual std::string NombreFuncion() {return "InvertirOrden";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class OrdenaPorReduccionErrorPonderado : public IntefazFuncion, Funcion
{
  public:
    OrdenaPorReduccionErrorPonderado(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    OrdenaPorReduccionErrorPonderado() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new OrdenaPorReduccionErrorPonderado();}
    virtual std::string NombreFuncion() {return "OrdenaPorReduccionErrorPonderado";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class PseudoBoosting : public IntefazFuncion, Funcion
{
  protected:
    int res;
  public:
    PseudoBoosting(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    PseudoBoosting() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new PseudoBoosting();}
    virtual std::string NombreFuncion() {return "PseudoBoosting";}

    virtual bool EsEntero(){return true;}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return ACad(res);}
    virtual double ComoNumero(){return (double)res;}
    virtual int ComoEntero(){return res;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class OrdenaPorReduccionErrorDoble : public IntefazFuncion, Funcion
{
  public:
    OrdenaPorReduccionErrorDoble(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    OrdenaPorReduccionErrorDoble() : Funcion() { CreateParams(); }


  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new OrdenaPorReduccionErrorDoble();}
    virtual std::string NombreFuncion() {return "OrdenaPorReduccionErrorDoble";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class SeleccionGenetica : public IntefazFuncion, Funcion
{
  protected:
    int res;

  public:
    SeleccionGenetica(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    SeleccionGenetica() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 6;}
    virtual Funcion* CrearFuncion() {return new SeleccionGenetica();}
    virtual std::string NombreFuncion() {return "SeleccionGenetica";}

    virtual bool EsEntero(){return true;}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return ACad(res);}
    virtual double ComoNumero(){return (double)res;}
    virtual int ComoEntero(){return res;}
    virtual bool ComoBooleano(){return res;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class AnalisisDeOrden : public IntefazFuncion, Funcion
{
  public:
    AnalisisDeOrden(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    AnalisisDeOrden() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 5;}
    virtual Funcion* CrearFuncion() {return new AnalisisDeOrden();}
    virtual std::string NombreFuncion() {return "AnalisisDeOrden";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class SeleccionCombinatoria : public IntefazFuncion, Funcion
{
  protected:
    Matriz *res;

  public:
    SeleccionCombinatoria(std::vector<Mandato *> mands) : Funcion(mands) {
      res = 0;
      CreateParams();
    }
    SeleccionCombinatoria() : Funcion() { res = 0; CreateParams();}
    ~SeleccionCombinatoria();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new SeleccionCombinatoria();}
    virtual std::string NombreFuncion() {return "SeleccionCombinatoria";}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)res;}
    virtual bool ComoBooleano(){return res;}
    virtual void* ComoDatos(){return res;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class Mejores : public IntefazFuncion, Funcion
{
  protected:
    Matriz *res;

  public:
    Mejores(std::vector<Mandato *> mands) : Funcion(mands) {
      res = 0;
      CreateParams();
    }
    Mejores() : Funcion() { res = 0; CreateParams();}
    ~Mejores();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new Mejores();}
    virtual std::string NombreFuncion() {return "Mejores";}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)res;}
    virtual bool ComoBooleano(){return res;}
    virtual void* ComoDatos(){return res;}

    static void DarDeAlta();

};
//----------------------------------------------------------------------
#ifdef COMP_SDP
class SeleccionSDP: public IntefazFuncion, Funcion
{
  protected:
    int res;

  public:
    SeleccionSDP(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    SeleccionSDP() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 6;}
    virtual Funcion* CrearFuncion() {return new SeleccionSDP();}
    virtual std::string NombreFuncion() {return "SeleccionSDP";}

    virtual bool EsEntero(){return true;}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return ACad(res);}
    virtual double ComoNumero(){return (double)res;}
    virtual int ComoEntero(){return res;}
    virtual bool ComoBooleano(){return res;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
#endif
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class OrdenaPorGGreedy : public IntefazFuncion, Funcion
{
  public:
    OrdenaPorGGreedy(std::vector<Mandato *> mands) : Funcion(mands)
                     {  CreateParams();}
    OrdenaPorGGreedy() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new OrdenaPorGGreedy();}
    virtual std::string NombreFuncion() {return "OrdenaPorGGreedy";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
double computeRow(double** G, int size, int* previousSelected, int previousSize, int trial);
};

//---------------------------------------------------------------------------

class OrdenaPorRERegularizadoPorDiversidad : public IntefazFuncion, Funcion
{
  public:
    OrdenaPorRERegularizadoPorDiversidad(std::vector<Mandato *> mands) : Funcion(mands)
                     {  CreateParams();}
    OrdenaPorRERegularizadoPorDiversidad() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new OrdenaPorRERegularizadoPorDiversidad();}
    virtual std::string NombreFuncion() {return "OrdenaPorRERegularizadoPorDiversidad";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
    double computeRow(double** G, int size, int* previousSelected, int previousSize, int trial);
};


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class FnsOrdClasUtils {
  public:
    // Prints an array
    void static print(int *data, int n, std::ostream &os);
    void static print(double **data, int n, std::ostream &os);
    void static swap(int &a, int &b);
    void static swap(int *array, int i, int j);
};


#endif
 
