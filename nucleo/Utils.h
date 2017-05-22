//---------------------------------------------------------------------------

#ifndef UtilsH
#define UtilsH

#include <string>
#include <vector>
//---------------------------------------------------------------------------
//--------------------------------------------  Definiciones parciales  -----
//---------------------------------------------------------------------------
class Ensemble;
class Data;
//---------------------------------------------------------------------------
//----------------------------------------------  Constantes generales  -----
//---------------------------------------------------------------------------
// Contantes generales
#ifndef WIN32
const char SEP = '/';
#else
const char SEP = '\\';
#endif

#ifndef WIN32
#define __int64 long int
#endif


//---------------------------------------------------------------------------
//-----------------------------------------------  Funciones generales  -----
//---------------------------------------------------------------------------
std::string GenerateNewFileName(std::string prev="", std::string ext="", 
                                                         std::string dir=".");
std::string GeneratePIDFileName(std::string prev="", std::string ext="", 
                                                         std::string dir=".");
std::string GenerateProcNameFileName(std::string prev="", std::string ext="", 
                                                         std::string dir=".");
void SetProcName(std::string prn);
std::string GetProcName();

int RandomInteger(int n);
double RandomDouble();
int RandomIntProportional(double *pop, int npop, double sumPop=1.0);
int RandomIntProportional(std::vector<double> pop, double sumPop=1.0);
bool Flip(double prob=0.5);
template<typename T> void RandomPermutation(T *elems, int size, int npermuts);

//---------------------------------------------------------------------------
//-----------------------------------------------------  Combinaciones  -----
//---------------------------------------------------------------------------
class Combinaciones
{
  protected:
    int nelementos;
    std::vector<int> indices;
    __int64 icomb;
    __int64 maxcomb;

  public:
    Combinaciones(int nelementos, int tomados_de);

  public:
    __int64 GetNComb();
    bool HayMasCombinaciones();
    std::vector<int> SiguienteCombinacion();
};

//---------------------------------------------------------------------------
//-----------------------------------------------------  FuncionObjeto  -----
//---------------------------------------------------------------------------
class FuncionObjeto
{
    public:
        FuncionObjeto(){;}
        virtual ~FuncionObjeto(){;}
    
    public:
        virtual bool operator()()=0;
};

//---------------------------------------------------------------------------
//-------------------------------------------------  FuncionDeProgreso  -----
//---------------------------------------------------------------------------
class FuncionDeProgreso
{
    protected:
        std::string texto;
        int numeroDeLlamadas;
        FuncionDeProgreso *f;

        bool next(int incr){ return f ? (*f)(incr) : true;}

    public:
        FuncionDeProgreso(FuncionDeProgreso *f=0){this->f=f; texto="call";}
        FuncionDeProgreso(std::string texto, FuncionDeProgreso *f=0) { 
            this->f=f;
            this->texto = texto;
        }
        virtual ~FuncionDeProgreso(){;}
    
    public:
        virtual void Start(int numeroDeLlamadas);
        virtual void End();
        virtual bool operator()(int incr=1)=0;
};

//---------------------------------------------------------------------------
//---------------------------------------------------  ProgresoConsola  -----
//---------------------------------------------------------------------------
class ProgresoConsola : public FuncionDeProgreso
{
    protected:
        double llamadasPorTic;
        int longDeBarra;
        int iLlamadas;
        int iTics;
        int t0;

    public:
        ProgresoConsola(FuncionDeProgreso *f=0, int longDeBarra=50) : 
                                                         FuncionDeProgreso(f) {
            this->longDeBarra = longDeBarra;
        }
        ProgresoConsola(std::string texto, FuncionDeProgreso *f=0, 
                            int longDeBarra=50) : FuncionDeProgreso(texto, f) {
            this->longDeBarra = longDeBarra;
        }
        virtual ~ProgresoConsola(){;}

    public:
        virtual void Start(int numeroDeLlamadas);
        virtual void End();
        virtual bool operator()(int incr=1);
};

//---------------------------------------------------------------------------
//--------------------------------------------------  ProgresoEnsemble  -----
//---------------------------------------------------------------------------
class ProgresoEnsemble : public FuncionDeProgreso
{
    protected:
        Ensemble *ens;
        Data *data;
        double **votos;
        int nClasf;
        std::string sal;

    public:
        ProgresoEnsemble(Ensemble *ens, Data *data, FuncionDeProgreso *f=0) 
                                                      : FuncionDeProgreso(f) {
            this->ens = ens;
            this->data = data;
            this->sal = "";
        }
        virtual ~ProgresoEnsemble(){;}

    public:
        virtual void Start(int numeroDeLlamadas);
        virtual void End();
        virtual bool operator()(int incr=1);
};
//---------------------------------------------------------------------------
#endif

