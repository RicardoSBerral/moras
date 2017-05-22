//---------------------------------------------------------------------------
// 
// 
//---------------------------------------------------------------------------
#ifndef ClassifierH
#define ClassifierH
//---------------------------------------------------------------------------
#include "data20.h"
#include <typeinfo>
//---------------------------------------------------------------------------
class FuncionDeProgreso;
class Evaluation;
class Matriz;
//---------------------------------------------------------------------------
typedef enum _clfFormato {clfCompleto, cflRapido} clfFormato;
//---------------------------------------------------------------------------
//--------------------------------------------------------  Classifier  -----
//---------------------------------------------------------------------------
class Classifier 
{
  private:
    static std::string LeeCabecera(std::istream &in, int &version);
    void GuardaCabecera(std::ostream &salida, int version);

  protected:
    Data *data;
    double estimError;

  public://Chapuza temporal, no sE cOmo hacer acceso trasversal
    std::vector<int> DatosUtilizados;
    std::vector<double> APrioriClas;
    std::vector<int> rc, pc;

  protected:
    //Datos del conjunto de datos con el que se construye el clasificador
    int NumVar;
    int NumVarOrd;
    int NumVarFuz;
    int NumVarNom;
    std::string nombre();

  public:
    Classifier();
    virtual ~Classifier();

    //Funciones virtuales
    virtual void Init(Data * data);

    //Una de las dos siguientes funciones debe ser implementada en las clases
    //que hereden de Classifier
    virtual int Classify(int ElementIndex);
    virtual std::vector<double> Distribution(int ElementIndex);

    virtual std::vector<double> UnnormalizedDistribution(int ElementIndex);
    virtual double Average(int ElementIndex) {return Classify(ElementIndex);}
    virtual std::vector<double> MultipleAverage(int ElementIndex)
    {
      throw std::runtime_error("Multi-objetivo no permitido");
    }

    virtual Evaluation* Evaluate(Data *data, std::string EvaluationClassName="", 
                                                           std::string params="");

    virtual int Classify(std::vector<double> distrb);

    virtual void Classify(double **dat, int NumeroDatos);
    double ClassificationCertainty(int ElementIndex, int &Class);

    virtual void Build(Data *data, FuncionDeProgreso *fp=0) = 0;

    virtual double Error(Data *d);
    virtual double Error(int first, int last);
    virtual double EstimError() {return estimError;}

    //
    Data* GetData(){return data;}
    virtual void SetData(Data *data);

    void ResetOriginalTrainingData();
    int OriginalTrainingDataCount(){return DatosUtilizados.size();}
    int OriginalTrainingDataPos(int Indice){
      return DatosUtilizados[Indice];
    }
    bool UsedInOriginalTrainingData(int value);
    virtual std::string Info(int value=0);

    static char nf_algo[260];
    static bool Leyendo;
    static FILE *f;
    virtual void ReadOrder(Data *data);
    virtual void WriteOrder(Data *data);

    //Funciones para guardar y leer de fichero
    void Guardar(std::string salida, int version=0);
    void Leer(std::string in);
    int LeerConCabecera(std::istream &in);
    virtual void Guardar(std::ostream &salida, int version=0);
    virtual void Leer(std::istream &in, int version=0);
    static Classifier *LeerClasificadorGeneral(std::istream &in);

    //Mas  
    static int WhichClass(double* distrb, int NumClass, int *winners=0, 
                                             std::vector<double> *APrioriClas=0);
    static int WhichClass(std::vector<double> distrb, int *winners=0, 
                                             std::vector<double> *APrioriClas=0);
    static double Margin(Instance& dis, int RealClass=-1, bool Normalized=true);
    static double Margin(std::vector<double> dis, int RealClass=-1, bool Normalized=true);
    static double Margin(double* distrb, int NumClass, int RealClass=-1, bool Normalized=true);

    //Funciones de medida
    //Estadistico kappa (ver margineantu97pruning)
    static double KappaStatistic(Classifier *c1, Classifier *c2, int ini=-1, int fin=-1);
    static double CommonErrors(int *c1, int *c2, int *trueClass, int size);
    static double KappaStatistic(int *c1, int *c2, int m, int L);
    friend Data* ReduceData(Data *data, Classifier *clsf, double Umbral);
    friend Data* ExpandData(Data *data, Classifier *clsf);
    friend Data* ExpandData2(Data *data, Classifier *clsf);
};
//---------------------------------------------------------------------------
//--------------------------------------------------  FilterClassifier  -----
//---------------------------------------------------------------------------
class FilterClassifier : public Classifier 
{

  protected:
    Classifier *c;

  public:
    FilterClassifier(Classifier *c);
    virtual ~FilterClassifier();

    //Funciones virtuales
    virtual void Init(Data * data);
    virtual int Classify(int ElementIndex);
    virtual int Classify(std::vector<double> distrb);
    virtual std::vector<double> Distribution(int ElementIndex);
    virtual void Build(Data *data, FuncionDeProgreso *fp=0) = 0;
    virtual double Error(int first, int last);
    virtual double EstimError() {return estimError;}
    virtual void SetData(Data *data);
    virtual std::string Info(int value=0);
    virtual void ReadOrder(Data *data);
    virtual void WriteOrder(Data *data);
    virtual void Guardar(std::ostream &salida);
    virtual void Leer(std::istream &in, bool LeerCabecera=true);
  
};
//---------------------------------------------------------------------------
//------------------------------------------------------  MIClassifier  -----
//---------------------------------------------------------------------------
class MIClassifier : public Classifier
{

  public:
    MIClassifier() : Classifier() { }
    virtual ~MIClassifier() { }

    //Funciones virtuales
    virtual void Init(MINomData * data);

    //Una de las dos siguientes funciones debe ser implementada en las clases
    //que hereden de Classifier
    virtual int Classify(int ElementIndex);
    virtual std::vector<double> Distribution(int ElementIndex);

    virtual void Build(MINomData *data, FuncionDeProgreso *fp=0) = 0;

    virtual double Error(Data *d);
    virtual double Error(int first, int last);

};
//---------------------------------------------------------------------------
//---------------------------------------------------------  ClsfUtils  -----
//---------------------------------------------------------------------------
int **MatClasif(std::vector<Classifier*> clsfs, Data *data);
double QStatistic(Classifier *c1, Classifier *c2, Data *data);
double QStatistic(int *c1, int *c2, int *clase_ok, int nDatos);
double QStatistic(int *er1, int *er2, int nClasf);
//---------------------------------------------------------------------------
//------------------------------------------------------  ClsfRegister  -----
//---------------------------------------------------------------------------
class ClsfReg;
class ClsfRegister
{
  private:
    std::vector<ClsfReg*> ClsfRegs;
    static ClsfRegister *_Register;
  public:
    Classifier *CreateClassifier(std::string clsfname);
    ClsfReg *GetRegByType(const std::type_info &type);
    unsigned size(){ return ClsfRegs.size();}
    ClsfReg *GetReg(int num) { return ClsfRegs[num];}
    bool Registrate(ClsfReg *value) {
      ClsfRegs.push_back(value);
      return true;
    }
    static ClsfRegister *Register();
};
//---------------------------------------------------------------------------
//-----------------------------------------------------------  ClsfReg  -----
//---------------------------------------------------------------------------
class ClsfReg
{
  public:
    ClsfReg(std::string _name, std::string _description, 
                                                       bool autoregister=true) {
      name = _name;
      description = _description;
      proto = 0;
      if (autoregister) ClsfRegister::Register()->Registrate(this);
    }
    virtual ~ClsfReg(){;}

  public:
    virtual Classifier *CreateClassifier()=0;
    bool ClassifierTypeIs(const std::type_info &type) {
      return typeid(*GetProto())==type;
    }
    virtual std::string ClassifierName(){ return name; }
    virtual std::string ClassifierDescription(){ return description; }

  protected:
    std::string name;
    std::string description;
    Classifier *proto;
  protected:
    Classifier *GetProto() {
      if(!proto) proto=CreateClassifier();
      return proto;
    }
};
#endif
