//---------------------------------------------------------------------------

#ifndef FnsClasificacionH
#define FnsClasificacionH

#include "Lector.h"
#include "Language.h"
#include <map> 
#include <vector> 

//---------------------------------------------------------------------------
class Data;
class Classifier;
class CARTTree;
class Ensemble;
class Matriz;

//---------------------------------------------------------------------------
void  DarDeAltaFnsClasificacion();
//---------------------------------------------------------------------------
void DoCreaConjuntosTrainAleat(Data *data, int numdatos,
                            int numconjuntos, std::string salida, int opts=1);
/*void DoMapaDeClasif2D(Classifier *c, int x1, int x2,
                                                double x1min, double x2min,
                                                double x1max, double x2max,
                                                int x1div, int x2div,
                                                std::vector<double> std::vector_const,
                                                std::string nomfichsal);
void DoMapaDeAlturas2D(Classifier *c, int x1, int x2, double x1min, double x2min,
                                                double x1max, double x2max,
                                                int x1div, int x2div,
                                                std::vector<double> std::vector_const,
                                                std::string nomfichsal,
                                                double zmin, double zmax);*/
//---------------------------------------------------------------------------
class TipoData : public TipoDatos
{
  protected:
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoData*>(tipo);
    }

  public:
    virtual std::string nombre(){return "data";}

  public:
    static TipoData *TData(){return _TData;}
  private:
    static TipoData *_TData;
};
class ParametroData : public Parametro
{
  public:
    ParametroData() : Parametro() {
      PonPropiedades(false, 0, "Ejemplos", "");
    }

    virtual Tipo* DameTipo(){return TipoData::TData();}
};

//---------------------------------------------------------------------------
/*class TipoMultiData : public TipoData
{
  protected:
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoMultiData*>(tipo);
    }

  public:
    virtual std::string nombre(){return "multidata";}

  public:
    static TipoMultiData *TMultiData(){return _TMultiData;}
  private:
    static TipoMultiData *_TMultiData;
};
class ParametroMultiData : public Parametro
{
  public:
    ParametroMultiData() : ParametroData() {
      PonPropiedades(false, 0, "Multiple Instances", "");
    }

    virtual Tipo* DameTipo(){return TipoMultiData::TMultiData();}
};
*/
//---------------------------------------------------------------------------
class TipoClasificador : public TipoDatos
{
  protected:
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoClasificador*>(tipo);
    }

  public:
    virtual std::string nombre(){return "clasificador";}

  public:
    static TipoClasificador* TClasificador() {return _TClasificador;}
  private:
    static TipoClasificador* _TClasificador;
};
class ParametroClasificador : public Parametro
{
  public:
    ParametroClasificador() : Parametro() {
      PonPropiedades(false, 0, "Clasificador", "");
    }

    virtual Tipo* DameTipo(){return TipoClasificador::TClasificador();}
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class TipoEnsemble : public TipoClasificador
{
  protected:
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoEnsemble*>(tipo);
    }

  public:
    virtual Variable *NuevaVariable(std::string nombre, Mandato *inic=0);
    virtual std::string nombre(){return "conjunto_de_clasificadores";}

  public:
    static TipoEnsemble* TEnsemble() {return _TEnsemble;}
  private:
    static TipoEnsemble* _TEnsemble;
};
class VariableEnsemble : public VariableDatos, public IVariableDesglosable
{
    std::map<Classifier *, int> pos;
    std::vector<Variable*> clsfs;

  private:
    void init(Ensemble *);
    void clear();

  public:
    VariableEnsemble(std::string nombre, void* valor=0) :
                        VariableDatos(nombre, valor) {init((Ensemble*)valor);}
    VariableEnsemble(void* valor) : VariableDatos(valor) 
                                                     {init((Ensemble*)valor);}
    VariableEnsemble() : VariableDatos() {}
    ~VariableEnsemble();

  public:
    virtual Tipo* DameTipo(){return TipoEnsemble::TEnsemble();}
    virtual void Set(Mandato *m){
      SetDatos(m->ComoDatos());
      init((Ensemble*)FValor);
    }

    int NumElementos();
    Variable* Elemento(int i);

};
class ParametroEnsemble : public ParametroClasificador
{
  public:
    ParametroEnsemble() : ParametroClasificador() {
      PonPropiedades(false, 0, "Conjunto", "");
    }

    virtual Tipo* DameTipo(){return TipoEnsemble::TEnsemble();}
};
//---------------------------------------------------------------------------
class LoadClassifier : public IntefazFuncion, Funcion
{
   Classifier *Clasif;
  public:
    LoadClassifier(std::vector<Mandato *> mands) : Funcion(mands) { 
      Clasif = 0;
      CreateParams();
    }
    LoadClassifier() : Funcion() { Clasif = 0; CreateParams(); }

    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new LoadClassifier();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("LoadClassifier");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo();
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)Clasif;}
    virtual bool ComoBooleano(){return Clasif;}
    virtual void* ComoDatos(){return Clasif;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class LoadDataset : public IntefazFuncion, Funcion
{
   Data *nd;

  public:
    LoadDataset(std::vector<Mandato *> mands) : Funcion(mands) 
                                                   { nd = 0; CreateParams();}
    LoadDataset() : Funcion() { nd = 0; CreateParams();}
    virtual ~LoadDataset();

    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new LoadDataset();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("LoadDataset");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoData::TData();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)nd;}
    virtual bool ComoBooleano(){return nd;}
    virtual void* ComoDatos(){return nd;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class GenerateWrongLabels : public IntefazFuncion, Funcion
{
   Matriz *mat;

  public:
    GenerateWrongLabels(std::vector<Mandato *> mands) : Funcion(mands) 
                                                   { mat = 0; CreateParams();}
    GenerateWrongLabels() : Funcion() { mat = 0; CreateParams();}
    virtual ~GenerateWrongLabels();

    void CreateParams();
    virtual unsigned NumParams(){return 1;}
    virtual Funcion* CrearFuncion() {return new GenerateWrongLabels();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("GenerateWrongLabels");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoData::TData();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)mat;}
    virtual bool ComoBooleano(){return mat;}
    virtual void* ComoDatos(){return mat;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class Corrected : public IntefazFuncion, Funcion
{

  double corr;

  public:
    Corrected(std::vector<Mandato *> mands) : Funcion(mands) 
                                                   { CreateParams();}
    Corrected() : Funcion() { CreateParams();}
    virtual ~Corrected(){;}

    void CreateParams();
    virtual unsigned NumParams(){return 5;}
    virtual Funcion* CrearFuncion() {return new Corrected();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("Corrected");}

    virtual Mandato* Ejecutar();

    virtual bool EsNumero(){return true;}

    virtual std::string ComoCadena(){return ACad(corr);}
    virtual double ComoNumero(){return corr;}
    virtual int ComoEntero(){return (int)corr;}
    virtual bool ComoBooleano(){return (int)corr;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class LoadLabels : public IntefazFuncion, Funcion
{

  public:
    LoadLabels(std::vector<Mandato *> mands) : Funcion(mands) 
                                                   { CreateParams();}
    LoadLabels() : Funcion() { CreateParams();}
    virtual ~LoadLabels(){;}

    void CreateParams();
    virtual unsigned NumParams(){return 4;}
    virtual Funcion* CrearFuncion() {return new LoadLabels();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("LoadLabels");}

    virtual Mandato* Ejecutar();

    virtual std::string ComoCadena(){ return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class FilterData : public IntefazFuncion, Funcion
{
   Data *nd;

  public:
    FilterData(std::vector<Mandato *> mands) : Funcion(mands) 
                                                   { nd = 0; CreateParams();}
    FilterData() : Funcion() { nd = 0; CreateParams();}
    virtual ~FilterData();

    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new FilterData();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("FilterData");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoData::TData();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)nd;}
    virtual bool ComoBooleano(){return nd;}
    virtual void* ComoDatos(){return nd;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class SetSplitCriterium : public IntefazFuncion, Funcion
{
   Data *nd;

  public:
    SetSplitCriterium(std::vector<Mandato *> mands) : Funcion(mands) 
                                                   { nd = 0; CreateParams();}
    SetSplitCriterium() : Funcion() { nd = 0; CreateParams();}
    virtual ~SetSplitCriterium() { }

    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new SetSplitCriterium();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("SetSplitCriterium");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoData::TData();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)nd;}
    virtual bool ComoBooleano(){return nd;}
    virtual void* ComoDatos(){return nd;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class Classify : public IntefazFuncion, Funcion
{
   Matriz *mat;

  public:
    Classify(std::vector<Mandato *> mands) : Funcion(mands) { 
      mat=0;
      CreateParams();
    }
    Classify() : Funcion() { mat=0; CreateParams();}
    ~Classify();

    void CreateParams();
    virtual unsigned NumParams(){return 5;}
    virtual Funcion* CrearFuncion() {return new Classify();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("Classify");}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return mat;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class Error : public IntefazFuncion, Funcion
{
   double err;

  public:
    Error(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    Error() : Funcion() { CreateParams(); }
    ~Error() {  }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new Error();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("Error");}

    virtual Mandato* Ejecutar();

    virtual bool EsNumero(){return true;}

    virtual std::string ComoCadena(){return ACad(err);}
    virtual double ComoNumero(){return err;}
    virtual int ComoEntero(){return (int)err;}
    virtual bool ComoBooleano(){return (int)err;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class ResubstituteNodeStats : public IntefazFuncion, Funcion
{

  public:
    ResubstituteNodeStats(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    ResubstituteNodeStats() : Funcion() { CreateParams(); }
    ~ResubstituteNodeStats() {  }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 4;}
    virtual Funcion* CrearFuncion() {return new ResubstituteNodeStats();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("ResubstituteNodeStats");}

    virtual Mandato* Ejecutar();

    virtual std::string ComoCadena(){return "true";}
    virtual double ComoNumero(){return 1;}
    virtual int ComoEntero(){return 1;}
    virtual bool ComoBooleano(){return true;}
    virtual void* ComoDatos(){return NULL;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class EvaluateClassifier : public IntefazFuncion, Funcion
{
   std::string output;

  public:
    EvaluateClassifier(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    EvaluateClassifier() : Funcion() { CreateParams(); }
    ~EvaluateClassifier() {  }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 4;}
    virtual Funcion* CrearFuncion() {return new EvaluateClassifier();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("EvaluateClassifier");}

    virtual Mandato* Ejecutar();

    virtual bool EsCadena(){return true;}

    virtual std::string ComoCadena(){return output;}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class Mrse : public IntefazFuncion, Funcion
{
   double tot;

  public:
    Mrse(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    Mrse() : Funcion() { CreateParams(); }
    ~Mrse() {  }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new Mrse();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("Mrse");}

    virtual Mandato* Ejecutar();

    virtual bool EsNumero(){return true;}

    virtual std::string ComoCadena(){return ACad(tot);}
    virtual double ComoNumero(){return tot;}
    virtual int ComoEntero(){return (int)tot;}
    virtual bool ComoBooleano(){return (int)tot;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class SequentialError : public IntefazFuncion, Funcion
{
  protected:
    Matriz *err;

  public:
    SequentialError(std::vector<Mandato *> mands) : Funcion(mands) {
      err = 0;
      CreateParams();
    }
    SequentialError() : Funcion() { err = 0; CreateParams();}
    ~SequentialError();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new SequentialError();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("SequentialError");}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)err;}
    virtual bool ComoBooleano(){return err;}
    virtual void* ComoDatos(){return err;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class ClassificationMatrix : public IntefazFuncion, Funcion
{
  protected:
    Matriz *mat;

  public:
    ClassificationMatrix(std::vector<Mandato *> mands) : Funcion(mands) {
      mat = 0;
      CreateParams();
    }
    ClassificationMatrix() : Funcion() { mat = 0; CreateParams();}
    ~ClassificationMatrix();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 4;}
    virtual Funcion* CrearFuncion() {return new ClassificationMatrix();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("ClassificationMatrix");}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)mat;}
    virtual bool ComoBooleano(){return mat;}
    virtual void* ComoDatos(){return mat;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class ConfusionMatrix : public IntefazFuncion, Funcion
{
   Matriz *mat;

  public:
    ConfusionMatrix(std::vector<Mandato *> mands) : Funcion(mands) { 
      mat=0;
      CreateParams();
    }
    ConfusionMatrix() : Funcion() { mat=0; CreateParams();}
    ~ConfusionMatrix();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new ConfusionMatrix();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("ConfusionMatrix");}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return mat;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class Margin : public IntefazFuncion, Funcion
{
   Matriz *mat;

  public:
    Margin(std::vector<Mandato *> mands) : Funcion(mands) {
      mat=0;
      CreateParams();
    }
    Margin() : Funcion() { mat=0;CreateParams(); }
    ~Margin();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new Margin();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("Margin");}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return mat;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class MarginOfInstances : public IntefazFuncion, Funcion
{
   Matriz *mat;

  public:
    MarginOfInstances(std::vector<Mandato *> mands) : Funcion(mands) {
      mat=0;
      CreateParams();
    }
    MarginOfInstances() : Funcion() { mat=0;CreateParams(); }
    ~MarginOfInstances();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new MarginOfInstances();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("MarginOfInstances");}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return mat;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class Certainty : public IntefazFuncion, Funcion
{
   Matriz *mat;

  public:
    Certainty(std::vector<Mandato *> mands) : Funcion(mands) {
      mat=0;
      CreateParams();
    }
    Certainty() : Funcion() { mat=0;CreateParams(); }
    ~Certainty();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new Certainty();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("Certainty");}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return mat;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class MarginSum : public IntefazFuncion, Funcion
{
   Matriz *mat;

  public:
    MarginSum(std::vector<Mandato *> mands) : Funcion(mands)
                                { mat=0;CreateParams(); }
    MarginSum() : Funcion() { mat=0;CreateParams(); }
    ~MarginSum();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new MarginSum();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("MarginSum");}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return mat;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class EstimErrorVal : public IntefazFuncion, Funcion
{
   Matriz *mat;

  public:
    EstimErrorVal(std::vector<Mandato *> mands) : Funcion(mands)
                                { mat=0;CreateParams(); }
    EstimErrorVal() : Funcion() { mat=0;CreateParams(); }
    ~EstimErrorVal();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new EstimErrorVal();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("EstimErrorVal");}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return mat;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class EstimErrorTest : public IntefazFuncion, Funcion
{
   Matriz *mat;

  public:
    EstimErrorTest(std::vector<Mandato *> mands) : Funcion(mands)
                                { mat=0;CreateParams(); }
    EstimErrorTest() : Funcion() { mat=0;CreateParams(); }
    ~EstimErrorTest();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 4;}
    virtual Funcion* CrearFuncion() {return new EstimErrorTest();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("EstimErrorTest");}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return mat;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class DiversityMeasures : public IntefazFuncion, Funcion
{
  protected:
    Matriz *div;

  public:
    DiversityMeasures(std::vector<Mandato *> mands) : Funcion(mands) {
      div = 0;
      CreateParams();
    }
    DiversityMeasures() : Funcion() { div = 0; CreateParams();}
    ~DiversityMeasures();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new DiversityMeasures();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("DiversityMeasures");}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)div;}
    virtual bool ComoBooleano(){return div;}
    virtual void* ComoDatos(){return div;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
class EnsembleMeasures : public IntefazFuncion, Funcion
{
  protected:
    Matriz *mat;

  public:
    EnsembleMeasures(std::vector<Mandato *> mands) : Funcion(mands) {
      mat = 0;
      CreateParams();
    }
    EnsembleMeasures() : Funcion() { mat = 0; CreateParams();}
    ~EnsembleMeasures();

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new EnsembleMeasures();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("EnsembleMeasures");}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)mat;}
    virtual bool ComoBooleano(){return mat;}
    virtual void* ComoDatos(){return mat;}

    static void DarDeAlta();

};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class SetPropertyValue : public IntefazFuncion, Funcion
{
  public:
    SetPropertyValue(std::vector<Mandato *> mands) : Funcion(mands) { 
      CreateParams(); }
    SetPropertyValue() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new SetPropertyValue();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("SetPropertyValue");}

    virtual Mandato* Ejecutar();

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class SetNTrain : public IntefazFuncion, Funcion
{
  public:
    SetNTrain(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    SetNTrain() : Funcion() { CreateParams(); }

    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new SetNTrain();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("SetNTrain");}

    virtual Mandato* Ejecutar();

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}
  
    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class GeneratePartitions : public IntefazFuncion, Funcion
{
  public:
    GeneratePartitions(std::vector<Mandato *> mands) : Funcion(mands) { 
      CreateParams(); }
    GeneratePartitions() : Funcion() { CreateParams(); }

    void CreateParams();
    virtual unsigned NumParams(){return 4;}
    virtual Funcion* CrearFuncion() {return new GeneratePartitions();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("GeneratePartitions");}

    virtual Mandato* Ejecutar();

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class SaveDataset : public IntefazFuncion, Funcion
{
  public:
    SaveDataset(std::vector<Mandato *> mands) : Funcion(mands) { 
      CreateParams(); }
    SaveDataset() : Funcion() { CreateParams(); }

    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new SaveDataset();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("SaveDataset");}

    virtual Mandato* Ejecutar();

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class SaveClassifier : public IntefazFuncion, Funcion
{
  public:
    SaveClassifier(std::vector<Mandato *> mands) : Funcion(mands) { 
      CreateParams(); }
    SaveClassifier() : Funcion() { CreateParams(); }

    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new SaveClassifier();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("SaveClassifier");}

    virtual Mandato* Ejecutar();

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class ClassifierInfo : public IntefazFuncion, Funcion
{
  std::string kk;
  public:
    ClassifierInfo(std::vector<Mandato *> mands) : Funcion(mands) { 
      CreateParams(); }
    ClassifierInfo() : Funcion() { CreateParams(); }

    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new ClassifierInfo();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("ClassifierInfo");}

    virtual Mandato* Ejecutar();

    virtual bool EsCadena() {return true;}

    virtual std::string ComoCadena() {return kk;}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class BuildIGPEnsemble : public IntefazFuncion, Funcion
{
   Classifier *Clasif;
  public:
    BuildIGPEnsemble(std::vector<Mandato *> mands) : Funcion(mands) {
      Clasif = 0;
      CreateParams();
    }
    BuildIGPEnsemble() : Funcion() { Clasif = 0; CreateParams();}
    virtual ~BuildIGPEnsemble();

    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new BuildIGPEnsemble();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("BuildIGPEnsemble");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoEnsemble::TEnsemble();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)Clasif;}
    virtual bool ComoBooleano(){return Clasif;}
    virtual void* ComoDatos(){return Clasif;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class BuildBaggingCART : public IntefazFuncion, Funcion
{
   Classifier *Clasif;
  public:
    BuildBaggingCART(std::vector<Mandato *> mands) : Funcion(mands) {
      Clasif = 0;
      CreateParams();
    }
    BuildBaggingCART() : Funcion() { Clasif = 0; CreateParams();}
//    virtual ~BuildBaggingCART() { if (Clasif) delete Clasif; }

    void CreateParams();
    virtual unsigned NumParams(){return 8;}
    virtual Funcion* CrearFuncion() {return new BuildBaggingCART();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("BuildBaggingCART");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoEnsemble::TEnsemble();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)Clasif;}
    virtual bool ComoBooleano(){return Clasif;}
    virtual void* ComoDatos(){return Clasif;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class BuildBaggingC45 : public IntefazFuncion, Funcion
{
   Classifier *Clasif;
  public:
    BuildBaggingC45(std::vector<Mandato *> mands) : Funcion(mands) {
      Clasif = 0;
      CreateParams();
    }
    BuildBaggingC45() : Funcion() { Clasif = 0; CreateParams();}
    virtual ~BuildBaggingC45();

    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new BuildBaggingC45();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("BuildBaggingC45");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoEnsemble::TEnsemble();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)Clasif;}
    virtual bool ComoBooleano(){return Clasif;}
    virtual void* ComoDatos(){return Clasif;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class BuildBaggingNNet : public IntefazFuncion, Funcion
{
  protected:
    Classifier *Clasif;

  public:
    BuildBaggingNNet(std::vector<Mandato *> mands) : Funcion(mands) {
      Clasif = 0;
      CreateParams();
    }
    BuildBaggingNNet() : Funcion() { Clasif = 0; CreateParams();}
    virtual ~BuildBaggingNNet();

    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new BuildBaggingNNet();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("BuildBaggingNNet");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoEnsemble::TEnsemble();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)Clasif;}
    virtual bool ComoBooleano(){return Clasif;}
    virtual void* ComoDatos(){return Clasif;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class BuildRandomForest : public IntefazFuncion, Funcion
{
   Classifier *Clasif;
  public:
    BuildRandomForest(std::vector<Mandato *> mands) : Funcion(mands) {
      Clasif = 0;
      CreateParams();
    }
    BuildRandomForest() : Funcion() { Clasif = 0; CreateParams();}
//    virtual ~BuildRandomForest() { if (Clasif) delete Clasif; }

    void CreateParams();
    virtual unsigned NumParams(){return 9;}
    virtual Funcion* CrearFuncion() {return new BuildRandomForest();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("BuildRandomForest");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoEnsemble::TEnsemble();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)Clasif;}
    virtual bool ComoBooleano(){return Clasif;}
    virtual void* ComoDatos(){return Clasif;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
/*class BuildBaggingComMod : public IntefazFuncion, Funcion
{
   Classifier *Clasif;
  public:
    BuildBaggingComMod(std::vector<Mandato *> mands) : Funcion(mands) { Clasif = 0; }
    BuildBaggingComMod() : Funcion() { Clasif = 0; }

  virtual unsigned NumParams(){return 2;}
  virtual Funcion* CrearFuncion() {return new BuildBaggingComMod();}
  virtual std::string NombreFuncion() {return StringRepository::GetString("BuildBaggingComMod");}

    virtual Mandato* Ejecutar();
    virtual bool EsDatos(){return true;}
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (int)Clasif;}
    virtual bool ComoBooleano(){return Clasif;}
    virtual void* ComoDatos(){return Clasif;}

    static void DarDeAlta();
};*/
//---------------------------------------------------------------------------
class BuildBoosting : public IntefazFuncion, Funcion
{
   Classifier *Clasif;
  public:
    BuildBoosting(std::vector<Mandato *> mands) : Funcion(mands) { Clasif = 0; CreateParams();}
    BuildBoosting() : Funcion() { Clasif = 0; CreateParams();}
  //  virtual ~BuildBoosting() { if (Clasif) delete Clasif; }

    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new BuildBoosting();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("BuildBoosting");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoEnsemble::TEnsemble();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)Clasif;}
    virtual bool ComoBooleano(){return Clasif;}
    virtual void* ComoDatos(){return Clasif;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class BuildBoostingC45 : public IntefazFuncion, Funcion
{
  protected:
   Classifier *Clasif;

  public:
    BuildBoostingC45(std::vector<Mandato *> mands):Funcion(mands) {Clasif = 0; CreateParams();}
    BuildBoostingC45() : Funcion() { Clasif = 0; CreateParams();}
    //virtual ~BuildBoostingC45() { if (Clasif) delete Clasif; }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new BuildBoostingC45();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("BuildBoostingC45");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoEnsemble::TEnsemble();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)Clasif;}
    virtual bool ComoBooleano(){return Clasif;}
    virtual void* ComoDatos(){return Clasif;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class BuildClassSwitching : public IntefazFuncion, Funcion
{
  protected:
   Classifier *Clasif;

  public:
    BuildClassSwitching(std::vector<Mandato *> mands):Funcion(mands) {Clasif = 0; CreateParams();}
    BuildClassSwitching() : Funcion() { Clasif = 0; CreateParams();}

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 5;}
    virtual Funcion* CrearFuncion() {return new BuildClassSwitching();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("BuildClassSwitching");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoEnsemble::TEnsemble();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)Clasif;}
    virtual bool ComoBooleano(){return Clasif;}
    virtual void* ComoDatos(){return Clasif;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
/*class BuildBoostingWithBaggingInfo : public IntefazFuncion, Funcion
{
   Classifier *Clasif;
  public:
    BuildBoostingWithBaggingInfo(std::vector<Mandato *> mands) : Funcion(mands) { Clasif = 0; CreateParams();}
    BuildBoostingWithBaggingInfo() : Funcion() { Clasif = 0; CreateParams();}

    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new BuildBoostingWithBaggingInfo();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("BuildBoostingWithBaggingInfo");}

    virtual Mandato* Ejecutar();
    virtual bool EsDatos(){return true;}
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (int)Clasif;}
    virtual bool ComoBooleano(){return Clasif;}
    virtual void* ComoDatos(){return Clasif;}

    static void DarDeAlta();
};                                                                    */
//---------------------------------------------------------------------------
class BuildCART : public IntefazFuncion, Funcion
{
   CARTTree *Clasif;
  public:
    BuildCART(std::vector<Mandato *> mands) : Funcion(mands) { Clasif = 0; CreateParams();}
    BuildCART() : Funcion() { Clasif = 0; CreateParams();}

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 4;}
    virtual Funcion* CrearFuncion() {return new BuildCART();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("BuildCART");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoClasificador::TClasificador();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)Clasif;}
    virtual bool ComoBooleano(){return Clasif;}
    virtual void* ComoDatos(){return Clasif;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class BuildC45 : public IntefazFuncion, Funcion
{
   Classifier *c45tree;
  public:
    BuildC45(std::vector<Mandato *> mands) : Funcion(mands) { c45tree = 0; CreateParams();}
    BuildC45() : Funcion() { c45tree = 0; CreateParams();}

    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new BuildC45();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("BuildC45");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoClasificador::TClasificador();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)c45tree;}
    virtual bool ComoBooleano(){return c45tree;}
    virtual void* ComoDatos(){return c45tree;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class BuildNNet : public IntefazFuncion, Funcion
{
   Classifier *nnet;
  public:
    BuildNNet(std::vector<Mandato *> mands) : Funcion(mands) { nnet = 0; CreateParams();}
    BuildNNet() : Funcion() { nnet = 0; CreateParams();}

    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new BuildNNet();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("BuildNNet");}

    virtual Mandato* Ejecutar();

    virtual Tipo* DameTipo(){return TipoClasificador::TClasificador();}
    virtual bool EsDatos(){return true;}

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)nnet;}
    virtual bool ComoBooleano(){return nnet;}
    virtual void* ComoDatos(){return nnet;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
/*class EnsembleReport : public IntefazFuncion, Funcion
{
  public:
    EnsembleReport(std::vector<Mandato *> mands) : Funcion(mands) { }
    EnsembleReport() : Funcion() {  }

  virtual unsigned NumParams(){return 3;}
  virtual Funcion* CrearFuncion() {return new EnsembleReport();}
  virtual std::string NombreFuncion() {return StringRepository::GetString("EnsembleReport");}

    virtual Mandato* Ejecutar();
    virtual bool EsCadena(){return true;}
    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
};*/
//---------------------------------------------------------------------------
class ClassificationMap2D : public IntefazFuncion, Funcion
{
  public:
    ClassificationMap2D(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    ClassificationMap2D() : Funcion() { CreateParams(); }

    void CreateParams();
    virtual unsigned NumParams(){return 10;}
    virtual Funcion* CrearFuncion() {return new ClassificationMap2D();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("ClassificationMap2D");}

    virtual Mandato* Ejecutar();

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class TestAll : public IntefazFuncion, Funcion
{
  void *algo;

  public:
    TestAll(std::vector<Mandato *> mands) : Funcion(mands) { algo = 0; CreateParams(); }
    TestAll() : Funcion() { algo = 0; CreateParams();}

    void CreateParams();
    virtual unsigned NumParams(){return 1;}
    virtual Funcion* CrearFuncion() {return new TestAll();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("TestAll");}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return (double)((long int)algo);}
    virtual int ComoEntero(){return (long int)algo;}
    virtual bool ComoBooleano(){return algo;}
    virtual void* ComoDatos(){return algo;}

    virtual bool EsDatos(){return true;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
class Clear : public IntefazFuncion, Funcion
{
  protected:
    Mandato *mand;

  public:
    Clear(std::vector<Mandato *> mands) : Funcion(mands) {CreateParams();  }
    Clear() : Funcion() {CreateParams();  }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 1;}
    virtual Funcion* CrearFuncion() {return new Clear();}
    virtual std::string NombreFuncion() {return StringRepository::GetString("Clear");}

    virtual Mandato* Ejecutar();

    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------

/**
 * Computes the G matrix of an ensemble a 
 * sum all the values of it.
 */
class GValue : public IntefazFuncion, Funcion
{
  double gValue;

  public:
    GValue(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    GValue() : Funcion() { CreateParams(); }
    ~GValue() {  }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new GValue();}
    virtual std::string NombreFuncion() {return "GValue";}

    virtual Mandato* Ejecutar();

    virtual bool EsNumero(){return true;}

    virtual std::string ComoCadena(){return ACad(gValue);}
    virtual double ComoNumero(){return gValue;}
    virtual int ComoEntero(){return (int)gValue;}
    virtual bool ComoBooleano(){return (int)gValue;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();

};

#endif
 
