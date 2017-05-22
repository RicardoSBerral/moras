//---------------------------------------------------------------------------
// C45Tree: wrapper for the adapted version of c4.5  
//---------------------------------------------------------------------------
#ifndef C45TreeH
#define C45TreeH
//---------------------------------------------------------------------------
#include "Classifier.h"
#include "data20.h"
#include "types.i"
#include "extern.i"
#include <vector>
#include <map>
//---------------------------------------------------------------------------
class Data;
//---------------------------------------------------------------------------
class C45Tree : public Classifier
{
protected:
  c45::Tree c45tree;
  c45::Tree c45tree_unpruned;
  bool UsePruned;
  bool GetUsePruned();
  c45::Tree GetTree();
  static Data *lastdataset;
public:
  C45Tree(bool UsePruned=true); 
  virtual ~C45Tree();

  virtual void Init(Data * data);
  virtual int Classify(int ElementIndex);
  virtual std::vector<double> UnnormalizedDistribution(int ElementIndex);
  virtual void Build(Data *data, FuncionDeProgreso *fp=0);
  virtual void SetData(Data *data);
  virtual std::string Info(int value);

  static void SetGlobalOpt(int opt, char *arg);
  //
  //Funciones para guardar y leer de fichero
  virtual void Guardar(std::ostream &salida, int version=0);
  virtual void Leer(std::istream &in, int version=0);

  static int UsePrunedTrees;
};
//---------------------------------------------------------------------------
class C45Data
{
  protected:
    short MaxAtt;               /* max att number */
    short MaxClass;             /* max class number */
    short MaxDiscrVal;          /* max discrete values for any att */
    c45::ItemNo MaxItem;        /* max data item number */
    c45::Description *Item;     /* data items */
    c45::DiscrValue *MaxAttVal; /* number of values for each att */
    char *SpecialStatus;        /* special att treatment */
    c45::String *ClassName;     /* class names */
    c45::String *AttName;       /* att names */
    c45::String **AttValName;   /* att value names */
    c45::Description *TrainData;
    std::map<int, int> isin;

    Data *data;                 /* Datos de origen */

  public:
    C45Data(Data *data);
    ~C45Data();

    void Init(Data *data);
    void GetData();
    void SetData();
    void SetTrainData();
    c45::Description GetDescription(int ElementIndex);
};
class C45DataSets : public DataEvents
{
  protected:
    std::map<Data *, C45Data*> datos;

    static C45DataSets *def;
    virtual void OnDelete(Data *);
    
  public:
    C45DataSets(){}
    virtual ~C45DataSets(){}

    C45Data *GetC45Data(Data *data);
    void SetTrainData(Data *data);
    c45::Description GetDescription(Data *data, int ElementIndex);

    static C45DataSets *GetDefault();
};
C45DataSets *GetC45DataSets();
//---------------------------------------------------------------------------

#endif
