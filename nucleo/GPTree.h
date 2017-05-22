//---------------------------------------------------------------------------
// ClassExplorer Pro generated header file
// Created by  on 04/12/00, 18:48:28
//---------------------------------------------------------------------------
#ifndef GPTreeH
#define GPTreeH
#include "Classifier.h"

//---------------------------------------------------------------------------
class Tree;
//---------------------------------------------------------------------------
class GPTree : public Classifier
{
protected:
  Tree* tree;
  int CountGroup[2];
  int NMin;
  int NumIter;
  double GroupsError();
public:
  GPTree( void (*salida)(char *)=0, int _NMin=1);
  virtual ~GPTree();
  virtual void Init(Data * data);
  virtual double Error(int first, int last);
  virtual int Classify(int ElementIndex);
  virtual void Build(Data *data, FuncionDeProgreso *fp=0);
  void CreateGroups();
  void SetGroup(int i);
  //
  int GetNumIter(){return NumIter;}
  int GetGroupCount(int i){return CountGroup[i];}
  void SetGroupCount(int i, int Value){CountGroup[i] = Value;}
  Tree *GetTree();
  static bool ComoArticulo;
  //Funciones para guardar y leer de fichero
  virtual void Guardar(std::ostream &salida, int version=0);
  virtual void Leer(std::istream &in, int version=0);
};
//---------------------------------------------------------------------------
class DecisionTree : public Classifier
{
  protected:
    Tree* tree;
    Node *lastClassifyingNode;
    bool PruneTree;
    bool SE_0;
    bool MultiSplits;
    int NMin;

  public:
    DecisionTree( void (*salida)(char *)=0, int _NMin=1, bool PruneTree=true, 
               bool MultiSplits=false, bool SE_0=false);
    virtual ~DecisionTree();

    virtual void Init(Data * data);
    virtual double Error(int first, int last);
    virtual int Classify(int ElementIndex);
    virtual double Average(int ElementIndex);
    virtual std::vector<double> MultipleAverage(int ElementIndex);
    virtual std::vector<double> UnnormalizedDistribution(int ElementIndex);
    virtual std::string Info(int value=0);

    //
    Tree *GetTree();
    Node *GetLastClassifyingNode(){return lastClassifyingNode;}

    //Funciones para guardar y leer de fichero
    virtual void Guardar(std::ostream &salida, int version=0);
    virtual void Leer(std::istream &in, int version=0);

    int GetK();

    static int K;
};
//---------------------------------------------------------------------------
class CARTTree : public DecisionTree
{
  protected:
 // double GroupsError();

  public:
    CARTTree( void (*salida)(char *)=0, int _NMin=1, bool PruneTree=true, 
               bool MultiSplits=false, bool SE_0=false);
    virtual ~CARTTree();

    virtual void Build(Data *data, FuncionDeProgreso *fp=0);

};
//---------------------------------------------------------------------------
class RandomForestTree : public DecisionTree
{
  protected:
    int NumberOfAttributes;
    std::vector<int> *atts;

  public:
    RandomForestTree(void (*salida)(char *)=0, int _NMin=1, 
         int NumberOfAttributes=2, bool MultiSplits=false, int AttsToRandomize = -1);
    virtual ~RandomForestTree();

    virtual void Build(Data *data, FuncionDeProgreso *fp=0);

};
//---------------------------------------------------------------------------
class DecisionStump : public DecisionTree
{
  public:
    DecisionStump(void (*s)(char *)=0, int _NMin=1, bool MultiSplits=false);
    virtual ~DecisionStump();

    virtual void Build(Data *data, FuncionDeProgreso *fp=0);
};
//---------------------------------------------------------------------------
class CARTTreeReg : public ClsfReg
{
  private:
    static CARTTreeReg *reg;
    static CARTTreeReg *autoreg(){return new CARTTreeReg();}
    CARTTreeReg() : ClsfReg("CARTTree", "CARTTree", true){}

  public:
    Classifier *CreateClassifier(){ return new CARTTree();}

    static CARTTreeReg *Reg(){return reg;}
};
//---------------------------------------------------------------------------
class RandomForestTreeReg : public ClsfReg
{
  private:
    static RandomForestTreeReg *reg;
    static RandomForestTreeReg *autoreg(){return new RandomForestTreeReg();}
    RandomForestTreeReg() : ClsfReg("RandomForestTree", "RandomForestTree", true){}

  public:
    Classifier *CreateClassifier(){ return new RandomForestTree();}

    static RandomForestTreeReg *Reg(){return reg;}
};
//---------------------------------------------------------------------------
class DecisionStumpReg : public ClsfReg
{
  private:
    static DecisionStumpReg *reg;
    static DecisionStumpReg *autoreg(){return new DecisionStumpReg();}
    DecisionStumpReg() : ClsfReg("DecisionStump", "DecisionStump", true){}

  public:
    Classifier *CreateClassifier(){ return new DecisionStump();}

    static DecisionStumpReg *Reg(){return reg;}
};
//---------------------------------------------------------------------------
#endif
