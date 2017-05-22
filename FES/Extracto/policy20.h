#ifndef _POLICY_H_
#define _POLICY_H_


#include "variable20.h"
#include <vector>
#include <string>
#define _THEN 0


// Rule types
#define _CONDITIONAL 0
#define _DECLARATIVE 1

// Reasoning Methods
#define _FLC 0
#define _FES 1

// Defuz methods
#define _CENTROIDS 0
#define _COMPOS_MOM 1

//1class CConsole;
class Hedge
{
  public:
    Hedge(int c, const char* s)
      : Name(s),Code(c) {}
    std::string Name;
//2_    char Name[256];
    int Code;
};

class HedgedSet : public FuzzySet
{
public:
  HedgedSet(FuzzySet* fs, int* h, int n);
  virtual ~HedgedSet(){delete [] hedges;}
  int* hedges;
  int num;
  void Generate();
  FuzzySet* FS;
};


//namespace cagonto {
class Statement
{
public:
  Statement()
    {Prop = NULL; Sets = NULL;sib = parent = child = NULL;connector = -1;}
  Property* Prop;
  HedgedSet* Sets;
  void Evaluate(int);
  void AddSet(HedgedSet *h){h->SetNext(Sets); Sets = h;}
  Statement *sib, *parent, *child;
  int connector;
  double value,con_val;
};
//}

class Model;

class Rule
{
public:
  Rule(): Name("dummy") {root = NULL;Type = _CONDITIONAL;Consequent=NULL;Next = NULL;}
  ~Rule(){};
  void AddStatement(/*cagonto::*/class Statement* s);
  int Parse(Model*, const char*);
  void SetText(const char* buf){Text = buf;}
  void Read(FILE* is);
  void Write(FILE* os);
  Rule*  GetNext(){return Next;}
  void SetNext(Rule* r){Next = r;}
  void SetName(const char* s){Name = s;}
  std::string GetName(){return Name;}
  std::string GetText(){return Text;}
//2_  char* GetName(){return Name;}
//2_  char* GetText(){return Text;}
  void Evaluate(double temp[VECMAX], int infer);

protected:
  FuzzySet* FindTerm(const char*, Property*);
  /*cagonto::*/class Statement* root;
  HedgedSet *Consequent;
  int Type;
  std::string Text,Name;
//2_  char Text[256],Name[256];
  Rule* Next;
};

class CSSView;
class CCartParmsDlg;
class DataSet
{

  public:
  void GenerateLED();
  int ExtractSet(int label);
  void LabelData(int nlabels);

    bool GetPonder(){ return ponder;}
    bool GetScramble(){ return scramble;}
    void SetPonder(bool p){ ponder = p;}
    void SetScramble(bool s){ scramble = s;}
//    void AddFuzzVar(Property* tempp);
    void SetFuzzVars(Property* prop){FuzzVars = prop;}
    void InitFuzzVars();
    void InitializeNomTerms(int n);
    VarType VarTypeOf(const std::string cs, int* loc);
//2_    VarType VarTypeOf(const char* cs, int* loc);
    VarType GetTypeDepVar(){return DepVarType;}
    void SetDepVarType(VarType vt){DepVarType = vt;}
    void InitialiseValues(double** values);
    DataSet(){label =""; NumNoms =0;VarNames = NULL; VarValue =NULL; NumVars =0;NumData =0;}
    DataSet(int nvars, int ndata);
    DataSet(const char* tag,  char** names, double** values,
      int nvars, int ndata);
    DataSet(char** names, int nvars, int ndata=0);
    ~DataSet();
  //  void SetValueVar(int i, int j,double value){ VarValue[j][i] =value;}
  //  void SetNameVar(int i, char* name){ delete[] VarNames[i]; VarNames[i] = new char[strlen(name)+1]; strcpy(VarNames[i],name);}
    std::string GetLabel(){return label;}
//2_    char* GetLabel(){return label;}
    int GetNumVars(){return NumVars;}
    void SetNumVars(int i){NumVars =i;}
    int GetNumData(){return NumData;}
    void SetNumData(int i){NumData =i;}
    int GetNumNoms(){return NumNoms;}
    void SetNumVarsOrdSplit(int n){ NumVarsOrdSplit =n;}
    int GetNumVarsOrdSplit(){return NumVarsOrdSplit;}
    void SetNumNoms(int n){ NumNoms =n;}
    int GetNumOrds(){return NumOrds;}
    void SetNumOrds(int n){ NumOrds =n;}
    int GetNumFuzz(){return NumFuzz;}
    void SetNumFuzz(int n){ NumFuzz =n;}
    std::string GetNomTerm(int i,int j){/*3return NomTerms[i].GetAt(j);*/
    return NomTerms[i].at(j);}
//2_    char* GetNomTerm(int i,int j){return NomTerms[i].GetAt(j);}
    int GetSizeNomTerm(int i){/*3return NomTerms[i].GetSize();*/
     return NomTerms[i].size();}
    std::string GetNameVar(int i){ return VarNames[i];}
//2_    char* GetNameVar(int i){ return VarNames[i];}
    void SetNameVar(int i,std::string cs){ VarNames[i] = cs;}
    Property* GetFuzzVars(){ return FuzzVars;}
    std::string GetDepVar(){ return DepVar;}
//2_    char* GetDepVar(){ return DepVar;}
    void  SetDepVar(std::string cs){VarNames[NumVars-1] = DepVar = cs;}
    double GetValueVar(int i, int j){return VarValue[j][i];}
    void Exchange(int col1,int col2);
    void AddTerm(int i, std::string newterm);

  protected:
    bool ponder;
    bool scramble;
    std::string label;
    std::string DepVar;
    VarType DepVarType;
    int NumVars;
    int NumNoms;
    int NumOrds;
    int NumFuzz;
    int NumVarsOrdSplit;
    int NumData;
    Property*  FuzzVars;
    std::string* VarNames;
    double** VarValue;
    std::vector<std::string> *NomTerms;// Translation table for nominals
    void Defaults(int ndata);
public:
    int NPrune;
    int ncv;
    int nmin;
    int ntest;
    int NSel;
  // Parameters Multiple splits
    bool multsplits;
    double beta;
    double epsilon;


  // Parameters PLM
    bool plm;
    bool deladd;
    double MinFToAdd;
    double MaxFToDelete;
    double tol;
  // Flag for fuzzification of ordinal splits
    bool fuzzify;

  // Parameters for optimisation
    int degree;
  // Flag for multiple trees (only for multivalued nominal dep var)
    bool multree;
};

class Policy
{
public:
  Policy(const char* name, const char* con,
    int defuz = _CENTROIDS, int infer = _FLC);
  ~Policy();
   friend class Model;
   std::string GetName(){return Name;}
   std::string GetConsequent(){return Consequent;}
   std::string GetInferString()
    {if(Infer == _FLC) return std::string("FLC"); return std::string("FES");}
  std::string GetDefuzString()
  {if(Infer == _CENTROIDS) return std::string("Centroids");
      return std::string("Comp. Mom.");}
/*2_   char* GetName(){return Name;}
   char* GetConsequent(){return Consequent;}
   char* GetInferString()
    {if(Infer == _FLC) return "FLC"; return "FES";}
  char* GetDefuzString()
  {if(Infer == _CENTROIDS) return "Centroids";
      return "Comp. Mom.";}*/
  int GetInfer(){return Infer;}
    int GetDefuz(){return DeFuz;}
  Policy* GetNext(){return Next;}
  void AddRule(Rule* r){r->SetNext(Rules); Rules = r;}
  Rule* GetRules(){return Rules;}
    void DeleteRule(Rule* r);
  void Read(FILE* is);
  void Write(FILE* os);
//  double Evaluate(Property*);
  double DeFuzzify(int defuz = -1);
  void EvaluateRule(Rule* r,int infer = -1);
  FuzzySet Work;
protected:
  std::string Name;
//2_  char Name[256];
  Rule* Rules;
  Policy* Next;
  std::string Consequent;
//2_  char Consequent[256];
  double value;
  int DeFuz;
  int Infer;
};

class Tree;
class Data;
class ParmsCART;
class DTPolicy;

void BuildTreeThread(DTPolicy* dtp);
class DTPolicy : public Policy
{
  public:
    double OptimiseParams();

    DTPolicy(/*1CConsole* pt*/void (*_console)(char *));
    DTPolicy(/*1CConsole* pt,*/void (*_console)(char *), DataSet* ds, int defuz = _CENTROIDS, int infer = _FLC);
    ~DTPolicy();
    void SetNmin(int i){Nmin = i;}
    Tree* GetTree(){return tree;}
    Data* GetData(){ return data;}
    int GetNmin(){return Nmin;}
    void BuildTree();
    friend class Model;

  protected:
    std::string* VarNames;
//2_    char* VarNames[256];
    Data* data;
    Tree* tree;
    int Nmin;
//1    CConsole* console;
void (*console)(char *);
};


class Model
{
  public:
  void OptimiseParams();
  Property* FindVariable(const char* buf);

    Model(const char *name);
    ~Model();
    void AddPolicy(const char* name, const char* con, int defuz, int infer);
    Property* AddProperty(const char* ,int type,double domain[2], int dec=0);
    void Read(FILE* is);
    void Write(FILE* os);
    Property* GetProps(){return Props;}
    Policy* GetPolicies(){return Policies;}
    int FindHedge(const char*);
    void GenerateCART(DataSet* ds, /*1CConsole* pt,*/
      int defuz = _CENTROIDS, int infer = _FLC);

  protected:

    Policy* Policies;
    std::string Name;
//2_    char* Name[256];
    int NumHedges;
    Property* Props;

};




#endif
