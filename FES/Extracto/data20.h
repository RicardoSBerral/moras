#ifndef __DATA_H
#define __DATA_H


#define MIN_BELONG_LOW 1e-16
#define MIN_BELONG_HIGH 0.5

#include "variable20.h"
#include "Matriz.h"
#include <vector>
#include <iterator>
#include <string>
#include<bitset>
#include<map>

class Node;
class Data;
class NomData;

class DataEvents
{
  friend class Data;

  public:
    virtual ~DataEvents() {}

  protected:
    virtual void OnDelete(Data *data)=0;
};

class Instance
{
  friend class Data;
  friend class NomData;
  friend class MINomData;

  protected:
    int nvars;
    int Class;
    int IniPos;
    int Group;
    double CvTreeLabel;
    double Membership;
    double Weight;

  public:
    Instance(int nvars){this->nvars = nvars;}
    virtual ~Instance(){ }

  public:
    static Instance *NewInstanceByName(std::string name);

  public:
    virtual Instance *NewInstance(int nvars) = 0;

  public:
    virtual double operator[](int index) = 0;
    virtual void operator()(int index, double value) = 0;
    virtual Instance& operator = (const Instance & other) = 0;
    virtual operator std::vector<double>() = 0;
    virtual void DeleteColumn(int icol) = 0;

  public:
    virtual inline int GetNVars() {return nvars;}

    virtual inline int GetDependentVariable() {return GetClass();}
    virtual inline int GetClass() {return Class;}
    virtual inline void SetClass(int cls) {Class = cls;}

    inline int GetIniPos(){return IniPos;}
    inline int GetGroup() {return Group;}
    inline double GetCvTreeLabel() {return CvTreeLabel;}
    
    virtual double GetWorkingVar(int i) = 0; 
    virtual void SetWorkingVar(int i, double value) = 0;
    
    inline double GetMembership() {return Membership;}
    inline double GetWeight() {return Weight;}
};
class BasicInstance : public Instance
{
  protected:
    double *data;

  public:
    BasicInstance(int nvars);
    BasicInstance(double *vals, int nvars);
    virtual ~BasicInstance();

  public:
    virtual Instance *NewInstance(int nvars) {return new BasicInstance(nvars);}

  public:
    inline virtual double operator[](int index);
    inline virtual void operator()(int index, double value);
    inline virtual Instance& operator = (const Instance & other);
    virtual operator std::vector<double>();
    virtual void DeleteColumn(int icol);

    virtual inline double GetWorkingVar(int i) {return data[nvars+1+i];}
    virtual inline void SetWorkingVar(int i, double value) {data[nvars+1+i] = value;}
};
class UInt8Instance : public Instance  // Tipico para hacer con templates....
{
  protected:
    unsigned char *data;
    double wv[2];

  public:
    UInt8Instance(int nvars);
    UInt8Instance(double *vals, int nvars);
    virtual ~UInt8Instance();

  public:
    virtual Instance *NewInstance(int nvars) {return new UInt8Instance(nvars);}

  public:
    inline virtual double operator[](int index);
    inline virtual void operator()(int index, double value);
    inline virtual Instance& operator = (const Instance & other);
    virtual operator std::vector<double>();
    virtual void DeleteColumn(int icol);

    virtual inline double GetWorkingVar(int i) {return wv[i];}
    virtual inline void SetWorkingVar(int i, double value) {wv[i] = value;}
};
class ShortInstance : public Instance
{
  protected:
    short *data;
    double wv[2];

  public:
    ShortInstance(int nvars);
    ShortInstance(double *vals, int nvars);
    virtual ~ShortInstance();

  public:
    virtual Instance *NewInstance(int nvars) {return new ShortInstance(nvars);}

  public:
    inline virtual double operator[](int index);
    inline virtual void operator()(int index, double value);
    inline virtual Instance& operator = (const Instance & other);
    virtual operator std::vector<double>();
    virtual void DeleteColumn(int icol);

    virtual inline double GetWorkingVar(int i) {return wv[i];}
    virtual inline void SetWorkingVar(int i, double value) {wv[i] = value;}
};
class Int32Instance : public Instance
{
  protected:
    int *data;
    double wv[2];

  public:
    Int32Instance(int nvars);
    Int32Instance(double *vals, int nvars);
    virtual ~Int32Instance();

  public:
    virtual Instance *NewInstance(int nvars) {return new Int32Instance(nvars);}

  public:
    inline virtual double operator[](int index);
    inline virtual void operator()(int index, double value);
    inline virtual Instance& operator = (const Instance & other);
    virtual operator std::vector<double>();
    virtual void DeleteColumn(int icol);

    virtual inline double GetWorkingVar(int i) {return wv[i];}
    virtual inline void SetWorkingVar(int i, double value) {wv[i] = value;}
};
class FloatInstance : public Instance
{
  protected:
    float *data;
    double wv[2];

  public:
    FloatInstance(int nvars);
    FloatInstance(double *vals, int nvars);
    virtual ~FloatInstance();

  public:
    virtual Instance *NewInstance(int nvars) {return new FloatInstance(nvars);}

  public:
    inline virtual double operator[](int index);
    inline virtual void operator()(int index, double value);
    inline virtual Instance& operator = (const Instance & other);
    virtual operator std::vector<double>();
    virtual void DeleteColumn(int icol);

    virtual inline double GetWorkingVar(int i) {return wv[i];}
    virtual inline void SetWorkingVar(int i, double value) {wv[i] = value;}
};
class ComboInstance : public Instance
{
  protected:
    unsigned char *data_description;
    float *data_location_bb;
    double wv[2];

  public:
    ComboInstance(int nvars);
    ComboInstance(double *vals, int nvars);
    virtual ~ComboInstance();

  public:
    virtual Instance *NewInstance(int nvars) {return new ComboInstance(nvars);}

  public:
    inline virtual double operator[](int index);
    inline virtual void operator()(int index, double value);
    inline virtual Instance& operator = (const Instance & other);
    virtual operator std::vector<double>();
    virtual void DeleteColumn(int icol);

    virtual inline double GetWorkingVar(int i) {return wv[i];}
    virtual inline void SetWorkingVar(int i, double value) {wv[i] = value;}
};

class MultilabelInstance : public BasicInstance
{
  protected:
    std::bitset<1024> labels;
    static int num_multilabels;
    static int label_index;

  public:
    MultilabelInstance(int nvars);
    MultilabelInstance(double *vals, int nvars);

  public:
    virtual Instance *NewInstance(int nvars) {return new MultilabelInstance(nvars);}

  public:
    inline virtual double operator[](int index);
    inline virtual void operator()(int index, double value);
    inline virtual Instance& operator = (const Instance & other);
    virtual operator std::vector<double>();
    virtual inline int GetClass() {return label_index < 0 ? Class : GetMultilabel(label_index);}

  public:
    static int GetNumMultilabels() { return num_multilabels; }
    static void SetNumMultilabels(int valor) { num_multilabels = valor; }
    static void SetSplitlabelIndex(int index) {label_index = index;}
    static int GetSplitlabelIndex() { return label_index;}

  public:
    bool GetMultilabel(int j) { return labels[j];}
    void SetMultilabel(int j, bool valor) {labels[j] = valor;}
};

class MultiobjectiveInstance : public BasicInstance
{
  protected:
    std::vector<double> objectives;
    static int num_multiobjectives;
    static int objective_index;

  public:
    MultiobjectiveInstance(int nvars);
    MultiobjectiveInstance(double *vals, int nvars);

  public:
    virtual Instance *NewInstance(int nvars) {return new MultiobjectiveInstance(nvars);}

  public:
    inline virtual double operator[](int index);
    inline virtual void operator()(int index, double value);
    inline virtual Instance& operator = (const Instance & other);
    virtual operator std::vector<double>();
    virtual inline int GetClass() {return objective_index < 0 ? Class : GetMultiobjective(objective_index);}

  public:
    static int GetNumMultiobjectives() { return num_multiobjectives; }
    static void SetNumMultiobjectives(int valor) { num_multiobjectives = valor; }
    static void SetSplitobjectiveIndex(int index) {objective_index = index;}
    static int GetSplitobjectiveIndex() { return objective_index;}

  public:
    double GetMultiobjective(int j) { return objectives[j];}
    void SetMultiobjective(int j, double valor) {objectives.insert(objectives.begin() + j, valor);}
};

class Instances
{
  std::vector<Instance*> instances;

  public:
    Instances(){};

  inline Instance& operator[](int index) {
    Instance* ins = instances[index];
    return *ins;
  }

  void push_back(Instance* ins) {
    instances.push_back(ins);
  }

  void pop_back() {
    delete instances.back();
    instances.pop_back();
  }

  void clear() {
    if (instances.size()==0) return;
    for(unsigned int i = 0; i < instances.size(); i++) {
      delete instances[i];
    }
    instances.clear();
  }

  std::vector<Instance*>::iterator begin() {
    return instances.begin();
  }
};
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
class SplitCriterium
{
  protected:
    Data *data;

  public:
    virtual ~SplitCriterium() {}

  public:
    void SetData(Data *data) {this->data = data;}

  public:
    virtual double GetLastSplitImpurity()=0;
    virtual double Initialise(int first, int last, int curSplit, double gamma = 100, int iX= -1, bool FuzzFlag = false)=0;
    virtual double MovePoint(int curSplit, double gamma, int iX= -1,int first=0,int last =0)=0;
    virtual void FixLeaf(Node* n, bool add=false)=0;
};

class GiniCriterium : public SplitCriterium
{
  protected:
    double R[2];
    double nL, nR;
    double *PopR, *PopL;

  public:
    GiniCriterium(Data *data) {
      SetData(data);
      PopR = PopL = 0;
    }
    virtual ~GiniCriterium() {
      if (PopL) delete []PopL;
      if (PopR) delete []PopR;
    }


  public:
    static double Impurity(double *Pop, int NumClass);

  public:
    virtual double GetLastSplitImpurity() {return R[0]+R[1];}
    virtual double Initialise(int first, int last, int curSplit, double gamma = 100, int iX= -1, bool FuzzFlag = false);
    virtual double MovePoint(int curSplit, double gamma, int iX= -1,int first=0,int last =0);
    virtual void FixLeaf(Node* n, bool add=false);
};
class MultilabelGiniCriterium : public GiniCriterium
{

  public:
    MultilabelGiniCriterium(Data *data) : GiniCriterium(data){}

  public:
    virtual double Initialise(int first, int last, int curSplit, double gamma = 100, int iX= -1, bool FuzzFlag = false);
    virtual void FixLeaf(Node* n, bool add=false);
};
class BoostedGiniCriterium : public GiniCriterium
{
  protected:
    Node *root;

  public:
    BoostedGiniCriterium(Data *data) : GiniCriterium(data) {
    }

  public:
    virtual void FixLeaf(Node* n, bool add=false);
};
class MSECriterium : public SplitCriterium
{
  protected:
    double R[2];
    double nL, nR;
		double s1L, s1R, s2L, s2R;

  public:
    MSECriterium(Data *data) {
      SetData(data);
    }
    virtual ~MSECriterium() {
    }


  public:
    static double Impurity(double *Pop, int NumClass);

  public:
    virtual double GetLastSplitImpurity() {return R[0]+R[1];}
    virtual double Initialise(int first, int last, int curSplit, double gamma = 100, int iX= -1, bool FuzzFlag = false);
    virtual double MovePoint(int curSplit, double gamma, int iX= -1,int first=0,int last =0);
    virtual void FixLeaf(Node* n, bool add=false);
};
class MultipleMSECriterium : public SplitCriterium
{
  protected:
    double R[2];
    double nL, nR;
		std::vector<double> s1L, s1R, s2L, s2R;

  public:
    MultipleMSECriterium(Data *data) {
      SetData(data);
    }
    virtual ~MultipleMSECriterium() {
    }


  public:
    static double MultipleVariance (std::vector<double> sum, std::vector<double> sumSquares, double numberInstances, int numberObjectives);

  public:
    virtual double GetLastSplitImpurity() {return R[0]+R[1];}
    virtual double Initialise(int first, int last, int curSplit, double gamma = 100, int iX= -1, bool FuzzFlag = false);
    virtual double MovePoint(int curSplit, double gamma, int iX= -1,int first=0,int last =0);
    virtual void FixLeaf(Node* n, bool add=false);
};

typedef struct _ComboSplitNodeInfo {
  Matriz *hist_scale;
  Matriz *hist_delta;
  Matriz *instances_data;
//  int * indexes;
//  int count;
} ComboSplitNodeInfo;

class ComboCriterium : public SplitCriterium
{
  protected:
    double R[2];
    double nL, nR;
    double *PopR, *PopL;
    Matriz *hdeltaR, *hdeltaL;
    Matriz *hscaleR, *hscaleL;
    std::vector<double> *alphas;
    double *alpha;

  public:
    ComboCriterium(Data *data, std::vector<double> *alphas=0);
    virtual ~ComboCriterium();


  public:
    static double Impurity(double *Pop, int NumClass, Matriz &delta, 
                                     Matriz &scale, double *alpha, int NumBins);
    static double Entropy(double *Pop, int NumClass);

  public:
    virtual double GetLastSplitImpurity() {return R[0]+R[1];}
    virtual double Initialise(int first, int last, int curSplit, double gamma = 100, int iX= -1, bool FuzzFlag = false);
    virtual double MovePoint(int curSplit, double gamma, int iX= -1,int first=0,int last =0);
    virtual void FixLeaf(Node* n, bool add=false);
};
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
class Data
{
  public:
    static int CvTreeLabelIndex;
    static int MembershipIndex;
    static int WeightIndex;
    static int IniPosIndex;
    static int GroupIndex;
    static int MultipleSplitWs1;
    static int MultipleSplitWs2;

  public:
    void redim(int ndatos);

protected:
    void GetNames(char *nf);
    void GetNames(FILE *f);
    void GetData(char *nf, int colIni=-1, int colFin=-1);
    void GetData(FILE *f, int colIni=-1, int colFin=-1);
    std::string GetNamesAsString();
    std::string GetExampleAsString(int i, std::string sep="\t");
    void DefaultInitValues();
    void SetDefaultNames(int nvar);
    double DameRepresentacion(std::string val, int ivarnom, bool force=true);

  public:
    Data(char *filename, Instance *builder=0);
    Data(int NTotalData, int NVar, Instance *builder=0);
    Data(int NTotalData, int NVar, double **data, Instance *builder=0);
    Data(Instance *builder=0) { 
      instanceBuilder =  builder ? builder->NewInstance(0) : new BasicInstance(0);
      DefaultInitValues();
    }

    virtual ~Data();

    virtual Data* NewData()=0;
    virtual Data* Clone(int ini=-1, int fin=-1, int factor=1);
    virtual Data* SelectDataOfClass(int iClass)=0;
    virtual void AddData(Data *d2, int ini=-1, int fin=-1);
//    virtual void AddFeature(Data *d2, int ifeature=-1);
    virtual void DeleteColumn(int i_col);
    virtual void init(){;}

    static int Proportional;
    static int GetCvTreeLabelIndex(){return CvTreeLabelIndex;}
    static int GetMembershipIndex(){return MembershipIndex;}
    static int GetWeightIndex(){return WeightIndex;}
    static int GetIniPosIndex(){return IniPosIndex;}
    static int GetGroupIndex() {return GroupIndex;}

    // Scaling data
    void ScaleVariableIQ(int iVar,int scfirst, int first, int last);
    double UnScaleValueIQ(int iVar, double x);
    void UnScaleVariableIQ(int iVar,int first, int last);
    // Sorting data
    int PreSort(int first, int last);
    void CreateGroups(int nGroups, int first=-1, int last=-1, 
                                                          bool doCVFolds=false);
    void SortByClass(int first=-1, int last=-1, bool _InvertSort=false);
    void SortByGroup(int first=-1, int last=-1, bool _InvertSort=false);
    void MarkOrder();
    void ResetOrder();
    void OriginalOrder(bool _InvertSort=false);
    void SortOn(int iVar, int first=0, int last=-1, bool _InvertSort=false);
    int SortNomOn(int Membesrship, int iX, int first, int last);
    void SortCV(int cv, int begin=0, int end = -1);
    virtual int SetCV(int cv);
    int FindFirstTest(int cv, int first, int last);
    void Scramble(int begin, int end);
    void PermuteAttributeValues(int iatt);

    // Building tree
    double FindMembership(const Node* n);
    void FindMembershipTest(const Node* n);
    virtual double FindNomSplit(int first, int last, int iX, int Nmin,
                                               double NMembers, int* pNomSplit);
    virtual double FindFuzzSplit(int first, int last, int iX, int Nmin,
                                               double NMembers, FuzzySet** pfz);
    virtual double FindOrdSplit(int first, int last, int iX, int Nmin,
                              double NMembers, int *iSplit, double gamma = 100);
    double FuzzifySplit(Node* n, double chiBest, int first, int Nmin);
    virtual int FindBestSplit(Node* n, int Nmin, int* leftLast,double *chiBest,
                                      std::vector<bool> *attributes_to_use = 0);
    virtual double FindMultiSplit(int first, int last, int Nmin, 
                  double NMembers, int *use, double* a, double* c, double* min);
    virtual double EliminateVariables(int first, int last, int Nmin, 
                                  double NMembers, double* a, int* use, 
                                  double chiOld, double* fSplit, double *min);

    // Abstract methods
    virtual double Impurity(int first, int last,double tot)=0;
//    virtual double Initialise(int first, int last, int curSplit, double gamma = 100, int iX= -1, bool FuzzFlag = false)=0;
//    virtual double MovePoint(int curSplit, double gamma, int iX= -1,int first=0,int last =0)=0;
    virtual void FixLeaf(Node* n, bool add=false)=0;
    virtual void FixCVError(Node* root)=0;
    virtual double Error(Node* root, int K, int begin,int end)=0;

    virtual bool GoLeft(Node * n, int idat);
    virtual double Error2(Node* root, int K, int begin,int end);

    // Returns the degree of fuzzification
    virtual double DegreeOfFuzz(Node* root, int K); 

    virtual double Defuzzify(Node* root, int K, int begin,int end, 
                double* errQuartile = NULL ){return Error(root, K, begin, end);}
    virtual void printout(FILE * fp, Node* root,int K, int begin, int end) {
      return;
    }
//    virtual void FixNodeCVError(Node* cur)=0;
    virtual void InitialiseErrorCV()=0;
    virtual void AccumulateErrorCV(Node* root)=0;
    virtual void GetErrorCV(double* RCV, double* SE)=0;
    virtual void WriteNode(Node*,char*,int K =-1)=0;
    virtual double Score(int begin, int end,Node* root, int k);
    virtual double GetYGvscale(){return 1.0;}
    virtual double Derivs(Node* root, int K,int begin,int end){return 1.0;}
    virtual void Ponderate(int){;}
    virtual void EqualizeClassWeight(){;}

    // Access functions
    const std::string GetProblemName(){ return name;}
    int GetTotalData(){return NTotal;}
    int GetNumData(){return Ndata;}
    int GetNumVar(){ return Nvar;}
    int GetNumVarOrd(){ return Nord;}
    int GetNumVarFuz(){ return Nfuzz;}
    int GetNumVarNom(){ return Nnom;}
    int GetNumCV(){return Ncv;}
//    void SetNudmCV(int Value){ 
    void ResetCV(int Value=-1); 
    int GetNumVarsOrdSplit(){ return NVarsOrdSplit;}
    int GetNPrune(){return NPrune;}
    int GetNTest(){return NTest;}
    int GetNTrain(){return NTrain;}
    virtual void SetNTrain(int first, int last);
    virtual void SetNTrain(int n){NTrain = n; NTest = NTotal-NTrain-NSel;
      SetNGrow(NTrain/2);}
    int GetNGrow(){return NGrow;}
    void SetNGrow(int n){NGrow=n; Ndata = NGrow; NPrune=NTrain-NGrow;}
    int GetNSel(){return NSel;}
    void SetNSel(int n){NTrain = n;}
    int GetNTotal(){return NTotal;}
/*    crepo::*/VarType  GetDepVarType(){ return DepVarType;}
    void GetScale(int first, int last, double& xmin, double& xmax, double *coeff, int att);
    int GetNordfuzz(){return Nordfuzz;}

    inline Instances& GetM()            { return instances;    } //Pasamos de la encapsulacion
    inline Instance &GetInstance(int i) { return instances[i]; } //Pasamos de la encapsulacion
    double GetValueVar(int i, int j);

    inline int    GetDatClass(int i)  { return instances[i].GetClass();  }
    inline double GetDatWeight(int i) { return instances[i].GetWeight(); }
    inline int    GetDatIniPos(int i) { return instances[i].GetIniPos(); }
    inline int    GetDatGroup(int i)  { return instances[i].GetGroup();  }

    inline void SetDatClass(int i, int clase)           { instances[i].Class  = clase; }
    inline void SetDatWeight(int i, double val)         { instances[i].Weight = val;   }
    inline void SetDatGroup(int i, int igrp)            { instances[i].Group  = igrp;  }
    void SetValueVar(int i, int j, double value);

    std::string GetNameVar(int i){ return VarNames[i];}
    virtual void CalculateAndArrangeTestData(Node * root) {;}
    virtual int GetLeftLast(Node * n, int begin, int end, int K);
    virtual int Classify(int dat, double **Class, Node *root, int K, Node**finalnode=0){return -1;}
    virtual void ResetWeights(int begin=-1, int end=-1);


    bool GroupsCreated;//Indicate if data groups has been created using CreateGroups
    int *GroupCount;   //Indicate the number of elements in each group
    void SetMultSplits(bool Value) {multsplits=Value;}
    bool GetMultSplits() {return multsplits;}
    void SelectNewTrainingSet();
    void SelectNewTrainingSetEqDis();
    virtual Data *GetBootstrapSample(int N=0, bool Weighted=false, bool IncludeOobAsTest=true, std::vector<int> *idxs=0);
    virtual Data* GetBootstrapSampleStratified(int N=0)=0;
    void SaveToFile(char *FileName, int ini=-1, int fin=-1);
    void SaveToStream(std::ostream &out, int ini=-1, int fin=-1);
    static Data* DataFromFile(char *filename, Instance* builder = 0);
    static Data* LoadFromCre(char *filename, Instance* builder = 0);
    static Data* LoadFromAsc(char *filename, Instance* builder = 0);

  protected:
    Instances instances;

    Instance *instanceBuilder;
    SplitCriterium *splitCriterium;

    int Ndata,NTotal,NTest,NGrow,NTrain, NPrune, NSel;
    int n_reservados;
      // NTrain = NGrow +NPrune ;   NTotal = NTrain + NSel + NTest
      // NSel is used to select the best fuzzification scheme after optimization
      //    or to determine best cut for degree of certainty / fuzzification
    int Nvar;
    double min_in_child; // Stores the smaller number of points after split
    int Ncv;
//    double nL, nR;
    std::string name;
    bool multsplits;  // if TRUE, allow  splits on lc's of variables

    std::vector<VarType> OriginalVarType;   // the independent var can be ORD, FUZZ or NOM
    VarType DepVarType;         // the dependent varis can be ORD, FUZZ or NOM
    bool plm;       // if TRUE and dep var is ordinal, allow linear fit on leaf
    bool fuzzify;   // if TRUE, allows fuzzification of ordinal splits
    double *Gvshift, *Gvscale;  // Only used for RPLMData
    double *vscale, *vshift;
    double cartBeta, cartEps;
//    double R[2];

    // output messages
    void (*console)(char *);

    std::string* VarNames;  // names of variables
    int Nord;    // number of ordinal variables
    int Nordfuzz;    // number of ordinal + fuzzy indep variables
                     // ordinal indep variables are labelled in instances[j][i] from 
                     // i=0 to i = Nord-1
    int NVarsOrdSplit; // Number of variables to be used in ord
                       // splits(either Nord or Nordfuzz)
    int Nfuzz; // fuzzy indep variables
               // fuzzy indep variables are labelled in instances[j][i] from 
               // i=Nord to i= Nordfuzz-1
    Property* FuzzVars; // contains decorations of the fuzzy variables

    int Nnom; // number of nominal indep variables;
              // nominal indep variables are labelled in instances[j][i] from 
              // i=Nordfuzz to Nordfuzz+Nnom-1
    std::vector<std::string>* NomTerms;    // terms of nominal variables
    std::vector<DataEvents*> reg_eventos;

  public:
    void SaveAsC45(char *FileName);
    void SaveC45Data(char *DF, int ini, int fin);
    void SaveAsArff(char *fichero, int ini=-1, int fin=-1);

    void AddListener(DataEvents *de){reg_eventos.push_back(de);}
    void DeleteListener(DataEvents *de){;}

    std::string GetTerm(int i, int ivar) {return NomTerms[ivar][i];}
    std::vector<std::string> GetTermVector(int ivar) {return NomTerms[ivar];}
    int GetNumTerms(int ivar) {return (int)NomTerms[ivar].size();}

    void SetSplitCriterium(SplitCriterium *splitCriterium);
    void SetSplitCriterium(std::string splitClassName);

    double Min(int icol);
    double Max(int icol);

    Data *SinteticDataNaiveBayes(int ntot, int ini, int fin);
    Data *SinteticDataNN(int ntot, int ini, int fin, bool ByClass=true);
    Data *BorderData(int ini, int fin, int knn=5);

    double DistanciaEjem(int d1, int d2);
};


class NomData:public Data
{
  public:
    void init();

  public:
    NomData(char *filename, Instance *builder=0);
    NomData(int NTotalData, int NVar, Instance *builder=0);
    NomData(int NTotalData, int NVar, double **data, Instance *builder=0);
    NomData(Instance *builder=0) : Data(builder) { }

    virtual ~NomData();

    virtual Data* NewData(){return new NomData(instanceBuilder);}
    virtual Data* Clone(int ini=-1, int fin=-1, int factor=1);
    virtual Data* SelectDataOfClass(int iClass);

    void FindPops(double* HoldPop, int first, int last);
    int  WhichClass(double* HoldPop, int NumClass=-1);
    inline double FindImpurity(double* HoldPop,int first =0,int last=0);
    virtual double Impurity(int first, int last,double tot);
//    virtual double Initialise(int first, int last, int curSplit, double gamma, int iX, bool FuzzFlag);
//    inline virtual double MovePoint(int curSplit, double gamma, int iX,int first,int last);
    virtual void FixLeaf(Node* n, bool add=false);
    virtual void FixCVError(Node* root);
    virtual void InitialiseErrorCV();
    virtual void AccumulateErrorCV(Node* root);
    virtual double Derivs(Node* root, int K, int begin,int end);
    virtual double Defuzzify(Node* root, int K, int begin, int end, double* errQuartile = NULL);
    virtual double Error(Node* root, int K, int begin,int end);
    double Fuzziness(int i, double fuzz);
    double CrispError(Node* root, int K, int begin,int end);
    virtual double DegreeOfFuzz(Node* root, int K); // Reassigns labels to leafs after fuzzification and returns degree of fuzzification
    void Filter(int i, double* Class);
    void AssignClass(int labelClass, Node* root,int K, double* Item,int begin,int end);
    double OrdCost(Node* root, int K, int begin, int end);
    virtual void GetErrorCV(double* RCV, double* SE);
    virtual void WriteNode(Node*,char*,int K=-1);
    virtual double Score(int begin, int end,Node* root, int K);
    virtual void Ponderate(int);
    virtual void EqualizeClassWeight();
    int Classify(int idat, double **Class, Node *root, int K, Node**finalnode=0);
    void CalculateAndArrangeTestData(Node * root);
    bool Getatdrawchosefathers() {return atdrawchosefathers;};
    void Setatdrawchosefathers(bool value) { atdrawchosefathers = value; }
    double GetPopClass(int j){return(Pop[j]);};

  protected:
    double NextDifferentValue(int var, int index);

  public:
    NomData *GenerateSinteticData(int order, int ini, int fin);
    NomData *GenerateSinteticData2(int ntot, int ini, int fin);
    NomData *ClassNoise(double prob);
    virtual void SubstituteClassLabels(Data *data, int attribute_idx=-1);

    virtual Data* GetBootstrapSampleStratified(int N=0);
    static Data* GenerateMIDataSet(NomData *data, int micolumn);

//    virtual void DeleteColumn(int i_col);
    int AddClass(std::string name);
    void RemoveInstancesOfClass(int class_index);

  public:
    int NumClass;
    int degree; // For optimization: exponent of the leaf membership
          // in error function.
    bool multree; // Flag to build several binary trees (only multivalued nominal)

  protected:
    double* Pop;
    int iClass;
    double rcv;
    int iClassL, iClassR;
//    double *PopR, *PopL;
    bool ponder;
    //If there are more then one class with the same number of samples in a node
    bool atdrawchosefathers;//then the class of the father, if possible is chosen
};

class MINomData:public NomData
{
  protected:
    int NMITotal;
    int NMITrain;
    std::vector<double> Weights;
    
  public:
    MINomData(char *filename, int mi_attrib, Instance *builder=0);
    MINomData(int NTotalData, int NVar, Instance *builder=0) : 
                                                  NomData(NTotalData, NVar, builder) { }
    MINomData(int NTotalData, int NVar, double **data, Instance *builder=0) : 
                                            NomData(NTotalData, NVar, data, builder) { }
    MINomData(Instance *builder=0)          : NomData(builder)                       { }

    virtual ~MINomData();

    virtual Data* NewData(){return new MINomData(instanceBuilder);}
    virtual Data *GetBootstrapSample(int N=0, bool Weighted=false, bool IncludeOobAsTest=true, std::vector<int> *idxs=0);
    virtual void AddData(Data *d2, int ini=-1, int fin=-1);
    virtual void ResetWeights(int begin=-1, int end=-1);
    virtual void SubstituteClassLabels(Data *data, int attribute_idx=-1);
    virtual Data* Clone(int ini=-1, int fin=-1, int factor=1);
    
    virtual void SetNTrain(int n);
    int GetMITrain() {return NMITrain;}
    int GetMITotal() {return NMITotal;}
    int CountMIInstances(int first, int last);
    double GetWeight(int i) {return Weights[i];}
    void SetWeight(int i, double value) {Weights[i] = value;}

};


#endif //__DATA_H

