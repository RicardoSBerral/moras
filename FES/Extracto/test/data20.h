#ifndef __DATA_H
#define __DATA_H

#define MIN_BELONG_LOW 1e-16
#define MIN_BELONG_HIGH 0.5

#include "policy20.h"

//1class CConsole;
class Node;


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
    Data(DataSet* dataset/*1, CConsole**/, void (*console)(char *));
    Data(std::istream &in, void (*console)(char *));
    Data(int NTotalData, int NVar);
    Data(int NTotalData, int NVar, double **data);
    virtual ~Data();

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
    void CreateGroups(int nGroups, int first, int last);
    void SortOn(int iVar, int first, int last, int _InvertSort=0);
    int SortNomOn(int Membesrship, int iX, int first, int last);
    void SortCV(int cv, int begin=0, int end = -1);
    virtual int SetCV(int cv);
    int FindFirstTest(int cv, int first, int last);
    void Scramble(int begin, int end);
    // Building tree
    double FindMembership(const Node* n);
    void FindMembershipTest(const Node* n);
    virtual double FindNomSplit(int first, int last, int iX, int Nmin,double NMembers, int* pNomSplit);
    virtual double FindFuzzSplit(int first, int last, int iX, int Nmin,double NMembers, FuzzySet** pfz);
    virtual double FindOrdSplit(int first, int last, int iX, int Nmin,double NMembers, int *iSplit, double gamma = 100);
    double FuzzifySplit(Node* n, double chiBest, int first, int Nmin);
    virtual int FindBestSplit(Node* n, int Nmin, int* leftLast,double *chiBest);
    virtual double FindMultiSplit(int first, int last, int Nmin, double NMembers, int *use, double* a, double* c, double* min);
    virtual double EliminateVariables(int first, int last, int Nmin, double NMembers, double* a,
      int* use, double chiOld, double* fSplit,double *min);
    // Abstract methods
    virtual double Impurity(int first, int last,double tot)=0;
    virtual double Initialise(int first, int last, int curSplit, double gamma = 100, int iX= -1, bool FuzzFlag = false)=0;
    virtual double MovePoint(int curSplit, double gamma, int iX= -1,int first=0,int last =0)=0;
    virtual void FixLeaf(Node* n)=0;
    virtual void FixCVError(Node* root)=0;
    virtual double Error(Node* root, int K, int begin,int end)=0;
    virtual double Error2(Node* root, int K, int begin,int end);
    virtual double DegreeOfFuzz(Node* root, int K); // Returns the degree of fuzzification
    virtual double Defuzzify(Node* root, int K, int begin,int end, double* errQuartile = NULL ){return Error(root, K, begin, end);}
    virtual void printout(FILE * fp, Node* root,int K, int begin, int end) {return;}
//    virtual void FixNodeCVError(Node* cur)=0;
    virtual void InitialiseErrorCV()=0;
    virtual void AccumulateErrorCV(Node* root)=0;
    virtual void GetErrorCV(double* RCV, double* SE)=0;
    virtual void WriteNode(Node*,char*,int K =-1)=0;
    virtual double Score(int begin, int end,Node* root, int k);
    virtual double GetYGvscale(){return 1.0;}
    virtual double Derivs(Node* root, int K,int begin,int end){return 1.0;}
    virtual void Ponderate(int){;}
    // Access functions
    int GetTotalData(){return NTotal;}
    int GetNumData(){return Ndata;}
    int GetNumVar(){ return Nvar;}
    int GetNumVarOrd(){ return Nord;}
    int GetNumVarFuz(){ return Nfuzz;}
    int GetNumVarNom(){ return Nnom;}
    int GetNumCV(){return Ncv;}
    void SetNumCV(int Value){Ncv = Value;}
    int GetNumVarsOrdSplit(){ return NVarsOrdSplit;}
    int GetNPrune(){return NPrune;}
    int GetNTest(){return NTest;}
    int GetNTrain(){return NTrain;}
    void SetNTrain(int n){NTrain = n; NTest = NTotal-NTrain-NSel;
      SetNGrow(NTrain-NTrain/2);}
    int GetNGrow(){return NGrow;}
    void SetNGrow(int n){NGrow=n; Ndata = NGrow; NPrune=NTrain-NGrow;}
    int GetNSel(){return NSel;}
    void SetNSel(int n){NTrain = n;}
    int GetNTotal(){return NTotal;}
    crepo::VarType  GetDepVarType(){ return DepVarType;}
    void GetScale(int first, int last, double& xmin, double& xmax, double *coeff, int att);
    int GetNordfuzz(){return Nordfuzz;}
    inline double GetValueVar(int i, int j){return m[i][j];}
    inline int GetDatClass(int i){return (int)m[i][Nvar-1];}
    inline double GetDatWeight(int i){return m[i][Nvar+WeightIndex];}
    inline int GetDatIniPos(int i){return (int)m[i][Nvar+IniPosIndex];}
    inline void SetValueVar(int i, int j, double value){m[i][j] = value;}
    CString GetNameVar(int i){ return VarNames[i];}
    virtual void CalculateAndArrangeTestData(Node * root) {;}
    virtual int GetLeftLast(Node * n, int begin, int end, int K);
    virtual int Classify(int dat, double **Class, Node *root, int K){return -1;}
    void ResetWeights(int begin=-1, int end=-1);
    bool GroupsCreated;//Indicate if data groups has been created using CreateGroups
    int *GroupCount;   //Indicate the number of elements in each group
    void SetMultSplits(bool Value) {multsplits=Value;}
    bool GetMultSplits() {return multsplits;}
    void SelectNewTrainingSet();
    void SaveToFile(char *FileName, int ini=-1, int fin=-1);
    void SaveToStream(std::ostream &out, int ini=-1, int fin=-1);
    static Data* DataFromFile(char *filename);
    static Data* LoadFromCre(char *filename);
    static Data* LoadFromAsc(char *filename);
  protected:
    double** m;
    int Ndata,NTotal,NTest,NGrow,NTrain, NPrune, NSel;
      // NTrain = NGrow +NPrune ;   NTotal = NTrain + NSel + NTest
      // NSel is used to select the best fuzzification scheme after optimization
      //    or to determine best cut for degree of certainty / fuzzification
    int Nvar;
    double min_in_child; // Stores the smaller number of points after split
    int Ncv;
    double nL, nR;

    bool multsplits;  // if TRUE, allow  splits on lc's of variables

    crepo::VarType DepVarType;   // the dependent variable can be ORD, FUZZ or NOM
    bool plm;      // if TRUE and dep var is ordinal, allow linear fit on leaf
    bool fuzzify;   // if TRUE, allows fuzzification of ordinal splits
    double *Gvshift, *Gvscale; // Only used for RPLMData
    double *vscale, *vshift;
    double cartBeta, cartEps;
    double R[2];
//1    CConsole* console;
    void (*console)(char *);

    string* VarNames;  // names of variables
    int Nord;      // number of ordinal variables
    int Nordfuzz;    // number of ordinal + fuzzy indep variables
              // ordinal indep variables are labelled in m[j][i] fr0m i=0 to i = Nord-1
    int NVarsOrdSplit;      // Number of variables to be used in ord splits(either Nord or Nordfuzz)
    int Nfuzz;      // fuzzy indep variables
              // fuzzy indep variables are labelled in m[j][i] from i=Nord to i= Nordfuzz-1
    Property* FuzzVars; // contains decorations of the fuzzy variables


    int Nnom;      // number of nominal indep variables;
              // nominal indep variables are labelled in m[j][i] from i=Nordfuzz to Nordfuzz+Nnom-1
    CStringArray* NomTerms;    // terms of nominal variables

    public:
    void SaveAsC45(char *FileName);
};

class NomData:public Data
{
  public:
    NomData(DataSet* ds/*1, CConsole* pt*/,void(*_console)(char*));
    NomData(std::istream &in, void (*console)(char *));
    NomData(int NTotalData, int NVar);
    NomData(int NTotalData, int NVar, double **data);
    virtual ~NomData();
    int NumClass;
    int degree; // For optimization: exponent of the leaf membership
          // in error function.
    bool multree; // Flag to build several binary trees (only multivalued nominal)
    double GetPopClass(int j){return(Pop[j]);};
    void FindPops(double* HoldPop, int first, int last);
    int  WhichClass(double* HoldPop);
    double FindImpurity(double* HoldPop,int first =0,int last=0);
    virtual double Impurity(int first, int last,double tot);
    virtual double Initialise(int first, int last, int curSplit, double gamma, int iX, bool FuzzFlag);
    virtual double MovePoint(int curSplit, double gamma, int iX,int first,int last);
    virtual void FixLeaf(Node* n);
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
    int Classify(int idat, double **Class, Node *root, int K);
    void CalculateAndArrangeTestData(Node * root);
    bool Getatdrawchosefathers() {return atdrawchosefathers;};
    void Setatdrawchosefathers(bool value) { atdrawchosefathers = value; }
  protected:

    double* Pop;

    int iClass;
    double rcv;
    int iClassL, iClassR;
    double *PopR, *PopL;
    bool ponder;
    //If there are more then one class with the same number of samples in a node
    bool atdrawchosefathers;//then the class of the father, if possible is chosen
    double NextDifferentValue(int var, int index);
  public:
    NomData *GenerateSinteticData(int order, int ini, int fin);
    NomData *GenerateSinteticData2(int ntot, int ini, int fin);
};
#endif //__DATA_H
