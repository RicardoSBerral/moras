//---------------------------------------------------------------------------

#ifndef Tree20H
#define Tree20H
#include<ostream>
#include<istream>
#include<map>
#include<vector>

//---------------------------------------------------------------------------
class Data;
class Node;
//---------------------------------------------------------------------------
class PruneInfo
{
  public:
    Node *node;
    std::vector<int> min_cost_left;
    std::vector<int> min_cost_right;
    std::vector<double> min_cost;
    std::vector<double> min_cost_alpha;

  public:
    PruneInfo(Node *n);

    double error(int size);
    void merge(int left, int right, double error);
};
//---------------------------------------------------------------------------
class Tree
{
  public:
    double sbest;  // Test error
    double RMSres; // resubstitution error
    int SizeBest;
    double bestSlope;

  public:
    Tree(/*CConsole *cc*/void (*_console)(char *)=0);
    Tree(Tree *t);
    ~Tree();

  public:
    void BuildRandomForestMember(Data* data,int Nmin, int num_random_vars,
                                             std::vector<int> *atts_to_use=0);
    void BuildDecisionStump(Data* data,int Nmin);
    void Build(Data* data,int Nmin, int GroupNumber=-1);
    void CreateAndFillChildren(Node* n, int last, Data* data, int Nmin);
    void FreeSubTree(Node* n, int K=-1);
    void FreeAuxData();

    void RecomputeNodeStats(Data *data, int first, int last, bool Add=false);

      // For CART
    void CART(Data* data,int Nmin, bool PruneTree=true, bool SE_0=false);
    double FindAlpha(int Ndata,int K=-1, bool UseSqrt=false);
    void IndepPrune();
    void Prune(double alpha, int K=-1);
    void DeleteNodesAndResetK(Node *n, int K);
    void Prune2(double alpha, int ndata, int K=-1, bool ComoArticulo=true);
  
      // For Clayton:
    std::map<Node*, PruneInfo*> MinimumCostTrees();
    std::vector<double> SelectMinimumCostTrees(std::map<Node*, PruneInfo*> min_cost_trees);
    void UpdateMinimumCostTrees(std::map<Node*, PruneInfo*> min_cost_trees, 
                                          Node* cur, int index, double alpha);
    void PruneMCT(Data* data, int Nmin, bool SE_0);
    void PruneMCT2(double alpha);
    void FreeSubTreeMCT(Node* n, double alpha);
  
    void Write(Data* data,int K = -1, char *buf=0);

    // For parameter optimisation
    double OptimiseParams(Data* data);
//1    CConsole* GetConsole(){return console;}
    Node* GetRoot(){return root;}
    int GetK(){return _Kfinal;}
    void PrepareOpt(double* a, int K, Data*);
    void InitOpt(double*,int,int);
    void CollectOpt(double*,int);
    int CountParams(int K);
    int CountNodes(int K=-1);
    void IndepPrune2();
    double ThrowData(Data * data, int begin, int end, int K, Node *n=0);
    void Assign(Tree * t);

  protected:
    void frprmn(double p[], int n, double ftol, int *iter, double *fret,
               double (*func)(double []), void (*dfunc)(double [], double []));

    int maxdepth; // maximum depth of tree
    Node* root;
//1    CConsole* console;
    void (*console)(char *);
    int _Kfinal;
    int Ndata;
    //If there are more then one class with the same number of 
    //samples in a node then the class of the father, if possible is chosen
    bool atdrawchosefathers;

  public:
    int Tmax;

    //Funciones para guardar y leer de fichero
    void Guardar(std::ostream &salida, int version=0);
    void Leer(std::istream &in, int version=0);
};

#endif
