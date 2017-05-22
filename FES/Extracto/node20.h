#ifndef __NODE_H
#define __NODE_H

#include "definitions20.h"
#include <string>
#include <vector>

class FuzzySet;
class ostream;
class istream;
class Instance;

class sigmoid{
  public:
    sigmoid(double v1, double v2){ a = v1; b = v2;}
    double y(double x);
    double dy_da(double x);
    double dy_db(double x);
    double a,b;
};

class Data;

//TH Trying to reduce the size of class Node by removing attributes
typedef struct _NodeData
{
    //
    //  Points assigned to node
    //
    int firstTest, lastTest;   // From the test set

    double min_in_child;    // Minimum number of points in child node


    double normcoeff;  // Norm of coefficients vector
    double *plm;    // Coefficients for piecewise linear models

    FuzzySet* FuzzSplit; // Pointer to FuzzySet that determines fuzzy split

    double RPrune,RPruneSubTree; // Error in the pruning set when  node is considered as leaf
    double Rleaf,RsubTree, RleafCV, RsubTreeCV;
    double s1leaf,s12leaf, s1,s12;
    double fRTest;
    int fuzzindex; // Indicates the degree of fuzzification of given ordinal variable
    bool CrispFlag;  // TRUE if all previous splits are crisp FALSE otherwise

    sigmoid* sig;      // Fuzzy split

    double tempY;
    double tempY2;
    double tempMU;
    double tempMU2;
    double cumMU;
    double derL;   // Derivative wrt fitting constant (regression only)
    double dera;   // Derivative wrt inverse width of sigmoid.
    double derb;   // Derivative wrt center of sigmoid.
    double* dercoeff;  // Derivative wrt parameters for multivariate splits

    double** storage;    // For Optimization of nominals
    double MembershipTotal;
    double entropy;
    double norm;
    bool IsLeafAfterPrune;  //Indicates if a Node has been a leaf after a call to Tree::Prune2
    bool DoesntWantSplitting;//If a Node has been a Leaf twice it doesnt need to be grow again
    int GeneratedByGroup;

} NodeData;

class Node 
{
  public:
    Node(Node* par, int M, int max);
    Node(Node* n, Node* _parent=0);  // Copy constructor
    ~Node();
    //
    // Navigation through tree
    //
    Node* nextUp();
    Node* nextDown(int K=-1);
    //
    //  Test to check whether node is leaf
    //
    inline bool IsLeaf(int K){ 
      if (!child) return true;
      return ( !(child->K < 0 || child->K > K) || depth==max_depth );
    }
    static int max_depth;
    //
    //  Used for pruning
    //
    Node *GetRoot();

    void ComputeAlpha(int Ndata, Node *root=0, bool UseSqrt=false);

    // 
    //  Fuctions for the computation of derivatives
    //
    void FixDers(Instance& m, double fac,int K,int NVarsOrdSplit);  // For Ordinal variables
    void FixDers(Instance& m,int K,int NVarsOrdSplit); // For Nominal variables

    double ParamMem(Instance& m, int NVarsOrdSplit);
    
    
    int CountParams(int K);
    void InitParams(double** p,int K, Data*);
    void FixParams(double** p,int K,int);
    void ColParams(double** p,int K);
    void FreeAuxData();

  public:
    NodeData *d;

    //  Points assigned to node
    int first,last;        // From the training set

    int iClass;            // Classification givenb by node
    double* NodePop;       // Contains population of classes for Nominal dependent variables
    double fClass;

    int att;        // Attribute for univariate split
    double *coeff;  // Coefficients for multivariate splits
    int sizecoeff;  // Size of coefficients vector 
    double fSplit;  // Threshold for crisp split
    int NomSplit;   // Contains subset that determines nominal split
                    // If = 0, then split is ordinal or no split
    VarType SplitType;    // Split can be ORD,FUZZ,NOM, or ERR

    Node *parent, *child, *sib;
    int depth; // Number of parents a node has
    double alpha;
    int T;
    int K;

    void *info;

  public:
    int GetGeneratedByGroup(){ return d ? d->GeneratedByGroup : -1;}
    void SetGeneratedByGroup(int Value){ if (d) d->GeneratedByGroup = Value;}
    void Assign(Node * n);

  public:
    //Funciones para guardar y leer de fichero
    void Guardar(std::ostream &salida, int version=0);
    void Leer(std::istream &in, int version=0);
    bool IsUnivariateSplit();

};

#endif // __NODE_H sentinal


