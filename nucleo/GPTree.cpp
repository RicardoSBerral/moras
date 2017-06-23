//---------------------------------------------------------------------------
// 
// 
//---------------------------------------------------------------------------
#include "GPTree.h"
#include "tree20.h"
#include "node20.h"
//#include "UnHilo.h"

using namespace std;
//---------------------------------------------------------------------------
CARTTreeReg *CARTTreeReg::reg = CARTTreeReg::autoreg();
DecisionStumpReg *DecisionStumpReg::reg = DecisionStumpReg::autoreg();
RandomForestTreeReg *RandomForestTreeReg::reg = RandomForestTreeReg::autoreg();
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void salypimienta2(char *texto)
{
//  if (texto[strlen(texto)-1]=='\n') printf(texto);
//  else printf("%s\n", texto);
  return;
}
//---------------------------------------------------------------------------
bool GPTree::ComoArticulo = true;

GPTree::GPTree(void (*salida)(char *), int _NMin):Classifier()
{
  tree = new Tree((salida) ? salida : salypimienta2);
  NMin = _NMin;
}

GPTree::~GPTree()
{
  delete tree;
}
//---------------------------------------------------------------------------

void GPTree::Init(Data * data)
{
  Classifier::Init(data);
  if (!data->GroupsCreated)
    CreateGroups();
    
  CountGroup[0] = data->GroupCount[0];
  CountGroup[1] = data->GroupCount[1];
  ((NomData*)data)->Setatdrawchosefathers(true);
}
//---------------------------------------------------------------------------

double GPTree::Error(int first, int last)
{
  if (first > last) return -1.0;
  return data->Error2(tree->GetRoot(), -2, first, last )/(last-first+1);
}
//---------------------------------------------------------------------------

int GPTree::Classify(int ElementIndex)
{
  return data->Classify(ElementIndex, 0, tree->GetRoot(), -2);
}
//---------------------------------------------------------------------------

void GPTree::Build(Data *data, FuncionDeProgreso *fp)
{
  FILE *kk=fopen("podilla.txt", "wt");
  fclose(kk);
  int iGroup=0, prevN;
  //unsigned i;
  NumIter = 0;
  Init(data);
//  hilillo->LimpiarArbolExt();
  if (data) {
    prevN = 0;
    vector<int> nnodos;
    SetGroup(0);
    do {
      int K=-2;
//      double alpha, error;
      prevN = tree->CountNodes();
/*      for(i=0;i<nnodos.size();i++) {
        if (prevN==nnodos.at(i)) {
          WriteOrder(data);
          fclose(this->f);
          this->f = 0;
          throw exception();
        }
      }*/
      nnodos.push_back(prevN);
      tree->Build(data, NMin, iGroup);
      if (tree->Tmax<tree->GetRoot()->T) tree->Tmax = tree->GetRoot()->T;
//      hilillo->MostrarArbolExt((NomData*)data, tree);
      iGroup++;
      SetGroup(iGroup%2);
      NumIter++;
      tree->ThrowData(data, 0, data->GetNGrow()-1, -2);
//      hilillo->MostrarArbolExt((NomData*)data, tree);
      tree->Prune2(0, data->GetNumData(), K, ComoArticulo);
//      hilillo->MostrarArbolExt((NomData*)data, tree);
    }
    while(prevN!=tree->CountNodes());
  }
  estimError = GroupsError();
}
//---------------------------------------------------------------------------

Tree* GPTree::GetTree()
{
  return tree;
}
//---------------------------------------------------------------------------

double GPTree::GroupsError()
{
  double error[2] = {0.0, 0.0};
  double TotWeights[2] = {0.0, 0.0};
  double EstimError;
  int iClass = data->GetNumVar()-1;
//  int iMemership = iClass + 4;
  int iWeight = iClass + 5;
  int clase;
  int K = -2;
  for (int iGrupo=0;iGrupo<2;iGrupo++) {
    SetGroup(iGrupo);
    for(int j=0;j<CountGroup[iGrupo];j++)
      TotWeights[iGrupo]+=data->GetValueVar(j, iWeight);
    tree->ThrowData(data, 0, CountGroup[iGrupo]-1, K);
    Node *cur=tree->GetRoot();
    while (cur->child) cur = cur->child;
    while(cur) {
      if (cur->IsLeaf(K) && iGrupo!=cur->GetGeneratedByGroup()) {
        for(int j=cur->first;j<=cur->last;j++) {
          clase = data->GetDatClass(j);//data->GetValueVar(j, iClass);
          if ( cur->iClass!=clase )
//            error[iGrupo] += data->GetValueVar(j, iMemership);
            error[iGrupo] += data->GetValueVar(j, iWeight);
        }

      }
      cur = cur->nextDown(K);
    }
  }
  EstimError =( error[0]/TotWeights[0]) + (error[1]/TotWeights[1]);
  if (EstimError>1) throw exception();
  return EstimError;
}
//---------------------------------------------------------------------------

void GPTree::CreateGroups()
{
  data->CreateGroups(2,0,data->GetNTrain()-1);
  CountGroup[0] = data->GroupCount[0];
  CountGroup[1] = data->GroupCount[1];
}
//---------------------------------------------------------------------------

void GPTree::SetGroup(int i)
{
  if (i==0) {
    data->SortOn(data->GetNumVar()+Data::GetGroupIndex(), 0, data->GetNTrain()-1, false);
    data->SetNGrow(CountGroup[0]);
  }
  else {
    data->SortOn(data->GetNumVar()+Data::GetGroupIndex(), 0, data->GetNTrain()-1, true);
    data->SetNGrow(CountGroup[1]);
  }
}
//---------------------------------------------------------------------------
//Funciones para guardar y leer de fichero
void GPTree::Guardar(std::ostream &sal, int version)
{
  Classifier::Guardar(sal, version);
  tree->Guardar(sal, version);
}
//---------------------------------------------------------------------------
void GPTree::Leer(std::istream &in, int version)
{
  Classifier::Leer(in, version);
  if (!tree) tree = new Tree(salypimienta2);
  tree->Leer(in, version);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//-1=cada arbol decide, -2=DecisionStump. otro=Valor de K, 0 sin poda
int DecisionTree::K = -1; 

//---------------------------------------------------------------------------
DecisionTree::DecisionTree(void (*salida)(char *), int _NMin, bool PruneTree, 
              bool MultiSplits, bool SE_0) : Classifier()
{
  tree = new Tree((salida) ? salida : salypimienta2);
  NMin = _NMin;
  this->PruneTree = PruneTree;
  this->MultiSplits = MultiSplits;
  this->SE_0 = SE_0;
}

DecisionTree::~DecisionTree()
{
  delete tree;
}
//---------------------------------------------------------------------------

int DecisionTree::GetK()
{
  //return UsePrunedTrees==2 ? 0 : tree->GetK();
  return K==-1 ? tree->GetK() : (K==-2 ? tree->GetRoot()->child->K-1 : K);
}
//---------------------------------------------------------------------------

void DecisionTree::Init(Data * data)
{
  Classifier::Init(data);
  data->SetMultSplits(MultiSplits);
}
//---------------------------------------------------------------------------

double DecisionTree::Error(int first, int last)
{
  double tot_peso = 0.0;
  for(int i=first; i<=last;i++) tot_peso += data->GetDatWeight(i);
  double error2 = data->Error2(tree->GetRoot(), GetK(), first, last );
  return error2/tot_peso;
}
//---------------------------------------------------------------------------

int DecisionTree::Classify(int ElementIndex)
{
  return data->Classify(ElementIndex, 0, tree->GetRoot(), GetK(),
                                                          &lastClassifyingNode);
}
//---------------------------------------------------------------------------

double DecisionTree::Average(int ElementIndex)
{
  data->Classify(ElementIndex, 0, tree->GetRoot(), GetK(), &lastClassifyingNode);
  return lastClassifyingNode->fClass;
}

double DecisionTree::Average(int ElementIndex, int AttributeIndex)
{
  if (MultiobjectiveInstance::GetNumMultiobjectives() <= 0) {
    throw std::logic_error("Don't specify the attribute index unless dealing with multiobjective instances");
  }
  data->Classify(ElementIndex, 0, tree->GetRoot(), GetK(), &lastClassifyingNode);
  return lastClassifyingNode->NodePop[AttributeIndex - (data->GetNumVar() - MultiobjectiveInstance::GetNumMultiobjectives())];
}
std::vector<double> DecisionTree::MultipleAverage(int ElementIndex)
{
  data->Classify(ElementIndex, 0, tree->GetRoot(), GetK(), &lastClassifyingNode);
  double* array = lastClassifyingNode->NodePop;
  std::vector<double> ret;
  ret.assign(array, array + MultiobjectiveInstance::GetNumMultiobjectives());
  return ret;
}
//---------------------------------------------------------------------------
vector<double> DecisionTree::UnnormalizedDistribution(int ElementIndex)
{
  vector<double> distrib( ((NomData*)data)->NumClass, 0.0 );
  
  data->Classify(ElementIndex, 0, tree->GetRoot(), GetK(), &lastClassifyingNode);

  for(int i = 0; i < ((NomData*)data)->NumClass; i++) { 
    distrib[i] = lastClassifyingNode->NodePop[i];
  }

  return distrib;
}
//---------------------------------------------------------------------------
string DecisionTree::Info(int value)
{
  char *kk;

  kk = value > 2 ? new char[2000000] : new char[4096];

  if ( value == 33 ) {
    sprintf(kk, "%d ", tree->CountNodes(GetK()));
    return string(kk);
  }

  sprintf(kk, "%d", tree->CountNodes(GetK()));
  string res = Classifier::Info(value);
  res = res + " -Num nodos: " + kk; 
  sprintf(kk, "%d", GetK());
  res = res + " -K: " + kk; 
  if (value>2) {
    kk[0] = '\0';
    tree->Write(data, GetK(), kk);
    res = res + "\n" + kk;
  }

  delete []kk;

  return res;
}
//---------------------------------------------------------------------------

Tree* DecisionTree::GetTree()
{
  return tree;
}
//---------------------------------------------------------------------------
//Funciones para guardar y leer de fichero
void DecisionTree::Guardar(std::ostream &sal, int version)
{
  Classifier::Guardar(sal, version);
  tree->Guardar(sal, version);
}
//---------------------------------------------------------------------------
void DecisionTree::Leer(std::istream &in, int version)
{
  Classifier::Leer(in, version);
  if (!tree) tree = new Tree(salypimienta2);
  tree->Leer(in, version);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
CARTTree::CARTTree(void (*salida)(char *), int _NMin, bool PruneTree, 
              bool MultiSplits, bool SE_0) : 
                    DecisionTree(salida, _NMin, PruneTree, MultiSplits, SE_0)
{
}

CARTTree::~CARTTree()
{
}
//---------------------------------------------------------------------------

void CARTTree::Build(Data *data, FuncionDeProgreso *fp)
{
  if (tree->GetRoot()) {
    delete tree;
    tree = new Tree(salypimienta2);
  }
  Init(data);
  data->SetNGrow(data->GetNTrain());
  tree->CART(data, NMin, PruneTree, SE_0);
  tree->GetRoot()->T = tree->SizeBest;
  estimError = tree->RMSres/data->GetNGrow();
  data->SetMultSplits(false);
  tree->FreeAuxData();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
RandomForestTree::RandomForestTree(void (*salida)(char *), int _NMin,
               int NumberOfAttributes, bool MultiSplits, int AttsToRandomize) :
                    DecisionTree(salida, _NMin, PruneTree, MultiSplits, SE_0)
{
  if (AttsToRandomize > 0) {
    atts = new vector<int>;
    for (int i = 0; i < AttsToRandomize; i++) {
      atts->push_back(i);
    }
  }
  else {
    atts = 0;
  }

  this->NumberOfAttributes = NumberOfAttributes;
}

RandomForestTree::~RandomForestTree()
{
  if (atts) delete atts;
}
//---------------------------------------------------------------------------

void RandomForestTree::Build(Data *data, FuncionDeProgreso *fp)
{
  if (tree->GetRoot()) {
    delete tree;
    tree = new Tree(salypimienta2);
  }
  Init(data);
  data->SetNGrow(data->GetNTrain());
  tree->BuildRandomForestMember(data, NMin, NumberOfAttributes, atts);
//  tree->GetRoot()->T = tree->SizeBest;
  estimError = tree->RMSres/data->GetNGrow();
  data->SetMultSplits(false);
  tree->FreeAuxData();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
DecisionStump::DecisionStump(void (*s)(char *), int _NMin, bool MultiSplits) :
                              DecisionTree(s, _NMin, false, MultiSplits, false)
{

}

DecisionStump::~DecisionStump()
{

}
//---------------------------------------------------------------------------

void DecisionStump::Build(Data *data, FuncionDeProgreso *fp)
{
  if (tree->GetRoot()) {
    delete tree;
    tree = new Tree(salypimienta2);
  }
  Init(data);
  data->SetNGrow(data->GetNTrain());
  tree->BuildDecisionStump(data, NMin);
//  tree->GetRoot()->T = tree->SizeBest;
  estimError = tree->RMSres/data->GetNGrow();
  data->SetMultSplits(false);
}
//---------------------------------------------------------------------------

