//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#include "C45Tree.h"
#include "extern.i"
#include "buildex.i"
#include "classify.h"
//#include "besttree.h"
#include "c4.5.h"
#include "getnames.h"
#include "getdata.h"
#include "trees.h"
#include "build.h"
#include "prune.h"
#include "besttree.h"
#include "Utils.h"

using namespace std;

//InicializaciOn de miembros estAticos
C45DataSets *C45DataSets::def = 0;
Data * C45Tree::lastdataset = 0;

//0=cada arbol decide, 1=Todos arboles podados, 2=Todos arboles sin podar
int C45Tree::UsePrunedTrees = 0; 
//---------------------------------------------------------------------------
//-----------------------------------------------------------  C45Tree  -----
//---------------------------------------------------------------------------
C45Tree::C45Tree(bool UsePruned):Classifier()
{
  c45tree = 0;
  c45tree_unpruned = 0;
  this->UsePruned = UsePruned;
}

C45Tree::~C45Tree()
{
  if (c45tree) c45::ReleaseTree(c45tree);
  if (c45tree_unpruned) c45::ReleaseTree(c45tree_unpruned);
}
//---------------------------------------------------------------------------

bool C45Tree::GetUsePruned()
{
  return UsePrunedTrees == 0 ? UsePruned : UsePrunedTrees==1;
}
//---------------------------------------------------------------------------

c45::Tree C45Tree::GetTree()
{
  return GetUsePruned() ? c45tree : c45tree_unpruned;
}
//---------------------------------------------------------------------------

void C45Tree::Init(Data * data)
{
  Classifier::Init(data);
  SetData(data);
}
//---------------------------------------------------------------------------

int C45Tree::Classify(int ElementIndex)
{
  c45::Description ejm = GetC45DataSets()->GetDescription(data, ElementIndex);  
  int c = c45::Category(ejm, GetTree());
  return c;
/*  int clase = c45::Category(c45::Item[ElementIndex], c45tree);
  if (ElementIndex==0) {
    SetData(data);*/
  /*  printf("%g %g %d\n", 
      c45::Item[ElementIndex][0]._cont_val,
      c45::Item[ElementIndex][1]._cont_val,
      c45::Item[ElementIndex][2]._discr_val);*/
/*  }
  return clase;//c45::Category(c45::Item[ElementIndex], c45::Pruned[0]);*/
}
//---------------------------------------------------------------------------
std::vector<double> C45Tree::UnnormalizedDistribution(int ElementIndex)
{
  vector<double> distrib( ((NomData*)data)->NumClass, 0.0 );

  c45::Description ejm = GetC45DataSets()->GetDescription(data, ElementIndex);  
  c45::Category(ejm, GetTree());

  for(int i = 0; i < c45::MaxClass + 1; i++) {
    distrib[i] = c45::FinalNode->ClassDist[i];
  }

  return distrib;
}
//---------------------------------------------------------------------------
/*namespace c45 {
extern Tree            *Raw;
extern Tree     *Pruned;
}*/
void C45Tree::Build(Data *data, FuncionDeProgreso *fp)
{
  //Si ya estaba creado liberamos memoria
  if (c45tree) c45::ReleaseTree(c45tree);
  if (c45tree_unpruned) c45::ReleaseTree(c45tree_unpruned);
  c45tree = 0;
  c45tree_unpruned = 0;

  Init(data);
  c45::InitialiseTreeData();
//  c45::InitialiseWeights();the weightd are init in SetTrainData
  GetC45DataSets()->SetTrainData(data);
  c45::TRIALS = 1;
  c45tree = c45::FormTree(0, c45::MaxItem);
  c45tree_unpruned = c45::CopyTree(c45tree);
  c45::Prune(c45tree);

  
/*  	c45::PrintTree(c45tree);

  printf("\nNumber of Leaves  : \t\n");
  printf("\nSize of the tree : \t%d\n\n", c45::TreeSize(GetTree()));*/
/*  if (!c45::Raw){printf("INITTREE<---\n");
    c45::Raw = (c45::Tree*)calloc(1, sizeof(c45::Tree));
    c45::Pruned = (c45::Tree*)calloc(1, sizeof(c45::Tree));
}
c45::Raw[0]=c45tree_unpruned;
c45::Pruned[0]=c45tree;
  c45::InitialiseTreeData();
  GetC45DataSets()->SetTrainData(data);
  c45::Evaluate(false, 0);*/
//  printf("ERROR: %f", Error(0, data->GetNTrain()-1));
}
//---------------------------------------------------------------------------
void C45Tree::SetData(Data *data)
{
  //Si ya estA creado el Arbol cargamos todos los datos,
  //si no cargamos sOlo los datos de train
/*  int ntrain;
  if (c45tree) {
    ntrain = data->GetNTrain();
    data->SetNTrain(data->GetNTotal());
  }*/

  Classifier::SetData(data);
  //SOlo se cargan los datos otra vez si han cambiado
  GetC45DataSets()->GetC45Data(data);
  if (!c45tree/* || c45::MaxItem+1!=data->GetNTrain() || lastdataset!=data*/) {
/*    data->SaveAsC45("ABCDEF");
    c45::SetGlobalOpt('f', "ABCDEF");
    c45::GetNames();
    c45::GetData(".data");
    lastdataset = data;*/
  }

  //Recuperamos los NTrain
/*  if (c45tree) 
    data->SetNTrain(ntrain);*/
}
//---------------------------------------------------------------------------
string C45Tree::Info(int value)
{
  SetData(data);
  string res = Classifier::Info(value);
  char kk[9000];

  sprintf(kk, "%d", c45::TreeSize(GetTree()));
  res = res + " -Pruned tree: " + (GetUsePruned() ? "YES" : "NO");
  res = res + " -Tree size: " + kk; 
  if (value>2) {
        //se imprime directamente en vez de guardarlo en la cadena
     if (!GetUsePruned() || value>3) {
        printf("SIN PODAR ------------------------------------");
  	c45::PrintTree(c45tree_unpruned);
     }
     if (GetUsePruned() || value>3) {
        printf("PODADO  --------------------------------------");
  	c45::PrintTree(c45tree);
     }
  }
  return res;
}
//---------------------------------------------------------------------------
void C45Tree::SetGlobalOpt(int opt, char *arg)
{
  c45::SetGlobalOpt(opt, arg); 
}
//---------------------------------------------------------------------------
//Funciones para guardar y leer de fichero
void C45Tree::Guardar(std::ostream &sal, int version)
{
//  Classifier::Guardar(sal, version);
//  tree->Guardar(sal);
}
//---------------------------------------------------------------------------
void C45Tree::Leer(std::istream &in, int version)
{
//  Classifier::Leer(in, version);
//  if (!tree) tree = new Tree(salypimienta);
//  tree->Leer(in);
}
//---------------------------------------------------------------------------
//-----------------------------------------------------------  C45Data  -----
//---------------------------------------------------------------------------
C45Data::C45Data(Data *data) 
{
  Init(data);
  TrainData = 0;
}
C45Data::~C45Data() 
{
  SetData();
  c45::ReleaseData();
  c45::ReleaseNames();
  if (TrainData) delete[]TrainData;
}
//---------------------------------------------------------------------------
void C45Data::Init(Data *data) 
{
  this ->data = data;

  if (!data) {
    MaxAtt = 0;
    MaxClass = 0;
    MaxDiscrVal = 2;
    MaxItem = 0;
    Item = 0;
    MaxAttVal = 0;
    SpecialStatus = 0;
    ClassName = 0;
    AttName = 0;
    AttValName = 0;
  }
  else {
    int prevtrain = data->GetNTrain();
    int ntotal = data->GetTotalData();
    int caux = data->GetNumVar();
    int cpos = data->GetNumVar()+Data::GetIniPosIndex();

    // Ponemos totdos los datos como train para poder cargarlos todos
    data->SetNTrain(ntotal);

    // Marcamos el orden en el que vienen los datos y ordenamos por posicion ini
    for(int i=0;i<ntotal;i++) {
      data->SetValueVar(i, caux, i);
    }
    data->SortOn(cpos, 0, ntotal-1);

    // Guardamos los datos como c45 y los recuperamos
    string dn = GeneratePIDFileName("ABC", ""); 
    data->SaveAsC45((char*)dn.c_str());
    c45::SetGlobalOpt('f', (char*)dn.c_str());
    c45::GetNames();
    c45::GetData((char*)".data");
    GetData();

    // Guardamos donde estA, en items, cada dato de data 
    for(int i=0;i<ntotal;i++) isin[data->GetDatIniPos(i)] = i;

    // Y lo dejamos todo como estaba
    data->SortOn(caux, 0, ntotal-1);
    data->SetNTrain(prevtrain);
  }
}
//---------------------------------------------------------------------------
void C45Data::GetData() 
{
  MaxAtt        = c45::MaxAtt;
  MaxClass      = c45::MaxClass;
  MaxDiscrVal   = c45::MaxDiscrVal;
  MaxItem       = c45::MaxItem;
  Item          = c45::Item;
  MaxAttVal     = c45::MaxAttVal;
  SpecialStatus = c45::SpecialStatus;
  ClassName     = c45::ClassName;
  AttName       = c45::AttName;
  AttValName    = c45::AttValName;
} 

//---------------------------------------------------------------------------
void C45Data::SetData() 
{
  c45::MaxAtt        = MaxAtt;
  c45::MaxClass      = MaxClass;
  c45::MaxDiscrVal   = MaxDiscrVal;
  c45::MaxItem       = MaxItem;
  c45::Item          = Item;
  c45::MaxAttVal     = MaxAttVal;
  c45::SpecialStatus = SpecialStatus;
  c45::ClassName     = ClassName;
  c45::AttName       = AttName;
  c45::AttValName    = AttValName;
}
//---------------------------------------------------------------------------
void C45Data::SetTrainData() 
{
  SetData();
  if (TrainData) delete[]TrainData;
  TrainData = new c45::Description[data->GetNTrain()];
  for(int i=0;i<data->GetNTrain();i++) {
    TrainData[i] = GetDescription(i);
    c45::Weight[i] = data->GetDatWeight(i);
    c45::InitWeight[i] = data->GetDatWeight(i);
  }
  c45::MaxItem = data->GetNTrain()-1;
  c45::Item    = TrainData;
}
//---------------------------------------------------------------------------
c45::Description C45Data::GetDescription(int ElementIndex) 
{
  return Item[isin[data->GetDatIniPos(ElementIndex)]];
}
//---------------------------------------------------------------------------
//-------------------------------------------------------  C45DataSets  -----
//---------------------------------------------------------------------------
C45DataSets *C45DataSets::GetDefault()
{
  if (!def) def = new C45DataSets();
  return def;
}
//---------------------------------------------------------------------------
void C45DataSets::SetTrainData(Data *data) 
{
  GetC45Data(data)->SetTrainData();
}
//---------------------------------------------------------------------------
c45::Description C45DataSets::GetDescription(Data *data, int ElementIndex) 
{
  return GetC45Data(data)->GetDescription(ElementIndex);
}
//---------------------------------------------------------------------------
C45Data *C45DataSets::GetC45Data(Data *data)
{
  C45Data* c45data = datos[data];
  if (!c45data) {
    c45data = new C45Data(data);
    datos[data] = c45data;
    data->AddListener(this);
  }
  c45data->SetData();
  return c45data;
}
//---------------------------------------------------------------------------
void C45DataSets::OnDelete(Data *data)
{
  C45Data* c45data = datos[data];
  if (c45data) {
    delete c45data;
    datos[data] = 0;
  } 
  return;
}
//---------------------------------------------------------------------------
C45DataSets *GetC45DataSets()
{
  return C45DataSets::GetDefault();
}
//---------------------------------------------------------------------------

