//#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "node20.h"
#include "data20.h"
#include<string>
#include<iostream>

/*UNIXextern*/double theSlope;
int Node::max_depth = -100;

void _printf(char*,...);

Node::Node(Node* n, int M, int max)
{
  d = new NodeData;

  info = 0;
  depth = 0;
  d->fuzzindex = 0;
  first = 0;
  last = 0;
  parent = n;
  child = sib = NULL;
  coeff = d->plm = NULL;
  NomSplit = 0;
  NodePop = NULL;
  d->FuzzSplit = NULL;

  SplitType = ERR;
  att = -1;
  iClass = -1;
  fClass = 0.0;
  d->Rleaf = d->RsubTree = d->RleafCV = d->RsubTreeCV = 0.0;
  d->s1 = d->s12 = d->s1leaf = d->s12leaf = 0.0;
  d->fRTest = 0.0;
  K = -1;
  d->normcoeff = 1.0;
  coeff = d->dercoeff = 0;
  sizecoeff = 0;
  if(M)
  {
    sizecoeff = M;
    coeff = 0;//new double[M]; // only if using lin combo of atts
    d->dercoeff = 0;//new double[M];
   // for(int i=0;i<M;i++) {
    //  coeff[i] = 0.0;
     // d->dercoeff[i] = 0.0;
   // }
  }
  // Size stuff
  T = 1;

  if(n)
  {
    if(n->child)
      n->child->sib = this;
    else
      n->child = this;
  }

  d->CrispFlag = true;
  first=0;
  last=max-1;
  d->sig = new sigmoid(1,1);
  d->GeneratedByGroup = -1;
  d->IsLeafAfterPrune = false;
  d->DoesntWantSplitting = false;
  //univariate_split = -1;

  alpha = 0.0; 
}


Node::Node(Node * n, Node *_parent)
{
  *this = *n;

  info  = 0; //Not copied

  parent = _parent;
  if (n->child) {
    this->child = new Node(n->child, this);
    this->child->sib = new Node(n->child->sib, this);
  }

  if (n->coeff) {
    coeff = new double[sizecoeff];
    for(int i=0;i<sizecoeff;i++) coeff[i] = n->coeff[i];
  }
  if (n->d->dercoeff) {
    d->dercoeff = new double[sizecoeff];
    for(int i=0;i<sizecoeff;i++) d->dercoeff[i] = n->d->dercoeff[i];
  }
  d = new NodeData;
  //Mal copiado
  NodePop = d->plm = 0;
  d->storage = 0;
  d->FuzzSplit = 0;
  d->sig = new sigmoid(1,1);
  d->IsLeafAfterPrune = false;
  d->DoesntWantSplitting = false;
  d->GeneratedByGroup = -1;
}

Node::~Node()
{
  FreeAuxData();
  if(coeff) delete[] coeff;
  if(NodePop) delete[] NodePop;
}
void Node::FreeAuxData()
{
  if (d) {
    if(d->dercoeff) delete [] d->dercoeff;
    if(d->plm) delete[] d->plm;
    if(d->FuzzSplit) delete[] d->FuzzSplit;
    delete d->sig;
    delete d;
    d = 0;
  }
  if (coeff && (IsUnivariateSplit() || child == 0)) {
    delete []coeff;
    coeff = 0;
  }
}
Node* Node::nextDown(int K)
{
  Node* cur = this;
  if(cur->sib)
  {
    cur = cur->sib;
    while(cur->child && (K < 0 || cur->child->K < 0)) cur = cur->child;
  } else cur = cur->parent;
  return cur;
}

Node* Node::nextUp()
{
  Node* cur = this;
  if(cur->child) cur = cur->child;
  else {
    while(cur && !cur->sib) cur = cur->parent;
    if(cur) cur = cur->sib;
  }
  return cur;
}
Node *Node::GetRoot() 
{
  Node *aux;

  aux=this;
  while (aux->parent) 
    aux = aux->parent;

  return aux;
}

void Node::ComputeAlpha(int Ndata, Node *root, bool UseSqrt) 
{
  Node *raiz = root ? root : GetRoot();

  double denominador = UseSqrt ? sqrt(raiz->T) - sqrt(raiz->T - T + 1) : T - 1;

  alpha = T>1 ? (d->Rleaf - d->RsubTree)/(Ndata*denominador) : 0.0; 
}

void Node::FixDers(Instance& m, int K,int NVarsOrdSplit)
{
   this->FixDers(m,1.0,K,NVarsOrdSplit);
}

void Node::FixDers(Instance& m, double fac, int K,int NVarsOrdSplit)
{

  if (IsLeaf(K))
  {
    if(iClass == -1)
      d->derL += fac * d->cumMU;   // Only Ordinal dep var
  }
  else if(SplitType == ORD)           // Ordinal + Nominal
  {  
    double value = 0.0;

    if(coeff) {
      for(int j=0;j<NVarsOrdSplit;j++) {
        value += coeff[j] * m[j];
      }
    }
    else {
      value = m[att];
    }

    value = value/d->normcoeff;  

    fac *= d->cumMU*(child->d->tempY - child->sib->d->tempY);
    d->dera += 0.0*fac*(d->sig->dy_da(value));
    d->derb += fac*(d->sig->dy_db(value));
    if(coeff)
    {
      for(int j=0;j<NVarsOrdSplit;j++)
      {
        // if(coeff[j]/normcoeff < 1.0e-6)
        //  dercoeff[j] = 0.0;
        // else
        // {
      if(iClass == -1)
        d->dercoeff[j] -= fac*m[j]*(d->sig->dy_db(value));   // for regression
      else 
        d->dercoeff[j] -= fac*((m[j]-value*coeff[j]/d->normcoeff)/d->normcoeff)*(d->sig->dy_db(value));  // for classification
        
      }
    }
  }
}


int Node::CountParams(int K)
{

  if(IsLeaf(K))
  {
    if(iClass == -1) return 1; // Only Ordinal
    return 0;
  }
  if(SplitType != ORD)
    return 0;
  if(coeff)
    return (2 + sizecoeff);
  return 2;  // Ordinal +Nominal

}

  // Called once to fix initial parameters
void Node::InitParams(double** p,int K, Data* data)
{
  double xmin,xmax;

  data->GetScale(first, last, xmin, xmax, coeff, att);


  if(IsLeaf(K) )
  {
    if(iClass == -1) // Only Ordinal
    {
      if(!d->plm) *(*p) = fClass; //y
      else *(*p) = d->plm[data->GetNordfuzz()];
      (*p) ++;
    }
  }
  else    // Ordinal and Nominal
  {
    double range = fSplit-xmin;
    if(xmax-fSplit < range) range = xmax-fSplit;

  ////////////////////////////////
  //
  //  Normalization for ordinals: Better splits are made stiffer
  //
    if(data->GetDepVarType() == ORD)
    {
      double aux1 = child->d->Rleaf + child->sib->d->Rleaf;
      double aux2; 
      if(aux1 < 1e-6)
        aux1 = 1e-6;
      else if((aux2 =  (sqrt(d->Rleaf/aux1) -1.0)) < 1.e-6)
        range = 1.0e6;
      else range /= aux2; 
    }
  //*/
    //
    // 
    //
    if(iClass != -1)   // For classification problems
    {
      d->normcoeff = 0.0;
      if(coeff){ 
        for(int j=0;j<sizecoeff;j++) {
          d->normcoeff += coeff[j]*coeff[j];
        }
      }
      d->normcoeff = sqrt(d->normcoeff);
    }
    else d->normcoeff = 1.0;
    //
    //
    
    if(SplitType ==ORD)
    {
      d->sig->a = *(*p) = theSlope/(range);  //a =150 Large for Crisp split
      (*p) ++;
      d->sig->b = *(*p) = fSplit/d->normcoeff; //b
      (*p) ++;
      if(coeff)
      {
        for(int j=0;j<sizecoeff;j++) 
        {
          *(*p) = coeff[j];
          (*p)++;
        }
      }
    }
  }
}
  // Called often to translate variables
  // from array to tree
  // and initalise accumulators
void Node::FixParams(double** p,int K, int Nordfuzz)
{
  if(IsLeaf(K))
  {
    if(iClass == -1)  // Only Ordinal dep vars
    {
      d->derL = 0.0;
      if(!d->plm) fClass = *(*p); //y
      else d->plm[Nordfuzz] = *(*p);
      (*p) ++;
    }
  } 
  else if (SplitType == ORD)      // Ordinal and Nominal dep vars
  {
    d->dera = d->derb = 0.0;
  
    d->sig->a = *(*p); // inverse width
    (*p) ++;
    d->sig->b = *(*p); // center of split
    (*p) ++;
  
  
    if(coeff)
    {
      for(int j=0;j<sizecoeff;j++)
      { 
        d->dercoeff[j] = 0.0;
        coeff[j] = *(*p);
        (*p)++;
      }
      //
      //  For classification problems, find normalization of 
      //  the vector of multivariate split coefficients 
      //
      if(iClass != -1)
      {
        d->normcoeff = 0.0;
        for(int j=0;j<sizecoeff;j++)
          d->normcoeff += coeff[j]*coeff[j];
        d->normcoeff = sqrt(d->normcoeff);  
      }
    }

  }
}
  // Called often to translate derivs
  // FROM TREE TO ARRAY
void Node::ColParams(double** p,int K)
{
  if(IsLeaf(K))
  {
    if(iClass == -1)  // Only Ordinal
    {
      *(*p) = d->derL; //y
      (*p) ++;
    }
  } 
  else  if (SplitType ==ORD)  // Ordinal + Nominal
  {
  
    *(*p) = d->dera; //a
    (*p) ++;
    *(*p) = d->derb; //b
    (*p) ++;
    
    if(coeff)
    {
      for(int j=0;j<sizecoeff;j++)
      {
        *(*p) = d->dercoeff[j];  
        (*p)++;
      }
    }
  }
}

double sigmoid::y(double x)
{
  double val = a*(x-b);
  if(val > 0.0)  return exp(-val)/(1.0+exp(-val));
    //return 0.0;
  return (1.0/(1.0+ exp(val)));
}

double sigmoid::dy_da(double x)
{
  return ((x-b)*y(x)*(y(x)-1.0));
  // return 0.0;   // ASG: To remove (Fixed width)
}

double sigmoid::dy_db(double x)
{
  return -a*(y(x)*(y(x)-1.0));
}

double Node::ParamMem(Instance& m, int NVarsOrdSplit)
{
  double val=0.0;

  switch(SplitType)
  {
    case NOM:
      if((1<<(int)m[att]) & NomSplit) // Could be completely wrong!
        val = 1;
      else val = 0;
      break;
    case ORD:
      if(coeff)
        for(int j=0;j<NVarsOrdSplit;j++) 
          val += coeff[j] * m[j]; 
      else val = m[att];
      val = d->sig->y(val/d->normcoeff);
      break;
        
    case FUZZ:
      // TO Be implemented
      break;
    default:
      break;
  };
  return val;
}


void Node::Assign(Node * n)
{
  //Copia los datos de un nodo en otro. NO TODO DE MOMENTO.
  iClass = n->iClass;
//  sizecoeff = n->sizecoeff;
//  coeff = new double[sizecoeff];
  for(int i=0;i<sizecoeff;i++) coeff[i] = n->coeff[i];
  fSplit = n->fSplit;
  if (child->sib) {delete child->sib; child->sib = 0;}
  if (n->child->sib) {
    child->sib = new Node(this, n->child->sib->sizecoeff, 0);
    child->sib->Assign(n->child->sib);
  }
  if (child) {delete child; child = 0;}
  if (n->child) {
    child = new Node(this, n->child->sizecoeff, 0);
    child->Assign(n->child);
  }
}
//Funciones para guardar y leer de fichero
void Node::Guardar(std::ostream &salida, int version)
{
/*  salida.write((char*)first, sizeof(first));
  salida.write((char*)last, sizeof(last));
  salida.write((char*)firstTest,
  salida.write((char*)lastTest
  salida.write((char*)iClass
  salida.write((char*)fClass
  salida.write((char*)min_in_child
  salida.write((char*)att
  salida.write((char*)normcoeff
  salida.write((char*)sizecoeff
  for(int i=0;i<sizecoeff;i++)
    salida.write((char*)coeff[i];
  salida.write((char*)(void*)plm
  salida.write((char*)fSplit
  salida.write((char*)NomSplit
  salida.write((char*)(void*)FuzzSplit
  salida.write((char*)SplitType
  salida.write((char*)depth
  salida.write((char*)RPrune
  salida.write((char*)RPruneSubTree
  salida.write((char*)Rleaf
  salida.write((char*)RsubTree
  salida.write((char*)RleafCV
  salida.write((char*)RsubTreeCV
  salida.write((char*)alpha
  salida.write((char*)s1leaf
  salida.write((char*)s12leaf
  salida.write((char*)s1
  salida.write((char*)s12
  salida.write((char*)fRTest
  salida.write((char*)fuzzindex
  salida.write((char*)CrispFlag
  salida.write((char*)T
  salida.write((char*)K
  salida.write((char*)sig->a
  salida.write((char*)sig->b
  salida.write((char*)tempY
  salida.write((char*)tempY2
  salida.write((char*)tempMU
  salida.write((char*)tempMU2
  salida.write((char*)cumMU
  salida.write((char*)derL
  salida.write((char*)dera
  salida.write((char*)derb
  salida.write((char*)dercoeff;
  if (dercoeff) for(int i=0;i<sizecoeff;i++) salida.write((char*)dercoeff[i];

  //Paso de estos: double* NodePop y double** storage;

  salida << MembershipTotal << entropy << norm;
  salida << IsLeafAfterPrune << DoesntWantSplitting;
  salida << (void*)child << (void*)sib;*/
  if (version==0) {
    salida.write((char*)this, sizeof(Node));

    if (this->d) salida.write((char*)this->d, sizeof(NodeData));
    else         salida.write("\0\0\0\0", sizeof(int));

    salida.write((char*)d->sig, sizeof(d->sig));

    for(int i=0;i<sizecoeff;i++)
      salida.write((char*)&coeff[i], sizeof(coeff[i]));

    if (d->dercoeff) {
      for(int i=0;i<sizecoeff;i++) {
        salida.write((char*)&d->dercoeff[i], sizeof(d->dercoeff[i]));
      }
    }
  }
  else if(version==1) {
    //Clase
    salida.write((char*)&iClass,    sizeof(int));
    //Split Info
    salida.write((char*)&att,       sizeof(int));
    VarType SaveSplitType = SplitType;
    if (SplitType==ORD && !IsUnivariateSplit()) SaveSplitType = ORDM;
    salida.write((char*)&SaveSplitType, 1);//sizeof(int));
    if (SplitType==NOM) {
      salida.write((char*)&NomSplit,  sizeof(int));
    }
    else {
      salida.write((char*)&fSplit,    sizeof(double));
      if (SaveSplitType==ORDM) {
        for(int i=0;i<sizecoeff;i++)
          salida.write((char*)&coeff[i],  sizeof(double));
      }
    }

    //K
    salida.write((char*)&K ,        2);//sizeof(int));
    

    //Familia
    int des = 0;
    if (child) des |= 0x1;
    if (sib) des |= 0x2;
    salida.write((char*)&des,    1);//sizeof(int));
  }

  if (child) child->Guardar(salida, version);
  if (sib)   sib->Guardar(salida, version);
}
void Node::Leer(std::istream &in, int version)
{
/*  in >> first >> last >> firstTest >> lastTest >> iClass;
  in >> fClass >> min_in_child >> att;
  in >> normcoeff >> sizecoeff;
  coeff = (sizecoeff>0) ? new double[sizecoeff] : 0;
  for(int i=0;i<sizecoeff;i++)
    in >> coeff[i];
  in >> (void*)plm >> fSplit >> NomSplit >> (void*)FuzzSplit;
  in >> ((int)SplitType) >> depth >> RPrune >> RPruneSubTree;
  in >> Rleaf >> RsubTree >> RleafCV >> RsubTreeCV >> alpha;
  in >> s1leaf >> s12leaf >> s1 >> s12 >> fRTest >> fuzzindex;
  sig = new sigmoid(1,1);
  in >> CrispFlag >> T >> K >> sig->a >> sig->b;
  in >> tempY >> tempY2 >> tempMU >> tempMU2 >> cumMU;
  in >> derL >> dera >> derb >> (void*)dercoeff;
  dercoeff = (dercoeff) ? new double[sizecoeff] : 0;
  if (dercoeff) for(int i=0;i<sizecoeff;i++) in >> dercoeff[i];

  //Paso de estos: double* NodePop y double** storage;

  in >> MembershipTotal >> entropy >> norm;
  in >> IsLeafAfterPrune >> DoesntWantSplitting;
  in >> (void*)child >> (void*)sib;*/

  if (version==0) { 
//ESTO HA DEJADO DE FUNCIONAR AL SACAR  LOS DATOS DEL NODO FUERA
    if(coeff) delete[] coeff;
    if(d->dercoeff) delete [] d->dercoeff;
//ESTO HA DEJADO DE FUNCIONAR AL SACAR  LOS DATOS DEL NODO FUERA
    delete d->sig;
    if (d) delete d;
    d = new NodeData;
//ESTO HA DEJADO DE FUNCIONAR AL SACAR  LOS DATOS DEL NODO FUERA
    //Leyendo objeto
    in.read((char*)&first, sizeof(int));
    in.read((char*)&last, sizeof(int));
//ESTO HA DEJADO DE FUNCIONAR AL SACAR  LOS DATOS DEL NODO FUERA
    in.read((char*)&d->firstTest, sizeof(int));
    in.read((char*)&d->lastTest, sizeof(int));
//ESTO HA DEJADO DE FUNCIONAR AL SACAR  LOS DATOS DEL NODO FUERA
    in.read((char*)&iClass, sizeof(int));
    in.read((char*)&fClass, sizeof(double));

//ESTO HA DEJADO DE FUNCIONAR AL SACAR  LOS DATOS DEL NODO FUERA
    in.read(((char*)this)+sizeof(int)*5+sizeof(double), 
                   sizeof(Node)-sizeof(int)*5-sizeof(double));

////  in >> fClass >> min_in_child >> att;
////  printf("%g:%d\n", norm, GeneratedByGroup);

    //Reconstruir punteros
    coeff = (sizecoeff>0) ? new double[sizecoeff] : 0;
    d->dercoeff = (d->dercoeff) ? new double[sizecoeff] : 0;
    d->plm = 0;
    d->FuzzSplit = 0;
    NodePop = 0;
    d->sig = new sigmoid(1, 1);
    in.read((char*)d->sig, sizeof(d->sig));
    for(int i=0;i<sizecoeff;i++)
      in.read((char*)&coeff[i], sizeof(coeff[i]));
    if (d->dercoeff)
      for(int i=0;i<sizecoeff;i++)
        in.read((char*)&d->dercoeff[i], sizeof(d->dercoeff[i]));;

    //Leyendo nodos hijo y hermano
    if (child) {
      child = new Node(0, 0, 0);
      child->Leer(in, version);
      child->parent = this;
      child->sib->parent = this;
    }
    if (sib) {
      sib = new Node(0, 0, 0);
      sib->Leer(in, version);
    }
//ESTO HA DEJADO DE FUNCIONAR AL SACAR  LOS DATOS DEL NODO FUERA
  }
  else if (version==1) {
//ESTO SIGUE FUNCIONANDO
    in.read((char*)&iClass,    sizeof(int));
    in.read((char*)&att,       sizeof(int));
    SplitType = (VarType)0;
    in.read((char*)&SplitType, 1);//sizeof(int));
    if (SplitType==NOM) {
      in.read((char*)&NomSplit,  sizeof(int));
    }
    else {
      in.read((char*)&fSplit,    sizeof(double));
      if (SplitType==ORDM) {//Multivariate split
        SplitType=ORD;
        for(int i=0;i<sizecoeff;i++)
          in.read((char*)&coeff[i],  sizeof(double));
      }
      else if (coeff) { //Univariate split
        for(int i=0;i<sizecoeff;i++)
          coeff[i] = i==att ? 1.0 : 0.0;
        //univariate_split = 1;
      }
    }
    K = 0;
    in.read((char*)&K,         2);//sizeof(int));

    //Reconstruimos algunas vars en pre orden
    if (parent) depth = parent->depth + 1;

    //Leyendo nodos hijo y hermano
    int des = 0;
    in.read((char*)&des,    1);//sizeof(int));
    if (des & 0x1) { //child pointer
      child = new Node(this, sizecoeff, 0);
      child->Leer(in, version);
    }
    if (des & 0x2) { //sib pointer
      sib = new Node(parent, sizecoeff, 0);
      sib->Leer(in, version);
    }

    //Reconstruimos algunas vars en post orden
    T = 1;
    if (child) T += child->T + child->sib->T;
  }
}
bool Node::IsUnivariateSplit()
{
//  if (univariate_split!=-1) return univariate_split;
  if (SplitType != ORD) return false;
  if (!coeff) return true;
  for(int i=0;i<sizecoeff;i++)
    if (i!=att && coeff[i]!=0.0) {
      //univariate_split = 0;
      return false;
    }

//  univariate_split = coeff[att] == 1.0 ? 1 : 0;
  return coeff[att]==1.0;
}
//---------------------------------------------------------------------------
#ifdef _DEBUG_NODE_CPP
int main(int argc, char* argv[])
{
  printf("%d\n", sizeof(Node));
}
#endif


