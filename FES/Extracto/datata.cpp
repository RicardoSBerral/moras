//------------------------------------------------------------------


#include "data20.h"
#include <time.h>
#include <values.h>
#include <sstream>
//------------------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "util20.h"
#include "node20.h"

#include<string>
#include<iostream>
#include<istream>
#include<fstream>
#include<limits>
using namespace std;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/*
NTotalOrdfuzz = Nordfuzz+1
The indexing scheme:
The ordinal independent varibles are 0->(Nord-1) in m
The fuzzy independent variables are Nord->(Nordfuzz-1) in m
The nominal independent variables are Nordfuzz->(Nordfuzz+Nnom-1) in m
The dependent variable  is (Nvar-1) in m
Nvar = Nordfuzz+Nnom+1

The independent variables are 0-(Nordfuzz) in St
  with one of these (Nordfuzz) the constant.
The dep       "         is (NtotalOrdfuzz) in St

Index[b][i] = index in m of var i in St if Nordfuzz>i>=0
Index[b][i] < 0 if i refers to constant
Index[b][i] = Nordfuzz refers to dependent variable
*/

// The following macro gives index in m
//    of var j on branch b
#define VAR(b,j) ((Index[b][j]<0 ? Nordfuzz : Index[b][j]))
// The following macro gives the actual value of
// the variable, in the n-th case, of var j on
// branch b.
#define VAL(n,b,j) (Index[b][j]<0 ? 1.0 : m[n][(Index[b][j] == Nordfuzz) ? Nvar-1 : Index[b][j]])

// 0->Nvar-2: Independent variables
// Nvar-1: Dependent variable
// Nvar,Nvar+1: Multiple split working space
// Nvar+2: CV tree label
// Nvar+3: Membership
// Nvar+4: Weight
// Nvar+5: Initial position
// Nvar+6: Group number for Iterative growing pruning method

int Data::MultipleSplitWs1 = 0;
int Data::MultipleSplitWs2 = 1;
int Data::CvTreeLabelIndex = 2;
int Data::MembershipIndex  = 3;
int Data::WeightIndex      = 4;
int Data::IniPosIndex      = 5;
int Data::GroupIndex       = 6;
int Data::Proportional     = true;

//---------------------------------------------------------------------------
int rand2()
{
/*  int r;
  static char nf[100] = "d:\\temp.rnd";*/
/*  FILE *f = fopen(nf, "a");
  r=rand();
  fwrite(&r,1,sizeof(int),f);*/
/*  static int pos=0;
  FILE *f = fopen(nf, "r");
  fseek(f, sizeof(int)*pos, 0);
  fread(&r,1,sizeof(int),f);
  fclose(f);
  pos++;
  return TUtil::rand();*/
  return -123;
}
double drand(unsigned int seed =0)
{
  if(seed) srand(seed);
  double dd = TDebugRand::Rand();
  return dd/RAND_MAX;
}


void salypimienta(char *texto)
{
  //printf(texto);
  return;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//-----------------------------------------------------------  NomData ------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


NomData::NomData(std::istream &in, void (*_console)(char *))
  : Data(in,  _console)
{
  init();
}
NomData::NomData(char *filename)
  : Data(filename)
{
  init();
}
NomData::NomData(int NTotalData, int NVar)
  : Data(NTotalData, NVar)
{
  init();
}
NomData::NomData(int NTotalData, int NVar, double**data)
  : Data(NTotalData, NVar, data)
{
  init();
}

NomData::~NomData()
{
  delete[] Pop;
  delete[] PopR;
  delete[] PopL;
}

void NomData::init()
{
  atdrawchosefathers = false;
  NumClass = NomTerms[Nnom].size();

//  printf("Class variables to %d\n", NumClass);
  Pop = new double[NumClass];
  PopR = new double[NumClass];
  PopL = new double[NumClass];
  ponder = false;
  multree = false;

  if (NumClass) {
    Ponderate(Ndata);
    iClass = WhichClass(Pop);
  }

  degree = 1;
}


void NomData::Ponderate(int Num)
{
  for(int i = 0; i < NTotal ; ++i)
    m[i][Nvar+3] = /*m[i][Nvar+4] =*/1.0;//???????????????


  FindPops(Pop, 0,Num-1);

  double total = 0.0;
  if(ponder)
  {
    double* weight = new double[NumClass];
    int i;
    for(i =0; i< NumClass; ++i)
      weight[i] = Num/(NumClass*Pop[i]);

    for(i = 0; i < NTotal ; ++i)
    {
      m[i][Nvar+4] = weight[(int)m[i][Nvar-1]];
      if(i<Num)
        total += m[i][Nvar+4];
    }
    if(fabs(total-Num) > 1e-6)
      throw new exception();//console->printf("Warning: Weights are not normalized");
    delete[] weight;
  }
}
void NomData::EqualizeClassWeight()
{
  double* class_weights = new double[NumClass];

  for(int i = 0; i < NTotal ; ++i)
    m[i][Nvar+3] = 1.0;


  FindPops(class_weights, 0, NTotal-1);

  for(int i = 0; i< NumClass; ++i)
    class_weights[i] = NTotal/(NumClass*class_weights[i]);

  for(int i = 0; i < NTotal ; ++i) {
    m[i][Nvar+4] = class_weights[ (int)m[i][Nvar-1] ];
  }

  delete[] class_weights;

}

void NomData::FindPops(double* HoldPop, int first, int last)
{
  

  for (int j=0; j < NumClass; j++)
    HoldPop[j]=0.0;

  for(int i=first; i<= last; i++)
    HoldPop[(int)m[i][Nvar-1]] += m[i][Nvar+3];


}


int NomData::WhichClass(double* pt)
{

  double max=0.0;
  int maxj = -1;
  for(int j=0; j<NumClass; j++,pt++)
  {
    if(j == 0 || *pt>max)
    {
      maxj = j;
      max = *pt;
    }
  }


  return maxj;
}



double NomData::Impurity(int first, int last,double tot) 
{

  double aux = 0.0;
  double TotPop = 0.0;

  FindPops(Pop,first,last);

  // Gini criterion

  for(int j=0; j< NumClass; j++)
  {
    TotPop += Pop[j];
    aux -=  Pop[j]*Pop[j];
  }
  return (TotPop + aux/TotPop);

/*/
  // Ordinal-like distance

  for(int j=0; j< NumClass; j++)
    TotPop += Pop[j];

  for(int i= first; i<=last; ++i)
  {
    double aux1 = 0.0;
    for(int j=0; j< NumClass; ++j)
    {
      double aux2 = ((i==(int)m[i][Nvar-1]) - Pop[i]/TotPop);
      aux1 += aux2*aux2;

    }
    aux += aux1*m[i][Nvar+3];
  }
  return aux;

//  

  // Entropy (gain criterion)

    for(int j=0; j< NumClass; j++)
    TotPop += Pop[j];

  for(j=0; j< NumClass; j++)
  {
    
    if(Pop[j] > 1e-16)
      aux -=  Pop[j]*log(Pop[j]/TotPop);
  }
  return aux;
*//*

  //      Resubstitution criterion

  for(int j=0; j< NumClass; j++)
  {
    TotPop += Pop[j];

    if(j ==0 || Pop[j] > aux)
      aux = Pop[j];
  }

  return TotPop - aux;
*/


}

double NomData::FindImpurity(double* HoldPop, int first, int last)
{
  double aux = 0.0;
  double TotPop = 0.0;

//
  // Gini criterion multiplied by node membership

  for(int j=0; j< NumClass; j++)
  {

     aux -=  HoldPop[j]*HoldPop[j];
     TotPop += HoldPop[j];
  }

  if (TotPop==0) return 1E308;
  return (TotPop + aux/TotPop);// = n*i(t) = n*(1 - sum_j(p(j|t)^2)
/*/

  // Ordinal-like distance

  for(int j=0; j< NumClass; j++)
    TotPop += Pop[j];

  for(int i= first; i<=last; ++i)
  {
    double aux1 = 0.0;
    for(int j=0; j< NumClass; ++j)
    {
      double aux2 = ((i==(int)m[i][Nvar-1]) - Pop[i]/TotPop);
      aux1 += aux2*aux2;

    }
    aux += aux1*m[i][Nvar+3];
  }
  return aux;

  //
  
  // Entropy (Gain) criterion

  for(int j=0; j< NumClass; j++)
    TotPop += HoldPop[j];
  

  for(j=0; j< NumClass; j++)
  {
    if(HoldPop[j] > 1e-16)
      aux -=  HoldPop[j]*log(HoldPop[j]/TotPop);
  }
  return  aux;

*//*
  //      Resubstitution criterion

  for(int j=0; j< NumClass; j++)
  {
    TotPop += HoldPop[j];
    if(j ==0 || HoldPop[j] > aux)
      aux = HoldPop[j];
  }

  return TotPop - aux/TotPop;
*/   
  
}


double NomData::Initialise(int first, int last, int curSplit, double gamma, int iX, bool FuzzFlag)
{

  bool flagleft = true;
  if(FuzzFlag)
  {
    if(curSplit == last)
      nL = R[0] = 0.0;
    else
    {  
      flagleft = false; 
      nR = R[1] = 0.0;
    } 
  }
  else R[0] = R[1] = nL = nR = 0.0;                 // Initialise
  
  for (int j=0; j < NumClass; j++)
    PopL[j] = PopR[j] = 0.0;

  for(int i=first;i<=last;i++)
  {
    if(m[i][Nvar+3] > 0.0)
    {
      if( (gamma > 1 && i <= curSplit) // normal
        || (gamma < 1 && m[i][iX]+gamma < 0))
      {
        PopL[(int)m[i][Nvar-1]] += m[i][Nvar+3];
        nL += m[i][Nvar+3];
      } else {
        PopR[(int)m[i][Nvar-1]] += m[i][Nvar+3];
        nR += m[i][Nvar+3];
      }
    }
  }
  

  if(nL >0.0)
  {
    if(!FuzzFlag   || (FuzzFlag && flagleft))
      R[0] = FindImpurity(PopL,first,last);
  }
  if(nR > 0.0)
  {
    if(!FuzzFlag  || (FuzzFlag && !flagleft))
      R[1] = FindImpurity(PopR,first,last);
  }

  if(!FuzzFlag)
    return (nR < nL ? nR : nL);
  
  if(flagleft) return nL;
  
  return nR; 
}

double NomData::MovePoint(int curSplit, double gamma, int iX,int first,int last)
{

  if(m[curSplit][Nvar+3] > 0.0)
  {
    if(gamma > 1 || (gamma<1 & m[curSplit][iX]+gamma <0))
    {
      PopL[(int)m[curSplit][Nvar-1]] -= m[curSplit][Nvar+3];
      PopR[(int)m[curSplit][Nvar-1]] += m[curSplit][Nvar+3];

      nL -= m[curSplit][Nvar+3];
      nR += m[curSplit][Nvar+3];

    } else {
      PopL[(int)m[curSplit][Nvar-1]] += m[curSplit][Nvar+3];
      PopR[(int)m[curSplit][Nvar-1]] -= m[curSplit][Nvar+3];  

      nL += m[curSplit][Nvar+3];
      nR -= m[curSplit][Nvar+3];

    } 

    
    R[0] = R[1] = 0.0;
  
    if(nL > 0.0) R[0] = FindImpurity(PopL,first,last);
    if(nR > 0.0) R[1] = FindImpurity(PopR,first,last);
  }

  return (nR < nL ? nR : nL);
}

void NomData::FixLeaf(Node* n)
{

  if(NumClass) {
    if (n->NodePop) delete []n->NodePop;
    n->NodePop = new double[NumClass];
  }
  n->d->Rleaf = 0.0;
  n->iClass = 0;

  FindPops(n->NodePop, n->first, n->last);
  n->iClass = WhichClass(n->NodePop);

  if (atdrawchosefathers && n->parent && n->iClass != n->parent->iClass &&
      n->NodePop[n->iClass]>0 && n->NodePop[n->iClass]==n->NodePop[n->parent->iClass]) {
    n->iClass = n->parent->iClass;
/*    char puf[128];
    sprintf(puf, "JAJAJAJAJALLL%f", (float)n->parent->fSplit);
    if (console) (*console)(puf);
  */}

  for (int j=0; j<NumClass; ++j)
  {
    if(j != n->iClass)
      n->d->Rleaf += n->NodePop[j];
  }

  /*// ASG: Gain-ratio//////////
  n->entropy = FindImpurity(n->NodePop);
  // Gain-ratio //////////////*/


  if(NPrune)
  {
    n->d->RPrune = 0;
    n->d->firstTest = Ndata; n->d->lastTest = Ndata+ NPrune - 1;
    FindMembershipTest(n);
    for(int i = n->d->firstTest; i <= n->d->lastTest; ++i)
    {
      if((int)m[i][Nvar-1] != n->iClass)
        n->d->RPrune += m[i][Nvar+3];
    }
  }
}

void NomData::FixCVError(Node* root)
{
  Node* cur = root;  // To Check ASG

  cur->d->RleafCV = 0.0;
  cur->d->firstTest = Ndata; cur->d->lastTest = NGrow-1;

  while(cur)
  {
    
    if(cur->d->CrispFlag)
    {
    
      for(int i=cur->d->firstTest;i<=cur->d->lastTest;i++)
      {  

        if((int)m[i][Nvar-1] != cur->iClass)
          cur->d->RleafCV += m[i][Nvar+4];  

        if(cur->child && cur->SplitType == ORD)
        {
          m[i][Nvar] = -cur->fSplit;
          if(cur->coeff)
          {
            for(int j=0;j<NVarsOrdSplit;j++) 
              m[i][Nvar] += cur->coeff[j]*m[i][j];
          } else 
            m[i][Nvar] += m[i][cur->att];
          
        }
      }

      if(cur->child)
      {
      
        if(cur->d->firstTest <= cur->d->lastTest)
        {
          if(cur->SplitType==NOM)
          {
            cur->child->d->firstTest =  cur->d->firstTest;
            cur->child->d->lastTest =  cur->d->firstTest - 1 + 
              SortNomOn(cur->NomSplit, cur->att,cur->d->firstTest,cur->d->lastTest);
          
          } else {
            SortOn(Nvar,cur->d->firstTest,cur->d->lastTest);
            cur->child->d->firstTest = cur->child->d->lastTest = cur->d->firstTest;
            cur->child->d->lastTest--;
            while(cur->child->d->lastTest < cur->d->lastTest && m[cur->child->d->lastTest+1][Nvar] < 0.0)
              cur->child->d->lastTest++;
          }
        } else {
          cur->child->d->firstTest = cur->d->firstTest;
          cur->child->d->lastTest = cur->d->lastTest;
        }

        cur->child->sib->d->firstTest = cur->child->d->lastTest+1;
        cur->child->sib->d->lastTest = cur->d->lastTest;
      }
    }
    else
    {
      FindMembershipTest(cur);
      for(int i=cur->d->firstTest;i<=cur->d->lastTest;i++)
      {  
        if((int)m[i][Nvar-1] != cur->iClass && m[i][Nvar+3] > 0.0)
          cur->d->RleafCV += m[i][Nvar+3];  
      }
      if(cur->child)
      {
        cur->child->d->firstTest = cur->child->sib->d->firstTest = cur->d->firstTest;
        cur->child->d->lastTest = cur->child->sib->d->lastTest = cur->d->lastTest;
      }
    }

    cur = cur->nextUp();
  }
    // Now fix the sub-tree RCV
  cur = root; while(cur->child) cur = cur->child;
  while(cur)
  {
    cur->d->RsubTreeCV = 0.0;
    cur = cur->nextDown();
  }
  cur = root; while(cur->child) cur = cur->child;
  while(cur)
  {
    if(cur->child)
    {
      cur->d->RsubTreeCV = cur->child->d->RsubTreeCV+cur->child->sib->d->RsubTreeCV;
    } else {
      cur->d->RsubTreeCV = cur->d->RleafCV;
    }
    cur = cur->nextDown();
  }
}

void NomData::InitialiseErrorCV()
{
  rcv = 0.0;
}

void NomData::AccumulateErrorCV(Node* root)
{
  rcv += root->d->RsubTreeCV;
}

void NomData::GetErrorCV(double* RCV, double* SE)
{
  rcv /= NGrow;

  *SE = sqrt(rcv*fabs(1.0-rcv)/NGrow);
  *RCV = rcv;
}



void NomData::WriteNode(Node* n, char *buf, int K)
{
//  if(!n->child || (n->child->K <= K && K>=0)){
  if(n->IsLeaf(K)){
    if (n->parent&& n->parent->child == n)
      sprintf(buf,"%s is %s (K=%d) (%g(%g)/%d) (a=%g) T=%d\n",
                         VarNames[Nvar-1].c_str(),
                         NomTerms[Nnom].at(n->iClass).c_str(), 
                         n->K, n->d ? n->d->Rleaf : 0, n->d ? n->d->RsubTree : 0, (n->last - n->first + 1),
                         n->alpha, n->T);
    else
      sprintf(buf,"ELSE var %s is %s (K=%d) (%g(%g)/%d) (a=%g) T=%d\n",
                         VarNames[Nvar-1].c_str(),
                         NomTerms[Nnom].at(n->iClass).c_str(), 
                         n->K, n->d ? n->d->Rleaf : 0, n->d ? n->d->RsubTree : 0, (n->last - n->first + 1),
                         n->alpha, n->T);
  } 
  else if (n->SplitType==NOM)  {
    if (n->parent) {
      if(n->parent->child == n) 
        sprintf(buf,"IF the value of (%s) is in ",VarNames[n->att].c_str());
      else 
        sprintf(buf,"ELSE IF the value of (%s) is in ",
                                                     VarNames[n->att].c_str());
    }
    else
      sprintf(buf,"IF the value of %s is in ",VarNames[n->att].c_str());

    int flag =0;
    for(unsigned i=0;i<NomTerms[n->att-Nordfuzz].size(); ++i)  {
      if(1<<i & n->NomSplit) {
        if(flag==0) {
          sprintf(buf+strlen(buf),"{ %s",
                                      NomTerms[n->att-Nordfuzz].at(i).c_str());
          flag=1;
        }
        else 
          sprintf(buf+strlen(buf),", %s", 
                                      NomTerms[n->att-Nordfuzz].at(i).c_str());
      }
    }

    sprintf(buf+strlen(buf)," } (K=%d) (%g(%g)/%d) (a=%g) T=%d\n", 
        n->K, n->d ? n->d->Rleaf : 0, n->d ? n->d->RsubTree : 0, (n->last - n->first + 1), n->alpha, n->T);
  }
  else if (n->SplitType == ORD) {
    if(n->parent)  {
      if(n->parent->child == n) sprintf(buf,"IF (%g) ",-(n->fSplit));
      else sprintf(buf,"ELSE IF (%g) ",-(n->fSplit));
    }
    else 
      sprintf(buf,"IF (%g) ",-(n->fSplit));
    for(int i=0;i<NVarsOrdSplit;i++) {
      if(fabs(n->coeff[i]) > 0.0) {
        sprintf(buf+strlen(buf)," + (%g) %s",n->coeff[i],VarNames[i].c_str());
      }
    }
    sprintf(buf+strlen(buf)," < 0 (K=%d) (%g(%g)/%d) (a=%g) T=%d\n",
        n->K, n->d ? n->d->Rleaf : 0, n->d ? n->d->RsubTree : 0, (n->last - n->first + 1), n->alpha, n->T);
  }
  else if(n->SplitType == FUZZ)  {
    if(n->d && n->d->fuzzindex)  {
      if(n->parent)  {
        if(n->parent->child == n)
          sprintf(buf,"IF fuzzification(%d) of (%g) ", n->d->fuzzindex,
                                                                 -(n->fSplit));
        else
          sprintf(buf,"ELSE IF fuzzification(%d) of (%g) ", n->d->fuzzindex,
                                                                 -(n->fSplit));
      }
    else sprintf(buf,"IF fuzzification(%d) of (%g) ",n->d->fuzzindex,-(n->fSplit));
    for(int i=0;i<NVarsOrdSplit;i++) {
          if(fabs(n->coeff[i]) > 0.0) {
              sprintf(buf+strlen(buf)," + (%g) %s", n->coeff[i],
                                                          VarNames[i].c_str());
          }
    }
    sprintf(buf+strlen(buf)," < 0 (K=%d) (%g(%g)/%d) (a=%g) T=%d\n",
        n->K, n->d->Rleaf, n->d->RsubTree, (n->last - n->first + 1), n->alpha, n->T);
    }
    else  {
      if(n->parent)  {
        if(n->parent->child == n)
          sprintf(buf,"IF the value of (%s) is %s", VarNames[n->att].c_str(), 
                                              n->d->FuzzSplit->GetName().c_str());
        else
          sprintf(buf,"ELSE IF the value of (%s) is  %s", 
                    VarNames[n->att].c_str(), n->d->FuzzSplit->GetName().c_str());
      }
      else
        sprintf(buf,"IF the value of %s is %s", VarNames[n->att].c_str(), 
                                              n->d->FuzzSplit->GetName().c_str());
    }
  }
}


double NomData::Score(int begin, int end, Node* root, int K)
{
  double ress = 0.0;
  if(begin<0)
  {
    begin = NTotal-NTest;
    end = NTotal-1;
  }
  double Npoints = end-begin+1; 
  for(int i=begin;i<= end;i++)
  {
    double unity =0.0;
    Node* cur = root;
    int y= (int)m[i][Nvar-1];
    int ypred; 
    
    while(cur)
    {
     //  Node* temp = cur;

      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {
        cur->d->firstTest = cur->d->lastTest = i;
        FindMembershipTest(cur);

        unity += m[i][Nvar+3];

        ypred = cur->iClass;
        if(ypred != y)
          ress += m[i][Nvar+3];
      }

      if(!curIsLeaf) cur = cur->child;
      else {
        while(cur && !cur->sib) cur = cur->parent;
        if(cur) cur = cur->sib;
      }


    }

    if(fabs(unity-1.0) > 1e-6)
    {
;//1      console->printf("Warning: error in memberships");
    }

  }

  if (Npoints==0) return 1.0;//km
  ress /= Npoints;
  return ress;
}




/////////////////////////////////////////////////////
//                          //
//   Optimization based on IMPURITY          //
//                          //
//////////////////////////////////////////////////////

double NomData::Derivs(Node* root, int K,int begin,int end)
{
  double ss = 0.0;

  int i;
  for(i=begin;i<=end;i++)
  {

    // Going up : fix memberships && partial predictions
    Node* cur = root;
    root->d->cumMU =1.0;
    root->d->tempMU = 1.0;

    while(cur)
    {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {
        double aux = cur->d->tempMU*m[i][Nvar+4];
        if(cur->parent) aux *= cur->parent->d->cumMU;
        int Class = (int)m[i][Nvar-1];
        if(i==begin)
        {
          for(int j=0; j<NumClass; ++j)
            cur->NodePop[j]  = 0.0;
          cur->d->MembershipTotal =0.0;
        }

        cur->NodePop[Class] += aux;
        cur->d->MembershipTotal += aux;

        if(i== end)
        {
          //  Calculate class of leaf node and impurity

          double aux =0.0;
          for(int j=0; j< NumClass; ++j)
            aux += cur->NodePop[j]*cur->NodePop[j];
          cur->iClass = WhichClass(cur->NodePop);
          if(cur->d->MembershipTotal > 1.0e-10)
            aux /= cur->d->MembershipTotal;
          ss -= aux;
        }

      } else {

        double mu =  cur->ParamMem(m[i],NVarsOrdSplit);
        cur->child->d->tempMU = mu;
        cur->child->d->cumMU = cur->d->cumMU*mu;
        cur->child->sib->d->tempMU = (1-mu);
        cur->child->sib->d->cumMU = cur->d->cumMU*(1-mu);
      }

      if(!curIsLeaf) cur = cur->child;
      else {
        while(cur && !cur->sib) cur = cur->parent;
        if(cur) cur = cur->sib;
      }
    }
  }

  for(i=begin;i<=end;i++)
  {
    // Going up : fix memberships && partial predictions
    Node* cur = root;
    root->d->cumMU =1.0;
    root->d->tempMU = 1.0;
    int Class = (int)m[i][Nvar-1];
    while(cur)
    {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {
        if(cur->d->MembershipTotal > 1.0e-10)
          cur->d->tempY  = -2.0*cur->NodePop[Class]/cur->d->MembershipTotal;
        double aux =0.0; 
        for(int j=0; j<NumClass; ++j)
          aux += cur->NodePop[j] * cur->NodePop[j];
        if(cur->d->MembershipTotal > 1.0e-10)
          cur->d->tempY +=  aux/(cur->d->MembershipTotal * cur->d->MembershipTotal);
      } else {

        double mu = cur->ParamMem(m[i],NVarsOrdSplit);
        cur->child->d->tempMU = mu;
        cur->child->d->cumMU = cur->d->cumMU*mu;
        cur->child->sib->d->tempMU = (1-mu);
        cur->child->sib->d->cumMU = cur->d->cumMU*(1-mu);
      }

      if(!curIsLeaf) cur = cur->child;
      else {
        while(cur && !cur->sib) cur = cur->parent;
        if(cur) cur = cur->sib;        
      }  
    }
      
    // Going down : accumulate predictions && derivatives
    cur = root; while(!(cur->IsLeaf(K))) cur = cur->child;
        
    while(cur)
    {
      cur->FixDers(m[i],K,NVarsOrdSplit);
      bool curIsLeaf = cur->IsLeaf(K);
      if(!curIsLeaf)
      {
        Node* temp = cur->child;
        cur->d->tempY = temp->d->tempY*temp->d->tempMU;
        temp = temp->sib;
        cur->d->tempY += temp->d->tempY*temp->d->tempMU;
      } 
      

      if(cur->sib)
      {
        cur = cur->sib;
        while(!(cur->IsLeaf(K))) cur = cur->child;
      } else cur = cur->parent;
    }
  }

  return ss;
}



double NomData::Error(Node* root, int K, int begin,int end)
{
  double ss = 0.0;
  double total =0.0;

  for(int i=begin;i<=end;i++) {
    // Going up : fix memberships && partial predictions
    Node* cur = root;
    root->d->cumMU =1.0;
    root->d->tempMU = 1.0;
    int Class = (int)m[i][Nvar-1];
    while(cur)  {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf) {
        if(Class != cur->iClass) {
          double aux = m[i][Nvar+4]*cur->d->tempMU;
          if(cur->parent) aux *= cur->parent->d->cumMU;
          ss += aux;
        }
      }
      else {
        double mu =  cur->ParamMem(m[i],NVarsOrdSplit);
        cur->child->d->tempMU = mu;
        cur->child->d->cumMU = cur->d->cumMU*mu;
        cur->child->sib->d->tempMU = (1-mu);
        cur->child->sib->d->cumMU = cur->d->cumMU*(1-mu);
      }
  
      if(!curIsLeaf) cur = cur->child;
      else {
        while(cur && !cur->sib) cur = cur->parent;
        if(cur) cur = cur->sib;
      }
    }
    total += m[i][Nvar+4];
  }
  return ss/total;  // if all data points have weight (m[i][Nvar+4] =1),
            // total = end - begin +1
}

double NomData::CrispError(Node* root, int K, int begin,int end)
{
  double ss = 0.0;
  double total =0.0;
  
  for(int i=begin;i<=end;i++)
  {
    // Going up : fix memberships && partial predictions
    Node* cur = root;
    root->d->cumMU =1.0;
    root->d->tempMU = 1.0;
    int Class = (int)m[i][Nvar-1];    
    m[i][Nvar+2] = 0.0;
    while(cur)
    {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {    
        if(Class != cur->iClass)
        {  
          
          double aux = m[i][Nvar+4]*cur->d->tempMU;
          if(cur->parent) aux *= cur->parent->d->cumMU;
          ss += aux;
          if(aux>0.99) m[i][Nvar+2] = 1.0;
        }
      } else {

        double mu =  cur->ParamMem(m[i],NVarsOrdSplit);
        cur->child->d->tempMU = mu;
        cur->child->d->cumMU = cur->d->cumMU*mu;
        cur->child->sib->d->tempMU = (1-mu);
        cur->child->sib->d->cumMU = cur->d->cumMU*(1-mu);
      }

      if(!curIsLeaf) cur = cur->child;
      else {
        while(cur && !cur->sib) cur = cur->parent;
        if(cur) cur = cur->sib;
      }
    }
    total += m[i][Nvar+4];
  }  
  return ss/total;  // if all data points have weight (m[i][Nvar+4] =1),
            // total = end - begin +1
}

double NomData::Fuzziness(int i, double fuzz)
{

  int N = NTrain + NSel;

  if(i<N) 
    return fuzz;
  if(m[i][Nvar] > m[N-1][Nvar])
    return m[N-1][Nvar+1];
  if( m[i][Nvar] <  m[0][Nvar])
    return 0;

  double val = N/2.0;
  double delta = val;

  while((int)(delta /= 2.0) >= 1)
  {
    if(m[i][Nvar] > m[(int)val][Nvar])
      val += delta;
    else val -= delta;
  }
  
  return m[(int)val][Nvar+1];
}




///////////////////////////////////////////////////////
//
//     Defuzzify (new)
//
////////////////////////////////////////////////////////
double NomData::Defuzzify(Node* root,  int K, int begin,int end, double* errQuartile)
{
  return 1.0;
}


//*/

///////////////////////////////////////






double NomData::DegreeOfFuzz(Node* root, int K)
{
  double total =0.0;
  double fuzz =0.0;
  int begin =0;    
  int end = NTrain-1;
  for(int i= begin;i<= end ;i++)
  {
    // Going up : fix memberships && partial predictions
    Node* cur = root;

    root->d->cumMU =1.0;
    root->d->tempMU = 1.0;
      
    while(cur)
    {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {
        double aux = cur->d->tempMU*m[i][Nvar+4];
        if(cur->parent) aux *= cur->parent->d->cumMU;
      
        total += aux;
        fuzz += aux*(1.0-aux);

        if(i==begin)
        {
          for(int j=0; j<NumClass; ++j)
            cur->NodePop[j]  = 0.0;
        }
        
        cur->NodePop[(int)m[i][Nvar-1]] += aux;
        
        if (i == end)   
          cur->iClass = WhichClass(cur->NodePop);
          
        
        
      } else {

        double mu =  cur->ParamMem(m[i],NVarsOrdSplit);
        cur->child->d->tempMU = mu;
        cur->child->d->cumMU = cur->d->cumMU*mu;
        cur->child->sib->d->tempMU = (1-mu);
        cur->child->sib->d->cumMU = cur->d->cumMU*(1-mu);
      }

      if(!curIsLeaf) cur = cur->child;
      else {
        while(cur && !cur->sib) cur = cur->parent;
        if(cur) cur = cur->sib;
      }
    }
  }  
  return fuzz/total;
}






double NomData::OrdCost(Node* root, int K, int begin, int end)
{
  double ss = 0.0;

  double alpha = degree;

  int i;
  for(i=begin;i<=end;i++)
  {
  
    // Going up : fix memberships && partial predictions
    Node* cur = root;
    root->d->cumMU =1.0;
    root->d->tempMU = 1.0;
    int Class = (int)m[i][Nvar-1];    
    while(cur)
    {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {
        double aux = cur->d->tempMU*m[i][Nvar+4];
        if(cur->parent) aux *= cur->parent->d->cumMU;
        
        if(i==begin)
        {
          cur->d->storage = new double*[1];
          cur->d->storage[0] = new double[NumClass];
          for(int j=0; j<NumClass; ++j)
            cur->d->storage[0][j] = cur->NodePop[j]  = 0.0;
          cur->d->MembershipTotal =0.0;
        }

        cur->NodePop[Class] += aux;
        cur->d->MembershipTotal += aux;
        if(i == end)
        {  
          cur->d->norm =0.0;
          for(int j =0; j<NumClass;++j)
            cur->d->norm +=pow(cur->NodePop[j],alpha);
        }

      } else {

        double mu =  cur->ParamMem(m[i],NVarsOrdSplit);
        cur->child->d->tempMU = mu;
        cur->child->d->cumMU = cur->d->cumMU*mu;
        cur->child->sib->d->tempMU = (1-mu);
        cur->child->sib->d->cumMU = cur->d->cumMU*(1-mu);
      }

      if(!curIsLeaf) cur = cur->child;
      else {
        while(cur && !cur->sib) cur = cur->parent;
        if(cur) cur = cur->sib;
      }
    }
  }


  for(i=begin;i<=end;i++)
  {
  
    // Going up : fix memberships && partial predictions
    Node* cur = root;
    root->d->cumMU =1.0;
    root->d->tempMU = 1.0;
    int Class = (int)m[i][Nvar-1];    
    while(cur)
    {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {
        double aux = cur->d->tempMU*m[i][Nvar+4];
        if(cur->parent) aux *= cur->parent->d->cumMU;
        if(i == begin)
        {
          if(aux > 1e-10)
          {
            for(int j =0; j<NumClass; ++j)
            {
              for(int k =0; k<NumClass; ++k)
              {
                double aux0 = ((k==j) - pow(cur->NodePop[k],alpha)/cur->d->norm);
                cur->d->storage[0][j] += aux0*aux0;
              }
  
            }
          
          } else {
            for(int j=0;j<NumClass; ++j)
            cur->d->storage[0][j] = 1.0;
          }
        }
          
        ss += aux*cur->d->storage[0][Class];
        if(i==end)
        {
          delete[] cur->d->storage[0];
          delete[] cur->d->storage;
        }
      } else {
        double mu = cur->ParamMem(m[i],NVarsOrdSplit);
        cur->child->d->tempMU = mu;
        cur->child->d->cumMU = cur->d->cumMU*mu;
        cur->child->sib->d->tempMU = (1-mu);
        cur->child->sib->d->cumMU = cur->d->cumMU*(1-mu);
      }

      if(!curIsLeaf) cur = cur->child;
      else {
        while(cur && !cur->sib) cur = cur->parent;
        if(cur) cur = cur->sib;        
      }  
    }
      
    
  }

  return ss;

}


void NomData::Filter(int index, double* Class)
{
  SortOn(0,0,NTotal);
  if(index <0)
  {
    for(int i=0; i<NTotal; ++i)
      Class[i] = m[i][Nvar-1];
    return;
  }

  for(int i=0; i < NTotal; ++i)
    m[i][Nvar-1] = 0.5 + ((int)Class[i] == index);
}

void NomData::AssignClass(int labelClass, Node* root, int K, double* Item, int begin,int end)
{
  Item = new double[end-begin+1];
  int i;
  for(i = begin; i<=end; ++i)
  {
    // Going up : fix memberships && partial predictions
    Node* cur = root;
    double* tempPop = new double[NumClass];
    
    root->d->cumMU =1.0;
    root->d->tempMU = 1.0;
    int j;
    for(j=0;j<NumClass; ++j)
      tempPop[j] =0.0;
    
    
    while(cur)
    {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {    
        double aux = cur->d->tempMU;
        if(cur->parent) aux *= cur->parent->d->cumMU;
      /*/
      //   Alternative 1: Assign weights in proportion to training set proportions in one node
      //    ASG: NOT GOOD
        for(j=0;j<NumClass; ++j)
        {
           tempPop[j] += aux*cur->NodePop[j]/cur->d->MembershipTotal;          
        }

      /*/
      //  Alternative 2: Assign all the weight to the node  class  
      //
        tempPop[cur->iClass] += aux;
      //*/

      } else {

        double mu =  cur->ParamMem(m[i],NVarsOrdSplit);
        cur->child->d->tempMU = mu;
        cur->child->d->cumMU = cur->d->cumMU*mu;
        cur->child->sib->d->tempMU = (1-mu);
        cur->child->sib->d->cumMU = cur->d->cumMU*(1-mu);
      }

      if(!curIsLeaf) cur = cur->child;
      else {
        while(cur && !cur->sib) cur = cur->parent;
        if(cur) cur = cur->sib;
      }
    }

    double max = -1.0;
    int Class = -1;
    for(j=0; j<NumClass; ++j)
    {
      if(max <0.0 || tempPop[j] > max)
      {
        max = tempPop[j];
        Class = j;
      }
    }

    delete[] tempPop;
    
    Item[i] = (Class == labelClass);
  
  }    
}
//Return the class of the element number idat, if it has many classes the percent
//of class membership is stored in Class
int NomData::Classify(int idat, double **Class, Node *root, int K, Node**finalnode)
{
  bool left;
  Node *cur = root; 
  double valor;

  while(cur) {

    if (cur->IsLeaf(K)) {
      if (finalnode) *finalnode = cur;
      break;
    }

    switch(cur->SplitType) {
      case NOM:
        if((1<<(int)m[idat][cur->att]) & cur->NomSplit) {
          left = true;
        }
        else {
          left = false;
        }
        break;
      case ORD:
        valor = -cur->fSplit;

        if(cur->coeff && !cur->IsUnivariateSplit()) {//!cur->univariate_split) { //
          for(int j=0;j<NVarsOrdSplit;j++) {
            valor += cur->coeff[j] * m[idat][j];
          }
        }
        else {
          valor += m[idat][cur->att];
        }

        if(valor <=0.0) {
          left = true;
        }
        else {
          left = false;
        }
        break;
      default:
        break;
    }
    cur = left ? cur->child : cur->child->sib;
  }

  return cur->iClass;



/*
  Node *cur = root;
  if (cur->IsLeaf(K)) {
    if (finalnode) *finalnode = cur;
    return cur->iClass;
  }
  else if (idat == GetLeftLast(cur, idat, idat, K))
    return Classify(idat, Class, cur->child, K, finalnode);
  else
    return Classify(idat, Class, cur->child->sib, K, finalnode);*/
}
//Calculate the prune errors and pops forthetree, it alse rearranges the data
void NomData::CalculateAndArrangeTestData(Node * root)
{
/*  n->RPrune = 0;
  n->d->firstTest = Ndata; n->d->lastTest = Ndata+ NPrune - 1;
  FindMembershipTest(n);
  for(int i = n->d->firstTest; i <= n->d->lastTest; ++i)
  {
    if((int)m[i][Nvar-1] != n->iClass)
      n->RPrune += m[i][Nvar+3];
  }*/
}
double NomData::NextDifferentValue(int var, int index)
{
  double val = m[index][var];
  for(int i=index+1;i<NTotal;i++) {
    if (val!=m[i][var]) return m[i][var];
  }
  return val;
}
NomData *NomData::GenerateSinteticData(int magnitude, int ini, int fin)
{
  int ntot = (fin-ini+1)*magnitude;
  NomData *nd = new NomData(ntot, Nvar);
  nd->Nord  = Nord;
  nd->Nfuzz = Nfuzz;
  nd->Nnom  = Nnom;
  nd->Nordfuzz = Nord+Nfuzz;
  nd->DepVarType = DepVarType;
  nd->NumClass = NumClass;
  int DepVarNom = DepVarType == NOM ? 1 : 0;
  nd->NomTerms = new std::vector<std::string>[Nnom+DepVarNom];
  for(int i=0; i<Nnom+DepVarNom; i++)
    for(unsigned j=0; j<NomTerms[i].size(); j++)
      nd->NomTerms[i].push_back(NomTerms[i][j]);

  //Copio los datos magnitude veces
  for(int i=ini;i<=fin;i++) {
    for(int j=0;j<magnitude;j++) {
      for(int k=0;k<Nvar;k++) {
        nd->m[i*magnitude+j][k] = m[i][k];
      }
    }
  }
  //Desvio los datos entre su valor superior e inferior
  double desde, rango, ant;
//  nd->SortOn(0, 0, ntot-1);
//  rango = nd->NextDifferentValue(0, magnitude) - nd->m[0][0];
//  desde = nd->m[0][0] - rango;
//  rango *= 2;
  int k;
  for(k=0;k<nd->Nordfuzz;k++) {
    nd->SortOn(k, 0, ntot-1);
    SortOn(k, ini, fin);
    rango = NextDifferentValue(k, ini) - m[ini][k];
    desde = m[ini][k] - rango;
    ant = m[ini][k];
    rango *= 2.0;
    for(int i=ini;i<=fin;i++) {
      for(int j=0;j<magnitude;j++) {
        nd->m[i*magnitude+j][k] = drand()*rango + desde;
      }
      desde = ant;
      if (i<fin-1) {
        rango = NextDifferentValue(k, i+1) - ant;
        ant = m[i+1][k];
      }
      else {
        rango = (m[fin][k] - ant)*2.0;
      }
    }
  }
  for(;k<Nvar-1;k++) {
  }


  nd->Pop  = new double[NumClass];
  nd->PopR = new double[NumClass];
  nd->PopL = new double[NumClass];
  nd->ponder = false;
  nd->Ponderate(Ndata);
  nd->iClass = WhichClass(nd->Pop);
  return nd;
}
NomData *NomData::GenerateSinteticData2(int ntot, int ini, int fin)
{
  NomData *nd = new NomData(ntot, Nvar);
  nd->Nord  = Nord;
  nd->Nfuzz = Nfuzz;
  nd->Nnom  = Nnom;
  nd->Nordfuzz = Nord+Nfuzz;
  nd->DepVarType = DepVarType;
  nd->NumClass = NumClass;
  int DepVarNom = DepVarType == NOM ? 1 : 0;
  nd->NomTerms = new std::vector<std::string>[Nnom+DepVarNom];
  for(int i=0; i<Nnom+DepVarNom; i++)
    for(unsigned j=0; j<NomTerms[i].size(); j++)
      nd->NomTerms[i].push_back(NomTerms[i][j]);

  //Desvio los datos entre su valor superior e inferior
  double desde, rango;
  int k;
  for(k=0;k<nd->Nordfuzz;k++) {
    SortOn(k, ini, fin);
    for(int j=0;j<ntot;j++) {
      int idx = ini + (int)(drand()*(fin-ini-1));
      desde = m[idx][k];
      rango = m[idx+1][k] - m[idx][k];
      nd->m[j][k] = drand()*rango + desde;
    }
  }
  for(;k<Nvar-1;k++) {
  }


  nd->Pop  = new double[NumClass];
  nd->PopR = new double[NumClass];
  nd->PopL = new double[NumClass];
  nd->ponder = false;
  nd->Ponderate(Ndata);
  nd->iClass = WhichClass(nd->Pop);
  return nd;
}

int NomData::AddClass(string name)
{
  int nnom = GetNumVarNom();
  unsigned ncl = NomTerms[nnom].size();
  int class_num = (int)DameRepresentacion(name, nnom, true);
  if (ncl < NomTerms[nnom].size()) NumClass++;
  return class_num; 
}

//-------------------------------------------------------------  Clone  -----
//---------------------------------------------------------------------------
Data* NomData::Clone(int ini, int fin, int factor)
{
  return Data::Clone(ini, fin, factor);
}

//--------------------------------------------------------  ClassNoise  -----
//---------------------------------------------------------------------------
NomData *NomData::ClassNoise(double prob)
{
  NomData *d = (NomData*)Clone();
  double nc = d->NumClass - 1;
//int k=0;

  for(int i=0;i<d->GetNTotal();i++) { //Para cada dato
    //Modificamos la etiqueta de clase con probabilidad "prob"
    double r = (1.0 * TDebugRand::Rand() / (RAND_MAX + 1.0));
    if (r < prob) {
      int c = ( nc * TDebugRand::Rand() / (RAND_MAX + 1.0));
      if (d->GetDatClass(i) <= c) c++;
//cout<<d->m[i][Nvar-1]<<"->" << c << endl;
      d->m[i][Nvar-1] = c; //Aqui se cambia la clase
//k++;
    }
  }
//cout<<"tot " << k << "/" << d->GetNTotal() << endl;
  return d;
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//----------------------------------------------------------  cartdata  -----
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------



Data::Data(std::istream &in, void (*_console)(char *))
{
//Constructor de Data a partir de un fichero con el siguiente formato:
//Fila1: "Número de variables(Nvar)" "Numero de ejemplos(NTotal)"
//Fila2: Nvar indicadores con el tipo de variable (0):Ordinal (1):Difusa (2): Nominal
//Fila3: Nombres de variables
//Datos un ejemplo por fila

  int *tipos;
  console = _console;

  name = ""; 
  in >> Nvar >> NTotal;
  Nord=Nfuzz=Nnom=0;
  tipos = new int[Nvar];
  for(int i=0;i<Nvar-1;i++) {
    in >> tipos[i];
    if (tipos[i]==0) Nord++;
    else if (tipos[i]==1) Nfuzz++;
    else Nnom++;
  }
  in >> tipos[Nvar-1];
  DepVarType = (tipos[Nvar-1]==0) ? ORD : (tipos[Nvar-1]==1) ? FUZZ : NOM;
  delete []tipos;
//  DepVarType = ds->GetTypeDepVar();

//  Nvar = ds->GetNumVars();
//  Nnom = ds->GetNumNoms();
//  Nord = ds->GetNumOrds();
//  Nfuzz = ds->GetNumFuzz();
  Nordfuzz = Nord+Nfuzz;
  NVarsOrdSplit = Nvar-1;//???????ds->GetNumVarsOrdSplit();
//  FuzzVars = ds->GetFuzzVars();
  NTest = NTotal/2;
  NSel = 0;
  NPrune = NTest/2;
//  NTotal = ds->GetNumData();
  NTrain = NTotal-NTest;
  NGrow = NTrain - NPrune;
  Ndata = NGrow;
  Ncv = 10;//ds->ncv;
  cartBeta = 0.1;//ds->beta;
  cartEps = 0.01;//ds->epsilon;
  multsplits = false;//ds->multsplits;
  plm = false;//ds->plm;
  fuzzify = false;//ds->fuzzify;
  n_reservados = NTotal;
  m = new double*[NTotal];
  int i;
  for(i=0;i<NTotal;i++)
    m[i] = new double[Nvar+7];  // Extra elements for working space later
    // 0->Nvar-2: Independent variables
    // Nvar-1: Dependent variable
    // Nvar,Nvar+1: Multiple split working space
    // Nvar+2: CV tree label
    // Nvar+3: Membership
    // Nvar+4: Weight
    // Nvar+5: Initial position
    // Nvar+6: Group number for Iterative growing pruning method
    // Nvar+7: Auxiliar

  //    IN NOMINAL OPTIMIZATION
    // Nvar Absolute Fuzziness
    // Nvar+1 Fuzziness by rank
    // Nvar+2: Error of point after defuzzification
    // Nvar+3: Predicted class
    // Nvar+4: Weight
    // Nvar+5: Initial position

  vscale = new double[NVarsOrdSplit+2];
  vshift = new double[NVarsOrdSplit+2];
  Gvshift = new double[Nvar + 1];
  Gvscale = new double[Nvar + 1];
  VarNames = new string[Nvar];
  int DepVarNom = (DepVarType == NOM ? 1:0);
  NomTerms = new std::vector<std::string>[Nnom + DepVarNom];
/*  for(i=0; i<Nnom+DepVarNom; ++i)
    for(int j=0; j < ds->GetSizeNomTerm(i); ++j)
      NomTerms[i].push_back((ds->GetNomTerm(i,j)));
  */
//  double amax,amin;
  for(i=0; i<=Nvar; ++i)
  {
    Gvshift[i] = 0.0;
    Gvscale[i] = 1.0;
  }

  char kk[256];
  for(i=0; i<Nvar; ++i)
  {
    in >> kk;
    VarNames[i] = kk;
  }

  for(i=0;i<NTotal;i++)
  {
    int j;
    string leo;
    m[i][Nvar+5] = i + 0.5;
    //  ordinal and fuzzy
    for(j = 0; j < Nordfuzz; ++j)
      in >> m[i][j];
    //Nominal
    for(; j < Nvar-1; ++j) {
      in >> leo;
      bool noesta = true;
      for(unsigned k=0; k < NomTerms[j-Nordfuzz].size() ; k++)
        if (0==NomTerms[j-Nordfuzz].at(k).compare(leo))
        {
          m[i][j] = k + 0.5;
          noesta = false;
          break;
        }
      if (noesta)
      {
        NomTerms[j-Nordfuzz].push_back(leo);
        m[i][j] = NomTerms[j-Nordfuzz].size() - 0.5;
      }
    }

    //Depvar
    if (DepVarNom) {
      vector<string>::iterator k;
      int ic=0;
      in >> leo;
      bool noesta = true;
      for(k=NomTerms[Nnom].begin(), ic=0; k!=NomTerms[Nnom].end() ; ic++, k++) {
        //Se ordenan las etiquetas de clase por nombre
        int comp = /*NomTerms[Nnom].at*/(*k).compare(leo);
        if (0==comp)
        {
          m[i][Nvar-1] = ic + 0.5;
          noesta = false;
          break;
        }
        else if (comp>0) {
          break;
        }
      }
      if (noesta)
      {
        if (k==NomTerms[Nnom].end()) {
          NomTerms[Nnom].push_back(leo);
        }
        else {
          NomTerms[Nnom].insert(k, leo);//push_back(leo);
          //Reasignamos las etiquetas
          for (int j=0;j<i;j++) {
            if (m[j][Nvar-1]>ic)
              m[j][Nvar-1] += 1;
          }

        }
        m[i][Nvar-1] = ic + 0.5;
      }
    }
    else //ordinal
      in >> m[i][Nvar-1];


    m[i][Nvar+2] = (int) (1+(Ncv)*drand()); // Label for CV group
    m[i][Nvar+1] = m[i][Nvar] = drand();
    m[i][Nvar+3] = 1.0; // Membership
    m[i][Nvar+4] = 1.0;  // Weight (only for NomData with non-homogeneous distribution of classes)
    m[i][Nvar+GetGroupIndex()] = 1.0;
  }

  // Uncomment to scramble test and training data
//  if(NTest && ds->GetScramble())
//    SortOn(Nvar,0,NTotal-1);  // Scramble test and training data
  GroupsCreated = false;
  GroupCount = 0;
}
Data::Data(char *filename)
{
  FILE *f =fopen(filename, "rt");
  if (!f) throw new exception();
  
  int ndatos;
  DefaultInitValues();
  fscanf(f, "%d", &ndatos);
//  printf("Getting names...\n");
  GetNames(f);
  redim(ndatos);
//  printf("Got names\nGetting data...\n");
  GetData(f);
  fclose(f);
}
Data::Data(int NTotalData, int _NVar)
{
  DefaultInitValues();
  SetDefaultNames(_NVar);
  redim(NTotalData);
//  init(NTotalData, _NVar);
}

Data::Data(int NTotalData, int _NVar, double **data)
{
  DefaultInitValues();
  SetDefaultNames(_NVar);
  redim(NTotalData);

  for(int i=0;i<NTotal;i++) {
    for (int j=0;j<Nvar-1;j++) m[i][j] = data[i][j];
    //Depvar
    char txt[256];
    sprintf(txt, "%g", data[i][Nvar-1]);
    m[i][Nvar-1] = DameRepresentacion(txt, 0);
  }
}


Data::~Data()
{
  for(unsigned i=0;i<reg_eventos.size();i++) reg_eventos[i]->OnDelete(this);

  for(int i = 0; i < NTotal; i++) {
    delete[] m[i];
  }

  delete[] m; m = NULL;
  delete[] vshift; vshift = NULL;
  delete[] vscale; vscale = NULL;
  delete[] Gvshift; Gvshift = NULL;
  delete[] Gvscale; Gvscale = NULL;
//  if(VarNames) delete[] VarNames; VarNames =NULL;
//  if(NomTerms) delete[] NomTerms;  NomTerms =NULL;
  if (GroupCount) delete []GroupCount;

  delete []NomTerms;
  delete []VarNames;
  delete []OriginalVarType;
}
void Data::SetDefaultNames(int nvar)
{
  //Variables "horizontales" (o las que tienen que ver con atributos...)
  Nvar = nvar;
  Nord = Nvar-1;
  Nfuzz=Nnom=0;
  Nordfuzz = Nord+Nfuzz;
  NVarsOrdSplit = Nord;
  DepVarType = NOM;
  OriginalVarType = new VarType[Nvar];
  NomTerms = new std::vector<std::string>[1];

  //Nombres y tipo de variables
  VarNames = new string[Nvar];
  char kk[256];
  for(int i=0; i<Nvar-1; ++i)  {
    sprintf(kk, "Att%d", i);
    VarNames[i] = kk;
    OriginalVarType[i] = ORD;
  }
  VarNames[Nvar-1] = "class";
  OriginalVarType[Nvar-1] = NOM;
  
  //A saber que son estas variables que vienen
  vscale = new double[NVarsOrdSplit+2];
  vshift = new double[NVarsOrdSplit+2];
  Gvshift = new double[Nvar + 1];
  Gvscale = new double[Nvar + 1];
  for(int i=0; i<=Nvar; ++i) {
      Gvshift[i] = 0.0;
      Gvscale[i] = 1.0;
  }
}
void Data::DefaultInitValues()
{
  //Variables "horizontales" (o las que tienen que ver con atributos...)
  Nvar = 0;
  Nord = 0;
  Nfuzz=Nnom=0;
  Nordfuzz = 0;
  NVarsOrdSplit = 0;
  DepVarType = NOM;

  //variables
  Ncv = 10;
  cartBeta = 0.1;
  cartEps = 0.01;
  multsplits = false;
  plm = false;
  fuzzify = false;
  console = salypimienta;
  min_in_child = 1;

  //Matriz de datos
  n_reservados = 0;
  m = 0;
  NTotal = 0;
  //Variables "verticales" (o las que tienen que ver con no. ejemplos de train, test...)
  NSel = 0;
  SetNTrain(NTotal/2);

  //Iterative growing and pruning mathod variables
  GroupsCreated = false;
  GroupCount = 0;
}
void Data::redim(int ndatos)
{

  if (NTotal == ndatos) return;

  if (NTotal > ndatos) {//Reducimos
    for(int i = ndatos; i < NTotal; i++) {
      delete []m[i];
      m[i] = 0;
    }
  }
  else {
    double **mm;

    if (n_reservados > ndatos) {//Ampliando 
      mm = m;
    }
    else {//Ampliando y cambiando de puntero
      n_reservados = ndatos;//n_reservados * 2 > ndatos ? n_reservados * 2 : ndatos;
      mm = new double*[n_reservados];
      //Los datos extras los marcamos a cero para que casque si alguien se pasa
      for(int i = ndatos; i < n_reservados; i++) mm[i] = 0;
      for(int i = 0; i < NTotal; i++) mm[i] = m[i];
    }

    for(int i = NTotal; i < ndatos; i++) {
      // Extra elements for working space later
      mm[i] = new double[Nvar+7];
      // (0->Nvar-2) Independent variables ; (Nvar-1) Dependent variable
      for (int j=0;j<Nvar+1;j++) mm[i][j] = 0.0;
      // (Nvar) and (Nvar+1) Multiple split working space
      mm[i][Nvar+1] = mm[i][Nvar] = drand();
      // (Nvar+2) Label for CV group
      mm[i][Nvar+2] = (int) (1+(Ncv)*drand());
      // (Nvar+3) Membership
      mm[i][Nvar+3] = 1.0;
      // (Nvar+4) Weight (only for NomData with non-homogeneous distribution of classes)
      mm[i][Nvar+4] = 1.0;
      // (Nvar+5) Initial position
      mm[i][Nvar+5] = i + 0.5;
      // (Nvar+6) Group number for Iterative growing pruning method
      mm[i][Nvar+6] = 1.0;
    }
    if (m && m != mm) delete []m;
    m = mm;
  }

  NTotal = ndatos;
//cout << NTotal << "x" << Nvar << endl;
  if (NTotal != NSel+NTest+NTrain) SetNTrain(NTotal);
}
void Data::ResetCV(int Value)
{
  if (Value<0) Value = Ncv;
  Ncv = Value;
  CreateGroups(Ncv, -1, -1, true);
//  for(int i=0;i<NTotal;i++)
//    m[i][Nvar+2] = (int) (1+(Ncv)*drand()); // Label for CV group
}

Data* Data::DataFromFile(char *filename)
{
  string nombre = filename;
  int pos = nombre.find_last_of('.');
  string ext = nombre.substr(pos);
  if (ext.compare(".cre")==0) {
    return LoadFromCre(filename);
  }
  else if (ext.compare(".asc")==0) {
    return LoadFromAsc(filename);
  }
  else {
    return LoadFromCre(filename);
  }
  return 0;
}
Data* Data::LoadFromCre(char *filename)
{
//  ifstream ifs(filename);
//  return new NomData(ifs, 0);
  return new NomData(filename);
}
Data* Data::LoadFromAsc(char *filename)
{
  string sfilename =filename;
  NomData *dat = new NomData();

  //Busco fichero con definición de los atributos
  int pos = sfilename.find_last_of('/');
//km  string nf = (pos<0 ? "" : sfilename.substr(pos+1)) + "default.names";
  string nf = (pos<0 ? string("") : sfilename.substr(pos+1)) + "default.names";
  FILE *f=fopen(nf.c_str(), "rt");
  if (f) {
    dat->GetNames(f);
    fclose(f);
  }
  else {
    int NVar = 1;
    char linea[4096];
    FILE *f=fopen(filename, "rt");
    if (!f) return 0;
    fgets(linea, 4096, f);
    char *tok = strtok(linea, " \r\t\n");
    while(tok) {
      NVar++;
      tok=strtok(0, " \r\t\n");
    }
//    char *aux = linea;
//    while (1==sscanf(aux, "%lf", &dato)) {
//    while (1==fscanf(f, "%lf", &dato)) {
//      NVar++;
//      aux+=16;
//    }
    fclose(f);

    dat->SetDefaultNames(NVar);
  }  
  
  string filelabels = filename;
  filelabels.replace(filelabels.find("data"), (unsigned)4, "labels");

  dat->GetData(filename);
  dat->GetData((char*)filelabels.c_str(), dat->GetNumVar()-1);
  dat->init();
  
  return dat;
/*  int NVar = 1;
  double dato;
  FILE *f=fopen(filename, "rt");
  if (!f) return 0;
  char linea[4096];
  int tam=0;
  int NTotal=0;
  double **m=0;

  fgets(linea, 4096, f);

  char *aux = linea;
  while (1==sscanf(aux, "%lf", &dato)) {
    NVar++;
    aux+=16;
  }

  do {
    if (tam<=NTotal) {  
      int hasta = tam==0 ? 100 : tam*2;
      m = (double**)realloc(m, sizeof(double*)*hasta);
      for(int i=tam;i<hasta;i++)
        m[i] = new double[NVar];
      tam = hasta;
    }
    aux = linea;
    for (int i=0;i<NVar-1;i++) {
      sscanf(aux, "%lf", &m[NTotal][i]);
      aux+=16;
    }
    NTotal++;
  } while(fgets(linea, 4095, f));

  fclose(f);
  string nf = filename;
  string hk = "labels";
  nf.replace(nf.find("data"), (unsigned)4, hk);
  f=fopen(nf.c_str(), "rt");
  for (int i=0;i<NTotal;i++) {
    fscanf(f, "%lf", &m[i][NVar-1]);
  }
  fclose(f);

  NomData *dat = new NomData(NTotal, NVar, m);

  for (int j=0;j<tam;j++) 
    delete []m[j];
  delete []m;
  
  return dat;
*/
}
void Data::GetNames(char *nf)
{
  FILE *f = fopen(nf, "rt");
  if (!f) return;
  GetNames(f);
  fclose(f);
}
void Data::GetNames(FILE *f)
{
  int tipo, iNom;
  char linea[1024];
  char varname[1024];
  
  fscanf(f, "%d", &Nvar);
//  init(NTotal?NTotal:1, Nvar);
//  printf("Nomtermsi, ... to %d\n", Nvar);
  NomTerms = new std::vector<std::string>[Nvar];
  VarNames = new string[Nvar];
  vector<string> varnames;
  OriginalVarType = new /*crepo::*/VarType[Nvar];
  
  Nord = Nfuzz = Nnom = iNom = 0;
  for(int i=0;i<Nvar;i++) {
    fscanf(f, "%s %d", varname, &tipo);
    if (tipo==0) {
      varnames.insert(varnames.begin()+Nord+Nfuzz, varname);
      Nord++; 
      OriginalVarType[i]=ORD;
    }
    else if (tipo==1) {
      varnames.insert(varnames.begin()+Nord+Nfuzz, varname);
      Nfuzz++; 
      OriginalVarType[i]=FUZZ;
    }
    else {
      varnames.push_back(varname);
      OriginalVarType[i]=NOM;
      char *tok;
      Nnom++;
      fgets(linea, 1024, f);
      tok = strtok(linea, " \r\t\n");
      while(tok) {
        DameRepresentacion(tok, iNom);
//        NomTerms[iNom++].push_back(tok);
        tok=strtok(0, " \r\t\n");
      }
      iNom++;
    }
  }
  Nordfuzz = Nord+Nfuzz;
  NVarsOrdSplit = Nord;
  for (int i = 0; i < Nvar; i++) VarNames[i] = varnames[i];
  DepVarType = OriginalVarType[Nvar-1];
  if (DepVarType==NOM) Nnom--;

  //A saber que son estas variables que vienen
//  printf("vscale... variables to %d\n", NVarsOrdSplit+2);
  vscale = new double[NVarsOrdSplit+2];
  vshift = new double[NVarsOrdSplit+2];
  Gvshift = new double[Nvar + 1];
  Gvscale = new double[Nvar + 1];
  for(int i=0; i<=Nvar; ++i) {
      Gvshift[i] = 0.0;
      Gvscale[i] = 1.0;
  }

  redim(1);
}
void Data::GetData(char *nf, int colIni, int colFin)
{
  FILE *f = fopen(nf, "rt");
  if (!f) return;
  GetData(f, colIni, colFin);
  fclose(f);
}
void Data::GetData(FILE *f, int colIni, int colFin)
{
  char linea[0x100000];
  char *tok;
  int ilinea = 0;
  int iOrd, iNom, iIgn, iCol;
  int iOrdIni=0, iNomIni=0, iIgnIni=0;

  //Initialising counters
  if (colIni==-1) colIni = 0;
  if (colFin==-1 || colFin>=Nvar) colFin = Nvar-1;
  for (int i = 0; i < colIni; i++) {
    if (OriginalVarType[i]==ORD) iOrdIni++;
    else if (OriginalVarType[i]==NOM) iNomIni++;
    else if (OriginalVarType[i]==IGN) iIgnIni++;
  }

  //loop for reading the data
  while(fgets(linea, 0x100000, f)) {
    if ((int)strlen(linea) < Nvar) break;
    if (NTotal==ilinea) {//Reservamos más espacio
      redim(NTotal?NTotal*2:1);
    }
//    if (NTotal>=70 && 0==ilinea%(NTotal/70)) printf("=");
    iOrd = iOrdIni;
    iNom = iNomIni;
    iIgn = iIgnIni;
    iCol = colIni;
    for(tok = strtok(linea, " \r\n\t"); tok && iCol<=colFin; 
                                            tok = strtok(0, " \r\n\t"),iCol++) {
      if (OriginalVarType[iCol]==ORD) {
        if (tok[0]=='?') m[ilinea][iOrd] = DBL_MAX;
        else sscanf(tok, "%lf",  &m[ilinea][iOrd]);
        iOrd++;
      }
      else if (OriginalVarType[iCol]==NOM) {
        m[ilinea][Nordfuzz + iNom] = DameRepresentacion(tok, iNom);
        iNom++;
      }
      else if (OriginalVarType[iCol]==IGN) {
      }
    }
    ilinea++;
  }

  //Ajustamos el tamaño
  if (NTotal!=ilinea) redim(ilinea);
}
double Data::DameRepresentacion(string val, int ivarnom, bool force)
{
  vector<string>::iterator k;
  int ic;
  bool noesta = true;
  vector<string> &noms = NomTerms[ivarnom];
  double ret;
  
  //Busco 'val' en noms (teniendo en cuenta que la lista 
  //  de etiquetas que esta ordenada por nombre)
  for(k=noms.begin(), ic=0; k!=noms.end() ; ic++, k++) {
    int comp = (*k).compare(val);
    if (0==comp) {
      noesta = false;
      break;
    }
    else if (comp>0) {
      break;
    }
  }

  //Si no se encuentra el nombre se crea otro (siempr que force==true)
  if (noesta && force)  {
    if (k==noms.end()) noms.push_back(val);
    else {
      int imcol = Nordfuzz + ivarnom;
      noms.insert(k, val);
      //Reasignamos las etiquetas
      for (int j=0;j<NTotal && m[j][imcol]>0.0;j++) {
        if (m[j][imcol]>ic)
          m[j][imcol] += 1;
      }
    }
  }

    ret = 0.5 + ic;
    return ret;
}
void Data::PermuteAttributeValues(int iatt)
{

  for(int i = 0; i < NTotal; i++) {
    int elem = i + (int) ( (NTotal - i) * ( (double)rand() / ( RAND_MAX + 1.0 )) );
    double hold = m[i][iatt];
    m[i][iatt] = m[elem][iatt] ;
    m[elem][iatt] = hold;
  }
}
void Data::Scramble(int begin,int end)
{

  for(int i = begin; i <= end ; ++i)
  {
    m[i][Nvar] = drand();
  }
  SortOn(Nvar,begin, end);
}


void Data::ScaleVariableIQ(int iVar, int scfirst,int first, int last)
{
  SortOn(iVar,scfirst,last);
  int N = last-scfirst+1;
  double val1 = m[scfirst+N/4][iVar];
  double val2 = m[scfirst+(3*N)/4][iVar];

  if(val2-val1 > 1e-6)
  {
    vscale[iVar] = 1.0/(val2-val1);
    vshift[iVar] = 0.5-vscale[iVar]*val2;
    for(int i=first;i<=last;i++)
      m[i][iVar] = vscale[iVar]*m[i][iVar]+vshift[iVar];
  } else {
    vscale[iVar] = 1.0;
    vshift[iVar] = 0.0;
  }
}

double Data::UnScaleValueIQ(int iVar, double x)
{
  return (x -(vshift[iVar]))/(vscale[iVar]);
}
/*
  Tira los datos de begin a end por el árbol y cuenta el número
  de errores que se comenten en las hojas
*/
double Data::Error2(Node* root, int K, int begin,int end)
{
  Node *cur=root;
  int i, leftLast=0;
  double nerror=0.0;

  if (end<begin) return 0.0;

  if (cur->IsLeaf(K)) {
    for(i=begin;i<=end;i++) {
      int Class = (int)m[i][Nvar-1];
      if(Class != cur->iClass) nerror+=m[i][Nvar+4];
    }
  }
  else /*if (begin<=(leftLast = GetLeftLast(cur, begin, end, K)))*/
//  else if (cur->SplitType == ORD)
  {
/*    for(i=begin;i<=end;i++) m[i][Nvar] = m[i][cur->att] - (cur->fSplit);
    // Sort
    SortOn(Nvar, begin, end);
    // and find boundary
    leftLast = begin-1;
    while(leftLast < end && m[leftLast+1][Nvar] < 0.0) leftLast++;*/
    leftLast = GetLeftLast(cur, begin, end, K);
    nerror = Error2(cur->child, K, begin, leftLast) +
                                    Error2(cur->child->sib, K, leftLast+1, end);
  }
/*  else if(cur->SplitType == NOM)
  {
    int MembersLeft = SortNomOn(cur->NomSplit, cur->att, begin, end);
    leftLast = begin + MembersLeft-1;
    nerror = Error2(cur->child, K, begin, leftLast) +
                                    Error2(cur->child->sib, K, leftLast+1, end);
  }*/

  return nerror;
}

void Data::UnScaleVariableIQ(int iVar,int first, int last)
{
  int ilabel;
  if(iVar == Nvar)
    ilabel = NVarsOrdSplit+1;
  else if(iVar == Nvar-1)
    ilabel = NVarsOrdSplit;
  else ilabel = iVar;

  for(int i=first;i<=last;i++)
    m[i][iVar] = (m[i][iVar] -vshift[ilabel])/vscale[ilabel];

  vscale[ilabel] = 1.0; vshift[ilabel] = 0.0;
}
/*******************
Divide los datos de first a last en nGroups grupos tal que cada grupo
quede con aprox. mismo no. de elementos de cada clase. Guarda el no.
de grupo en m[first..last][Nvar+Data::GetGroupIndex()]
*******************/
void Data::CreateGroups(int nGroups, int first, int last, bool doCVFolds)
{
  int j;
  unsigned i;
  int iclass;
  int nelem;
  int icol = doCVFolds ?  Nvar+GetCvTreeLabelIndex() : Nvar+GetGroupIndex();
  int igrp = doCVFolds ?  1 : 0;
  vector<int> nnclass(0);

  first = first<0 ? 0 : first;
  last  = last <0 ? NTotal-1 : last;

  GroupsCreated = true;
  if (GroupCount) delete []GroupCount;
  GroupCount = new int[nGroups];
  for(int i=0;i<nGroups;i++) GroupCount[i]=0;

  //Leave-one-out
  if (nGroups==last-first+1) {
    for(j=first;j<=last;j++) {
      m[j][icol] = j + igrp;
      GroupCount[j] = 1;
    }
    return;
  }

  //Sort by class
  SortByClass(first, last);

//  srand(time(0) + rand());//randomize();
  iclass=-1;
  for(i=first;i<=(unsigned)last;i++) {
    //A random number is thrown
    m[i][Nvar]=TDebugRand::Rand();//%10000;
    //The number of classes and number of elements of each class are obtained
    if (iclass!=GetDatClass(i)){ //((int)m[i][Nvar-1])) {
      iclass=GetDatClass(i);     //(int)m[i][Nvar-1];
      nnclass.resize(nnclass.size()+1);
    }
    nnclass.at(nnclass.size()-1)++;
  }

  //Each class is ordered by the random number just thrown
  int ClassBegin=first;
//  aux = 10000/nGroups;
//**************************
//PRUEBASSSSSSSSSSSSSSSSSSSS
if (!Proportional) {
  nnclass.clear();
  nnclass.push_back(last-first+1);
}
//**************************
  int gr = TDebugRand::Rand()%nGroups;
  for(i=0;i<nnclass.size();i++) {
    SortOn(Nvar, ClassBegin, ClassBegin+nnclass.at(i)-1);
    nelem = nnclass.at(i)/nGroups;
    nelem = nelem*nGroups + ClassBegin;
    //A group number is assigned
    for(j=ClassBegin;j<nelem;j++) {
      m[j][icol] = igrp + j%nGroups;
      GroupCount[(int)m[j][icol]-igrp]++;
    }
    //The elements from nGroups*nelem to nnclass.at(i) are added at random to
    // the groups
    for(;j<ClassBegin+nnclass.at(i);j++) {
      m[j][icol] = igrp + gr%nGroups;
      gr++;
      GroupCount[(int)m[j][icol] - igrp]++;
    }
    ClassBegin += nnclass.at(i);
  }

  //Sort by group
  SortOn(icol, first, last);
}

static int theSortAtt;
static bool InvertSort=false;
static int CV;

int compAtt(const void* AA, const void *BB)
{
  int ret;
  double* A = *((double**) AA);
  double* B = *((double**) BB);

  if(A[theSortAtt] > B[theSortAtt]) ret = 1;
  else if(A[theSortAtt] < B[theSortAtt]) ret = -1;
  else ret = 0;

  if (InvertSort) ret = -ret;

  return ret;
}
int compCV(const void* AA, const void *BB)
{
  int A  = (int)( (*((double**) AA))[theSortAtt]+0.5);
  int B = (int)((*((double**) BB))[theSortAtt]+0.5);
  if(A == CV && B != CV) return 1;
  if(A != CV && B == CV) return -1;
  return 0;
}
void Data::SortByClass(int first, int last, bool _InvertSort)
{
  first = first < 0 ? 0 : first;
  last = last < 0 ? NTotal - 1 : last;

  SortOn(Nvar-1, first, last, _InvertSort);
}
void Data::SortByGroup(int first, int last, bool _InvertSort)
{
  first = first < 0 ? 0 : first;
  last = last < 0 ? NTotal - 1 : last;

  SortOn(Nvar+Data::GroupIndex, first, last, _InvertSort);
}
void Data::OriginalOrder(bool _InvertSort)
{
  SortOn(Nvar+5, 0, NTotal-1, _InvertSort);
}
void Data::SortOn(int iVar, int first, int last, bool _InvertSort)
{
  InvertSort = _InvertSort;
  theSortAtt = iVar;
  if (last<0) last = GetNTotal()-1;
  qsort(m+first,last-first+1,sizeof(double**),compAtt);
}

int Data::PreSort(int first, int last)
{
  int i;
  SortOn(Nvar+3,first,last);
  for(i=first; i<=last; ++i)
  {
    if (m[i][Nvar+3] > MIN_BELONG_LOW)
      break;
  }
  return i;
}


void Data::SortCV(int cv, int begin, int end)
{
  CV = cv;
  theSortAtt = Nvar+2;
  if(end < 0) end = NGrow-1;
  qsort(m+begin,end-begin+1,sizeof(double**),compCV);
}

int Data::FindFirstTest(int cv, int first, int last)
{
  int ret = last+1;
  while(ret > first && fabs(m[ret-1][Nvar+2]-cv) < 0.5)
    ret--;
  return ret;
}

int Data::SetCV(int cv)
{
  Ndata = NGrow;

  if(cv != 0)
  {
    SortCV(cv);
    while(Ndata && fabs(m[Ndata-1][Nvar+2] - cv) < 0.1) Ndata--;
    if(DepVarType == NOM) Ponderate(Ndata);
  }


  return Ndata;
}

void Data::MarkOrder()
{ //It wont work if someone touches the values in column c
  // before reseting the order
  int c = Data::MultipleSplitWs2 + Nvar;
  for(int i = 0; i < GetNTotal(); i++) 
    m[i][c] = i;
}
void Data::ResetOrder()
{
  int c = Data::MultipleSplitWs2 + Nvar;
  SortOn(c);
}
void Data::SetNTrain(int first, int last)
{
  int c = Data::MultipleSplitWs2 + Nvar;
  for(int i=0; i<first;i++) m[i][c] = 1;
  for(int i=first; i<=last;i++) m[i][c] = 0;
  for(int i=last+1; i<GetNTotal();i++) m[i][c] = 2;
  SortOn(c, 0, GetNTotal()-1);
  SetNTrain(last-first+1);
}


int Data::FindBestSplit(Node* n, int Nmin, int* leftLast, double *chiBest,
                                              vector<bool> *attributes_to_use)
{
  double chi;
  int iSplit, i;

  n->att = -1; // Indicates an error
  *chiBest = -1;
  n->NomSplit=0;
  n->SplitType = ERR;

  // focus on data with significant membership (alpha-cut)
  int scfirst = PreSort(n->first, n->last); 

  /////////////////////////////////////////////////////////////////
  //
  //    Simple ordinal splits
  //
  for(i=0;i<NVarsOrdSplit;i++) {
    if (attributes_to_use && !((*attributes_to_use)[i])) continue;
    iSplit=-1;
    chi = FindOrdSplit(scfirst, n->last, i, Nmin,n->d->MembershipTotal, &iSplit);
    if(iSplit >=0 && (*chiBest < 0.0 || chi < *chiBest)) {
      n->att = i;
      int iS2 = iSplit+1;
      while(m[iS2][Nvar+3] < MIN_BELONG_LOW || m[iSplit][i] == m[iS2][i])
        iS2++;
      n->fSplit = 0.5*(m[iSplit][i]+m[iS2][i]);
      n->NomSplit = 0;
      n->SplitType = ORD;
      n->d->min_in_child = min_in_child;
      *chiBest = chi;

      if(fuzzify) {
        chi = FuzzifySplit(n,*chiBest,scfirst,Nmin);
        if(chi>0.0 && chi<*chiBest) {
          n->SplitType = FUZZ;
          *chiBest = chi;
        }
      }
    }
  }


  double oldSplit = n->fSplit;

  /////////////////////////////////////////////////////////////////
  //
  //    Fuzzy splits
  //
  for(i = Nord; i< Nordfuzz; ++i) {
    if (attributes_to_use && !((*attributes_to_use)[i])) continue;
    Property* tempprop = FuzzVars;
    while(tempprop->Name != VarNames[i])
      tempprop = tempprop->next;

    FuzzySet* FuzzSplit = tempprop->Types;
    chi = FindFuzzSplit(scfirst, n->last, i, Nmin, n->d->MembershipTotal, 
                                                                   &FuzzSplit);
    if(chi>=0.0  && (*chiBest < 0.0 || chi < *chiBest)) {
      n->att = i;
      if(n->d->fuzzindex) {
        delete n->d->FuzzSplit;
        n->d->fuzzindex = 0;
      }
      n->d->FuzzSplit = FuzzSplit;
      n->SplitType = FUZZ;
      n->d->min_in_child = min_in_child;
      *chiBest = chi;
    }
  }
  
  /////////////////////////////////////////////////////////////////
  //
  // Nominal splits
  //
  for(i=Nordfuzz;i<Nvar-1;i++) {
    if (attributes_to_use && !((*attributes_to_use)[i])) continue;
    int NomSplit = 0;
    chi = FindNomSplit(scfirst,n->last, i, Nmin, n->d->MembershipTotal,&NomSplit);
    if(NomSplit>0 && (*chiBest < 0.0 || chi < *chiBest)) {
      n->att = i;
      n->NomSplit = NomSplit;
      n->SplitType = NOM;
      n->d->min_in_child = min_in_child;
      *chiBest = chi;
    }
  }


  if(n->att < 0) return 0; // no good split

  double tempmin;// Do only if there is an ordinal split. If correct, always
  if(multsplits && iSplit >=0) {
    // multiple ordinal split
    if (!n->coeff) {
      n->coeff = new double[NVarsOrdSplit];
    }
    if(n->coeff) {
      int *use = new int[NVarsOrdSplit];
      for(i=0;i<NVarsOrdSplit;i++) {
        // Initialization of the multivariate split
        // ASG: Can be altered to avoid trapping in local min
        n->coeff[i] = 0.0;
        if(i == n->att) n->coeff[i] = 1.0;
        use[i] = attributes_to_use && !((*attributes_to_use)[i]) ? 0 : 1;
      }

      // Scale independent variables
      for(i=0;i<NVarsOrdSplit;i++) {
        ScaleVariableIQ(i, scfirst,n->first,n->last);
      }
      n->fSplit = vscale[n->att]*n->fSplit+vshift[n->att];

      for(i=n->first;i<=n->last;i++) {
        m[i][Nvar] = 0.0;
        for(int j=0;j<NVarsOrdSplit;j++) {
          m[i][Nvar] += m[i][j]*n->coeff[j];
        }
      }
      iSplit =-1;

      FindOrdSplit(scfirst,n->last,Nvar, Nmin,n->d->MembershipTotal, &iSplit);
      if(iSplit >= 0) {
        int iS2 = iSplit+1;
        while(m[iS2][Nvar+3] < MIN_BELONG_LOW)
          iS2++;
        n->fSplit = 0.5*(m[iSplit][Nvar]+m[iS2][Nvar]);
        tempmin = min_in_child;
        chi = FindMultiSplit(scfirst,n->last, Nmin,n->d->MembershipTotal, use,
                                             n->coeff, &(n->fSplit), &tempmin);
        EliminateVariables(scfirst, n->last, Nmin, n->d->MembershipTotal, 
                                   n->coeff,  use, chi,&(n->fSplit), &tempmin);
        for(i=n->first;i<=n->last;i++) m[i][Nvar] -= n->fSplit;
        for(i=0;i<NVarsOrdSplit;i++) {
          n->fSplit -= n->coeff[i]*vshift[i];
          n->coeff[i]*= vscale[i];
          UnScaleVariableIQ(i, n->first,n->last);
        }
      }
      else {
        for(i=0;i<NVarsOrdSplit;i++) {
          UnScaleVariableIQ(i, n->first,n->last);
        }
        chi = *chiBest+1;
      }
      delete[] use;
    } 
    else  { // DE if(n->coeff) {
      chi = *chiBest +1;
    }
  }
  else { //DE if(multsplits && iSplit >=0)
    chi = *chiBest +1;
  }

  // Finally, the decision
  if(chi < *chiBest) {
    // Note that in this case m[i][M] == m[i]*n->coeff
    *chiBest = chi;
    n->NomSplit=0;
    n->SplitType = ORD;
    n->d->min_in_child = tempmin;
    n->att = Nvar;

    if(fuzzify) {
      chi = FuzzifySplit(n,*chiBest,scfirst,Nmin);
      if(chi>0.0 && chi<*chiBest) {
        n->SplitType = FUZZ;
        *chiBest = chi;
      }
    }
  }
  else if (n->SplitType ==ORD || (n->SplitType == FUZZ && n->d->fuzzindex)) {
    n->fSplit = oldSplit;
    if(n->coeff) {
      for(i=0;i<NVarsOrdSplit;i++) {
        n->coeff[i] = (i == n->att ? 1.0 : 0.0);
      }
    }
    for(i=n->first;i<=n->last;i++) {
      m[i][Nvar] = m[i][n->att] - (n->fSplit);
    }
  }

  n->d->CrispFlag = (n->SplitType != FUZZ) && 
                            (n->parent == NULL ? true :  n->parent->d->CrispFlag);
  if(n->d->CrispFlag) {
    if(n->SplitType == ORD) {
      // Sort
      SortOn(Nvar,n->first,n->last);
      // and find boundary
      (*leftLast) = n->first - 1;
      while((*leftLast) < n->last && m[(*leftLast)+1][Nvar] < 0.0) 
        (*leftLast)++;
    }
    else if(n->SplitType == NOM) {
      int MembersLeft = SortNomOn(n->NomSplit, n->att, n->first, n->last);
      (*leftLast) = n->first + MembersLeft-1;
    }
  }
  else {
    *leftLast = -1;
  }

  return 1;
}

int Data::OLD_FindBestSplit(Node* n, int Nmin, int* leftLast, double *chiBest)
{
  n->att = -1; // Indicates an error
  *chiBest = -1;
  n->NomSplit=0;
  n->SplitType = ERR;
  double chi;
  int iSplit, i;

  int scfirst = PreSort(n->first, n->last); // focus on data with significant membership (alpha-cut)

  // Simple ordinal splits

  /*////  ASG: Gain-ratio

  int* SplitHolder = new int[NVarsOrdSplit];
  double* ChiHolder = new double[NVarsOrdSplit];
  double* MinHolder = new double[NVarsOrdSplit];
  double* GainRatio = new double[NVarsOrdSplit];
  double AvChi =0.0;
  int denom = NVarsOrdSplit;
  // Gain-ratio ////////*/


  for(i=0;i<NVarsOrdSplit;i++)
  {
    iSplit=-1;
    chi = FindOrdSplit(scfirst, n->last, i, Nmin,n->d->MembershipTotal, &iSplit);

    /*////  ASG: Gain-ratio
    *chiBest = -1;
    SplitHolder[i] = iSplit;
    if(iSplit<0) --denom;
    // Gain-ratio//////*/


    if(iSplit >=0 && (*chiBest < 0.0 || chi < *chiBest))
    {
      /*////  ASG: Gain-ratio

      ChiHolder[i] = chi;
      MinHolder[i] = min_in_child;
      double p1 = min_in_child/n->d->MembershipTotal;
      double p2 = 1.0 - p1;
      AvChi += chi;

      GainRatio[i] = (n->entropy - chi)/((n->d->MembershipTotal) * (- p1*log(p1) - p2*log(p2)));

      // Gain-ratio ///////*/


      n->att = i;
      int iS2 = iSplit+1;
      while(m[iS2][Nvar+3] < MIN_BELONG_LOW || m[iSplit][i] == m[iS2][i])
        iS2++;
      n->fSplit = 0.5*(m[iSplit][i]+m[iS2][i]);
      n->NomSplit = 0;
      n->SplitType = ORD;
      n->d->min_in_child = min_in_child;
      *chiBest = chi;

      if(fuzzify)
      {
        chi = FuzzifySplit(n,*chiBest,scfirst,Nmin);
        if(chi>0.0 && chi<*chiBest)
        {
          n->SplitType = FUZZ;
          *chiBest = chi;
        }
      }
    }

  }


  /*/// ASG; Gain-ratio /////////////

  *chiBest = -1;
  AvChi /= denom;
  for(i =0; i< NVarsOrdSplit; ++i)
  {

    iSplit = SplitHolder[i];
    if(iSplit>=0 && ChiHolder[i] <= AvChi && (GainRatio[i] > *chiBest))
    {
      SortOn(i,scfirst,n->last);
      n->att = i;
      int iS2 = iSplit+1;
      while(m[iS2][Nvar+3] < MIN_BELONG_LOW || m[iSplit][i] == m[iS2][i])
        iS2++;
      n->fSplit = 0.5*(m[iSplit][i]+m[iS2][i]);
      n->NomSplit = 0;
      n->SplitType = ORD;
      n->d->min_in_child = MinHolder[i];

      *chiBest = GainRatio[i];
    }
  }

  delete[] ChiHolder;
  delete[] GainRatio;
  delete[] SplitHolder;
  delete[] MinHolder;

  // Gain-ratio ////////////////*/



  double oldSplit = n->fSplit;

  /////////////////////////////////////////////////////////////////
  //
  //    Fuzzy splits
  //



/*  for(i = Nord; i< Nordfuzz; ++i)
  {
    Property* tempprop = FuzzVars;
    while(tempprop->Name != VarNames[i])
      tempprop = tempprop->next;

    FuzzySet* FuzzSplit = tempprop->Types;
    chi = FindFuzzSplit(scfirst, n->last, i, Nmin, n->d->MembershipTotal, &FuzzSplit);

    if(chi>=0.0  && (*chiBest < 0.0 || chi < *chiBest))
    {
      n->att = i;
      if(n->d->fuzzindex)
      {
        delete n->d->FuzzSplit;
        n->d->fuzzindex = 0;
      }
      n->d->FuzzSplit = FuzzSplit;
      n->SplitType = FUZZ;
      n->d->min_in_child = min_in_child;
      *chiBest = chi;

    }
  }
  */
  /////////////////////////////////////////////////////////////////
  //
  // Nominal splits
  //


  for(i=Nordfuzz;i<Nvar-1;i++)
  {
    int NomSplit = 0;
    chi = FindNomSplit(scfirst,n->last, i, Nmin, n->d->MembershipTotal,&NomSplit);

    if(NomSplit>0 && (*chiBest < 0.0 || chi < *chiBest))
    {
      n->att = i;
      n->NomSplit = NomSplit;
      n->SplitType = NOM;
      n->d->min_in_child = min_in_child;
      *chiBest = chi;
    }
  }


  if(n->att < 0) return 0; // no good split

  double tempmin;
  if(multsplits && iSplit >=0) // Do only if there is an ordinal split. If correct, always
  {

    // multiple ordinal split

    if(n->coeff)
    {
      int *use = new int[NVarsOrdSplit];
      for(i=0;i<NVarsOrdSplit;i++)
      {
        // Initialization of the multivariate split
        // ASG: Can be altered to avoid trapping in local min
        n->coeff[i] = 0.0;
        if(i == n->att) n->coeff[i] = 1.0;
        use[i] = 1;
      }


      // Scale independent variables
      for(i=0;i<NVarsOrdSplit;i++)
        ScaleVariableIQ(i, scfirst,n->first,n->last);
      n->fSplit = vscale[n->att]*n->fSplit+vshift[n->att];


      for(i=n->first;i<=n->last;i++)
      {
        m[i][Nvar] = 0.0;
        for(int j=0;j<NVarsOrdSplit;j++)
          m[i][Nvar] += m[i][j]*n->coeff[j];
      }
      iSplit =-1;

      FindOrdSplit(scfirst,n->last,Nvar, Nmin,n->d->MembershipTotal, &iSplit);
      if(iSplit >= 0)
      {
        int iS2 = iSplit+1;
        while(m[iS2][Nvar+3] < MIN_BELONG_LOW)
          iS2++;
        n->fSplit = 0.5*(m[iSplit][Nvar]+m[iS2][Nvar]);
        tempmin = min_in_child;
        chi = FindMultiSplit(scfirst,n->last, Nmin,n->d->MembershipTotal, use, n->coeff, &(n->fSplit), &tempmin);
        EliminateVariables(scfirst, n->last, Nmin, n->d->MembershipTotal, n->coeff,  use, chi,&(n->fSplit), &tempmin);
        for(i=n->first;i<=n->last;i++) m[i][Nvar] -= n->fSplit;
        for(i=0;i<NVarsOrdSplit;i++)
        {
          n->fSplit -= n->coeff[i]*vshift[i];
          n->coeff[i]*= vscale[i];
          UnScaleVariableIQ(i, n->first,n->last);
        }
      }
      else
      {
        for(i=0;i<NVarsOrdSplit;i++)
          UnScaleVariableIQ(i, n->first,n->last);
        chi = *chiBest+1;
      }
      delete[] use;
    } else  chi = *chiBest +1;
  }
  else chi = *chiBest +1;


  // Finally, the decision

  if(chi < *chiBest)
  {               // Note that in this case m[i][M] == m[i]*n->coeff
    *chiBest = chi;
    n->NomSplit=0;
    n->SplitType = ORD;
    n->d->min_in_child = tempmin;
    n->att = Nvar;

    if(fuzzify)
    {
      chi = FuzzifySplit(n,*chiBest,scfirst,Nmin);
      if(chi>0.0 && chi<*chiBest)
      {
        n->SplitType = FUZZ;
        *chiBest = chi;
      }
    }
  }
  else if (n->SplitType ==ORD || (n->SplitType == FUZZ && n->d->fuzzindex))
  {
    n->fSplit = oldSplit;
    if(n->coeff) for(i=0;i<NVarsOrdSplit;i++) n->coeff[i] = (i == n->att ? 1.0 : 0.0);
    for(i=n->first;i<=n->last;i++) m[i][Nvar] = m[i][n->att] - (n->fSplit);
  }


  n->d->CrispFlag = (n->SplitType != FUZZ) && (n->parent == NULL ? true :  n->parent->d->CrispFlag);
//  n->d->CrispFlag = false;
  if(n->d->CrispFlag)
  {

    if(n->SplitType == ORD)
    {
      // Sort
      SortOn(Nvar,n->first,n->last);
      // and find boundary
      (*leftLast) = n->first - 1;
      while((*leftLast) < n->last && m[(*leftLast)+1][Nvar] < 0.0) (*leftLast)++;
    }
    else if(n->SplitType == NOM)
    {
      int MembersLeft = SortNomOn(n->NomSplit, n->att, n->first, n->last);
      (*leftLast) = n->first + MembersLeft-1;
    }
  }
  else (*leftLast = -1);
  return 1;
}

double Data::FindFuzzSplit(int first, int last, int iX, int Nmin, double NMembers, FuzzySet** pfz)
{
/*  FuzzySet* tempfz = *pfz;
  *pfz = NULL;
  double chi = -1.0;

  if(Nmin < 1 || NMembers< 2*Nmin ) return chi; // Bad Nmin or Too few data or variable not nominal



  double* membership = new double[last-first+1];

  int i;
  for(i = first; i <= last; ++i)
    membership[i-first] = m[i][Nvar+3];



  while(tempfz)
  {
    int i;
    for(i= first; i <=last; ++i)
    {
      m[i][Nvar+3] = membership[i-first]*tempfz->Membership((m[i][iX] - Gvshift[iX])/Gvscale[iX]);
    }
    double min_in_node = Initialise(first,last,last,100,-1,true);

    for(i = first; i <= last; ++i)
      m[i][Nvar+3] = membership[i-first]*(1.0 - tempfz->Membership((m[i][iX] - Gvshift[iX])/Gvscale[iX]));

    min_in_node = min(Initialise(first,last,first-1,100,-1,true), min_in_node);

    if(min_in_node >= Nmin  && (chi<0.0 || (R[0]+R[1] < chi)))
    {
      chi = R[0]+R[1];
      *pfz = tempfz;
      min_in_child = min_in_node;
    }

    tempfz = tempfz->Next();

  }

  for(i = first; i <= last; ++i)
    m[i][Nvar+3] = membership[i-first];

  delete[] membership;

  return chi;*/                   return 0.0;
}


double Data::FindNomSplit(int first, int last, int iX, int Nmin, double NMembers, int* pNomSplit)
{
  *pNomSplit = -1; // there is an error
  if(NMembers< 2*Nmin || iX < Nordfuzz) return 0.0; // Bad Nmin or Too few data or variable not nominal
//  int dim = NomTerms[iX-Nordfuzz].GetSize();
  int dim = NomTerms[iX-Nordfuzz].size();

  double chi =-1.0;
  
  if(DepVarType == NOM && NomTerms[Nnom].size() == 2)
  {
    double** Pop = new double*[dim];
    
    double* Norm = new double[dim];
    int i;
    for(i =0; i < dim; ++i)
    {
      Pop [i] = new double[2];
      Pop[i][0] = i + 0.5;
      Pop[i][1] = 0.0;
      Norm[i] =0.0;
    }

    for (i = first; i<=last; ++i)
    {    
      if((int)m[i][Nvar-1] == 1)
        Pop[(int)m[i][iX]][1] += m[i][Nvar+3];
      Norm[(int)m[i][iX]] += m[i][Nvar+3];
    }
    


    for (i = 0; i< dim; ++i)
    {
      if(Norm[i] >= 1e-10)
        Pop[i][1] /= Norm[i];
      else Pop[i][1] = 1.0; 
    }

    theSortAtt = 1;
    qsort(Pop,dim,sizeof(double**),compAtt);
    
    int k = 0;
    for (i = 0; i< dim ; ++i)
    {
    
      k += (1 << (int)Pop[i][0]);

      R[0] = R[1] = 0.0;          
      int count = SortNomOn(k,iX, first,last);
      double min_in_node =Initialise(first,last, first+count-1);
      if(min_in_node >= Nmin  && (chi<0.0 || ((R[1]+R[0]) <= chi)))
      {
        chi = R[0]+R[1];
        *pNomSplit= k;
        min_in_child = min_in_node;
      }

    }

    for(i=0; i<dim; ++i)
      delete[] Pop[i];
  
    delete[] Norm;
    delete[] Pop;
  }
  else
  {

    for (int i=1; i < pow(2.0,dim-1); ++i)
    {
    //  membership of sets = i in binary form
      R[0] = R[1] = 0.0;          
      int count = SortNomOn(i,iX, first,last);
      double min_in_node =Initialise(first,last, first+count-1);
      if(min_in_node >= Nmin  && (chi<0.0 || ((R[1]+R[0]) <= chi)))
      {
        chi = R[0]+R[1];
        *pNomSplit= i;
        min_in_child = min_in_node;
      }
    }
  }
  
  return chi;
}

int Data::SortNomOn(int Membership, int iX,int first, int last)
{

  int count =0;
  for(int i = first; i<=last; ++i)
  {
    if((1<<(int)m[i][iX]) & Membership)
    {
      m[i][Nvar] = 0;
      ++count;
    } else {
      m[i][Nvar] = 1;
    }
  }
  theSortAtt = Nvar;
  qsort(m+first,last-first+1,sizeof(double**),compAtt);
  return count;
}



double Data::FindOrdSplit(int first, int last, int iX, int Nmin,double NMembers, int *iSplit, double gamma)
{
  *iSplit = -1; // Indicates an error
  if(NMembers < 2*Nmin) return 0.0; // Bad Nmin or Too few data or variable not ordinal
                                          
  int iVar = iX;
  if(gamma < 1) iVar = Nvar+1;

  SortOn(iVar,first,last);
    // We need a split between different values
  int curSplit = last;
  if(gamma > 1)
  {                                  // For single ordinal split
    curSplit--;//GMM
    while(curSplit >= (first) && fabs(m[curSplit][iVar] - m[curSplit+1][iVar]) < 1e-10 &&  NMembers >= Nmin)
    {
      curSplit--;
      NMembers -= m[curSplit+1][Nvar+3];
    }
//    if(curSplit == (first-1) || NMembers < Nmin) return 0.0;  // Error: no good split
  } else curSplit = last;
    // Begin main algorithm
  
  R[0] = 0.0, R[1] = 0.0;
  double chi= 0.0;
  double min_in_node = Initialise(first, last, curSplit,gamma, iX);  // ASG: ask James why iX and not iVar?
      
  while(curSplit >= first)
  {
    if((curSplit < last)
      && (fabs(m[curSplit][iVar] - m[curSplit+1][iVar]) >= 1e-10) // check for new best
      && (min_in_node >= Nmin)
      && ((R[0]+R[1] < chi) || *iSplit <0)
      && (m[curSplit][Nvar+3]>= MIN_BELONG_LOW))
    {
      *iSplit = curSplit; 
      chi = R[0]+R[1]; 
      min_in_child = min_in_node;
    }
    
    // Move current point from left to right
    if(m[curSplit][Nvar+3] > 0.0) min_in_node = MovePoint(curSplit,gamma,iX,first,last);
    curSplit--;
  }
  
  return chi; 
}

double Data::FindMultiSplit(int first, int last, int Nmin, double NMembers, 
                                   int *use, double* a,  double* c, double* min)
{

  double chiLoop = -1, chiPrev = -1;
  int i;
  for(i=first;i<=last;i++)
  {
    m[i][Nvar] = 0.0;  
    for(int j=0;j<NVarsOrdSplit;j++) m[i][Nvar] += m[i][j]*a[j];
  }

  do
  {
    for(int iX=0;iX<NVarsOrdSplit;iX++) // for each variable
    {
      double chi, chiBest = -1;
      double gammaBest=-1000, deltaBest=-1000;
      int iSplit;
      if(!use[iX]) continue;
      for(double gamma = -0.25001; gamma < 0.3; gamma += 0.25)
      {
        for(i=first;i<=last;i++)
          m[i][Nvar+1] = (m[i][Nvar]-(*c))/(m[i][iX]+gamma);
        chi = FindOrdSplit(first,last,iX,Nmin, NMembers, &iSplit, gamma);

        if(iSplit>=0 &&(chi < chiBest || chiBest < 0.0))
        {
          chiBest = chi;
        
          int iS2 = (iSplit == last) ? iSplit : (iSplit+1);
          // ASG: iSplit can be the last point in the node. 
          // In that case, it does not make sense to shift split
        /*  while(m[iS2][Nvar+3] <MIN_BELONG_LOW) 
          {
            if(iS2 == last) iS2 = iSplit;
            else ++iS2;
            }              //ASG: Unnecesary
          */
        
          

          deltaBest = 0.5*(m[iSplit][Nvar+1]+m[iS2][Nvar+1]);
          
          gammaBest = gamma;
        }
      } // end gamma-loop
      a[iX] -= deltaBest;
      *c += deltaBest*gammaBest;
      for(i=first;i<=last;i++) m[i][Nvar] -= m[i][iX]*deltaBest;
    } // end iX loop
    int iSplit;
    chiPrev = chiLoop;
    // Normalise
    double sum = 0.0;
    for(i=0;i<NVarsOrdSplit;i++) if(fabs(a[i]) > sum) sum = fabs(a[i]);
    for(i=0;i<NVarsOrdSplit;i++) a[i] /= sum;
    for(i=first;i<=last;i++) m[i][Nvar] /= sum;
    *c /= sum;
    chiLoop = FindOrdSplit(first,last,Nvar, Nmin  ,NMembers, &iSplit); // ASG: if problems change to Nmin - 1 avoids problems related to not finding a split
    if(iSplit >=0 /*&& (*chiBest < 0.0 || chi < *chiBest)*/)
    {
      int iS2 = iSplit+1;
      while(m[iS2][Nvar+3] < MIN_BELONG_LOW )
        iS2++;
      *c = 0.5*(m[iSplit][Nvar]+m[iS2][Nvar]);
      *min = min_in_child;
    }
    else 
    {
      chiLoop = -1; // Normally Never
    }
  } while (chiLoop > -0.1 && (chiPrev < -0.1 || chiLoop < chiPrev) && (chiPrev < -0.1 || fabs(chiLoop-chiPrev) > cartEps));
  return chiLoop;
}


double Data::EliminateVariables(int first, int last, int Nmin, double NMembers, double* a,
      int* use, double chiOld, double* fSplit, double* min)
{
  double I = Impurity(first,last,NMembers);
  double* old_a = new double[NVarsOrdSplit];
  
  int i;
  for(i = 0; i<NVarsOrdSplit; ++i)
  {
    old_a[i] = a[i];
  }
  
  double arem;
  
  for(i=first;i<=last;i++)
  {
    m[i][Nvar] = 0.0;
    for(int j=0;j<NVarsOrdSplit;j++)
      m[i][Nvar] += a[j]*m[i][j];
  }

  double delta = I-chiOld;
  double deltaMax=0,deltaMin=0;
  double chi,chimax;
  int imax=-1000;
  int iSplit;

  int stop=0;
  do
  {
    stop = 1;
    int flag = 0;
    for(int iX=0;iX<NVarsOrdSplit;iX++)
    {
      if(use[iX] == 0) continue;
      arem = a[iX]; a[iX] = 0.0;

      for(int i=first;i<=last;i++) m[i][Nvar] -= arem*m[i][iX];
      chi = FindOrdSplit(first,last,Nvar, Nmin, NMembers, &iSplit);

    //  if(iSplit >= 0) // To Review ASG
      {
        double deltaTemp = I-chi;
        if(flag == 0 || deltaTemp < deltaMin) {deltaMin = deltaTemp;}
        if(flag == 0 || deltaTemp > deltaMax) {deltaMax = deltaTemp;imax = iX;}
      }
      for(int i=first;i<=last;i++) m[i][Nvar] += arem*m[i][iX];
      a[iX] = arem; arem = 0.0; flag++;
    }
    if(flag>1)
      if(delta-deltaMax < cartBeta*(delta-deltaMin))
      {
        stop = 0;
        for(int i=first;i<=last;i++) m[i][Nvar] -= a[imax]*m[i][imax];
        a[imax] = 0.0;
        use[imax] = 0;
        chimax = I-deltaMax;
      } else stop = 1;
  } while(!stop);
  chimax = FindOrdSplit(first,last,Nvar, Nmin,NMembers,&iSplit);
  if(iSplit < 0) 
  {
    for(i = 0; i<NVarsOrdSplit; ++i)
    {
      a[i] = old_a[i];
    }
    delete old_a;
    return -1;   // ASG: Recover old coeffs
  }

  delete old_a;
  int iS2 = iSplit+1;
  while(m[iS2][Nvar+3] < MIN_BELONG_LOW)
    iS2++;
  *fSplit = 0.5*(m[iSplit][Nvar]+m[iS2][Nvar]);
  *min = min_in_child;


  return chimax;
}

double Data::Score(int begin, int end,Node* root, int K)
{
  //  console->printf("Scoring not implemented");
  return 0.0;
}
double Data::FindMembership(const Node* n)
{
  const Node* cur = n;

  int i;
  for(i= n->first; i<= n->last; ++i)
    m[i][Nvar+3] = 1.0*m[i][Nvar+4];



  while(cur->parent)
  {
    for(int i= n->first; i<= n->last; ++i)
    {
      switch(cur->parent->SplitType)
      {
      case NOM:

        if((1<<(int)m[i][cur->parent->att]) & cur->parent->NomSplit)
        {
          if(!cur->sib)
          {
            m[i][Nvar+3] = 0.0;
          }
        }
        else if(cur->sib)
        {
          m[i][Nvar+3] = 0.0;
        }
        break;
      case ORD:
        m[i][Nvar] = -cur->parent->fSplit;

        if (cur->parent->coeff && !cur->parent->IsUnivariateSplit()) // !cur->parent->univariate_split 
        {
          for(int j=0;j<NVarsOrdSplit;j++)
            m[i][Nvar] += cur->parent->coeff[j] * m[i][j];
if ((m[i][cur->parent->att]-cur->parent->fSplit)*m[i][Nvar]<=0.0) 
printf("MS<-------------%g--------%g-------------------------\n", m[i][Nvar], m[i][cur->parent->att]-cur->parent->fSplit);
        } else
{
          m[i][Nvar] += m[i][cur->parent->att];
/*if (m[i][Nvar]==0.0) 
printf("MS<-------------%g--------%g-------------------------\n", m[i][Nvar], m[i][cur->parent->att]-cur->parent->fSplit);*/
}


        if(m[i][Nvar] <=0.0) // ASG to Verify
        {
          if(!cur->sib)  // R Node
          {
            m[i][Nvar+3] = 0.0;
          }
        }
        else if(cur->sib)  // L Node
        {
          m[i][Nvar+3] = 0.0;
        }
        break;

      case FUZZ:

  /*      int iX = cur->parent->att;
        if(iX==Nvar)
        {
          m[i][Nvar] = -cur->parent->fSplit;
          for(int j=0;j<NVarsOrdSplit;j++)
            m[i][Nvar] += cur->parent->coeff[j]*m[i][j];
        }

        if(cur->sib)
        {
          m[i][Nvar+3] *= cur->parent->d->FuzzSplit->Membership((m[i][iX]- Gvshift[iX])/Gvscale[iX]);
        }
        else
        {
          m[i][Nvar+3] *= (1.0 - cur->parent->d->FuzzSplit->Membership((m[i][iX]- Gvshift[iX])/Gvscale[iX]));
        }
    */
        break;
        default:
        break;

      }
    }
    cur = cur->parent;
  }

  double total =0.0;

  for(i= n->first; i<= n->last; ++i)
    total += m[i][Nvar+3];

  return total;
}

void Data::FindMembershipTest(const Node* n)
{
  const Node* cur = n;

  for(int i= n->d->firstTest; i<= n->d->lastTest; ++i)
    m[i][Nvar+3] = 1.0*m[i][Nvar+4];

  
  while(cur->parent)
  {  
    for(int i= n->d->firstTest; i<= n->d->lastTest; ++i)
    {    
      switch(cur->parent->SplitType)
      {
      case NOM:
        
        if((1<<(int)m[i][cur->parent->att]) & cur->parent->NomSplit)
        {  
          if(!cur->sib)
            m[i][Nvar+3] = 0.0;
        }
        else if(cur->sib)
          m[i][Nvar+3] = 0.0;
        break;

      case ORD:
        m[i][Nvar] = -cur->parent->fSplit;
        
        if(cur->parent->coeff)
        {
          for(int j=0;j<NVarsOrdSplit;j++) 
            m[i][Nvar] += cur->parent->coeff[j] * m[i][j]; 
        } else   
          m[i][Nvar] += m[i][cur->parent->att];
      
        
                
        if(m[i][Nvar] <=0.0) // ASG to Verify 
        {
          if(!cur->sib)  // R Node
            m[i][Nvar+3] = 0.0;
        }
        else if(cur->sib)  // L Node
          m[i][Nvar+3] = 0.0;
        break;
        
      case FUZZ:
/*        int iX = cur->parent->att;

        if(iX == Nvar)
        {
          m[i][Nvar] = -cur->parent->fSplit;
          for(int j=0;j<NVarsOrdSplit;j++)
            m[i][Nvar] += cur->parent->coeff[j]*m[i][j];
        }

        if(cur->sib)
        {
          m[i][Nvar+3] *= cur->parent->d->FuzzSplit->Membership((m[i][iX]- Gvshift[iX])/Gvscale[iX]);
        }
        else
        {
          m[i][Nvar+3] *= (1.0 - cur->parent->d->FuzzSplit->Membership((m[i][iX]- Gvshift[iX])/Gvscale[iX]));
        }*/
        break;
        default:
        break;
      }  
    }
    cur = cur->parent;
  }  
}


double Data::FuzzifySplit(Node* n, double chiBest, int first, int Nmin)
{  
  
  double split = n->att == Nvar ?  0.0 : (n->fSplit-Gvshift[n->att])/Gvscale[n->att];


  double d1,d2;

  d1 = d2 = split;
  
  
/*  if(n->att != Nvar)
  {
    Property* tempp = FuzzVars;
    while(VarNames[n->att] != tempp->Name)
      tempp = tempp->next;

    d1 = tempp->Domain[0];
    d2 = tempp->Domain[1];
  }
  else */
  {
    for(int i = n->first; i <= n->last; ++i)
    {
      if(m[i][Nvar+3] > MIN_BELONG_HIGH)
      {
        double temp;
        if(n->att == Nvar)
        {
          m[i][Nvar] = -n->fSplit;
          for(int j=0;j<NVarsOrdSplit;j++)
            m[i][Nvar] += n->coeff[j]*m[i][j];
          temp = m[i][Nvar];
        }
        else temp = (m[i][n->att] - Gvshift[n->att])/Gvscale[n->att];

        if(i == n->first || temp > d2)
          d2 = temp;
        if(i==n->first || temp < d1)
          d1 = temp;
      }
    }

  }

/*  FuzzySet* newfz = new FuzzySet("new",SHOULDER_DEC);
  newfz->SetDomain(d1,d2);

  double interval = min(split-d1,d2-split)/10.0;

  for(int i=1; i<=10; ++i)
  {
    newfz->Generate(split-i*interval,split, split+i*interval);
    FuzzySet* fzsplit = newfz;
    double chi = FindFuzzSplit(first, n->last, n->att,  Nmin, n->d->MembershipTotal, &fzsplit);
    chi = 2.0*chi/2.0;
    if(chi>=0.0 && chi<chiBest)
    {
      chiBest = chi;
      if(n->d->fuzzindex) delete n->d->FuzzSplit;
      n->d->FuzzSplit = new FuzzySet(fzsplit); // ASG Danger: new objects that are not destroyed
      n->d->fuzzindex = i;
      n->d->min_in_child = min_in_child;
    }
  }
  delete newfz;*/
  return chiBest;
}

void Data::GetScale(int first, int last, double& xmin, double& xmax, double* coeff, int att)
{
  for(int i=first;i<=last;i++)
  {
    double val=0.0;

    if(coeff)
      for(int j=0;j<NVarsOrdSplit;j++)
          val += coeff[j] * m[i][j];
    else val = m[i][att];
    if(val > xmax || i == first) xmax = val;
    if(val < xmin || i == first) xmin = val;
  }

}


double Data::DegreeOfFuzz(Node* root, int K)
{
  double total =0.0;
  double fuzz =0.0;
  int begin =0;    
  int end = NTrain-1;
  for(int i= begin;i<= end ;i++)
  {
    // Going up : fix memberships && partial predictions
    Node* cur = root;
  
    root->d->cumMU =1.0;
    root->d->tempMU = 1.0;
      
    while(cur)
    {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {
        double aux = cur->d->tempMU*m[i][Nvar+4];
        if(cur->parent) aux *= cur->parent->d->cumMU;
        
        total += aux;
        fuzz += aux*(1.0-aux);

      } else {

        double mu =  cur->ParamMem(m[i],NVarsOrdSplit);
        cur->child->d->tempMU = mu;
        cur->child->d->cumMU = cur->d->cumMU*mu;
        cur->child->sib->d->tempMU = (1-mu);
        cur->child->sib->d->cumMU = cur->d->cumMU*(1-mu);
      }

      if(!curIsLeaf) cur = cur->child;
      else {
        while(cur && !cur->sib) cur = cur->parent;
        if(cur) cur = cur->sib;
      }
    }
  }
  return fuzz/total;
}
//Set the weights to 1.0 for the data starting at begin to end
void Data::ResetWeights(int begin, int end)
{
  begin = begin < 0 ? 0        : begin;
  end   =   end < 0 ? NTotal-1 : end;
 
  for (int i = begin; i <= end; i++)
    m[i][Nvar+4] = 1.0;
}
//Returns the last element from begin to end that fulfills the conditions of
//node n. If the node is a leaf then -1 is returned.
int Data::GetLeftLast(Node * n, int begin, int end, int K)
{
  int i, leftLast=-1;

  if (!n->IsLeaf(K)) {
    int holdFirst = n->child->first;
    int holdLast  = n->child->last;
    
    n->child->first = begin;
    n->child->last = end;
//    double acum=0;
/*    double mms = ???*/FindMembership(n->child);
    n->child->first = holdFirst;
    n->child->last = holdLast;
    SortOn(Nvar+3, begin, end, true);
  /*  int j
    for(j=begin;j<=end;j++) {
      acum += m[j][Nvar+4];
      if (acum>mms) {break;}
    }
    int leftLast2 = j - 1;*/
    for(i=begin;i<=end;i++) {
      if (m[i][Nvar+3]==0.0) {break;}
    }
    leftLast = i - 1;
/*    if (leftLast!=leftLast2){
  char nf[100] = "d:\\temp.txt";
  FILE *f = fopen(nf, "at");
  itoa(leftLast, nf, 10);
  strcat(nf,"\n");
  fwrite(nf,1,strlen(nf),f);
  fclose(f);
  char nf2[100] = "d:\\temp2.txt";
  f = fopen(nf2, "at");
  itoa(leftLast2, nf2, 10);
  strcat(nf2,"\n");
  fwrite(nf2,1,strlen(nf2),f);
  fclose(f);                 }*/
  }
/*  else if (n->SplitType == ORD)
  {
    for(i=begin;i<=end;i++) m[i][Nvar] = m[i][n->att] - (n->fSplit);
    // Sort
    SortOn(Nvar, begin, end);
    // and find boundary
    leftLast = begin-1;
    while(leftLast < end && m[leftLast+1][Nvar] <= 0.0) leftLast++;
  }
  else if(n->SplitType == NOM)
  {
    int MembersLeft = SortNomOn(n->NomSplit, n->att, begin, end);
    leftLast = begin + MembersLeft-1;
  }*/
  return leftLast;

}

//
void Data::SelectNewTrainingSetEqDis()
{
  int j;
  unsigned i;
  int iclass;
  vector<int> nnclass(0);

  //Sort by class
  SortOn(Nvar-1, 0, NTotal-1);

  GroupsCreated = false;
  iclass=-1;
  for(i=0;i<(unsigned)NTotal;i++) {
    //A random number is thrown
    m[i][Nvar]=TDebugRand::Rand();
    m[i][Nvar+2] = (int) (1+(Ncv)*drand()); // Label for CV group
    //The number of classes and number of elements of each class are obtained
    if (iclass!=((int)m[i][Nvar-1])) {
      iclass=(int)m[i][Nvar-1];
      nnclass.resize(nnclass.size()+1);
    }
    nnclass.at(nnclass.size()-1)++;
  }

  //Each class is ordered by the random number just thrown
  int ClassBegin=0;
  for(i=0;i<nnclass.size();i++) {
    SortOn(Nvar, ClassBegin, ClassBegin+nnclass.at(i)-1);
    //The elements in the class are numbered
    for(j=0;j<nnclass[i];j++)
      m[j + ClassBegin][Nvar]=j;
    ClassBegin += nnclass.at(i);
  }

  //Sort by order 
  SortOn(Nvar, 0, NTotal-1);
}
void Data::SelectNewTrainingSet()
{
  int j;
  unsigned i;
  int iclass;
  int nelem;
  vector<int> nnclass(0);

  //Sort by class
  SortOn(Nvar-1, 0, NTotal-1);

  GroupsCreated = false;
/*  GroupsCreated = true;
  if (GroupCount) delete []GroupCount;
  GroupCount = new int[nGroups];
  for(int i=0;i<nGroups;i++) GroupCount[i]=0;
*/
//  srand(time(0) + rand());//  randomize();
  iclass=-1;
  for(i=0;i<(unsigned)NTotal;i++) {
    //A random number is thrown
    m[i][Nvar]=iclass + 1 + TDebugRand::Rand();
    m[i][Nvar+2] = (int) (1+(Ncv)*drand()); // Label for CV group
    //The number of classes and number of elements of each class are obtained
    if (iclass!=((int)m[i][Nvar-1])) {
      iclass=(int)m[i][Nvar-1];
      nnclass.resize(nnclass.size()+1);
    }
    nnclass.at(nnclass.size()-1)++;
  }

  //Each class is ordered by the random number just thrown
  int ClassBegin=0;
  float moco = (float)NTrain/NTotal;
  for(i=0;i<nnclass.size();i++) {
    SortOn(Nvar, ClassBegin, ClassBegin+nnclass.at(i)-1);
    nelem = (int)(moco*nnclass.at(i));
    //The first nelem are selected
    for(j=0;j<nelem;j++)
      m[j + ClassBegin][Nvar]=-1;
    ClassBegin += nnclass.at(i);
  }

  //Sort by group
  SortOn(Nvar, 0, NTotal-1);
}
void Data::SaveToFile(char *FileName, int ini, int fin)
{
  ofstream ofs(FileName);
  SaveToStream(ofs, ini, fin);
  ofs.close();
}
string Data::GetExampleAsString(int iEx)
{
  int iNom, iOrd;
  std::ostringstream buf;
  
  buf.precision(14);
  iNom = iOrd = 0;
  for(int i=0; i<Nvar; i++) {
    if (OriginalVarType[i]==ORD) {
      if (m[iEx][iOrd]==DBL_MAX) buf << "?\t";
      else   buf << m[iEx][iOrd] << "\t";
      iOrd++;
    }
    else if (OriginalVarType[i]==NOM) {
      if (NomTerms[iNom].size() == 0) {
        buf << (int)m[iEx][Nordfuzz + iNom] << "\t";
      }
      else {
        buf << NomTerms[iNom][(int)m[iEx][Nordfuzz + iNom]] << "\t";
      }
      iNom++;
    }
  }

  return buf.str();
}
string Data::GetNamesAsString()
{
  std::ostringstream buf;
  int iNom, iOrd;
  
  buf << Nvar;
  iNom = iOrd = 0;
  for (int i = 0; i < Nvar; i++) {
    buf << "\n";
    if (OriginalVarType[i]==ORD)
      buf << VarNames[iOrd++] << "\t0";
    else if (OriginalVarType[i]==NOM) {
      buf << VarNames[Nordfuzz + iNom] << "\t2\t";
      for (unsigned iv = 0; iv < NomTerms[iNom].size(); iv++) {
        buf << NomTerms[iNom][iv] << " ";
      }
      iNom++;
    }
  }
  return buf.str();
}
void Data::SaveToStream(std::ostream &out, int ini, int fin)
{
  int i;
  int ntotal = fin<0 ? NTotal : fin - ini + 1;
  ini = ini<0 ? 0 : ini;
  out << ntotal << "\n";
  out << GetNamesAsString() << "\n";
  for(i=ini;i<ini+ntotal;i++)
    out << GetExampleAsString(i) << "\n";
  
/*  out << Nvar << "\t" << ntotal << "\n";
  //falta VarNames ????????????????????
  for(i=0;i<Nord;i++) out << "0\t";
  for(i=0;i<Nfuzz;i++) out << "1\t";
  for(i=0;i<Nnom;i++) out << "2\t";
  out << ((DepVarType==ORD) ? "0\n" : (DepVarType==FUZZ) ? "1\n" : "2\n");
  for(i=0;i<Nvar;i++) out << VarNames[i] << "\t";
  out << "\n";
  for(i=ini;i<ini+ntotal;i++)  
  {
    int j;
    //  ordinal and fuzzy
    for(j = 0; j < Nordfuzz; ++j)
      out << m[i][j] << "\t";
    //Nominal
    for(; j < Nvar-1; ++j)
      out << NomTerms[j-Nordfuzz][(int)m[i][j]] << "\t";

    //Depvar
    if (DepVarType==NOM)
      out << NomTerms[Nnom][(int)m[i][j]] << "\n";
    else //ordinal
      out << m[i][j] << "\n";
  }*/
}
void Data::SaveC45Data(char *fname, int ini, int fin)
{
  //
  //Se abre el fichero .data
  //
  int j;

  if (fin-ini+1<=0) return;

  FILE *f = fopen(fname, "wt");
  if (f==NULL) return;
  for(int i=ini;i<=fin;i++) {
    //  ordinal and fuzzy
    for(j=0;j<Nordfuzz;j++)
      if (m[i][j]==DBL_MAX) fprintf(f, "?,");
      else fprintf(f, "%f," , m[i][j]);

    //Nominal
    for(;j<Nvar-1;j++)
      fprintf(f, "%s,", NomTerms[j-Nordfuzz][(int)m[i][j]].c_str());

    //Depvar
    if (DepVarType==NOM)
      fprintf(f, "%s\n", NomTerms[Nnom][(int)m[i][j]].c_str());
    else //ordinal
      fprintf(f, "%f\n", m[i][j]);
  }
  fclose(f);

}
void Data::SaveAsC45(char *fichero)
{
  unsigned i, j;
  int lon = strlen(fichero);
  char fnames[512];//MAXPATH];

  //Se añade la extension .names si no la trae ya puesta
  strcpy(fnames, fichero);
  if (lon<7 || strcmpci(fichero+lon-6, ".names")!=0)
    strcat(fnames, ".names");

  //
  //Se abre el fichero .names
  //
  FILE *f = fopen(fnames, "wt");
  if (f==NULL) return;

  fprintf(f, "|%s\n", fnames);

  //Se imprimen las posibles clases
  for (i=0;i<NomTerms[Nnom].size()-1;i++) {
    fprintf(f, "%s, ", NomTerms[Nnom][i].c_str());
  }
  fprintf(f, "%s.\t\t\t\t|classes\n\n", NomTerms[Nnom][i].c_str());

  //  ordinal y borrosas
  for(i=0;i<(unsigned)Nordfuzz;i++)
    fprintf(f, "%s:\t\t\t\tcontinuous.\n", VarNames[i].c_str());

  //Se imprime el tipo de los de datos
  for(;i<(unsigned)Nvar-1;i++) {
    fprintf(f, "%s:   ", VarNames[i].c_str());
    for (j=0;j<NomTerms[i-Nordfuzz].size()-1;j++) {
      fprintf(f, "%s, ", NomTerms[i-Nordfuzz][j].c_str());
    }
    fprintf(f, "%s\n", NomTerms[i-Nordfuzz][j].c_str());
  }

  fclose(f);

  //
  // Guardo el fichero .data
  //
  fnames[strlen(fnames)-6] = '\0';//Se quita la extension .names
  strcat(fnames, ".data");
  SaveC45Data(fnames, 0, NTrain-1);
  //
  // Guardo el fichero .test
  //
  fnames[strlen(fnames)-5] = '\0';//Se quita la extension .names
  strcat(fnames, ".test");
  SaveC45Data(fnames, NTrain, NTotal-1);
}
//------------------------------------------------------  DeleteColumn  -----
//---------------------------------------------------------------------------
void Data::DeleteColumn(int icol)
{
  for(int i = 0; i < NTotal; i++) {
    memmove(m[i]+icol, m[i]+(icol+1), sizeof(double) * (7 + Nvar - icol - 1));
  }

  string *vn = new string[Nvar - 1];
  for(int i = 0; i < icol; i++) vn[i] = VarNames[i];
  for(int i = icol + 1; i < Nvar; i++) vn[i - 1] = VarNames[i];
  delete []VarNames;
  VarNames = vn;

  memmove(OriginalVarType+icol, OriginalVarType+(icol+1), 
                                           sizeof(VarType) * (Nvar - icol - 1));

  if (icol < Nord) {
    Nord--;
    Nordfuzz--;
    NVarsOrdSplit--;
  }
  else if (icol < Nordfuzz) {
    Nordfuzz--;
    Nfuzz--;
  }
  else if (icol < Nvar) {
    vector<string> *nt = new vector<string>[Nvar - 1];
    for(int i = 0; i < icol; i++) nt[i] = NomTerms[i];
    for(int i = icol + 1; i < Nvar; i++) nt[i - 1] = NomTerms[i];
    delete []NomTerms;
    NomTerms = nt;
    Nnom--;
  }

  Nvar--;

  return; 
}
//-----------------------------------------------------------  AddData  -----
//---------------------------------------------------------------------------
void Data::AddData(Data *d2, int ini, int fin)
{
  ini = ini < 0 ? 0 : ini;
  fin = fin < 0 ? d2->GetNTotal() - 1 : fin;

  d2->MarkOrder();
  d2->SortOn(Nvar+5, ini, fin);

  //Copiamos los datos
  int num_dif=fin-ini+1;
  int lastNTotal = NTotal;
  redim(NTotal+num_dif);
  for(int i=ini; i<=fin;i++) {
    for(int j=0; j<Nvar+5;j++) { // dont copy the original order as it makes no sense
      m[ lastNTotal + i - ini ][ j ] = d2->m[ i ][ j ];
    }
    m[ lastNTotal + i - ini ][ Nvar+6 ] = d2->m[ i ][ Nvar+6 ];
  }
  d2->ResetOrder();
  
}
//-----------------------------------------------------------  AddData  -----
//---------------------------------------------------------------------------
/*void Data::AddFeature(Data *d2, int ifeature)
{
  ini = ini < 0 ? 0 : ini;
  fin = fin < 0 ? d2->GetNTotal() - 1 : fin;

  //Copiamos los datos
  int num_dif=fin-ini+1;
  int lastNTotal = NTotal;
  redim(NTotal+num_dif);
  for(int i=ini; i<=fin;i++) {
    for(int j=0; j<Nvar+7;j++) {
      m[ lastNTotal + i - ini ][ j ] = d2->m[ i ][ j ];
    }
  }
  
}*/
//-------------------------------------------------------------  Clone  -----
//---------------------------------------------------------------------------
Data* Data::Clone(int ini, int fin, int factor)
{
  ini = ini<0 ? 0 : ini;
  fin = fin<0 ? GetNTotal()-1 : fin;

  Data *clon = NewData();

  // Cargamos los nombres
  FILE *f=tmpfile();
  fprintf(f, GetNamesAsString().c_str());
  rewind(f);
  clon->GetNames(f);
  fclose(f);

  //Copiamos los datos
  clon->init();
  int num_dif=fin-ini+1;
  clon->redim(num_dif*factor);
  for(int i=ini; i<=fin;i++) {
    for(int j=0; j<Nvar+7;j++) {
      for(int k=0; k<factor;k++) {
        clon->m[k*num_dif+i-ini][j] = m[i][j];
      }
    }
  }
  
  return clon;
}
//------------------------------------------------  GetBootstrapSample  -----
//---------------------------------------------------------------------------
Data* Data::GetBootstrapSample(int N, bool Weighted, bool IncludeOobAsTest, std::vector<int> *idxs)
{
  Data *bootstrap = NewData();

  //Por omision se extraen igual numero de ejemplos que datos de entrenamiento
  if (N<=0) N=NTrain;

  // Cargamos los nombres y redimensionamos los datos
  FILE *f=tmpfile();
  fprintf(f, GetNamesAsString().c_str());
  rewind(f);
  bootstrap->GetNames(f);
  fclose(f);
  bootstrap->init();

  if ( IncludeOobAsTest && !idxs ) {
    bootstrap->redim(N+NTrain);
  }
  else {
    bootstrap->redim(N);
  }

  //Ordeno por orden original
  //SortOn(Nvar+5);

  //Lo utilizo para marcar los datos que entran en la muestra de bootstrap
  bool *incluido = new bool[NTrain];
  double TotWeight = 0.0;
  for(int i=0; i<NTrain;i++) {
    incluido[i] = false;
    TotWeight += GetDatWeight(i); 
  }

  //Obtengo N elementos al azar con reposicion (muestra bootstrap)
  for(int i=0; i<N;i++) {
    int elem;
    if (Weighted) {
      double p = TotWeight*TDebugRand::Rand()/(RAND_MAX);
      double acum = 0.0;
      for(elem = 0; elem < NTrain; elem++) {
        acum += GetDatWeight(elem);
        if (p <= acum) break; 
      }
      if (elem == NTrain) elem = NTrain - 1;
    }
    else {
      elem = (int)((double)NTrain*TDebugRand::Rand()/(RAND_MAX+1.0));
    }
    incluido[elem] = true;
    for(int j=0; j<Nvar+7;j++) {
      bootstrap->m[i][j] = m[elem][j];
    }
  }

  if (Weighted)
    ResetWeights();

  int nNoUsados=0;
  if (IncludeOobAsTest) {
    //Los ejemplos no utilizados en el bootstrap los pongo como test del dataset
    if ( ! idxs ) {
      for(int i=0; i<NTrain;i++) {
        if (incluido[i]) continue;
        for(int j=0; j<Nvar+7;j++) {
          bootstrap->m[N+nNoUsados][j] =  m[i][j];
        }
        nNoUsados++;
      }
    }
    else {
      idxs->clear();
      for(int i = 0; i < NTrain; i++) {
        if (incluido[i] == 0) {
          idxs->push_back(i);
        }
      }
    }
  }

  delete[] incluido;
  
  //
  bootstrap->redim(N+nNoUsados);
  bootstrap->SetNTrain(N);

  return bootstrap;
}

//-------------------------------------------------  SelectDataOfClass  -----
//---------------------------------------------------------------------------
Data* NomData::SelectDataOfClass(int iClass)
{
  int nClass;

  nClass = 0;
  for(int i = 0; i < NTotal; i++) {
    if (GetDatClass(i) == iClass) {
      m[i][Nvar] = 0;
      nClass++;
    }
    else {
      m[i][Nvar] = 1;
    }
  }
  SortOn(Nvar);

  return Clone(0, nClass - 1);
}
//--------------------------------------  GetBootstrapSampleStratified  -----
//---------------------------------------------------------------------------
Data* NomData::GetBootstrapSampleStratified(int N)
{
  Data *bootstrap = Clone(0, 0);

  //Por omision se extraen igual numero de ejemplos que datos de entrenamiento
  if (N==0) N=NTrain;

  int n = (N + 0.5) / NumClass;

  for(int i = 0; i < NumClass; i++) {
    Data *d = SelectDataOfClass(i);
    d->Scramble(0, d->GetNTotal()-1);
//if (d->GetNTotal() < n+10){
//bootstrap->AddData(d, 0, d->GetNTotal()-1);
//}
//else{
    Data *d1 = d->Clone(0, n > d->GetNTotal() ? d->GetNTotal()-1 : n-1);
//int kk = n;
    Data *d2 = d1->GetBootstrapSample( n );
    bootstrap->AddData(d2, 0, n-1);
    delete d1;
    delete d2;
//}
    delete d;
  } 
  //
  bootstrap->redim(N);
  bootstrap->SetNTrain(N);

  return bootstrap;
}
//-------------------------------------------------  GenerateMIDataSet  -----
//---------------------------------------------------------------------------
Data* NomData::GenerateMIDataSet(NomData *data, int i_mi_column)
{
  MINomData *clon = new MINomData();

  // Cargamos los nombres
  FILE *f=tmpfile();
  fprintf(f, data->GetNamesAsString().c_str());
  rewind(f);
  clon->GetNames(f);
  fclose(f);

  //Copiamos los datos
  int Nvar = data->Nvar;
  clon->init();
  int n_total = data->GetNTotal();
  clon->redim(n_total);
  for(int i = 0; i < n_total; i++) {
    for(int j = 0; j < Nvar + 6; j++) {
      clon->m[ i ][ j ] = data->m[ i ][ j ];
    }
    clon->m[ i ][ Nvar + 6 ] = data->m[ i ][ i_mi_column ];
  }

  // Marcamos la columna de mi
  clon->MarkOrder();
  clon->SortByGroup();
  int i_mi = -1;
  double last_value = clon->m[ 0 ][ Nvar + 6 ] - 1.0;
  for(int i = 0; i < n_total; i++) {
    if (last_value !=  clon->m[ i ][ Nvar + 6 ] ) i_mi++;
    last_value = clon->m[ i ][ Nvar + 6 ];
    clon->m[ i ][ Nvar + 6 ] = i_mi;
  }
  clon->ResetOrder();

  clon->DeleteColumn( i_mi_column );
  clon->SetNTrain( data->GetNTrain() );
  clon->ResetWeights();

  return clon;
}

//---------------------------------------------------------------------------
double Data::DistanciaEjem(int d1, int d2)
{
  double dist = 0.0;
  for (int i=0;i<Nvar-1;i++)
    dist += pow(m[d1][i]-m[d2][i], 2.0);
  return sqrt(dist);
}
//---------------------------------------------------------------------------
double Data::Min(int icol)
{
  double min = numeric_limits<double>::infinity(); 

  for(int i=0; i<NTotal;i++)
    if (m[i][icol]<min) 
      min = m[i][icol];

  return min;
}
double Data::Max(int icol)
{
  double max = -numeric_limits<double>::infinity(); 

  for(int i=0; i<NTotal;i++)
    if (m[i][icol]>max) 
      max = m[i][icol];

  return max;
}
Data *Data::SinteticDataNaiveBayes(int ntot, int ini, int fin)
{
  Data *d = Clone(0, 0);
  d->redim(ntot);
  vector<int> tam_c_0, ini_c_0, tam_c_1;
  int tam_acum;
  double factor;

  factor = (double)ntot/(fin-ini+1);

  //Obtenemos las dimensiones de cada clase
  SortByClass(ini, fin);
  int clase = -1;
  tam_acum = 0;
  ini_c_0.push_back(ini);
  for(int i=ini;i<=fin;i++) {
    if (GetDatClass(i)!=clase) {
      clase = GetDatClass(i);
      tam_c_0.push_back(i-ini_c_0.back());
      tam_c_1.push_back((int)(factor*tam_c_0.back() + 0.5));
      tam_acum += tam_c_1.back();
      ini_c_0.push_back(i);
    }
  }
  tam_c_0.push_back(fin+1-ini_c_0.back());
  tam_c_1.push_back(ntot-tam_acum);

  //Muestreo dentro de cada clase cada atributo de un ejemplo
  //aleatorio
  tam_acum = 0;
  for(unsigned k=0;k<ini_c_0.size();k++) {
    for(int i=0;i<tam_c_1[k];i++) {
      for(int j=0;j<Nvar;j++) {
        int pos = ini_c_0[k] + 
              (int) ((double)tam_c_0[k]*TDebugRand::Rand()/(RAND_MAX+1.0));
        d->m[tam_acum+i][j] = m[pos][j];
      }
    }
    tam_acum += tam_c_1[k];
  }

  return d;
}

int comp(const void* AA, const void *BB)
{
  int ret;
  double* A = *((double**) AA);
  double* B = *((double**) BB);
  if(A[0] > B[0]) ret = 1;
  else if(A[0] < B[0]) ret = -1;
  else ret = 0;
  return ret;
}
Data *Data::SinteticDataNN(int ntot, int ini, int fin, bool ByClass)
{
  Data *d = Clone(0, 0);
  vector<int> tam_c_0, ini_c_0;
  int factor;
  int N = fin-ini+1;

  //solo multiplos enteros de los datos originales
  factor = (int)(0.5+(double)ntot/N);
  ntot = N*factor;
  d->redim(ntot);

  ini_c_0.push_back(ini);
  if (ByClass) {
    //Obtenemos las dimensiones de cada clase
    SortByClass(ini, fin);
    int clase = -1;
    for(int i=ini;i<=fin;i++) {
      if (GetDatClass(i)!=clase) {
        clase = GetDatClass(i);
        tam_c_0.push_back(i-ini_c_0.back());
        ini_c_0.push_back(i);
      }
    }
  }
  tam_c_0.push_back(fin+1-ini_c_0.back());

  //Se buscan los k vecinos mas proximos a cada dato
  int knn=5;

  double **distancias = new double*[fin-ini+1];
  for(int i=0;i<N;i++) {
    distancias[i] = new double[2];
  }

  int iSal = 0;
  //Para cada subconjunto (ByClass o todo)
  for(unsigned k=0;k<ini_c_0.size();k++) {
    for(int i=ini_c_0[k];i<tam_c_0[k]+ini_c_0[k];i++) { //Para cada dato
      //Calculamos distancias a todos los demas del grupo
      for(int j=ini_c_0[k], id=0;j<tam_c_0[k]+ini_c_0[k];j++, id++) {
        distancias[id][0] = j==i ? DBL_MAX : DistanciaEjem(j, i);
      //  printf("%5.3g ", distancias[id][0]); 
        distancias[id][1] = j;
      }
     // printf("\n"); 
      qsort(distancias, tam_c_0[k], sizeof(double**), comp);
     // for(int j=0;j<knn;j++) printf("%g ", distancias[j][1]); 
      //printf("\n"); 
      for(int j=0;j<factor;j++) {
        int pos = (int) ((double)knn*TDebugRand::Rand()/(RAND_MAX+1.0));
        pos = (int)distancias[pos][1];
        //printf("(%d-%d)", i, pos); 
        double gap = (double)TDebugRand::Rand()/RAND_MAX;
        for(int iAtt=0;iAtt<Nvar-1;iAtt++) {
          double dif = m[pos][iAtt]-m[i][iAtt];
          d->m[iSal][iAtt] = m[i][iAtt] + gap*dif;
        }
        int iclass = gap>drand() ? pos : i;     //Clase
        d->m[iSal][Nvar-1] = m[iclass][Nvar-1]; //Clase
        if (m[pos][Nvar-1]!=m[i][Nvar-1]) iSal++; //Generate only frontier examples
      }
     // printf("\n"); 
    }
  }

  d->redim(iSal);

  for(int i=0;i<N;i++) delete []distancias[i];
  delete []distancias;

  return d;
}

Data *Data::BorderData(int ini, int fin, int knn)
{
  Data *d = Clone(ini, fin);
  int N = fin-ini+1;

  double **distancias = new double*[fin-ini+1];
  for(int i=0;i<N;i++) {
    distancias[i] = new double[2];
  }

  int iSal = 0;
  for(int i=ini;i<fin+1;i++) { //Para cada dato
    //Calculamos distancias a todos los demas del grupo
    for(int j=ini, id=0;j<fin+1;j++, id++) {
      distancias[id][0] = j==i ? DBL_MAX : DistanciaEjem(j, i);
      distancias[id][1] = j;
    }
    qsort(distancias, fin-ini+1, sizeof(double**), comp);
    d->m[i][Nvar] = 10;
    for(int j=0;j<knn;j++) {
      int pos = (int)distancias[j][1];
      if (m[pos][Nvar-1]!=m[i][Nvar-1]) {//Distinta clase -> se aniade
        d->m[i][Nvar] = 1;
        iSal++;
        break; 
      }
    }
  }

  d->SortOn(Nvar);
  d->redim(iSal);

  for(int i=0;i<N;i++) delete []distancias[i];
  delete []distancias;

  return d;
}

//---------------------------------------------------------  MINomData  -----
//---------------------------------------------------------------------------
MINomData::MINomData(char *filename, int mi_attrib) : NomData(filename)
{
  NMITrain = 0;
  NMITotal = 0;
  int n_total = GetNTotal();
  for(int i = 0; i < n_total; i++) {
    m[ i ][ Nvar + 6 ] = m[ i ][ mi_attrib ];
  }
  // Marcamos la columna de mi
  MarkOrder();
  SortByGroup();
  int i_mi = -1;
  double last_value = m[ 0 ][ Nvar + 6 ] - 1.0;
  for(int i = 0; i < n_total; i++) {
    if (last_value !=  m[ i ][ Nvar + 6 ] ) i_mi++;
    last_value = m[ i ][ Nvar + 6 ];
    m[ i ][ Nvar + 6 ] = i_mi;
  }
  ResetOrder();

  DeleteColumn( mi_attrib );
  SetNTrain(n_total);
  ResetWeights();
}

MINomData::~MINomData()
{
}
//---------------------------------------------------------------------------
Data* MINomData::GetBootstrapSample(int N, bool Weighted, bool IncludeOobAsTest, std::vector<int> *idxs)
{
  Data *bootstrap = Clone(0, 0);

  //Por omision se extraen igual numero de ejemplos que datos de entrenamiento
  if ( N <= 0 ) N = NMITrain;

  int size = 20 * N;
  bootstrap->redim( size );

  //Ordena y busco donde empieza cada miistance
  SortByGroup(0, NTrain-1);
  int *posiciones = new int[NMITrain+1];
  posiciones[0] = 0;
  for(int i = 1; i < NTrain; i++) { 
    if ( GetDatGroup(i) != GetDatGroup(i-1) ) {
      posiciones[GetDatGroup(i)] = i;
    }
  }
  posiciones[NMITrain] = NTrain;

  //Lo utilizo para marcar los datos que entran en la muestra de bootstrap
  double TotWeight = 0.0;
  int *incluido = new int[NMITrain];
  for(int i = 0; i < NMITrain; i++) {
    incluido[i] = 0;
    TotWeight += Weights[i];
  }

  //Obtengo N elementos al azar con reposicion (muestra bootstrap)
  int imi = 0;
  int igrp = 0;
  for(int i = 0; i < N; i++) {
    int elem;
    if (Weighted) {
      double p = TotWeight*TDebugRand::Rand()/(RAND_MAX);
      double acum = 0.0;
      for(elem = 0; elem < NMITrain; elem++) {
        acum += Weights[elem];
        if (p <= acum) break; 
      }
      if (elem == NMITrain) elem = NMITrain - 1;
    }
    else {
      elem = (int)((double)NMITrain*TDebugRand::Rand()/(RAND_MAX+1.0));
    }
    incluido[elem]++;
    for(int k = posiciones[elem]; k < posiciones[elem+1]; k++) {
      for(int j = 0; j < Nvar + 7; j++) {
      	((MINomData*)bootstrap)->m[imi][j] = m[k][j];
      }
      bootstrap->SetDatGroup(imi, igrp); //Renumero las instancias
      //bootstrap->SetDatWeight(imi, (double)NTrain/NMITrain/(posiciones[elem+1]-posiciones[elem])); //Las doy un peso inv proporcional al no de instancias en el ejemplo
      imi++;
      if ( imi == size ) {
       // if (i < N/2 ) {
          size = 0.5 + (double)(N+5.0)*size/(i+1);
//cout << "estimating size to " << size << endl;
        //}
        //else {
//          size = 1000 + size;
//cout << "increasing to " << size << endl;
  //      }
        bootstrap->redim( size );
      }
    }
    igrp++;
  }
  int NTr = imi;

  if (Weighted)
    ResetWeights();

  if (IncludeOobAsTest) {
    if ( ! idxs ) {
      //Los ejemplos no utilizados en el bootstrap los pongo como test del dataset
      for(int i = 0; i < NMITrain; i++) {
        if (incluido[i] > 0) continue;
        for(int k = posiciones[i]; k < posiciones[i+1]; k++) {
          for(int j=0; j<Nvar+7;j++) {
          	((MINomData*)bootstrap)->m[imi][j] = m[k][j];
          }
          //bootstrap->SetDatGroup(imi, igrp); //Renumero las instancias
          imi++;
        }
        igrp++;
      }
    }
    else {
      idxs->clear();
      for(int i = 0; i < NMITrain; i++) {
        if (incluido[i] == 0) {
          idxs->push_back(i);
        }
      }
    }
  }

  delete[] incluido;
  delete[] posiciones;
  
  //
  bootstrap->redim(imi);
//cout << "reduced to " << imi << endl;
  bootstrap->SetNTrain(NTr); 
  ((MINomData*)bootstrap)->NMITotal = igrp;

  return bootstrap;
}
//---------------------------------------------------------------------------
void MINomData::ResetWeights(int ini, int end)
{
  NomData::ResetWeights(ini, end);
  Weights.clear();
//  cout << "ttt " << NMITotal << endl;
  Weights.resize(NMITotal, 1.0);
}
//---------------------------------------------------------------------------
void MINomData::SetNTrain(int n)
{
  NomData::SetNTrain(n); 
  NMITrain = CountMIInstances(0, n - 1);
  NMITotal = CountMIInstances(0, GetNTotal() - 1 );// No debe estar aqui...
  if ( Weights.size() != NMITotal ) {
    ResetWeights();
  }
}
//---------------------------------------------------------------------------
int MINomData::CountMIInstances(int first, int last)
{
  MarkOrder();
  SortByGroup(first, last);
  int NMI = 0;
  int last_grp = -1;
  for(int i = first; i <= last; i++) {
    if (last_grp != GetDatGroup(i)) {
      last_grp = GetDatGroup(i);
      NMI++;
    }
  }
  ResetOrder();

  return NMI;
}
//-----------------------------------------------------------  AddData  -----
//---------------------------------------------------------------------------
void MINomData::AddData(Data *d2, int ini, int fin)
{
  int ntotalprevio = GetNTotal();

  NomData::AddData(d2, ini, fin);

  //Ordenamos los nuevos datos por grupo
  SortByGroup(ntotalprevio, GetNTotal()-1); 

  //Marcamos los nuevos grupos
  int igrp = NMITotal;
  int last_value =  GetDatGroup(ntotalprevio);
  for(int i = ntotalprevio; i < NTotal; i++) {
    if ( last_value !=  GetDatGroup(i) ) igrp++;
    last_value = GetDatGroup(i);
    SetDatGroup(i, igrp);
  }

  NMITotal = igrp + 1;
  NMITrain = igrp + 1;

  ResetWeights();
  
}

//--------------------------------------------------  Main for testing  -----
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#ifdef _DEBUG_DATA_CPP
int main(int argc, char* argv[])
{
  NomData *data = new NomData("test.cre");

  MINomData *midata = (MINomData*)NomData::GenerateMIDataSet(data, 12);
//  Data *midata = data->Clone();
  delete midata;
  delete data;
}
#endif

