//---------------------------------------------------------------------------

#include "tree20.h"
#include "node20.h"
#include "data20.h"
#include "util20.h"
#include <iostream>
#include <string>
#include <string.h>
#include <limits>
#include <math.h>
#include <cstdlib>
using namespace std;
//---------------------------------------------------------------------------

extern double theSlope;

#define EPS 1e-4
#define EPS2 1e-6
#define EPS3 1e-10

//1 EliminaciOn de CConsole usado para salida
//2 CString de MFC: se sustituye por string de AnsiC++ con un typedef
//3 Anyadido por Gonzalo
//& cosas varias
void dummy_console(char *){}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
PruneInfo::PruneInfo(Node *n) 
{
   node = n;
   min_cost_left.push_back(-1);
   min_cost_right.push_back(-1);
   min_cost.push_back(node->d->Rleaf);
}
double PruneInfo::error(int size)
{
  return min_cost[size];
}
void PruneInfo::merge(int left, int right, double error)
{
  if (left+right+1>=(int)min_cost.size()) {
    min_cost_left.push_back(-1);
    min_cost_right.push_back(-1);
    min_cost.push_back(-1.);
    min_cost_alpha.push_back(-numeric_limits<double>::infinity());
  }
  min_cost_left[left+right+1] = left;
  min_cost_right[left+right+1] = right;
  min_cost[left+right+1] = error;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Tree::Tree(/*1 CConsole *pt*/void (*_console)(char *))
{
  root = NULL;
  maxdepth = 0;
  bestSlope = -1.0;
//1  console = pt;
  console = _console ? _console : dummy_console;
  _Kfinal = -1;
  Tmax = -1;
}
Tree::Tree(Tree*t)
{
  *this = *t;
  root  = new Node(t->root);
}
Tree::~Tree()
{
  Node* cur = root;
  Node* next = NULL;
  while(cur && cur->child) cur = cur->child;

  while(cur)
  {
    next = cur->nextDown();
    delete cur;
    cur = next;
  }

  root = NULL;
  maxdepth = 0;
  _Kfinal = -1;
  Tmax = -1;

}

/*void Tree::BuildRandom(Data* data)
{
  Node *cur,*ori;
  root = new Node(NULL,data->GetNumVarsOrdSplit(),data->GetNumData());
  Ndata = data->GetNumData();
  double chiBest;
  int last;
  int nmax = data->
    cur = root;
    while(cur && nmax>0) {
      int ret = 0;
      ret = data->FindRandomSplit(cur);
      if(ret > 0) {// successful split
        CreateAndFillChildren(cur,last,data,Nmin);
      }
      cur = cur->nextUp();
    }
}             */
void Tree::BuildRandomForestMember(Data* data,int Nmin, int num_random_vars, 
                                                      vector<int> *atts_to_use)
{
  double chiBest;
  int last;
  Node *cur;

  root = new Node(0, data->GetNumVarsOrdSplit(), data->GetNumData());
  cur = root;

  Ndata = data->GetNumData();

  vector<int> atts(0);

  if (atts_to_use) {
    for (unsigned i = 0; i < atts_to_use->size(); i++) {
      atts.push_back((*atts_to_use)[i]);
    }
  }
  else {
    for (int i = 0; i < data->GetNumVar()-1; i++) {
      atts.push_back(i);
    }
  }

  while(cur) {
    cur->d->MembershipTotal = data->FindMembership(cur);
    if(cur==root) data->FixLeaf(cur);
    if (cur->last-cur->first+1 > Nmin && cur->d->Rleaf > 0.0) {
       vector<bool> use(data->GetNumVar()-1, false);
       for (int i=0;i<num_random_vars;i++) {
         int j = (int)(((double)atts.size()*TDebugRand::Rand())/(RAND_MAX+1.0));
         int hold = atts[j];
         atts[j] = atts[i];
         atts[i] = hold;
       }
       for(int i=0;i<num_random_vars;i++) {
         use[atts[i]] = true;
       }
       int ret = data->FindBestSplit(cur,Nmin,&last, &chiBest, &use);
       if(ret > 0) {// successful split
         CreateAndFillChildren(cur,last,data,Nmin);
       }
    }
    cur = cur->nextUp();
  }

  _Kfinal = -1;
  SizeBest = root->T;
  Tmax = root->T;

  //Calculate family of K
  int K = -1; 
  if(root->T > 1) FindAlpha(Ndata, -1);
  double alpha = 0.0;
  K = 1;
  while(root->T > 1) {
    Prune(alpha,K);  // Construct T_K
    if(root->T > 1) alpha = FindAlpha(Ndata, K);
    K++;
  }

}
void Tree::BuildDecisionStump(Data* data,int Nmin)
{
  double chiBest;
  int last;

  root = new Node(NULL,data->GetNumVarsOrdSplit(),data->GetNumData());
  Ndata = data->GetNumData();

  root->d->MembershipTotal = data->FindMembership(root);
  data->FixLeaf(root);

  if((root->last - root->first + 1) > Nmin && root->d->Rleaf > 0.0) {
    if (data->FindBestSplit(root,Nmin,&last, &chiBest)>0) {
      CreateAndFillChildren(root,last,data,Nmin);
    }
  }
  _Kfinal = -1;
  SizeBest = root->T;
  Tmax = root->T;
}
void Tree::Build(Data* data,int Nmin, int GroupNumber)
{
  Node *cur;
  vector<Node*> nodos;//Added by GMM
  if (!root) {
    root = new Node(NULL,data->GetNumVarsOrdSplit(),data->GetNumData());
    root->SetGeneratedByGroup(GroupNumber);
    nodos.push_back(root);
  }
  else {                            //Added by GMM
    /*int y = */ThrowData(data, 0, data->GetNGrow()-1, -2, root);

    cur = root;
    while (cur->child) cur = cur->child;
    while (cur) {
      if (cur->IsLeaf(-2)) nodos.push_back(cur);
      cur = cur->nextDown(-2);
    }
  }
  Ndata = data->GetNumData();
  double chiBest;
  int last;
//  for (int j=0;j<nodos.size();j++) {
//    cur = nodos.at(j);
    cur = root;
    while (cur->child) cur = cur->child;          //Added by GMM

    while(cur) {
      int ret = 0;
      cur->d->MembershipTotal = data->FindMembership(cur);
      if(cur==root)
        data->FixLeaf(cur);

//      if( cur->d->MembershipTotal >= 2*Nmin && cur->d->Rleaf > 0.0) // Original
//      if( !cur->d->DoesntWantSplitting && (cur->last-cur->first+1) >= 2*Nmin &&  //Changed by GMM
//                                                               cur->d->Rleaf > 0.0)
      if( !cur->d->DoesntWantSplitting && ( cur->last - cur->first + 1 ) > Nmin &&  //Changed by GMM
                                                            cur->d->Rleaf > 0.0) {
        ret = data->FindBestSplit(cur,Nmin,&last, &chiBest);
      }
      if(ret > 0) {// successful split
        CreateAndFillChildren(cur,last,data,Nmin);
        cur->child->SetGeneratedByGroup(GroupNumber);
        cur->child->sib->SetGeneratedByGroup(GroupNumber);
      }
      cur = cur->nextUp();
      while (cur && cur->child)
        cur = cur->nextUp();
    }
//  }
}
void Tree::CreateAndFillChildren(Node* n, int last, Data* data, int Nmin)
{
  Node* n1 = new Node(n,data->GetNumVarsOrdSplit(),data->GetNumData()); // first child
  Node* n2 = new Node(n,data->GetNumVarsOrdSplit(), data->GetNumData()); // second child

  n1->depth = n2->depth = n->depth +1;
  if(n1->depth> maxdepth)
    maxdepth = n1->depth;
  n1->first = n->first;
  n2->last = n->last;
  if(n->d->CrispFlag)
  {
    n1->last = last;
    n2->first = last+1;
  }
  else
  {
    n1->last = n->last;
    n2->first = n->first;
  }


  n1->d->MembershipTotal = data->FindMembership(n1);
  data->FixLeaf(n1);
  if(n1->d->MembershipTotal < Nmin)
  {
//    (*console)("Error: CreateandFillChildren0");
//1    console->printf("Error: CreateandFillChildren");
  }
  n2->d->MembershipTotal = data->FindMembership(n2);
  data->FixLeaf(n2);
  if(n2->d->MembershipTotal <Nmin)
  {
//    (*console)("Error:CreateandFillChildren1");
//1    console->printf("Error:CreateandFillChildren");
  }

  if(fabs(n->d->MembershipTotal - n1->d->MembershipTotal - n2->d->MembershipTotal) > 1e-6)
  {
//    (*console)("Error: CreateandFillChildren2 ");
//1    console->printf("Error: CreateandFillChildren ");
  }

  if(fabs(min( n1->d->MembershipTotal, n2->d->MembershipTotal) - n->d->min_in_child) > 1e-1)
  {
//    (*console)("Error: CreateandFillChildren3 ");
//1    console->printf("Error: CreateandFillChildren ");
  }

  n1->d->RsubTree = n1->d->Rleaf;
  n2->d->RsubTree = n2->d->Rleaf;

  if(n)
  {
    double deltaR = n1->d->Rleaf + n2->d->Rleaf - n->d->RsubTree;

    while(n) {
      n->d->RsubTree += deltaR;
      n->T++;
      n = n->parent;
    }
  }
}

//MCT info
//int fam_mct[4096];

void Tree::CART(Data* data, int Nmin, bool PruneTree, bool SE_0)
{
//  using namespace crepo;
  int ncv = data->GetNumCV();
  char kk[1024];

  bool depvarnom = (data->GetDepVarType() == NOM);

  if(!ncv)
  {
    data->SetCV(0);
    Build(data,Nmin);
    FILE* fp = 0;//fopen("out.txt","a");
    if(fp) fprintf(fp,"%d*\t",root->T);
//    int K = 0;

//    sbest = data->Score(-1,-1,root,K);
    if(fp) fprintf(fp,"%g+\t",sbest);
    IndepPrune();
    SizeBest = root->T;
    RMSres = root->d->RsubTree;
  //  sbest = data->Score(-1,-1,root,K);
    if(depvarnom)
    {

      if(fp)
      {
        fprintf(fp,"%dç\t %g-",SizeBest,sbest);
      }
      sprintf(kk," Size = %d \t ErrorRS = %g \t ErrorTest = %g (%g) ", 
                            SizeBest, RMSres/data->GetNGrow(), sbest  ,1-sbest);
      (*console)(kk);
    }
    else   {
      sprintf(kk,"Size = %d \t RMSres = %g RMStest = %g",SizeBest,
                          RMSres/data->GetYGvscale(),sbest/data->GetYGvscale());
      (*console)(kk);
    }

    if (fp) fclose(fp);

    return;
  }

  sprintf(kk,"Building final tree ");
  (*console)(kk);
  data->SetCV(0);
  Build(data,Nmin);
  Tmax = root->T;
  FILE* fp = 0;//fopen("out.txt","a");
  if(fp) fprintf(fp,"%d#\t",root->T);

  bool UseSqrt = false;
//char arb[16000];
  // Just to fix the initial value of nodes' alpha
  if(root->T > 1) FindAlpha(Ndata, -1, UseSqrt);
 
  double alpha = 0.0, alpha_prev = 0.0,alpha_mean,RCV = 0.0;
  int K = 0;
  int KOld = 0;
  int SizeOld =0;
  _Kfinal = 0;
  double Rbest = -1.0,SE=0.0,SEbest = 0.0;;
  double ss = -1.0;
  sbest = 0.0 ;
  RMSres = 0.0;  // Resubstitution error

//FILE *tablas=fopen("tablas.txt", "a");//SG

//arb[0]='\0';
//Write(data, -1, arb);
//fprintf(tablas, "%s\n", arb);//SG

  if (PruneTree)  {
    Tree **aux = new Tree*[ncv];

    for(int i=0;i<ncv;i++) {
       sprintf(kk,"Building CV tree %d of %d ",i,ncv);
       (*console)(kk);
      data->SetCV(i+1);
      aux[i] = new Tree(console);
      aux[i]->Build(data,Nmin);
      data->FixCVError(aux[i]->root);
//S   Obtener familia MCT de arbol aux[i]
    }

//S Obtener familia MCT de arbol this
//MCT info
//double mat[5][4096];
//int col=0;
//int imct=0;
    while(root->T > 1) //S for i=->n_alphas
    {
      Prune(alpha,K);  // Construct T_K
//MCT info
//if (col >= 4096) {
//col = 4095;
//cout << "EEEEEEEEEERRRRRRRRRRRRRRRRRROOOOOOOOOOOOOOOOOOOOOOOOOOORRRRRRRRRRRRRRRRRRRRRRRR" << endl;
//}
//mat[0][col] = root->T;
//mat[1][col] = alpha;
//mat[2][col] = fam_mct[imct] == root->T ? 1: 0;
//if (fam_mct[imct] == root->T) imct++;
    //  ss = data->Score(-1,-1,root,K);
      if(K==0)
        if(fp)fprintf(fp,"%g/\t",ss); //GMM
      if(root->T > 1)
      {
        alpha_prev = alpha; alpha = FindAlpha(Ndata, K, UseSqrt);//S Sobra
        alpha_mean = sqrt(alpha_prev*alpha);
        data->InitialiseErrorCV();
        for(int i=0;i<data->GetNumCV();i++)
        {//GMM: Bug FIXED. No seleccionaba bien dentro de los arboles
         //     de cv las subfamilias, en ocasiones seleccionaba un arbol
         //     intermedio
          if(aux[i]->root->T > 1) {
            double a;
            do {
              a = aux[i]->FindAlpha(aux[i]->Ndata);
              if (a<=alpha_mean || fabs(alpha_mean - a) < EPS3) aux[i]->Prune(a);
              else break;
            } while (true);
          }
          data->AccumulateErrorCV(aux[i]->root);
//OLD          if(aux[i]->root->T > 1) aux[i]->FindAlpha(aux[i]->Ndata, -1, UseSqrt);
//OLD          aux[i]->Prune(alpha_mean);
//OLD          data->AccumulateErrorCV(aux[i]->root);

        }
        data->GetErrorCV(&RCV,&SE);
  
      } else RCV =  root->d->Rleaf/data->GetNGrow();
      if(Rbest < 0.0 || RCV < Rbest) {Rbest = RCV; SEbest = SE;}
  
//MCT info
//mat[3][col] = RCV;
//mat[4][col] = SE;
//col++;
      // The SE rule
      double SEaux = SE_0 ? 0.0 : 1.0;
      
      if(RCV <= Rbest+ SEaux*SEbest)
      {
        KOld = _Kfinal;
        _Kfinal = K;
        SizeOld = SizeBest;
        SizeBest = root->T;
        sbest = ss;
        if(depvarnom)
          RMSres = root->d->RsubTree;
        else RMSres= sqrt(root->d->RsubTree/data->GetNGrow())/data->GetYGvscale();
      }
//fprintf(tablas, "%5g  %5d, ", root->d->RsubTree, root->T);//SG
//arb[0]='\0';
//Write(data, K, arb);
//fprintf(tablas, "ARBOL INI:\n%s\nARBOL_FIN\n", arb);//SG
  //jfl
//printf("------------------------- K=%d \n", K);
      sprintf(kk,"K = %d Size = %d R = %g RCV = %g (%g) alpha = %g", K, 
                               root->T,root->d->RsubTree/root->d->Rleaf,RCV,SE,alpha);
      (*console)(kk);
//Write(data);
      double factor = data->GetYGvscale()*data->GetYGvscale();
  
      if(data->GetNTest())
      {
        if(depvarnom)
        {
  //GMM        sprintf(kk,"K = %d  Size = %d  R = %g  RCV = %g (%g)  alpha = %g  Rtest = %g ",
  //GMM          K,root->T,root->d->RsubTree/(data->GetNGrow()),RCV,SE,alpha, ss);
  //GMM        (*console)(kk);
  
        }
        else {
  
  //        sprintf(kk, "K = %d  Size = %d  R = %g  RCV = %g (%g)  alpha = %g  Rs = %g  RMS = %g \n  NormRMS = %g ",
    //        K,root->T,root->d->RsubTree/(data->GetNGrow()*factor),RCV/factor,SE/factor,alpha,ss*ss/factor,ss/sqrt(factor), ss/(sqrt(factor)* ((OrdData*)data)->GetRMSTest()));
      //    (*console)(kk);
        }
      }
      else
      {
        if(depvarnom)
        {
//      (*console)("Hola1");
          sprintf(kk, "K = %d  Size = %d  R = %g  RCV = %g (%g)   alpha = %g ",
            K,root->T,root->d->RsubTree/(data->GetNGrow()),RCV,SE,alpha);
          (*console)(kk);

        }
        else {
  
          sprintf(kk, "K = %d   Size = %d   R = %g  RCV = %g (%g)   "
"alpha = %g ", K,root->T,root->d->RsubTree/(data->GetNGrow()*factor),
RCV/factor,SE/factor,alpha);
          (*console)(kk);
        }
      }
    K++;
  }
//MCT info begin
//FILE *ff = fopen("datos_familia_cart.txt", "a");
//for(int i=0;i<5;i++) {
//for(int j=0;j<col;j++) {
//if (i==0 || i==2) fprintf(ff, "%d\t", (int)(0.5+mat[i][j]));
//else fprintf(ff, "%.3f\t",100.0* mat[i][j]);
//}
//fprintf(ff, "\n");
//}
//for(int j=0;j<col;j++){
//fprintf(ff, "%d\t", j==_Kfinal ? 1 : 0);
//}
//fprintf(ff, "FINFIN\n");
//fclose(ff);
//MCT info end

  //  if(_Kfinal == K-1) _Kfinal--;  // JFL version: Uncomment in order to avoid  tree = root (useless for further optimisation)
    if(_Kfinal == K-1)
    {
      _Kfinal = KOld;  // ASG version: Take the next best tree as determined by CV
      SizeBest = SizeOld;
    }
//    (*console)("Raw tree:");
//    Write(data);
  //  (*console)("Pruned tree:");
    //Write(data,_Kfinal);

  //6  for(int i=1;i<=ncv;i++) delete (aux[i]);
    for(int i=0;i<ncv;i++) delete (aux[i]);//6
    delete[] aux;

  //  double err = data->Score(-1,-1,root,_Kfinal);
  }// Fin if (PruneTree)
  else {
    _Kfinal = -2;
    //sbest = data->Score(-1,-1,root,_Kfinal);
    SizeBest = root->T;
    RMSres = root->d->RsubTree;
//    sbest = data->Score(-1,-1,root,_Kfinal);
}


//  if(depvarnom)
//  {

    if(fp)
    {
//GMM        fprintf(fp,"%d \t %g ",SizeBest,err);
        fclose(fp);
    }
//GMM      sprintf(kk, " ErrorRS = %g \t ErrorTest = %g (%g) ",RMSres/data->GetNGrow(), err  ,1-err);
//GMM      (*console)(kk);
//  }
//  else {
//GMM      sprintf(kk, "RMSres = %g RMStest = %g",RMSres/data->GetYGvscale(),err/data->GetYGvscale());
//GMM      (*console)(kk);
//  }

/*for(int i=-1;i<10;i++) {
arb[0]='\0';
Write(data, i, arb);
printf("------------------------- K=%d \n  %s \n------------------", i, arb);
}*/

//fprintf(tablas, "\n");//SG
//fclose(tablas);//SG
}

void Tree::PruneMCT(Data* data, int Nmin, bool SE_0)
{
  double alpha = 0.0, alpha_prev = 0.0,alpha_mean, RCV = 0.0;
  int K = 0;
  int KOld = 0;
  int SizeOld =0;
  _Kfinal = 0;
  double Rbest = -1.0,SE = 0.0,SEbest = 0.0;;
  double ss = -1.0;
  sbest = 0.0 ;
  RMSres = 0.0;  // Resubstitution error

  // Anadidas para que compile
  char kk[1024];
  FILE* fp = 0;
  bool depvarnom = (data->GetDepVarType() == NOM);

  int ncv = data->GetNumCV();
  Tree **aux = new Tree*[ncv];  

  for(int i=0;i<ncv;i++) {
    sprintf(kk,"Building CV tree %d of %d ",i,ncv);
    (*console)(kk);
    data->SetCV(i+1);
    aux[i] = new Tree(console);
    aux[i]->Build(data, Nmin);
    data->FixCVError(aux[i]->root);
    //S: Obtener familia MCT de arbol aux[i]
    map<Node*, PruneInfo*> trees = aux[i]->MinimumCostTrees();
    aux[i]->SelectMinimumCostTrees(trees);
  }
  
  //S: Obtener familia MCT de arbol this
  map<Node*, PruneInfo*> trees = MinimumCostTrees();
  vector<double> alphas = SelectMinimumCostTrees(trees); // Esto actualiza las K
  int ia = 0;
//Info MCT
//double mat[5][4096];
//int col=0;
//int imct=0;
  while(root->T > 1)
  {
    PruneMCT2(alpha);  // Construct T_K
//Info MCT
//fam_mct[imct] = root->T;
//imct++;
//mat[0][col] = root->T;
//mat[1][col] = alpha;
//mat[2][col] = 1;
    ss = data->Score(-1,-1,root,K);
    if(K==0)
      if(fp)fprintf(fp,"%g/\t",ss); //GMM
    if(root->T > 1)
    {
      // Tenemos todos los alphas en los PruneInfo, pero no sabemos con cuales
      // nos hemos quedado porque en SelectMinimumCostTrees() solamente actualizamos las K
      alpha_prev = alpha; alpha = alphas[ia];//FindAlpha(Ndata, K, true); // Calcula el alpha con la raiz
      ia++;
      alpha_mean = sqrt(alpha_prev*alpha);
      data->InitialiseErrorCV();
      for(int i=0;i<data->GetNumCV();i++)
      {
//        if(aux[i]->root->T > 1) aux[i]->FindAlpha(aux[i]->Ndata, -1, true); //S Sobra
        aux[i]->PruneMCT2(alpha_mean); //S Queda igual ??
        data->AccumulateErrorCV(aux[i]->root); //S Queda igual ??
      }
      data->GetErrorCV(&RCV,&SE);
    } 
    else RCV =  root->d->Rleaf/data->GetNGrow();

    if(Rbest < 0.0 || RCV < Rbest) {Rbest = RCV; SEbest = SE;}
  
//Info MCT
//mat[3][col] = RCV;
//mat[4][col] = SE;
//col++;
    // The SE rule
    double SEaux = SE_0 ? 0.0 : 1.0;
      
    if(RCV <= Rbest+ SEaux*SEbest)
    {
      KOld = _Kfinal;
      _Kfinal = K;
      SizeOld = SizeBest;
      SizeBest = root->T;
      sbest = ss;
      if(depvarnom)
        RMSres = root->d->RsubTree;
      else RMSres= sqrt(root->d->RsubTree/data->GetNGrow())/data->GetYGvscale();
    }

    K++;
  } // end while(root->T > 1)
//Info MCT
//FILE *ff = fopen("datos_familia_mct.txt", "a");
//for(int i=0;i<5;i++) {
//for(int j=0;j<col;j++) {
//if (i==0 || i==2) fprintf(ff, "%d\t", (int)(0.5+mat[i][j]));
//else fprintf(ff, "%.3f\t",100.0* mat[i][j]);
//}
//fprintf(ff, "\n");
//}
//for(int j=0;j<col;j++){
//fprintf(ff, "%d\t", j==_Kfinal ? 1 : 0);
//}
//fprintf(ff, "\n");
//fclose(ff);
  if(_Kfinal == K-1)
  {
    _Kfinal = KOld;  // ASG version: Take the next best tree as determined by CV
    SizeBest = SizeOld;
  }

  //6  for(int i=1;i<=ncv;i++) delete (aux[i]);
  for(int i=0;i<ncv;i++) delete (aux[i]);//6
  delete[] aux;
}
void Tree::PruneMCT2(double alpha)
{
  FreeSubTreeMCT(root, alpha);
}
void Tree::FreeSubTreeMCT(Node* n, double alpha)
{
//  Node* cur = n;

  if(n->child && (fabs(n->child->alpha - alpha) < EPS3  || n->child->alpha < alpha)) {
    double deltaR = n->d->Rleaf - n->d->RsubTree;
    double deltaRCV = n->d->RleafCV - n->d->RsubTreeCV;
    int deltaT = 1 - n->T;
    double deltas1 = n->d->s1leaf - n->d->s1;
    double deltas12 = n->d->s12leaf - n->d->s12;
    while(n)
    {
      n->d->RsubTree += deltaR; n->d->RsubTreeCV += deltaRCV;
      n->T += deltaT;
      n->d->s1 += deltas1; n->d->s12 += deltas12;
      n = n->parent;
    }
  }
  else {
    if (n->child) {
      FreeSubTreeMCT(n->child, alpha);
      FreeSubTreeMCT(n->child->sib, alpha);
    }
  }
}


void Tree::Prune(double alpha, int K)
{
  Node* cur = root;  // Construct T_K
  while(cur->child && (K < 0 || cur->child->K < 0) ) cur = cur->child;
  while(cur)
  {
    if(cur->child && (K < 0 || cur->child->K < 0))
      if(fabs(cur->alpha - alpha) < EPS3 || cur->alpha < alpha)
        FreeSubTree(cur, K);
    cur = cur->nextDown(K);
  }
}

void Tree::DeleteNodesAndResetK(Node *n, int K)
{
//Delete the nodes with K>=0 and K<(parameter K)
  vector<Node*> nodos;
  Node *cur;
  nodos.push_back(n);
  while(nodos.size()>0) {
    cur=nodos.at(0);
    nodos.erase(nodos.begin());
    cur->K = -1;
    if (!cur->child) continue;
    if (cur->child->K>=0 && cur->child->K<K) {
      FreeSubTree(cur, -2);
    }
    else {
      nodos.push_back(cur->child);
      nodos.push_back(cur->child->sib);
    }
  }
  cur = n;
  while (cur->child) cur=cur->child;
  while (cur && cur!=n) {
    cur->T = 0;
    cur = cur->nextDown();
  }
  n->T = 0;

  cur = n;
  while (cur->child) cur=cur->child;
  while (cur) {
    if (cur->IsLeaf(K)) {
      Node *cur2;
      cur2 = cur;
      while (cur2) {
        cur2->T++;
        cur2 = cur2->parent;
      }
    }
    cur = cur->nextDown();
  }
/*  n->K = -1;
  if (!n->child) {
    n->T =  1;
  }
  else if (n->child->K>=0 && n->child->K<K) {
    FreeSubTree(n, -2);
    n->T =  1;
  }
  else {
    DeleteNodesAndResetK(n->child, K);
    DeleteNodesAndResetK(n->child->sib, K);
  }*/
}
void Tree::Prune2(double alpha, int ndata, int K, bool ComoArticulo)
{
  FILE *kk=fopen("podilla.txt", "at");

  Node* cur = root;  // Construct T_K
  while(cur->child && (K < 0 || cur->child->K < 0)) cur = cur->child;
  while(cur)
  {
    cur->d->RPrune = cur->d->Rleaf + alpha*ndata;
    cur = cur->nextDown(K);
  }
  cur = root;
  while(cur->child && (K < 0 || cur->child->K < 0)) cur = cur->child;
  while(cur)
  {
    if (cur->IsLeaf(K)) {
      cur->d->RPruneSubTree = cur->d->RPrune;
      for(int jj=0;jj<cur->depth;jj++) fprintf(kk, "  ");
      fprintf(kk,"ESub=%g  RNodo=%g <-----\n", cur->d->RPruneSubTree, cur->d->RPrune);
    }
    else {
      bool HayQuePodar;
      cur->d->RPruneSubTree = cur->child->d->RPruneSubTree +
                                                 cur->child->sib->d->RPruneSubTree;
      if (ComoArticulo) {
        //the branch is pruned if the error of its leaves is equal or bigger
        HayQuePodar =fabs(cur->d->RPrune - cur->d->RPruneSubTree) < EPS2 ||
                                              cur->d->RPrune <= cur->d->RPruneSubTree;
      }
      else {
        //the branch is pruned if the error of its leaves is bigger
        HayQuePodar = fabs(cur->d->RPrune - cur->d->RPruneSubTree) > EPS3 &&
                                               cur->d->RPrune < cur->d->RPruneSubTree;
      }
      if(HayQuePodar) {
        fprintf(kk,"----------------CUT----------------\n");
        for(int jj=0;jj<cur->depth;jj++) fprintf(kk, "  ");
        fprintf(kk,"ESub=%g  RNodo=%g\n", cur->d->RPruneSubTree, cur->d->RPrune);
        FreeSubTree(cur, K);
        cur->d->RPruneSubTree = cur->d->RPrune;
      }
      else {
        for(int jj=0;jj<cur->depth;jj++) fprintf(kk, "  ");
        fprintf(kk,"ESub=%g  RNodo=%g\n", cur->d->RPruneSubTree, cur->d->RPrune);
      }
    }
    cur = cur->nextDown(K);
  }
  //Sets the labels of the leaves to true
  cur = root;
  while(cur->child && (K < 0 || cur->child->K < 0)) cur = cur->child;
  while(cur) {
    if (cur->IsLeaf(-2)) {
      if (cur->d->IsLeafAfterPrune)
        cur->d->DoesntWantSplitting = true;
      else
        cur->d->IsLeafAfterPrune = true;
    }
    cur = cur->nextDown(-2);
  }
 fclose(kk);
}
//--------------------------------------------------------------------------
std::map<Node*, PruneInfo*> Tree::MinimumCostTrees()
{
  std::map<Node*, PruneInfo*> min_cost_trees;
  int k=0, i=0, j=0, maximo=0, minimo=0;
  double mincost = std::numeric_limits<double>::infinity(), cost = 0;

  //Se recorre el arbol en profundidad
  Node* cur = root;
  while(cur->child) cur = cur->child;
  while(cur) {
    min_cost_trees[cur] = new PruneInfo(cur); // Liberar despues!!
    if (cur->child) {
      for (k=2; k<=cur->T; k++) {
        mincost = std::numeric_limits<double>::infinity();
        maximo = (1>k-(cur->child->sib->T)) ? 1:(k - cur->child->sib->T);
        minimo = (k-1<cur->child->T) ? (k-1):(cur->child->T);
        for (i = maximo; i <= minimo; i++) {
          j = k - i;
          cost = min_cost_trees[cur->child]->error(i-1) 
               + min_cost_trees[cur->child->sib]->error(j-1);
          if (cost < mincost) {
            mincost = cost;
            min_cost_trees[cur]->merge(i-1, j-1, cost);
          } // if cost end
        } // for i end    
      } // for k end
    } // if cur->child end
    cur = cur->nextDown(-1);
  }

  return min_cost_trees; 
}
//--------------------------------------------------------------------------
vector<double> Tree::SelectMinimumCostTrees(std::map<Node*, PruneInfo*> min_cost_trees) 
{
  PruneInfo *pi = min_cost_trees[root];
  int kl;  // Numero de hojas del arbol = indice en min_cost_trees + 1
  int new_kl = 1;
  double alpha = 0, new_alpha = 0;
  vector<double> alphas;

  // Empezamos por el arbol mas pequeño que tenga el mismo error que el completo
  kl = pi->min_cost.size()-1;
  for (int i = pi->min_cost.size()-2; i>=0; i--) {
    if (pi->min_cost[i] == pi->min_cost[pi->min_cost.size()-1])
      kl = i;
  }

  //Inicializamos la K
  Node* cur = root;
  while(cur->child) cur = cur->child;
  while(cur) {
    cur->K = 0;
    cur->alpha = 0.0;
    cur = cur->nextDown(-1);
  }

  // Mientras no lleguemos a la raiz 
  kl++;
  while (kl > 1) {
    alpha = std::numeric_limits<double>::infinity();
    for (int k = kl -1 ; k >= 1 ; k--) {
//      new_alpha = (pi->min_cost[k-1] - pi->min_cost[kl-1])/ ((kl - k)*Ndata);
      new_alpha = (pi->min_cost[k-1] - pi->min_cost[kl-1])/ (Ndata*(sqrt(kl) - sqrt(k)));
      if (new_alpha <= alpha) {
        alpha = new_alpha;
        new_kl = k;	
      }
      pi->min_cost_alpha[k-1] = new_alpha;
    }
    alphas.push_back(alpha);
//    cout << "[kl=" << kl << ":k=" << new_kl << "]Con alpha: " << alpha << "\n";

    // Recorrer arbol cambiando el numero de poda(solo se me ocurren recursivos)
    UpdateMinimumCostTrees(min_cost_trees, root, kl-1, alpha);
    kl = new_kl;
  }
  root->alpha = std::numeric_limits<double>::infinity();
  root->K = -1;

  return alphas;
}
//--------------------------------------------------------------------------
void Tree::UpdateMinimumCostTrees(std::map<Node*, PruneInfo*> min_cost_trees, 
                                            Node* cur, int index, double alpha)
{

  PruneInfo *pi = 0;

  cur->K++;
  cur->alpha = alpha;

  if (index == 0)
    return;

  pi = min_cost_trees[cur];
//cout << cur << ": " << pi->min_cost_alpha[index] << endl;

  UpdateMinimumCostTrees(min_cost_trees, cur->child, pi->min_cost_left[index], 
                                                                        alpha);
  UpdateMinimumCostTrees(min_cost_trees, cur->child->sib, 
                                             pi->min_cost_right[index], alpha);
}
//--------------------------------------------------------------------------
void Tree::IndepPrune()
{
  int K = 0;
  for(int i = maxdepth; i >=0; --i)
  {
    Node* cur =root;
    while(cur)
    {
      bool curIsLeaf = cur->IsLeaf(K);
      if(cur->depth == i)
      {
        if(curIsLeaf)
          cur->d->RPruneSubTree = cur->d->RPrune;
        else
        {
          Node* temp = cur->child;
          cur->d->RPruneSubTree = temp->d->RPruneSubTree + temp->sib->d->RPruneSubTree;
        }

        if(!curIsLeaf && cur->d->RPrune <= cur->d->RPruneSubTree)
        {
          FreeSubTree(cur,K);
          cur->d->RPruneSubTree = cur->d->RPrune;
          curIsLeaf = true;
        }
      }

      if(!curIsLeaf) cur = cur->child;
      else {
        while(cur && !cur->sib) cur = cur->parent;
        if(cur) cur = cur->sib;
      }
    }
  }

}



double Tree::FindAlpha(int Ndata, int K, bool UseSqrt)
{
  Node* cur = root;
  while(cur->child && (K < 0 || cur->child->K < 0)) cur = cur->child;
  double alpha = 10E300;//-1; (GMM)
  while(cur)  // Find alpha_K+1
  {
    cur->ComputeAlpha(Ndata, root, UseSqrt);
    if(cur->child && (K < 0 || cur->child->K < 0))
      if(cur->alpha < alpha/* || alpha < 0.0*/) //GMM
        alpha = cur->alpha;
    cur = cur->nextDown(K);
  }
  return alpha;
}

void Tree::RecomputeNodeStats(Data *data, int first, int last, bool Add)
{
  int Nvar = data->GetNumVar();
  Node *cur = root;

  cur->first = first;
  cur->last  = last;
  data->FixLeaf(cur);

  while (cur = cur->nextUp() ) {
    cur->first = first;
    cur->last  = last;
    data->FindMembership(cur);
    data->SortOn(Nvar+3, first, last, true);
    for(int i = first; i <= last; i++ ) {
      if (data->GetInstance(i).GetMembership()==0.0) {
         cur->first = first;
         cur->last = i - 1;
         data->FixLeaf(cur, Add); 
         break;
      }
    }
  }

}

void Tree::FreeAuxData()
{
  Node* cur = root;

  while (cur->child) { 
    cur = cur->child;
  }

  while (cur != root) {
    cur->FreeAuxData();
    cur = cur->nextDown();
  }
  root->FreeAuxData();
}

void Tree::FreeSubTree(Node* n, int K)
{
  Node* cur = n;
  while(cur->child && (K < 0 || cur->child->K < 0)) cur = cur->child;

  while(cur != n)
  {
    Node* next = cur->nextDown(K);
    if(K < 0) delete cur;
    else cur->K = K;
    cur = next;
  }
  if(K<0) n->child = NULL;
  double deltaR = n->d->Rleaf - n->d->RsubTree;
  double deltaRCV = n->d->RleafCV - n->d->RsubTreeCV;
  int deltaT = 1 - n->T;
  double deltas1 = n->d->s1leaf - n->d->s1;
  double deltas12 = n->d->s12leaf - n->d->s12;
  while(n)
  {
    n->d->RsubTree += deltaR; n->d->RsubTreeCV += deltaRCV;
    n->T += deltaT;
    n->d->s1 += deltas1; n->d->s12 += deltas12;
    n = n->parent;
  }
}

void Tree::Write(Data* data, int K, char *buff)
{
  char buf[32765];
  char palotes[1024];
  Node* cur;
  int level;

  for(int i=0;i<1024;i++) palotes[i] = ' ';

  cur = root;
  level = 0;
  while(cur) {
    if(level) {
      sprintf(buf,"   ");
      for(int i=0;i<level-1;i++) sprintf(buf+strlen(buf),"%c  ", palotes[i+1]);
    } 
    else { 
      buf[0] = '\0';
    }

    strcat(buf, "+--");
    data->WriteNode(cur,buf+strlen(buf),K);
    if (buff) {
      strcat(buff, buf);  
//      strcat(buff, "\n");  
    }
    else {
      (*console)(buf);
    }

//    if(cur->child && (K < 0 || cur->child->K > K)) {
    if(!cur->IsLeaf(K)) {
      cur = cur->child;
      level++;
      palotes[level] = '|';
    }
    else {
      while(cur->parent && !cur->sib) {
        cur = cur->parent;
        level--;
      }
      cur = cur->sib;
      palotes[level] = ' ';
    }
  }
}
int Tree::CountNodes(int K)
{
  /*Node* cur = root;
  int num=0;
  while(cur && cur->child && (K<0 || cur->child->K<0)) cur = cur->child;
  while(cur)
  {
    num++;
    cur = cur->nextDown(K);
  }
  return num;*/
  Node* cur = root;
  int num=0;
  while(cur)
  {
    num++;
    if(!cur->IsLeaf(K)) cur = cur->child;
    else {
      while(cur && !cur->sib) cur = cur->parent;
      if(cur) cur = cur->sib;
    }
  }
  return num;
}

int Tree::CountParams(int K)
{
  Node* cur = root;
  int num=0;
  while(cur)
  {
    num += cur->CountParams(K);
    if(!cur->IsLeaf(K)) cur = cur->child;
    else {
      while(cur && !cur->sib) cur = cur->parent;
      if(cur) cur = cur->sib;
    }
  }
  return num;
}

void Tree::PrepareOpt(double* a, int K, Data* data)
{
  Node* cur = root;
  double* p = a;
  while(cur)
  {
    cur->InitParams(&p,K, data);
    if(!cur->IsLeaf(K)) cur = cur->child;
    else {
      while(cur && !cur->sib) cur = cur->parent;
      if(cur) cur = cur->sib;
    }
  }
}

void Tree::InitOpt(double* a, int K, int Nordfuzz)
{
  Node* cur = root;
  double* p = a;
  while(cur)
  {
    cur->FixParams(&p,K, Nordfuzz);
    if(!cur->IsLeaf(K)) cur = cur->child;
    else {
      while(cur && !cur->sib) cur = cur->parent;
      if(cur) cur = cur->sib;
    }
  }
}

void Tree::CollectOpt(double* dyda, int K)
{
  Node* cur = root;
  double* p = dyda;
  while(cur)
  {
    cur->ColParams(&p,K);  
    if(!cur->IsLeaf(K)) cur = cur->child;
    else {
      while(cur && !cur->sib) cur = cur->parent;
      if(cur) cur = cur->sib;
    }
  }
}


static Data* _data;
static Tree* _tree;
static Node* _root;
static int _K;

double deriv1(double* a)
{
  _tree->InitOpt(a+1,_K, _data->GetNordfuzz()); //distribute params & set accumulators to 0
  int begin =  0;
  int end =_data->GetNTrain()- 1;
  return _data->Derivs(_root,_K,begin,end);
}


void derivs(double* a, double* da)
{
  _tree->InitOpt(a+1,_K,_data->GetNordfuzz()); //distribute params & set accumulators to 0
  int begin =  0;
  int end = _data->GetNTrain() - 1;
  /*GMM double deriv = */_data->Derivs(_root,_K,begin,end);
  if(da) _tree->CollectOpt(da+1,_K);
//  return deriv;
}



double Tree::OptimiseParams(Data* data)
{
  int n = CountParams(_Kfinal);
  double* a = new double[n];
  double* abest = new double[n];
  double* da = new double[n];
  /*UNIXextern*/ bool Gbflag;

  //
  //  Various Initializations
  //

  int NTotal = data->GetNTotal();
  int NTest = data->GetNTest();
  int NTrain = data->GetNTrain();
  int NSel = 0; //(int)(NTrain/3.0);  // Use NSel to select optimal fuzzification
  NTrain = NTrain - NSel;  // Optimize function only with NSel  
  // data->SetNTrain(NTrain);
  //
  //  Prepare data for Optimization
  //


  theSlope = -1000000.0;
  PrepareOpt(a,_Kfinal, data);

  _data = data;
  _root = root;
  _tree = this;
  _K = _Kfinal;

  memcpy(da,a,n*sizeof(double));


  /////////////////////////////////////////////////////////
  //  First calculate the Cost Function for a crisp tree
  //
  
  for(int i=0;i<n;i++)
  {
    if(da[i] < -100) 
      a[i] = -da[i];

//    (*console)("Hola1");
//1    // console->printf("%g\t",a[i]);
  }
//    (*console)("Hola1");
//1  console->printf("\n");
  memcpy(abest,a,n*sizeof(double));
  
  int begin,end;
  begin = NSel ? NTrain: 0;
  end = NSel ? NTrain+NSel-1: NTrain-1;
  _tree->InitOpt(a,_K,_data->GetNordfuzz()); 
    double fuzz = _data->DegreeOfFuzz(_root,_K) *  _tree->SizeBest  /  (_tree->SizeBest - 1); // Fix classes after reinitializing

  double bestSel = _data->Defuzzify(_root,_K,begin, end);
//  double bestSel = _data->Defuzzify(_root,_K,NTrain,NTrain+NSel-1); // Find Fuzzification scale
  double bestErrTest = _data->Defuzzify(_root,_K,NTotal-NTest,NTotal-1);
/*GMM  double factor =*/ _data->GetYGvscale() * _data->GetYGvscale();
//    (*console)("Hola1");
//1  console->printf("Crisp value of Cost Func = %g\n",deriv1(a-1)/(factor*Ndata));
//1  console->printf("bestSel= %g\n",bestSel);

  FILE * tempout = fopen("class_out.txt","a");
  fprintf(tempout,"%d \t %g \t", SizeBest,bestErrTest);

  bestSlope = -1.0;
  double B_estSlopeTest = -1.0;
  double B_estTest = bestErrTest;
  double bestFuzz=-1000, B_estFuzz=-1000;

  
  ///////////////////////////
  // ASG: FOR REGRESSION SEEMS TO WORK WELL WITH ONE
  // theSlope = 2.0/2.0;
  // int kmax = 1;

  theSlope = 0.0625/2.0;
  int kmax = 9;

  /*/ ASG:  Dialog Box for Optimization (TO REINTRODUCE)
  
  
  dbOptimize db;
  
  // Invoke the dialog box
  
  db.m_degree = ((NomData*)data)->degree;
  db.m_slope = 2.0* theSlope;
  db.m_kmax = kmax;
  db.DoModal();

  ((NomData*)data)->degree = db.m_degree;
  theSlope = db.m_slope/2.0;
  kmax = db.m_kmax;
  
  
  // */

  for(int k = 0; k < kmax; ++k)
  {
    //
    //    Fuzzification
    //
    
    theSlope *= 2.0;
    data->SetNTrain(NTrain);  // If NSel!= 0, NTrain is different for the tree 
                  // generation and for the optimization
    data->SortOn(data->GetNumVar()+5,0,NTrain+NSel-1);
    memcpy(a,da,n*sizeof(double));
    for(int i=0;i<n;i++)
    {
      if(da[i] < -100)
        a[i] = da[i]*(-theSlope/1000000.0);

 //   (*console)("Hola1");
//1      // console->printf("%g\t",a[i]);
    }
  //  (*console)("Hola1");
//1    console->printf("\n");

    _tree->InitOpt(a,_K,_data->GetNordfuzz()); //distribute params & set accumulators to 0
    fuzz = _data->DegreeOfFuzz(_root,_K) *  _tree->SizeBest  /  (_tree->SizeBest - 1); // Fix classes 
  
   // (*console)("Hola1");
//1    console->printf("Slope =  %g \t Degree of Fuzzification  = %g \n",theSlope,fuzz);
      

/*/ for each param, compare numeric && anayl derivs;
//  
  double*  dtemp= new double[n];
  derivs(a-1,dtemp-1);

  for( i=0;i<n;i++)
  {
    double eps =0.0000001;
    a[i] += eps;
    double ss = derivs(a-1,NULL);
    a[i] -= 2.0*eps;
    ss -= derivs(a-1,NULL);
    ss /= 2.0*eps;
    a[i] += eps;
//1    console->printf("i = %d a = %f da = %f da_num = %f\n",
      i,a[i],dtemp[i],ss);
  }
  delete[] dtemp;
// */  
    //
    //  PrintOuts
    //


//1    console->printf("Init Val Cost Func = %g\n",deriv1(a-1)/(factor*Ndata));
    
    //
    //  Optimization
    //
    
    /*double CostFunc = 0.0;
    int iter = 0;
    //UNIXfrprmn(a-1, n, 1e-6, &iter, &CostFunc,deriv1, derivs);*/

    _tree->InitOpt(a,_K,_data->GetNordfuzz()); //distribute params & set accumulators to 0
    fuzz = _data->DegreeOfFuzz(_root,_K) *  _tree->SizeBest  /  (_tree->SizeBest - 1); // Factor to normalize to 1 for completely fuzzy tree
    

    //
    //    Selection of best Fuzzification
    //

    int begin,end;
    begin = NSel ? NTrain: 0; 
    end = NSel ? NTrain+NSel-1: NTrain-1;
    double Sel = _data->Defuzzify(_root,_K,begin, end);
  //  double Sel = _data->Defuzzify(_root,_K,0,NTrain+NSel-1); // Find Fuzzification scale
    double errTest = _data->Defuzzify(_root,_K,NTotal-NTest,NTotal-1);
/*GMM    double noDefuzzErr =*/  _data->Error(_root,_K,NTotal-NTest,NTotal-1);
  
    if(bestSlope <0.0 || Sel < bestSel)  // Use test set error rate to determine best
    {
      bestFuzz = fuzz;
      bestSel = Sel;
      bestSlope = theSlope;
      bestErrTest = errTest;
      memcpy(abest,a,n*sizeof(double));
    }


    if(B_estSlopeTest < 0.0 || errTest< B_estTest)
    {
      B_estSlopeTest = theSlope;
      B_estTest = errTest;
      B_estFuzz = fuzz;

    }

  //
  //      VARIOUS PRINTOUTS
  //
  
    //(*console)("Hola1");
//1    console->printf("Iter= %d \t Cost Func= %g  \t SelFunc = %g  \t Fuzz"
//" = %g \n ErrTest(noDefuzz) = %g (%g) \t ErrTest (Defuzz) = %g (%g) \n\n",
//iter,CostFunc/(factor*Ndata), Sel, fuzz, noDefuzzErr, 1.0-noDefuzzErr, errTest, 
//1.0-errTest);

  }
 
   FILE * B_fp = fopen("Best.txt","a");
  fprintf(B_fp,"%g \t %g \t %g \t %g\t %g \t %g \n",bestSlope, bestFuzz,
   bestErrTest, B_estSlopeTest,B_estFuzz,B_estTest);
  fclose(B_fp);

  
  /*/////// ASG: TO DELETE  //////////////////////
  fprintf(out,"%g \t %g \t %g \t %g \n",bestSlope,bestTest,bestSlope);
  fflush(out);
  fclose(out);
  / */////////////////////////////////////////////

  //
  //  Optimize again with initial seed the best selection
  //  but with selection+training sets
  //
  _tree->InitOpt(abest,_K,_data->GetNordfuzz()); //distribute params & set accumulators to 0
  data->SetNTrain(NTrain+NSel);
//GMM  double CostFunc = 0.0;
//GMM  int iter = 0;
   // frprmn(abest-1, n, 1e-6, &iter, &CostFunc,deriv1, derivs);
  _tree->InitOpt(abest,_K,_data->GetNordfuzz()); //distribute params & set accumulators to 0

  //
  //
  //

  fuzz = _data->DegreeOfFuzz(_root,_K) *  _tree->SizeBest  /  (_tree->SizeBest - 1);

//    (*console)("Hola1");
//1  console->printf("Best Slope = %g \t BestSel = %g \n",bestSlope,bestSel);
  

  Gbflag = true;
  bestSel = _data->Defuzzify(_root,_K,0,NTrain+NSel-1);
  double tmpquartile[4];
  double errTest = _data->Defuzzify(_root,_K, NTotal - NTest, NTotal - 1, tmpquartile);
  Gbflag = false;
    //(*console)("Hola1");
//1  console->printf(" TestErr = %g (%g)\n\n",errTest,1-errTest);

  fprintf(tempout,"%g \t ",errTest);
  int countout;
  for(countout=0 ; countout<4 ; ++countout)
    fprintf(tempout,"%g \t ",tmpquartile[countout]);
  fprintf(tempout,"\n");
  fclose(tempout);

//    (*console)("Hola1");
//1  console->printf("Optimized tree:\n");
//  Write(data,_Kfinal);

  for(int jj=0;jj<n;jj++)
 //   (*console)("Hola1");
  //1    console->printf("%g\t",abest[jj]);

//    (*console)("Hola1");
//1  console->printf("\n");

    
  delete[] a;
  delete[] da;
  delete[] abest;

  return bestErrTest; 
  
}

void Tree::IndepPrune2()
{
  //TODO: Add your source code here

}

/*
  This function throws the data from begin to end throught the tree and
  calculates the node members first, last, Rleaf, RsubTree.
  returns the errors of the tree rooted at n.
*/
double Tree::ThrowData(Data * data, int begin, int end, int K, Node *n)
{
  Node *cur = n ? n : root;
  int leftLast=0;
  double nerror=0.0;
  int Nvar = data->GetNumVar();

  cur->first = begin;
  cur->last = end;
  cur->d->MembershipTotal = data->FindMembership(cur);
  int PrevClass = cur->iClass;
  data->FixLeaf(cur);
  cur->iClass = PrevClass;

  for(int i = begin; i <= end; i++ ) {
    int Class = (int)data->GetValueVar(i, Nvar-1);
    if (Class != cur->iClass) nerror += data->GetValueVar(i, Nvar+4);
  }
  cur->d->Rleaf = nerror;

  if (cur->IsLeaf(K)) {
/*    cur->first = begin;
    cur->last = end;
    cur->d->MembershipTotal = data->FindMembership(cur);
    int PrevClass = cur->iClass;
    data->FixLeaf(cur);
    cur->iClass = PrevClass;
    for(i=begin;i<=end;i++) {
      int Class = (int)data->GetValueVar(i, Nvar-1);
      if(Class != cur->iClass) nerror+=data->GetValueVar(i, Nvar+4);
    }
    cur->d->Rleaf = nerror;*/
    cur->d->RsubTree = nerror;
  }
  else /*if (begin<=(leftLast = data->GetLeftLast(cur, begin, end, K)))*/
  {
/*    cur->first = begin;
    cur->last = end;
    cur->d->MembershipTotal = data->FindMembership(cur);
    int PrevClass = cur->iClass;
    data->FixLeaf(cur);
    cur->iClass = PrevClass;
    for(i=begin;i<=end;i++) {
      int Class = (int)data->GetValueVar(i, Nvar-1);
      if(Class != cur->iClass) nerror+=data->GetValueVar(i, Nvar+4);
    }
    cur->d->Rleaf = nerror;*/
    leftLast = data->GetLeftLast(cur, begin, end, K);
    cur->d->RsubTree = nerror = ThrowData(data, begin, leftLast, K, cur->child) +
                           ThrowData(data, leftLast+1, end, K, cur->child->sib);
  }
  return nerror;
}

void Tree::Assign(Tree * t)
{
  //TODO: Add your source code here
  _Kfinal = t->_Kfinal;
  atdrawchosefathers = t->atdrawchosefathers;
  bestSlope = t->bestSlope;
  maxdepth = t->maxdepth;
  Ndata = t->Ndata;
  RMSres = t->RMSres;
  sbest = t->sbest;
  SizeBest = t->SizeBest;
  root = new Node(0, t->root->sizecoeff, 0);
  root->Assign(t->root);
}

void Tree::Guardar(std::ostream &salida, int version)
{
  if (version==0) {
    //salida.write((char*)this, 56);
    salida.write((char*)&sbest             , sizeof(double));
    salida.write((char*)&RMSres            , sizeof(double));
    salida.write((char*)&SizeBest          , sizeof(int));
    salida.write((char*)&bestSlope         , sizeof(double));
    salida.write((char*)&maxdepth          , sizeof(int));
    salida.write((char*)&root              , sizeof(int));
    salida.write((char*)&console           , sizeof(int));
    salida.write((char*)&_Kfinal           , sizeof(int));
    salida.write((char*)&Ndata             , sizeof(int));
    salida.write((char*)&atdrawchosefathers, sizeof(bool));
    salida.write("bool"                    , sizeof(int)-sizeof(bool));
    salida.write((char*)&Tmax              , sizeof(int));

    if (root) root->Guardar(salida, version);
  }
  else if (version==1) {
    salida.write((char*)&maxdepth          , sizeof(int));
    salida.write((char*)&root              , sizeof(int));
    salida.write((char*)&_Kfinal           , sizeof(int));
    
    if (root) {
      salida.write((char*)&root->sizecoeff , sizeof(int));
      root->Guardar(salida, version);
    }
  }

}
void Tree::Leer(std::istream &in, int version)
{
//  in >> sbest >> RMSres >> SizeBest >> bestSlope >> maxdepth >> _Kfinal;
//  in >> Ndata >> atdrawchosefathers >> (void*)root;
  if (version==0) {
    //in.read((char*)this, 56);//sizeof(Tree)/*-sizeof(Tmax)*/);
    int dummy;
    in.read((char*)&sbest             , sizeof(double));
    in.read((char*)&RMSres            , sizeof(double));
    in.read((char*)&SizeBest          , sizeof(int));
    in.read((char*)&bestSlope         , sizeof(double));
    in.read((char*)&maxdepth          , sizeof(int));
    in.read((char*)&root              , sizeof(int));
    in.read((char*)&dummy             , sizeof(int));
    in.read((char*)&_Kfinal           , sizeof(int));
    in.read((char*)&Ndata             , sizeof(int));
    in.read((char*)&atdrawchosefathers, sizeof(bool));
    in.read((char*)&dummy             , sizeof(int)-sizeof(bool));
    in.read((char*)&Tmax              , sizeof(int));

    if (root) {
      root = new Node(0, 0, 0);
      root->Leer(in, version);
    }
  }
  else if (version==1) {
    in.read((char*)&maxdepth          , sizeof(int));
    in.read((char*)&root              , sizeof(int));
    in.read((char*)&_Kfinal           , sizeof(int));

    if (root) {
      int aux;
      in.read((char*)&aux             , sizeof(int));
      root = new Node(0, aux, 0);
      root->Leer(in, version);
      Tmax = root->T;
    }
  }
}
//---------------------------------------------------------------------
#ifdef _DEBUG_TREE20_CPP

int main(int argc, char* argv[])
{
  NomData *dtr = new NomData((char*)"train.cre");

  NomData *dts = new NomData((char*)"test.cre");

  //MIBoosting *b = new MIBoosting();
  Tree *tree = new Tree();

  int ndatos = dtr->GetNTotal();
  dtr->SetNTrain(ndatos);

  tree->CART(dtr, 1);
  tree->FreeAuxData();

  double tot_peso = 0.0;
  for(int i=0; i<dts->GetNTotal();i++) tot_peso += dts->GetDatWeight(i);
  double error2 = dts->Error2(tree->GetRoot(), tree->GetK(), 0, dts->GetNTotal()-1 );

  cout << error2/tot_peso << endl;

  delete dtr;
  delete dts;
  delete tree;
}
#endif

