#include <math.h>
#include <string>
#include "policy20.h"
#include "data20.h"
#include "tree20.h"
//#include "..\CARTree\console.h"

extern "C" double gasdev();
extern "C" double ran3(long *idum);

static char* ConTable[20] = {  (char*)"and",
        (char*)"and.prod",
        (char*)"and.mean",
        (char*)"and.sum",
        (char*)"and.yager",
        (char*)"and.zimmerman",
        (char*)"or",
        (char*)"or.prod",
        (char*)"or.mean",
        (char*)"or.sum",
        (char*)"or.yager",
        (char*)"or.zimmerman"};
static char* HedgeTable[20] = {   (char*)"about",
          (char*)"above",
          (char*)"positively",
          (char*)"below",
          (char*)"around",
          (char*)"generally",
          (char*)"close",
          (char*)"not",
          (char*)"somewhat",
          (char*)"very",
          (char*)"extremely",
          (char*)"slightly"};



int GetConnector(char *buf, double& val);
double ConEval(double val1,double val2,int con,double con_val);



HedgedSet::HedgedSet(FuzzySet* fs, int* h, int n)
  : FuzzySet(fs)
{
  num = n;
  hedges = NULL;
  if(num) hedges = new int[num];
  for(int i=0;i<num;i++)
    hedges[i] = h[i];
  FS = fs;
  Generate();
}



// Allocate space 

DataSet::DataSet(int nvars, int ndata)
  : NumVars(nvars), NumData(ndata)
{
  
  VarNames = new std::string[NumVars];
  VarValue = new double*[NumVars];

  for (int i =0; i < NumVars; ++i)
    VarValue[i] = new double[NumData];    

  Defaults(ndata);
}

// Allocate space and initialize to specified values

DataSet::DataSet(const char* tag,  char** names, double** values,
         int nvars, int ndata)
  : label(tag)
{

  VarNames = new std::string[NumVars];
  VarValue = new double*[NumVars];

  for (int i =0; i < NumVars; ++i)
  {
    VarValue[i] = new double[NumData];  
    VarNames[i] = names[i];
    for(int j=0; j<NumData; ++j)
      VarValue[i][j] = values[i][j];
  }
  DepVar = VarNames[NumVars-1];

  Defaults(ndata);
}

// Allocates space and initializes only names
DataSet::DataSet(char** names, int nvars, int ndata)
  :NumVars(nvars), NumData(ndata)
{
  
  
  VarNames = new std::string[NumVars];
  VarValue = new double*[NumVars +1]; 
  VarValue[NumVars] = new double[NumData]; // ASG extra working space for glass.txt
  for (int i =0; i < NumVars; ++i)
  {  
    VarValue[i] = new double[NumData];  
    VarNames[i] = names[i];
  }
  
  DepVar = VarNames[NumVars-1];
  Defaults(ndata);
  
}

void DataSet::Defaults(int ndata)
{
  ponder = scramble = multree = false;
  FuzzVars = NULL;
  NumFuzz = NumNoms =0;
  NumVarsOrdSplit = NumOrds = NumVars-1;
  DepVarType = ORD;
  NPrune = NSel = 0;
  ncv = 10;
  if(ndata/2<100)
    nmin = (int)sqrt((double)(ndata/2));
  else nmin = 10* (int)(0.1*sqrt((double)(ndata/2)));

  ntest = ndata/2;
  
  multsplits = false;
  beta = 0.1;
  epsilon = 0.01;

  plm = deladd = false;
  MaxFToDelete = 0.1;
  MinFToAdd = 0.1;
  tol =0.01;

  fuzzify = false;
  
  degree = 1;

}

DataSet::~DataSet()
{  
  for (int i=0; i < NumVars; ++i)
    delete VarValue[i];

  delete[] VarValue; VarValue = NULL;
  delete[] VarNames; VarNames = NULL; 
}

void DataSet::AddTerm(int i, std::string newterm)
{
  NomTerms[i].push_back(newterm);
}

  
void DataSet::Exchange(int col1,int col2)
{
  if(col1 == col2) return;

  for(int i=0; i<NumData; ++i)  
  {
    double tempvalue = VarValue[col1][i];
    VarValue[col1][i] = VarValue[col2][i];
    VarValue[col2][i] = tempvalue;
  }

  std::string tempname;
  tempname = VarNames[col1]; 
  VarNames[col1] = VarNames[col2];
  VarNames[col2] =tempname;
  
}


Model::Model(const char *name): Name(name)
{
  Policies = NULL;
  Props = NULL;
  NumHedges = 0;
}

Model::~Model()
{
  Property *p, *pnext;
  p = Props; while(p){pnext = p->next;delete p; p=pnext;}
  Props = NULL;

  Policy *q, *qnext;
  q = Policies; while(q){qnext = q->Next;delete q; q=qnext;}
  Policies = NULL;
}

int Model::FindHedge(const char* buf)
{
  for(int i=0;i<_HEDGEMAX;i++)
    if(!strcmpci(HedgeTable[i],buf)) return i;
  return -1;
}

Property* Model::AddProperty(const char* name,int type,double domain[2], int dec)
{  

  Property* p = new Property(name,type,domain,Props);
  Props=p;
  if(dec > 0) p->Decorate(dec);
   return p;
}
void Model::AddPolicy(const char* name, const char* consequent, int defuz, int infer)
{
  Policy* temp = new Policy(name, consequent,defuz,infer);
  temp->Next = Policies;
  Policies = temp;
}

void Model::Read(FILE* is)
{
  char buf[255];
  fscanf(is,"Model: %[^\n]\n",buf);
  Name = buf;

  int count=0;
  fscanf(is,"Number of Policies: %d\n",&count);

  Policy *pcur=NULL, *pold=NULL;
  int i;
  for(i=0;i<count;i++)
  {
    pcur = new Policy("dummy","dummycon");
    pcur->Read(is);
    if(pold) pold->Next = pcur;
    else Policies = pcur;
    pold = pcur;
  }

  fscanf(is,"Number of Variables: %d",&count);

  Property *cur=NULL, *old=NULL;
  double d[2]; d[0] = 0; d[1] = 1;
  for(i=0;i<count;i++)
  {
    cur = new Property("dummy",0,d,NULL);
    cur->Read(is);
    if(old) old->next = cur;
    else Props = cur;
    old = cur;
  }

  pcur = Policies;
  while(pcur)
  {
    Rule* cur = pcur->Rules;
    while(cur){ cur->Parse(this, pcur->GetConsequent().c_str()); cur = cur->GetNext();}
    pcur = pcur->Next;
  }

}
void Model::Write(FILE* os)
{
  fprintf(os,"Model: %s\n",Name.c_str());

  int count = 0;
  Policy* pcur = Policies;
  while(pcur){count++; pcur = pcur->Next;}

  fprintf(os,"Number of Policies: %d\n",count);

  pcur = Policies;
  int i;
  for(i=0;i<count;i++,pcur = pcur->Next)
    pcur->Write(os);

  count=0;
  Property* cur = Props;
  while(cur){count++; cur = cur->next;}

  fprintf(os,"Number of Variables: %d\n",count);

  cur = Props;
  for(i=0;i<count;i++,cur = cur->next)
    cur->Write(os);
}


Policy::Policy(const char* name, const char* con, int defuz, int infer)
  : Work("work",-1), Name(name), Consequent(con)
{
  Next = NULL; Rules = NULL;
  DeFuz = defuz;
  Infer = infer;
}

Policy::~Policy()
{
  Rule *r, *rnext;
  r = Rules; while(r){rnext = r->GetNext();delete r; r=rnext;}
}

void Policy::Read(FILE* is)
{
  char buf[255];
  fscanf(is,"Policy: %[^\n]\n",buf);
  Name = buf;
  fscanf(is,"\tConsequent: %s\n",buf);
  Consequent = buf;
  fscanf(is,"\tReasoning: %d\n",&Infer);
  fscanf(is,"\tDeFuz: %d\n",&DeFuz);
  int count=0;
  fscanf(is,"Number of Rules: %d\n",&count);

  Rule *cur=NULL, *old = NULL;
  for(int i=0;i<count;i++)
  {
    cur = new Rule();
    cur->Read(is);
    if(old) old->SetNext(cur);
    else Rules = cur;
    old = cur;
  }

}
void Policy::Write(FILE* os)
{
  fprintf(os,"Policy: %s\n",Name.c_str());
  fprintf(os,"\tConsequent: %s\n",Consequent.c_str());
  fprintf(os,"\tReasoning: %d\n",Infer);
  fprintf(os,"\tDeFuz: %d\n",DeFuz);
  int count = 0;
  Rule* cur = Rules;
  while(cur){count++; cur = cur->GetNext();}

  fprintf(os,"Number of Rules: %d\n",count);

  cur = Rules;
  while(cur){cur->Write(os); cur = cur->GetNext();}

}

void Rule::Read(FILE* is)
{
  char buf1[64],buf[1024];
  fscanf(is,"Rule %[^:]: %[^$]$\n",buf1,buf);





  Name = buf1;
  Text = buf;
}
void Rule::Write(FILE* os)
{
  fprintf(os,"Rule %s: %s$\n",Name.c_str(),Text.c_str());
}



int Rule::Parse(Model* model, const char* consequent)
{
  char *work = (char*)malloc(strlen(Text.c_str())+1);
  strcpy(work, Text.c_str());
  // What type of rule?
  char *p = strtok(work," \r\n");
  if(!strcmpci(p,"if")) Type = _CONDITIONAL;
  else Type = _DECLARATIVE;

  p = strtok(NULL," \r\n"); // strip "IF"

  Statement* cur;
  root = new Statement();
  cur = root;
  Property* con = NULL;
  FuzzySet* Term = NULL;
  int h[20], pos,terms = 0;
  do
  {
    if( *p != '(')
    {
      cur->Prop = model->FindVariable(p); // This must be variable
      if(!cur->Prop) {goto _ERROR;} // variable not found

      p = strtok(NULL," \r\n"); // should be "is"
      if(strcmpci(p,"is")) {goto _ERROR;}

      Term = NULL;
      terms = 0;
      do {
        pos= -1;
        do {
          p = strtok(NULL," \r\n");  // This must be hedge or term
          pos++;
          h[pos] = model->FindHedge(p);
        } while(h[pos]>=0);

        Term = FindTerm(p,cur->Prop);
        if(Term) {cur->AddSet(new HedgedSet(Term,h,pos));terms++;}
      } while(Term);
      if(!terms || pos){goto _ERROR;}
    } else {
      cur->child = new Statement(); cur->child->parent = cur;
      cur = cur->child;
      p++;
      continue;
    }
    if(cur->parent)
    {
      do {
        if(*p == ')') {cur = cur->parent;p = strtok(NULL," \r\n");}  // parent finished
      } while(*p == ')');
    }
    if(cur->parent == NULL && !strcmpci(p,"then")) break; // get conclusion
    if(! *p)  break;//something wrong

    // get connector - only simple forms for now
    cur->connector = GetConnector(p,cur->con_val);
    if(cur->connector < 0) goto _ERROR;
    p = strtok(NULL," \r\n");
    // create sib and continue
    cur->sib = new Statement(); cur->sib->parent = cur->parent;
    cur = cur->sib;
  } while(cur && cur != root);

  p = strtok(NULL," \r\n");
  if( strcmpci(consequent,p))goto _ERROR; //error
  p = strtok(NULL," \r\n");
  if( strcmpci("is",p)) goto _ERROR; //error

  con  = model->FindVariable(consequent);
  Term = NULL;
 pos= -1;
  do {
    p = strtok(NULL," \r\n");  // This must be hedge or term
    pos++;
    h[pos] = model->FindHedge(p);
  } while(h[pos]>=0);
  Term = FindTerm(p,con);
  Consequent = new HedgedSet(Term,h,pos);

  if(work) delete work;
  return -1;
  _ERROR:
    if(root) delete root;
    if(work) delete work;
    return (int)(p-work);
}


int GetConnector(char *buf, double& val)
{
  char *p,*rem=NULL;
  int con = -1;

  p = buf;
  val = -1;
  while(*p)   // Has bracket - i.e., compesatory?
  {
    if(*p == '('){rem = p; *p = '\0'; break;}
    p++;
  }
  for(int i=0;i<_CONMAX_;i++)
    if(!strcmpci(buf,ConTable[i])){con = i; break;}
  if(rem) *rem = '(';
  if(con < 0) return -1;
  if(rem)
  {
    rem++;
    p = rem;
    while(*p && *p != ')')p++;
    if(!*p) return -1;
    *p = '\0';
    val = atof(rem);
    *p = ')';
  }
  return con;
}


FuzzySet* Rule::FindTerm(const char* buf, Property* cur)
{
  FuzzySet* term = cur->Types;
  while(term)
  {
    char st[64]; strcpy(st,term->GetName().c_str());
    if(!strcmpci(buf,term->GetName().c_str())) break;
    term = term->Next();
  }
  return term;
}



void Statement::Evaluate(int infer)
{
  value=0.0;

  if(Prop)  // simple statement
  {        // Note that for singleton input, necc = membership
    HedgedSet* hs = Sets;
    double val = Prop->GetVal();
    while(hs)    // Loop over Or'd possibilities
    {
      double f = hs->Membership(val);
      if(infer == _FES) f = 1.0-f;
      if(f > value) value = f;
      hs = (HedgedSet*) hs->Next();
    }
  } else {  // children are already evaluated
    Statement* cur = child;
    value = cur->value;
    while(cur->sib)
    {
      value = ConEval(value,cur->sib->value,cur->connector,cur->con_val);
      cur = cur->sib;
    }
  }
  return ;
}

double ConEval(double val1,double val2,int con,double con_val)
{

  switch(con)
  {
    case _AND:
      return std::min(val1,val2);
    case _AND_MEAN:
      return (val1+val2)/2.0;
    case _AND_PROD:
      return val1*val2;
    case _AND_SUM:
      return std::max((double) 0.0,val1+val2-1);
    case _AND_YAGER:
      return 1.0-std::min( 1.0,pow(pow(1-val1,con_val)+pow(1-val2,con_val),1.0/con_val));

    case _OR:
      return std::max(val1,val2);
    case _OR_MEAN:
      return (2*std::min(val1,val2)+4*std::max(val1,val2))/6.0;
    case _OR_PROD:
      return (val1+val2-val1*val2);
    case _OR_SUM:
      return std::min((double) 1.0,val1+val2);
    case _OR_YAGER:
      return std::min(1.0,pow(pow(val1,con_val)+pow(val2,con_val),1.0/con_val));
    case _AND_ZIM:
    case _OR_ZIM:
      return pow(val1*val2,1.0-con_val)*pow((1.0-(1.0-val1)*(1.0-val2)),con_val);
  };
  return 0.0;

}


void Rule::Evaluate(double temp[VECMAX], int infer)
{
  Statement *cur = root;

  do {
    while(cur->child) cur = cur->child;
    cur->Evaluate(infer);
    while(cur->parent && !cur->sib){cur = cur->parent; cur->Evaluate(infer);}
    cur = cur->sib;
  } while(cur);

  // Assume ands
  cur = root;
  double value = root->value;
  while(cur)
  {
    if(infer == _FLC && cur->value < value) value = cur->value;
    if(infer == _FES && cur->value > value) value = cur->value;
    cur = cur->sib;
  }

  double *vals = Consequent->GetValues();
  for(int i=0;i<VECMAX;i++)
  {
    if(infer == _FLC) temp[i] = (vals[i] > value ? value : vals[i]);
    else temp[i] = (vals[i]+value > 1) ? 1 : vals[i]+value;
  }
}
void Policy::EvaluateRule(Rule* r,int infer)
{
  if(infer < 0) infer = Infer;
  double temp[VECMAX];
  r->Evaluate(temp,infer); // evaluate rule
              // Composition
  double* vals = Work.GetValues();
  for(int i=0;i<VECMAX;i++)
  {
    if((infer == _FLC) && (temp[i] > vals[i]))
      vals[i] = temp[i];
    if((infer == _FES) && (temp[i] < vals[i]))
      vals[i] = temp[i];
  }
}

double Policy::DeFuzzify(int defuz)
{
  double d0 ,d1; Work.GetDomain(d0,d1);
  double* vals = Work.GetValues();
  if(defuz < 0) defuz = DeFuz;

  double value = 0.0, amin = 1.0;
  for(int i=0;i<VECMAX;i++) if(vals[i] < amin) amin = vals[i];
  if(defuz == _CENTROIDS)
  {
    double ret=0.0;
      d1 = (d1-d0)/VECMAX;
    for(int i=0;i<VECMAX;i++) {value += (d0+d1*i)*(vals[i]-amin); ret += (vals[i]-amin);}
    value /= ret;
  }
  if(defuz == _COMPOS_MOM)
    value = d0+d1/2;
  return value;
}



/*double Policy::Evaluate(Property* Props)
{
  Rule* r = Rules;
  double ret=0.0;
  double temp[VECMAX];
  double* vals = Work.GetValues();
  for(int i=0;i<VECMAX;i++) vals[i] = 0.0;

  while(r)
  {
    r->Evaluate(temp,Infer);
      // Composition
    for(i=0;i<VECMAX;i++) if(temp[i] > vals[i]) vals[i] = temp[i];
  }

  // de-fuzzify: Centroid for now
  ret = value = 0;
  Property* w = Rules->FindVariable(Consequent,Props);
  Work.SetDomain(w->Domain[0],w->Domain[1]);
  double d0 = w->Domain[0],d1= w->Domain[1]-w->Domain[0];
  for(i=0;i<VECMAX;i++) {value += (d0+d1*i)*vals[i]; ret += vals[i];}
  value /= ret;
  return value;
}
*/
void Policy::DeleteRule(Rule *r)
{
  Rule* cur = Rules;

  if(r == cur){ Rules = cur->GetNext(); delete cur;}
  else {
    while(r != cur->GetNext() && cur->GetNext()) cur = cur->GetNext();
    if(!cur) return;
    cur->SetNext(r->GetNext());
    delete r;
  }
}

void HedgedSet::Generate()
{
  int i;
  for(i=0;i<VECMAX;i++) _MemVec[i] = FS->GetValues()[i];


  for(i=0;i<num;i++)
  {
    double EXPVal = -1,x;
    int j,pos=-1000;

    switch(hedges[i])
    {
      case  _NOT:
        for(j=0;j<VECMAX;j++) _MemVec[j] =  1.0-_MemVec[j];
        break;
      case  _POSITIVE:
        for(j=0;j<VECMAX;j++)
        {
          if(_MemVec[j] > 0.5)
            _MemVec[j] =  std::min(2*_MemVec[j]*_MemVec[j],(double)1.0);
          else 
            _MemVec[j] =  std::max(1.0-2*(1-_MemVec[j]*_MemVec[j]),0.0);
        }
        break;

      case  _GENERALLY:
        for(j=0;j<VECMAX;j++)
        {
          if(_MemVec[j] > 0.5)_MemVec[j] =  0.8*sqrt(fabs(_MemVec[j]));
          else _MemVec[j] =  1.0-0.8*(1.0-sqrt(fabs(_MemVec[j])));
        }
        break;
      case  _SOMEWHAT:
        for(j=0;j<VECMAX;j++) _MemVec[j] =  sqrt(fabs(_MemVec[j]));
        break;
      case  _VERY:
        for(j=0;j<VECMAX;j++) _MemVec[j] =  _MemVec[j]*_MemVec[j];
        break;
      case  _EXTREMELY:
        for(j=0;j<VECMAX;j++) _MemVec[j] =  _MemVec[j]*_MemVec[j]*_MemVec[j];
        break;
      case  _SLIGHTLY:
        for(j=0;j<VECMAX;j++) _MemVec[j] =  sqrt(fabs(_MemVec[j]));
        break;
      case _ABOUT:
        EXPVal = 2;
        break;
      case  _CLOSE:
        EXPVal = 4;
        break;
      case  _AROUND:
        EXPVal = 1.2;
        break;

      case  _ABOVE:
        for(j=0,x=0;j<VECMAX;j++)
          if(_MemVec[j]  > x)
            {pos = j; x = _MemVec[j]; if(x == 1.0) j = VECMAX+1;}
        for(j=0;j<VECMAX;j++)
          { j <= pos ? _MemVec[j] = 0 : _MemVec[j] = 1.0-_MemVec[j];}
        break;
      case  _BELOW:
        for(j=0,x=0;j<VECMAX;j++)
          if(_MemVec[j]  > x)
            {pos = j; x = _MemVec[j]; if(x == 1.0) j = VECMAX+1;}
        for(j=0;j<VECMAX;j++)
          { j >= pos ? _MemVec[j] = 0 : _MemVec[j] = 1.0-_MemVec[j];}
        break;
    };

    if(EXPVal > 0.0)
    {
      double local[VECMAX];
      int l;
      for(l=0;l<VECMAX;l++) local[l] = _MemVec[l];
      for(l=0;l<VECMAX;l++)
        if(_MemVec[l] != 0)
        {
          double x = i;
          for(int k=0;k<VECMAX;k++)
          {
            double p1 = fabs(double(k)-x)*(_Domain[1]-_Domain[0])/VECMAX;
            double p2 = (1.0/(1.0+pow(p1,EXPVal)))*_MemVec[l];
            local[k] = std::max(local[k],p2);
          }
        }
      for(l=0;l<VECMAX;l++) _MemVec[l] = local[l];
    }
  }
  return ;
}

DTPolicy::DTPolicy(/*1CConsole* pt*/void (*_console)(char *))
  :Policy("","",0,0)
{
  data =NULL;
  tree = new Tree(/*1pt*/_console);
}

DTPolicy::DTPolicy(/*1CConsole* pt,*/void (*_console)(char *), DataSet *ds,int defuz , int infer )
  :Policy(ds->GetLabel().c_str(),ds->GetNameVar(ds->GetNumVars()-1).c_str(),defuz,infer)
{

  switch(ds->GetTypeDepVar())
  {

    case NOM: data = new NomData(ds/*1,pt*/,_console); break;
//    case ORD:  if (ds->plm) data = new RPLMData(ds/*1,pt*/,_console);
//          else  data = new OrdData(ds/*1,pt*/,_console);
//          break;
    case FUZZ: ;  // TODO ASG
    default:
      break;
  }

  tree = new Tree(/*1pt*/_console);

  (*_console)((char*)"Tree and data created \n");
//  pt->printf("Tree and data created \n");

}

DTPolicy::~DTPolicy()
{
  if(data) delete data;
  delete tree;
}

void DTPolicy::BuildTree()
{
  if(!(data->GetDepVarType() == NOM && ((NomData*)data)->multree))
    tree->CART(data,Nmin);
  else
  {
    NomData* dt = (NomData*)data;
    int Ntot = dt->GetTotalData();
    double* Class = new double[Ntot];
    dt->Filter(-1,Class); // store classes in Class
    Tree **aux = new Tree*[dt->NumClass];
    double** IsItemOfClass= new double*[dt->NumClass];

    FILE* fp = fopen("multinom.txt","a");
    int i;
    for (i=0; i < dt->NumClass; ++i)
    {

      dt->Filter(i,Class);
      aux[i] = new Tree(/*1console*/console);
      aux[i]->CART(dt,Nmin);
      dt->AssignClass(i,aux[i]->GetRoot(),aux[i]->GetK(),IsItemOfClass[i],0,Ntot);
      char name[256];
      sprintf(name,"X%d",i);
      fprintf(fp,"%s\t",name);
    }
    fprintf(fp,"Y\n");
    dt->SortOn(0,0,Ntot);
    for(i=0; i<dt->GetTotalData(); ++i)
    {
      for(int j =0; j<dt->NumClass; ++j)
        fprintf(fp,"%d\t",IsItemOfClass[i]);

      fprintf(fp,"%d\n",(int)Class[i]);
      delete IsItemOfClass[i];
    }
    fclose(fp);

    delete[] IsItemOfClass;

  }
}


void BuildTreeThread(DTPolicy* dtp)
{
  if(dtp->GetData())
  {
    dtp->BuildTree();
  }
  else
  {
//1    CConsole* pt = dtp->GetTree()->GetConsole();
    int nvars = 3;
    int ndata = 200;
    int ntest =100;
//    int NSel =0;

    char** names;
    double** values;
    values = new double*[nvars];
    names = new char*[nvars];
    int i;
    for(i =0; i<nvars; ++i)
    {
      values[i] = new double[ndata];
      names[i] = new char[255];
      sprintf(names[i],"x%d",i+1);
    }
    sprintf(names[nvars-1],"y");

    DataSet* SelectedSet = new DataSet(names,nvars,ndata);
    SelectedSet->SetDepVar("y");
    SelectedSet->SetDepVarType(ORD);
    SelectedSet->SetNumOrds(nvars-1);
    SelectedSet->SetNumNoms(0);
    SelectedSet->SetNumFuzz(0);
    SelectedSet->SetNumVarsOrdSplit(nvars-1);
    SelectedSet->SetScramble(false);
    SelectedSet->NPrune =0;
    SelectedSet->ncv = 10;
    SelectedSet->nmin = (int)sqrt((double)(ndata-ntest));
    SelectedSet->SetNumData(ndata);
    SelectedSet->ntest = ntest;
    SelectedSet->NSel =0;
    SelectedSet->plm = false;
    SelectedSet->multsplits = true;
    SelectedSet->beta = 0.1;
    SelectedSet->epsilon =0.01;
    SelectedSet->deladd = true;
    SelectedSet->MaxFToDelete =0.01;
    SelectedSet->MinFToAdd =0.01;
    SelectedSet->fuzzify  = false;
    FILE* out = fopen("1_PLM.txt","a");
    if(out){
    //  fprintf(out,"multiple splits  %d \t  PLM %d\n",SelectedSet->multsplits,SelectedSet->plm);
    //  if(SelectedSet->plm)
    //    fprintf(out,"MaxFToDelete %g \t MinFToAdd  %g \n",SelectedSet->MaxFToDelete,SelectedSet->MinFToAdd);
    //  if(SelectedSet->multsplits)
    //    fprintf(out,"Beta %g \n",SelectedSet->beta);
    //  fprintf(out, "Size \t bestSlope  \t  RMS/VarY  \t  RMS/VarNoise \n");
    }
    for(int ii=0;ii<1; ++ii)
    {
//1       pt->printf("ITERATION # %d\n",ii+1);
      double sumNoise=0.0;
      double sumNoise2=0.0;
      double sumY=0.0;
      double sumY2=0.0;
      FILE* out2 = fopen("temp.txt","w");
      for(i = 0; i<ndata; ++i)
      {
        static long dum = -58;

    //    double a = ran3(&dum);
        double a = 2.0*ran3(&dum)-1.0;
    //    double b = ran3(&dum);
        double Noise =0.05*gasdev();

          // pi with 20 digits = 3.14159265358979323846
//        double pi2=2.0* 3.14159265358979323846;
        values[0][i] = a;
        values[1][i] = a*a;
      //  values[2][i] = cos(a*a+b*b);
      //  values[3][i] = a*a;
      //  values[4][i] = a*a*a;
        values[nvars-1][i] = sqrt(1-0.5*a*a*(1+a*a))+ Noise;

        for(int kk=0; kk < nvars; ++kk)
          if(out2) fprintf(out2,"%g\t",values[kk][i]);
        if(out2) fprintf(out2,"%g\n",values[nvars-1][i] -Noise);

        if (i>= ndata-ntest)
        {
          sumNoise += Noise;
          sumNoise2 += Noise*Noise;
          sumY += values[nvars-1][i];
          sumY2 += values[nvars-1][i]*values[nvars-1][i];
        }

      }
      if(out2) fclose(out2);
      

      double RMSNoise = sqrt((sumNoise2 - sumNoise*sumNoise/ntest)/(ntest-1));
      double RMSY = sqrt((sumY2-sumY*sumY/ntest)/(ntest-1));
      FILE* B_fp = fopen("Best.txt","a");
      fprintf(B_fp,"%g \t %g\t",RMSNoise,RMSY);
      fclose(B_fp);

      
    
      SelectedSet->InitialiseValues(values);  

      dtp = new DTPolicy(/*1pt,*/0, SelectedSet);
      dtp->SetNmin(SelectedSet->nmin);
      dtp->BuildTree();
      double err = dtp->GetTree()->sbest/dtp->GetData()->GetYGvscale();
//1      pt->printf("err/VarY = %g \t err/VarNoise = %g \n", err/RMSY, err/RMSNoise);
      if(out)  fprintf(out," %d \t %g \t %g \t", dtp->GetTree()->SizeBest,err/RMSY, err/RMSNoise);
      err = dtp->OptimiseParams();
    
      if(out)  fprintf(out," %g \t %g \n",err/RMSY, err/RMSNoise);
      fflush(out);
//1      pt->printf("Optimization leads to: err/VarY = %g \t err/VarNoise = %g \n", err/RMSY, err/RMSNoise);
      delete dtp;
    }

    if(out) fclose(out);
    delete SelectedSet;
    for(i =0; i<nvars; ++i)
    {
      delete[] values[i];
      delete[] names[i];
    }
    delete names;
    delete values;
    
  }
}


void Model::GenerateCART(DataSet* ds/*1, CConsole* pt*/, int defuz,int infer)
{  
  // ASG Comment for average
  DTPolicy* temp;
//1  pt->printf("");
  if(ds)
  {  ds->SetFuzzVars(Props); 
    temp = new DTPolicy(/*1pt,*/0, ds);
    temp->SetNmin(ds->nmin); 
    temp->Next = Policies;  
    Policies = temp;
  }
  else temp = new DTPolicy(/*1pt*/ 0);
/*
  /////////////////////////////////////////
  //
  //    Begin a Thread
  //

  DWORD dwThreadId; 
    HANDLE hThread; 
    hThread = CreateThread( 
        NULL,                        // no security attributes        
        0,                           // use default stack size         
        (LPTHREAD_START_ROUTINE) BuildTreeThread, // thread function        
        temp,                // argument to thread function   
        0,                           // use default creation flags     
        &dwThreadId);                // returns the thread identifier  
 
    // Check the return value for success. 
 
    if (hThread == NULL) 
        printf("CreateThread error\n");
  

    


  //
  //    End Thread
  //
  ///////////////////////////////////////
*/
//  temp->BuildTree();
  BuildTreeThread(temp);

}

Property* Model::FindVariable(const char* buf)
{
  Property* cur = Props;
  while(cur)
  {
    if(!strcmpci(buf,cur->Name.c_str())) break;
    cur = cur->next;
  }

  return cur;
}

void DataSet::InitialiseValues(double** values)
{
  for (int i=0; i<NumVars; ++i)
    for(int j=0; j<NumData; ++j)
      VarValue[i][j] = values[i][j];
}



VarType DataSet::VarTypeOf(const std::string cs, int* loc)
{

  int i;
  for(i=0; i< NumVars; ++i)
  {
    if(cs == VarNames[i])  break;

  }

  *loc = i;

  if(i == NumVars -1) return DepVarType;
  if(i<NumOrds) return ORD;
  i -= NumOrds;
  if(i< NumFuzz) return FUZZ;
  i -= NumFuzz;
  if(i< NumNoms) return NOM;
  return ERR;

}

void DataSet::InitializeNomTerms(int n)
{
  NomTerms = new std::vector<std::string>[n];

}

/*
void DataSet::AddFuzzVar(Property* tempp)
{
  Property* p = new Property(tempp->Name,tempp->Type,tempp->Domain,FuzzVars);
  p->Types = tempp->Types;
  FuzzVars = p;
}
*/

void DataSet::InitFuzzVars()
{
  while(FuzzVars)
  {
    Property* temp = FuzzVars;
    FuzzVars = temp->next;
    delete temp;
  }
}

void Model::OptimiseParams()
{
  DTPolicy* temp = (DTPolicy*) Policies;
/*  double err = */temp->OptimiseParams();
}

double DTPolicy::OptimiseParams()
{
  if(this) return  tree->OptimiseParams(data);
  return -1;
}


double gasdev()
{
  static long idum = -10;
  double ran3(long *idum);
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if  (iset == 0) {
    do {
      v1=2.0*ran3(&idum)-1.0;
      v2=2.0*ran3(&idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software *!^3,z3+1. */
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *idum)
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;

  if (*idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software *!^3,z3+1. */




void DataSet::LabelData(int nlabels)
{

  for(int i=0; i<NumData ; ++i)
  {
    static long dum = -51;
    VarValue[NumVars][i] = nlabels*ran3(&dum);
  }
}

int DataSet::ExtractSet(int label)
{
  double** m = new double*[NumData];
    
    
  int begin = 0;
  int end = NumData;
  int i;
  for(i =0; i < NumData; ++i)
  {  
    
    int index = (((int)VarValue[NumVars][i] == label) ? --end : begin++);
  
    m[index] = new double[NumVars];
    for(int j =0; j<NumVars; ++j)
      m[index][j] = VarValue[j][i];
  }

  if(end != begin)
    printf("Warning: Wrong count in ExtractSet");

  for(i =0; i < NumData; ++i)
  {  
    for(int j =0; j<NumVars; ++j)
      VarValue[j][i] = m[i][j];
    delete[] m[i];
  }

  delete[] m;

  return NumData-end;



}

void DataSet::GenerateLED()
{
  ;
  /*
  int nvars = 8;
  int ndata = 5200;
  int nord = 7;
  int nnom = 0;
  int nfuzz = 0;
  double** values;
  VarValue = new double*[nvars];
  VarNames =  new std::string[nvars];
  int* aux[10];
  aux[0][] = {1,1,1,0,1,1,1};
  aux[1][] = {0,0,1,0,0,1,0};
  aux[2][] = {1,0,1,1,1,0,1};
  aux[3][] = {1,0,1,1,0,1,1};
  aux[4][] = {0,1,1,1,0,1,0};
  aux[5][] = {1,1,0,1,0,1,1};
  aux[6][] = {1,1,0,1,1,1,1};
  aux[7][] = {1,0,1,0,0,1,0};
  aux[8][] = {1,1,1,1,1,1,1};
  aux[9][] = {1,1,1,0,1,1,1};
  
  for(int i =0; i < nvars; ++i)
  {
    char* str;
    sprintf(str,"X%d",i);
    VarNames[i] = str;
    VarValue[i] = new double[ndata];

    for(int j =0;j<ndata; ++j)
    {
      i
      if(i ==nvars-1)
        VarValue[i][j] = j/10;
      else
      {
        BOOL temp = (BOOL)(aux[i][j])
        Var[i][j] = drand()>0.9 ?  temp : !temp;
      }
  }
  names[nvars-1] ="class";

  */
}

