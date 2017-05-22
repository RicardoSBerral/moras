#include "SpaceCubes.h"
#include "GPTree.h"
#include "tree20.h"
#include "node20.h"
#include "data20.h"
#include <limits>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;

//---------------------------------------------------------------------------
class compare_tag {
public:
  int operator()(HyperCube* const& hc1, HyperCube* const& hc2) {
    return hc1->tag<hc2->tag;
  }
};
class compare_var_min {
  int ivar;
public:
  compare_var_min(int ivar){this->ivar=ivar;}

  int operator()(HyperCube* const& hc1, HyperCube* const& hc2) {
    return hc1->min(ivar)<hc2->min(ivar);
  }
};
class compare_var_max {
  int ivar;
public:
  compare_var_max(int ivar){this->ivar=ivar;}

  int operator()(HyperCube* const& hc1, HyperCube* const& hc2) {
    return hc1->max(ivar)<hc2->max(ivar);
  }
};
std::ostream& operator << (std::ostream &out, ClassTag tag)
{

  for(int i=0;i<tag.size();i++) {
    out << tag[i] << " ";
  }

  return out;
}
//---------------------------------------------------------------------------
//-------------------------------------------------------  HyperCube  -----
//---------------------------------------------------------------------------
HyperCube::HyperCube(std::vector<double> mins, std::vector<double> maxs, 
                                                                  ClassTag tag)
{
  for(unsigned i=0;i<mins.size();i++) {
    this->mins.push_back(mins[i]);
    this->maxs.push_back(maxs[i]);
  }
//  this->maxs = maxs;
//  this->mins = mins;
  this->intouch.resize(2*dimension());
  this->tag = tag;
  volume = vol();
}
HyperCube::HyperCube(int d, int tag)
{
  this->maxs.resize(d, numeric_limits<double>::infinity());
  this->mins.resize(d, -numeric_limits<double>::infinity());
  this->intouch.resize(2*d);
  this->tag = tag;
  volume = numeric_limits<double>::infinity();
}
//---------------------------------------------------------------------------
HyperCube* HyperCube::split(int d, double val)
{
  if (mins[d]>val ||maxs[d]<val) return 0;

  std::vector<double> mins2 = this->mins;
  std::vector<double> maxs2 = this->maxs;
 
  volume = volume*(val-mins[d])/(maxs[d]-mins[d]);

  mins2[d] = val;
  this->maxs[d] = val;

  HyperCube *ret = new HyperCube(mins2, maxs2, tag);
/*  ret->intouch = intouch;
  ret->intouch[d].resize(1);
  ret->intouch[d][0] = this;
  this->intouch[d+dimension()].resize(1);
  this->intouch[d+dimension()][0] = ret;*/

  return ret; 
}
//---------------------------------------------------------------------------
bool HyperCube::Merge(HyperCube *hc)
{
  if (tag != hc->tag) return false;

  int  v = -1;
  for(int i=0;i<dimension();i++) {
    if (mins[i]!=hc->mins[i] || maxs[i]!=hc->maxs[i]) {
      if (v!=-1) return false;
      else v = i;
    }
  }

  if (v==-1) return true;

  if (mins[v]==hc->maxs[v])
    mins[v] = hc->mins[v];
  else if (maxs[v]==hc->mins[v])
    maxs[v] = hc->maxs[v];
  else 
    return false;

  volume = volume + hc->volume;

  return true;
}
//---------------------------------------------------------------------------
/*bool HyperCube::touching(HyperCube *hc)
{

  for(int i=0;i<dimension();i++) {
    double min = mins[i] > hc->mins[i] ? mins[i] : hc->mins[i];
    double max = maxs[i] < hc->maxs[i] ? maxs[i] : hc->maxs[i];
    if (min>max) {
      return false;
    }
  }

  return true;
}*/
//---------------------------------------------------------------------------
HyperCube *HyperCube::intersec(HyperCube *hc)
{
  HyperCube *r = new HyperCube(dimension());

  for(int i=0;i<dimension();i++) {
    r->mins[i] = mins[i] > hc->mins[i] ? mins[i] : hc->mins[i];
    r->maxs[i] = maxs[i] < hc->maxs[i] ? maxs[i] : hc->maxs[i];
    if (r->mins[i]>=r->maxs[i]) {
      delete r;
      return 0;
    }
  }

  r->volume = r->vol();

  return r;
}
//---------------------------------------------------------------------------
string HyperCube::Print()
{
  std::ostringstream puf;
  int i;

  puf.precision(15);
  puf << "<" ;
  for(i=0;i<dimension()-1;i++) {
    puf << setw(10) << mins[i];
    puf << ", ";
  }
  puf << setw(10) << mins[i];
  puf << "> <";
  for(i=0;i<dimension()-1;i++) {
    puf << setw(10) << maxs[i];
    puf << ", ";
  }
  puf << setw(10) << maxs[i];
  puf  << "> " <<  tag;

  return puf.str();
}
//---------------------------------------------------------------------------
double HyperCube::vol()
{
  double v = 1.0;

  for(int i=0;i<dimension();i++)
    v *=arista(i);

  return v;
}
//---------------------------------------------------------------------------
//------------------------------------------------------  SpaceCubes  -----
//---------------------------------------------------------------------------
SpaceCubes::SpaceCubes(int d)
{
  cubos.push_back(new HyperCube(d));
}
SpaceCubes::SpaceCubes(DecisionTree *tree, vector<double>mins, 
                                                           vector<double> maxs)
{
  for(unsigned i=0;i<mins.size();i++) {
    this->mins.push_back(mins[i]);
    this->maxs.push_back(maxs[i]);
  }
  cubos = TreeToHyperCubes(tree, new HyperCube(mins, maxs));
}
SpaceCubes::SpaceCubes(Tree *tree, vector<double>mins, vector<double> maxs)
{
  for(unsigned i=0;i<mins.size();i++) {
    this->mins.push_back(mins[i]);
    this->maxs.push_back(maxs[i]);
  }
  cubos = TreeToHyperCubes(tree, -1, new HyperCube(mins, maxs));
}
//---------------------------------------------------------------------------
SpaceCubes SpaceCubes::operator &(SpaceCubes other)
{
  SpaceCubes sc(dimension());

  for(int i=0;i<dimension();i++) {
    sc.mins.push_back(mins[i]<other.mins[i] ? mins[i] : other.mins[i]);
    sc.maxs.push_back(maxs[i]>other.maxs[i] ? maxs[i] : other.maxs[i]);
  }

  sc.cubos.pop_back();
  for(int i=0;i<size();i++) {
    for(int j=0;j<other.size();j++) {
      HyperCube *hc = cubos[i]->intersec(other.cubos[j]);
      if (hc) {
        hc->tag = cubos[i]->tag + other.cubos[j]->tag;
        sc.cubos.push_back(hc);
      }
    }
  }

  return sc;
}
//---------------------------------------------------------------------------
string SpaceCubes::Print()
{
  std::ostringstream puf;

  puf << "-----------------------------------------------------" << endl;
  for(unsigned i=0;i<cubos.size();i++) {
    puf << setw(3) << i+1;
    puf << cubos[i]->Print() << endl;
  }

  return puf.str();
}
//---------------------------------------------------------------------------
vector<HyperCube*>::iterator  SpaceCubes::Merge(
            vector<HyperCube*>::iterator ini, vector<HyperCube*>::iterator fin)
{
  vector<HyperCube*>::iterator i, j; 

  for(i=ini; i != fin; i++) {
    for(j=i+1; j != fin ;j++) {
      if ((*i)->Merge(*j)) {
        j = cubos.erase(j);
        fin--;
        if (j==fin) break;
      }
    }
  }

  return fin;
}
//---------------------------------------------------------------------------
SpaceCubes* SpaceCubes::split(int d, double val)
{
  if (mins[d]>val ||maxs[d]<val) return 0;

  SpaceCubes *ret = new SpaceCubes(dimension());
  ret->mins = this->mins;
  ret->maxs = this->maxs;
  ret->mins[d] = val;
  this->maxs[d] = val;
  ret->cubos.clear();

  vector<HyperCube*> remaining;
  vector<HyperCube*>::iterator ihc;
  for(ihc=cubos.begin();ihc!=cubos.end();ihc++) {
    if ((*ihc)->min(d)>=val) {
      ret->cubos.push_back(*ihc);
    }
    else if ((*ihc)->max(d)<=val) {
      remaining.push_back(*ihc);
    }
    else {
      ret->cubos.push_back((*ihc)->split(d, val));
      remaining.push_back(*ihc);
    }
  }
  this->cubos = remaining;

  return ret;
}
//---------------------------------------------------------------------------
int SpaceCubes::Merge()
{
  vector<HyperCube*>::iterator i, j;
  int elems = size();
 
  stable_sort(cubos.begin(), cubos.end(), compare_tag());
  for(i=cubos.begin(); i != cubos.end(); i++) {
    vector<HyperCube*>::iterator ini = i;
    while (i+1!=cubos.end() && (*i)->tag==(*(i+1))->tag) i++;
    i = Merge(ini, i+1);
    i--;
  }

  return elems - size();
}
//---------------------------------------------------------------------------
Data* SpaceCubes::ToData()
{
  NomData *d;
  int id, dummy_id, dummy_total;
  bool valid_example;
  //double inf_neg = 0.0;//-numeric_limits<double>::infinity(); 
  //double inf_pos = 1.0;//numeric_limits<double>::infinity();

  int d2 = 1<<dimension();
  d = new NomData(2, dimension()+1);

  //Busco el cubo con menor arista
  double larista = numeric_limits<double>::infinity(); 
  for(int i=0;i<size();i++) { 
    for(int k=0;k<dimension();k++) {
      if (cubos[i]->arista(k)<larista) {
        larista = cubos[i]->arista(k);
      }
    }
  }
  larista /= 2.0;
 
  double cl;
  char ncl[100];
  dummy_id = id = 0;
  dummy_total = 2;
  ncl[0] = '0'; ncl[1] = '\0';  d->AddClass(string(ncl));
  ncl[0] = '1'; ncl[1] = '\0';  d->AddClass(string(ncl));
  for(int i=0;i<size();i++) {
    //int previd = id;
    HyperCube *c = cubos[i];
    sprintf(ncl, "%d", int(c->tag));
    cl = 0.5 + d->AddClass(string(ncl));
/*    for(int j=0;j<d2;j++) { //Aristas
      valid_example = true;
      for(int k=0;k<dimension();k++) {
        double v = j & (1<<k) ? c->max(k)-larista : c->min(k)+larista;
        if (v==inf_neg || v==inf_pos) {
//          valid_example = false;
          v = v==inf_neg ? 0.0 : 1.0;
//          break;
        }
        d->SetValueVar(id, k, v);
      }
      d->SetValueVar(id, dimension(), cl); //clase
      if (valid_example) id++; 
    }*/
 
    for(int k=0;k<dimension();k++) { //Bolardos en los lados para inundaciones
      double holdMin = c->min(k);
      double holdMax = c->max(k);
      c->mins[k] = mins[k];
      c->maxs[k] = maxs[k];
      for(int ii=i+1;ii<size();ii++) {
        HyperCube *in;
        if (c->tag==cubos[ii]->tag ) continue;
        if (!(in=c->intersec(cubos[ii]))) continue;
        for(int j=0;j<d2;j++) {
          for(int kk=0;kk<dimension();kk++) {
            double v, v2;
            v = v2 = j & (1<<kk) ? in->max(kk)-larista : in->min(kk)+larista;
            //  double v = in->centre(kk);
            if (k==kk) {
              v = v > holdMax ? holdMax - larista : holdMin + larista;
              v2 = cubos[ii]->centre(k) > holdMax ? 
                    cubos[ii]->min(k) + larista : cubos[ii]->max(k) - larista;
            }
            if (v==mins[kk] || v==maxs[kk]) {
              valid_example = false;
              break;
            }
            d->SetValueVar(id, kk, v);
            d->SetValueVar(id+1, kk, v2);
          }
          if (valid_example) {
            d->SetValueVar(id, dimension(), cl); //clase
            d->SetValueVar(id+1, dimension(), 0.5 + int(cubos[ii]->tag));//clase
//            id+=2; 
            dummy_id+=2;
            if (dummy_id>=dummy_total/*d->GetNTotal()*/) {
              dummy_total *= 2;//d->redim(d->GetNTotal()*2);
              printf("Resizing to %d\n", dummy_total); //d->GetNTotal());
            }
          }
        }
      }
      c->mins[k] = holdMin;
      c->maxs[k] = holdMax;
    }
//    if (previd<id) id = uniq(d, 0, id-1);
  }

  d->redim(id);

  return d;
}
//---------------------------------------------------------------------------
Tree* SpaceCubes::ToTree()
{
  SpaceCubesTree *sct = new SpaceCubesTree(this);

  return sct;
}
//---------------------------------------------------------------------------
//--------------------------------------------------  SpaceCubesTree  -----
//---------------------------------------------------------------------------
SpaceCubesTree::SpaceCubesTree(SpaceCubes *sc)
{
  vector<HyperCube*> c = sc->cubos;
  vector<Split*> splits;

  //Obtain all possible splits and number of classes
  for(unsigned i=0;i<sc->cubos.size();i++) {
    classes[c[i]->tag] = i; 
  }
  map<ClassTag, int>::iterator it=classes.begin();
  for(int i=0;it!=classes.end();it++, i++)
    classes[it->first] = i;

  for(int i=0;i<sc->dimension();i++) {
    sort(c.begin(), c.end(), compare_var_min(i));
    double v = c[0]->min(i);
    for(unsigned j=0;j<c.size();j++) {
      while(j<c.size() && c[j]->min(i)==v) j++;
      if (j<c.size()) {
        splits.push_back(new Split(i, c[j]->min(i)));
        v = c[j]->min(i);
//cout << i << "\t" << c[j]->min(i) << endl;
      }
    }
  }

  cout << "Splits: " << splits.size();
  cout << " --------------- Classes: " << classes.size() << endl;
/*  for(unsigned i=0;i<sc->cubos.size();i++) 
    cout << classes[c[i]->tag] << "|" << c[i]->tag  <<  endl;*/
 
  root = new Node(0, sc->dimension(), sc->size());

  BuildTree(root, sc, splits);

}
//---------------------------------------------------------------------------
void SpaceCubesTree::BuildTree(Node *n, SpaceCubes *sc, vector<Split*> splits)
{
  vector<Split*>::iterator isplit = GetBestSplit(sc, splits);
  if (isplit!=splits.end()) {
    SpaceCubes *cR = sc->split((*isplit)->iSplit, (*isplit)->fSplit);
    n->fSplit = (*isplit)->fSplit;
    n->att = (*isplit)->iSplit;
    n->coeff[n->att] = 1.0;
    n->SplitType = ORD;
    splits.erase(isplit);
    Node *cur = n;
    while (cur) {
      cur->T++;
      cur = cur->parent;
    }

    //Hijo 1
    n->child = new Node(n, n->sizecoeff, sc->size());
    if (!IsPure(sc)) BuildTree(n->child, sc, splits);
    else n->child->iClass = classes[sc->cubos[0]->tag];

    //Hijo 2
    n->child->sib = new Node(n, n->sizecoeff, cR->size());
    if (!IsPure(cR)) BuildTree(n->child->sib, cR, splits);
    else n->child->sib->iClass = classes[cR->cubos[0]->tag];
  }
}
//---------------------------------------------------------------------------
double FindImpurity(vector<double> d)
{
  double aux = 0.0;
  double TotPop = 0.0;

//
  // Gini criterion multiplied by node membership

  for(unsigned j=0; j<d.size() ; j++)  {
     aux -=  d[j]*d[j];
     TotPop += d[j];
  }

  return (TotPop + aux/TotPop);// = n*i(t) = n*(1 - sum_j(p(j|t)^2)
}
vector<Split*>::iterator SpaceCubesTree::GetBestSplit(SpaceCubes *sc, 
                                                       vector<Split*> & splits)
{
  unsigned n_classes = classes.size();
  int d = sc->dimension();
  vector<Split*>::iterator is, ibest;
  vector<HyperCube*>::iterator ihc;
  double imp_min = numeric_limits<double>::infinity();

  for(int i=0;i<d;i++) {
    sort(sc->cubos.begin(), sc->cubos.end(), compare_var_min(i));

    for(is=splits.begin();is!=splits.end();is++) {
      if ((*is)->iSplit!=i) continue;
      vector<double> dL(n_classes, 0.0);
      vector<double> dR(n_classes, 0.0);
      bool algoR = false;
      bool algoL = false;
      for(ihc=sc->cubos.begin();ihc!=sc->cubos.end();ihc++) {
        int iii = classes[(*ihc)->tag];
        if ((*ihc)->min(i)>=(*is)->fSplit) {
          dR[iii] += 1.0;
          algoR = true;
        }
        else if ((*ihc)->max(i)<=(*is)->fSplit) {
          dL[iii] += 1.0;
          algoL = true;
        }
        else {
          dR[iii] += 0.5;
          dL[iii] += 0.5;
          algoR = true;
          algoL = true;
        }
      }
      if (algoR && algoL) {
        double imp = FindImpurity(dL) + FindImpurity(dR);
        if (imp<imp_min) {
          ibest = is;
          imp_min = imp;
          (*ibest)->impR = FindImpurity(dR);//Por no declararme una var mas
          (*ibest)->impL = imp_min - (*ibest)->impR;
        }
      }
    }
  }

  return ibest;
}
//---------------------------------------------------------------------------
bool SpaceCubesTree::IsPure(SpaceCubes *sc)
{
  for(unsigned i=1;i<sc->cubos.size();i++) {
    if (sc->cubos[i-1]->tag!=sc->cubos[i]->tag) return false;
  }

  return true;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
string printall(vector<HyperCube*> lst)
{
  std::ostringstream puf;

  for(unsigned i=0;i<lst.size();i++) {
    puf << lst[i]->Print() << endl;
  }

  return puf.str();
}
//---------------------------------------------------------------------------
void banana_split(Node *banana, HyperCube *choco, vector<HyperCube*> &late,
                                                                        int K)
{
//  static int i = 0;
//  char fn[10];

  HyperCube *fresa;

  if (!banana->IsLeaf(K)) {

    choco->tag = banana->iClass;
//    sprintf (fn, "recs%d", i++);
//    FILE *f=fopen(fn , "w");
//    fprintf(f, "%s", choco->Print().c_str());
//    fclose(f);

    fresa = choco->split(banana->att, banana->fSplit);
    banana_split(banana->child->sib, fresa, late, K);
    banana_split(banana->child, choco, late, K);
  }
  else {
    choco->tag = banana->iClass;
    late.push_back(choco);
  }
}
//---------------------------------------------------------------------------
vector<HyperCube*> TreeToHyperCubes(DecisionTree *tree, HyperCube *hc)
{
  return TreeToHyperCubes(tree->GetTree(), tree->GetK(), hc);
}
//---------------------------------------------------------------------------
vector<HyperCube*> TreeToHyperCubes(Tree *tree, int K, HyperCube *hc)
{
  vector<HyperCube*> lst;
  HyperCube *root = hc ? hc : new HyperCube(tree->GetRoot()->sizecoeff);

  banana_split(tree->GetRoot(), root, lst, K);
  
  return lst;
}
//---------------------------------------------------------------------------
int uniq(Data *d, int first, int last)
{

  d->SortOn(0, first, last);

  for(int i=1;i<d->GetNumVar();i++) {
    for(int j=first+1;j<=last;j++) {
      int prev = j-1;
      while(j<=last && d->GetValueVar(j-1, i-1)==d->GetValueVar(j, i-1)) j++;
      d->SortOn(i, prev, j-1);
    }
  }

  int ndistintos = 1;
  d->SetValueVar(first, d->GetNumVar(), ndistintos);
  for(int j=first;j<last;j++) {
    bool distintos = false;
    for(int i=0;i<d->GetNumVar();i++) {
      if (d->GetValueVar(j, i)!=d->GetValueVar(j+1, i)) {
        distintos = true;
        ndistintos++;
        break;
      }
    }
    double v = ndistintos + (distintos ? 0.0 : last+1);
    d->SetValueVar(j+1, d->GetNumVar(), v);
  }

  d->SortOn(d->GetNumVar(), first, last);

  return ndistintos;
}

