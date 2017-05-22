#include <vector>
#include <map>
#include "tree20.h"
#include "GPTree.h"

#ifndef _SPACE_CUBES_H
#define _SPACE_CUBES_H

class Node;
class Data;

//---------------------------------------------------------------------------
class ClassTag
{
  int imax;
  std::vector<double> dis;

  public:
    ClassTag() {imax=-1;}
    ClassTag(int iclass, double w = 1.0) {
      imax = iclass;
      if (iclass<0) return; 
      dis.resize(iclass+1, 0.0);
      dis[iclass] += w;
    }

  public:
    int size() const {return dis.size();}

    double operator () (int i) const {return dis[i];}

    double & operator [] (int i) 
    {
      if (i>=size()) dis.resize(i+1, 0.0);
      return dis[i];
    }

    ClassTag& operator = (const ClassTag & other)
    {
      imax = other.imax;
  dis.clear();
  for(unsigned i=0;i<other.dis.size();i++) {
    dis.push_back(other.dis[i]);
  }
//      dis = other.dis;
      return *this;
    }

    bool operator == (const ClassTag & other) const {return dis==other.dis;}
    bool operator != (const ClassTag & other) const {return !((*this)==other);}
    bool operator <(ClassTag const& other) const 
    {
      if (this->size()!=other.size())
        return this->size() < other.size();

      for (int i=0;i<this->size();i++)
        if ((*this)(i)!=other(i)) 
          return (*this)(i) > other(i);

      return false;
    }
    bool operator >(ClassTag const& other) const 
    {
      return (*this)!=other && !((*this)<other);
    }

    void operator += (const ClassTag & other)
    {
      if (other.dis.size()>dis.size())
        dis.resize(other.dis.size(), 0.0);

      double max = -1; //dis[imax] + other.dis[imax];
      for(unsigned i=0;i<other.dis.size();i++) {
        dis[i] += other.dis[i];
        if (dis[i]>max) {
          imax = i;
          max = dis[i];
        }
      }
    }
    operator int() {return imax;}

    friend ClassTag operator + (const ClassTag & one, const ClassTag & another)
    {
      ClassTag res;
      res += one;
      res += another;
      return res;
    }

};
std::ostream& operator << (std::ostream &out, ClassTag tag);
//---------------------------------------------------------------------------
class HyperCube
{
  public:
    ClassTag tag;

  protected:
    std::vector<double> mins;
    std::vector<double> maxs;
    double volume;

    std::vector<std::vector<HyperCube*> > intouch;

  public: //Constructores
    HyperCube(std::vector<double> mins, std::vector<double> maxs, 
                                                     ClassTag tag=ClassTag());
    HyperCube(int d, int tag=-1);

  public: //funciones
    HyperCube* split(int d, double val); //Division
    bool Merge(HyperCube *hc); //Union (si podemos)
    HyperCube *intersec(HyperCube *hc); //Interseccion
    std::string Print();

  public:
    int dimension() {return mins.size();}
    double min(int i) {return mins[i];}
    double max(int i) {return maxs[i];}
    double arista(int i) {return maxs[i]-mins[i];}
    double centre(int i) {return mins[i]+(maxs[i]-mins[i])/2.0;}
    double vol();

  friend class SpaceCubes;
};

//---------------------------------------------------------------------------
class SpaceCubes
{
  protected:
    std::vector<HyperCube*> cubos;
    std::vector<double> mins;
    std::vector<double> maxs;

  protected:
    std::vector<HyperCube*>::iterator Merge(
                                        std::vector<HyperCube*>::iterator ini, 
                                        std::vector<HyperCube*>::iterator fin);

  public:
    SpaceCubes(int d);
    SpaceCubes(DecisionTree *tree, std::vector<double> mins, 
                                                     std::vector<double> maxs);
    SpaceCubes(Tree *tree, std::vector<double> mins, std::vector<double> maxs);

  public:
    SpaceCubes operator &(SpaceCubes other);
    SpaceCubes* split(int d, double val); //Division
    int Merge();
    std::string Print();

  public:
    int dimension() {return cubos[0]->dimension();}
    int size() {return cubos.size();}
    HyperCube* GetHyperCube(int i) {return cubos[i];}
    Data* ToData();
    Tree* ToTree();

  friend class SpaceCubesTree;
};
//---------------------------------------------------------------------------
class Split
{
  public:
    int iSplit;
    double fSplit;
    double impL;
    double impR;

  public:
    Split(int i, double f) {iSplit=i; fSplit=f;}
};
//---------------------------------------------------------------------------
class SpaceCubesTree : public Tree
{
  protected:
    std::map<ClassTag, int> classes;

  public:
    SpaceCubesTree(SpaceCubes *sc);

  public:
    void BuildTree(Node *n, SpaceCubes *sc, std::vector<Split*> splits);
    std::vector<Split*>::iterator GetBestSplit(SpaceCubes *sc, 
                                                 std::vector<Split*> & splits);
    bool IsPure(SpaceCubes *sc);
};
//---------------------------------------------------------------------------
std::string printall(std::vector<HyperCube*> lst);
void banana_split(Node *banana, HyperCube *choco, 
                                        std::vector<HyperCube*> &late, int K);
std::vector<HyperCube*> TreeToHyperCubes(DecisionTree *tree, HyperCube *hc=0);
std::vector<HyperCube*> TreeToHyperCubes(Tree *tree, int K, HyperCube *hc=0);
int uniq(Data *d, int first, int last);

#endif

