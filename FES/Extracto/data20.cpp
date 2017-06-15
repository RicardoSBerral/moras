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
#include<algorithm>
#include<map>
#include<numeric>

using namespace std;

//#define ALLINONE

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/*
NTotalOrdfuzz = Nordfuzz+1
The indexing scheme:
The ordinal independent varibles are 0->(Nord-1) in instances
The fuzzy independent variables are Nord->(Nordfuzz-1) in instances
The nominal independent variables are Nordfuzz->(Nordfuzz+Nnom-1) in instances
The dependent variable  is (Nvar-1) in instances
Nvar = Nordfuzz+Nnom+1

The independent variables are 0-(Nordfuzz) in St
  with one of these (Nordfuzz) the constant.
The dep       "         is (NtotalOrdfuzz) in St

Index[b][i] = index in instances of var i in St if Nordfuzz>i>=0
Index[b][i] < 0 if i refers to constant
Index[b][i] = Nordfuzz refers to dependent variable
*/

// The following macro gives index in instances
//    of var j on branch b
//#define VAR(b,j) ((Index[b][j]<0 ? Nordfuzz : Index[b][j]))
// The following macro gives the actual value of
// the variable, in the n-th case, of var j on
// branch b.
//#define VAL(n,b,j) (Index[b][j]<0 ? 1.0 : m[n][(Index[b][j] == Nordfuzz) ? Nvar-1 : Index[b][j]])

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

int MultilabelInstance::num_multilabels = -1;
int MultilabelInstance::label_index     = -1;

int MultiobjectiveInstance::num_multiobjectives = -1;
int MultiobjectiveInstance::objective_index     = -1;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Instance *Instance::NewInstanceByName(string name)
{
  if (name.compare("BasicInstance")==0) {
    return new BasicInstance(0);
  }
  else if (name.compare("UInt8Instance")==0) {
    return new UInt8Instance(0);
  }
  else if (name.compare("FloatInstance")==0) {
    return new FloatInstance(0);
  }
  else if (name.compare("Int32Instance")==0) {
    return new Int32Instance(0);
  }
  else if (name.compare("ShortInstance")==0) {
    return new ShortInstance(0);
  }
  else if (name.compare("ComboInstance")==0) {
    return new ComboInstance(0);
  }
  else if (name.compare("MultilabelInstance")==0) {
    return new MultilabelInstance(0);
  }
  else if (name.compare("MultiobjectiveInstance")==0) {
    return new MultiobjectiveInstance(0);
  }

  return 0;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
BasicInstance::BasicInstance(int nvars) : Instance (nvars)
{
  this->nvars = nvars < 0 ? 0 : nvars;
  data = new double[this->nvars + 3];
  for ( int i = 0; i < this->nvars + 3; i++) {
    data[i] = 0.0;
  }
}
BasicInstance::BasicInstance(double *vals, int nvars) : Instance (nvars)
{
  this->nvars = nvars;
  data = new double[nvars+3];
  for ( int i = 0; i < nvars+3; i++) {
    data[i] = vals[i];
  }
}
BasicInstance::~BasicInstance()
{
  delete []data;
}
inline double BasicInstance::operator[](int index)
{
/*  if (index==nvars) {
    cout << "WARNING!!! accesing class through [] operator" << endl << endl;
  }
  else if (index > nvars) {
    cout << "WARNING!!! accesing WV through [] operator" << endl << endl;
  }*/

  return data[index];
}
inline void BasicInstance::operator()(int index, double valor)
{
/*  if (index==nvars) {
    cout << "WARNING!!! setting class through () operator" << endl << endl;
  }
  else if (index > nvars) {
    cout << "WARNING!!! setting WV through () operator" << endl << endl;
  }*/
  data[index]=valor;
}
BasicInstance::operator std::vector<double>()
{
  vector<double> ins;
  for (int i = 0; i < nvars; i++) {
    ins.push_back(data[i]);
  }
  ins.push_back(Class+0.5);
  ins.push_back(data[nvars+1]);   //MultipleSplitWs1);
  ins.push_back(data[nvars+2]); //MultipleSplitWs2);
  ins.push_back(CvTreeLabel);
  ins.push_back(Membership);
  ins.push_back(Weight);
  ins.push_back(IniPos);
  ins.push_back(Group);

  return ins;
}
void BasicInstance::DeleteColumn(int icol)
{
  nvars--;
  memmove(data+icol, data+(icol+1), sizeof(double) * (nvars+3 - icol));
}
inline Instance& BasicInstance::operator = (const Instance & other)
{
  const BasicInstance *pother = dynamic_cast<const BasicInstance*>(&other);
  if (!pother || pother->nvars != nvars) cout << "Error assigning instance" << endl;

  for (int i = 0; i < nvars+3; i++) data[i] = pother->data[i];

  Class             = pother->Class;
  IniPos            = pother->IniPos;
  Group             = pother->Group;
//  MultipleSplitWs1  = pother->MultipleSplitWs1;
//  MultipleSplitWs2  = pother->MultipleSplitWs2;
  CvTreeLabel       = pother->CvTreeLabel;
  Membership        = pother->Membership;
  Weight            = pother->Weight;

  return *this;
}
void BasicInstance::ChangeDependentColumn(int newDependentColumnIndex) {
  double hold = (*this)[newDependentColumnIndex];
  (*this)(newDependentColumnIndex, (*this)[nvars]);
  (*this)(nvars, hold);
}
void BasicInstance::CopyDependentColumn() {
  // We construct the new data, same as the old one, except for the duplicated column
  double* newData = new double[nvars + 4];
  std::copy(data, data + nvars + 1, newData);
  newData[nvars + 1] = newData[nvars];
  std::copy(data + nvars + 1, data + nvars + 3, newData + nvars + 2);
  
  // Substitution of the relevant fields
  delete [] data;
  data = newData;
  nvars += 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
UInt8Instance::UInt8Instance(int nvars) : Instance (nvars)
{
  this->nvars = nvars;
  data = new unsigned char[nvars];
  for ( int i = 0; i < nvars; i++) {
    data[i] = 0;
  }
}
UInt8Instance::UInt8Instance(double *vals, int nvars) : Instance (nvars)
{
  this->nvars = nvars;
  data = new unsigned char[nvars];
  for ( int i = 0; i < nvars; i++) {
    data[i] = (unsigned char)vals[i];
  }
}
UInt8Instance::~UInt8Instance()
{
  delete []data;
}
inline double UInt8Instance::operator[](int index)
{
  return index < nvars ? data[index] : wv[index-nvars-1];
}
inline void UInt8Instance::operator()(int index, double valor)
{
  if (index < nvars) data[index] = (unsigned char) (valor + 0.5);
  else wv[index-nvars-1] = valor;
}
UInt8Instance::operator std::vector<double>()
{
  vector<double> ins;
  for (int i = 0; i < nvars; i++) {
    ins.push_back(data[i]);
  }
  ins.push_back(Class+0.5);
  ins.push_back(wv[0]);   //MultipleSplitWs1);
  ins.push_back(wv[1]); //MultipleSplitWs2);
  ins.push_back(CvTreeLabel);
  ins.push_back(Membership);
  ins.push_back(Weight);
  ins.push_back(IniPos);
  ins.push_back(Group);

  return ins;
}
void UInt8Instance::DeleteColumn(int icol)
{
  nvars--;
  memmove(data+icol, data+(icol+1), sizeof(unsigned char) * (nvars - icol));
}
inline Instance& UInt8Instance::operator = (const Instance & other)
{
  const UInt8Instance *pother = dynamic_cast<const UInt8Instance*>(&other);
  if (!pother || pother->nvars != nvars) cout << "Error assigning instance" << endl;

  for (int i = 0; i < nvars; i++) data[i] = pother->data[i];

  Class             = pother->Class;
  IniPos            = pother->IniPos;
  Group             = pother->Group;
  wv[0]             = pother->wv[0];
  wv[1]             = pother->wv[1];
  CvTreeLabel       = pother->CvTreeLabel;
  Membership        = pother->Membership;
  Weight            = pother->Weight;

  return *this;
}
void UInt8Instance::ChangeDependentColumn(int newDependentColumnIndex) {
  throw std::logic_error("UInt8Instance::ChangeDependentColumn is not implemented yet.");
}
void UInt8Instance::CopyDependentColumn() {
  throw std::logic_error("UInt8Instance::CopyDependentColumn is not implemented yet.");
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
ShortInstance::ShortInstance(int nvars) : Instance (nvars)
{
  this->nvars = nvars;
  data = new short[nvars];
  for ( int i = 0; i < nvars; i++) {
    data[i] = 0;
  }
}
ShortInstance::ShortInstance(double *vals, int nvars) : Instance (nvars)
{
  this->nvars = nvars;
  data = new short[nvars];
  for ( int i = 0; i < nvars; i++) {
    data[i] = (short)vals[i];
  }
}
ShortInstance::~ShortInstance()
{
  delete []data;
}
inline double ShortInstance::operator[](int index)
{
  return index < nvars ? data[index] : wv[index-nvars-1];
}
inline void ShortInstance::operator()(int index, double valor)
{
  if (index < nvars) data[index] = (short)(valor + 0.5);
  else wv[index-nvars-1] = valor;
}
ShortInstance::operator std::vector<double>()
{
  vector<double> ins;
  for (int i = 0; i < nvars; i++) {
    ins.push_back(data[i]);
  }
  ins.push_back(Class+0.5);
  ins.push_back(wv[0]);   //MultipleSplitWs1);
  ins.push_back(wv[1]); //MultipleSplitWs2);
  ins.push_back(CvTreeLabel);
  ins.push_back(Membership);
  ins.push_back(Weight);
  ins.push_back(IniPos);
  ins.push_back(Group);

  return ins;
}
void ShortInstance::DeleteColumn(int icol)
{
  nvars--;
  memmove(data+icol, data+(icol+1), sizeof(short) * (nvars - icol));
}
inline Instance& ShortInstance::operator = (const Instance & other)
{
  const ShortInstance *pother = dynamic_cast<const ShortInstance*>(&other);
  if (!pother || pother->nvars != nvars) cout << "Error assigning instance" << endl;

  for (int i = 0; i < nvars; i++) data[i] = pother->data[i];

  Class             = pother->Class;
  IniPos            = pother->IniPos;
  Group             = pother->Group;
  wv[0]             = pother->wv[0];
  wv[1]             = pother->wv[1];
  CvTreeLabel       = pother->CvTreeLabel;
  Membership        = pother->Membership;
  Weight            = pother->Weight;

  return *this;
}
void ShortInstance::ChangeDependentColumn(int newDependentColumnIndex) {
  throw std::logic_error("UInt8Instance::ShortInstance is not implemented yet.");
}
void ShortInstance::CopyDependentColumn() {
  throw std::logic_error("ShortInstance::CopyDependentColumn is not implemented yet.");
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Int32Instance::Int32Instance(int nvars) : Instance (nvars)
{
  this->nvars = nvars;
  data = new int[nvars];
  for ( int i = 0; i < nvars; i++) {
    data[i] = 0;
  }
}
Int32Instance::Int32Instance(double *vals, int nvars) : Instance (nvars)
{
  this->nvars = nvars;
  data = new int[nvars];
  for ( int i = 0; i < nvars; i++) {
    data[i] = (int)vals[i];
  }
}
Int32Instance::~Int32Instance()
{
  delete []data;
}
inline double Int32Instance::operator[](int index)
{
  return index < nvars ? data[index] : wv[index-nvars-1];
}
inline void Int32Instance::operator()(int index, double valor)
{
  if (index < nvars) data[index] = (int)(valor + 0.5);
  else wv[index-nvars-1] = valor;
}
Int32Instance::operator std::vector<double>()
{
  vector<double> ins;
  for (int i = 0; i < nvars; i++) {
    ins.push_back(data[i]);
  }
  ins.push_back(Class+0.5);
  ins.push_back(wv[0]);   //MultipleSplitWs1);
  ins.push_back(wv[1]); //MultipleSplitWs2);
  ins.push_back(CvTreeLabel);
  ins.push_back(Membership);
  ins.push_back(Weight);
  ins.push_back(IniPos);
  ins.push_back(Group);

  return ins;
}
void Int32Instance::DeleteColumn(int icol)
{
  nvars--;
  memmove(data+icol, data+(icol+1), sizeof(int) * (nvars - icol));
}
inline Instance& Int32Instance::operator = (const Instance & other)
{
  const Int32Instance *pother = dynamic_cast<const Int32Instance*>(&other);
  if (!pother || pother->nvars != nvars) cout << "Error assigning instance" << endl;

  for (int i = 0; i < nvars; i++) data[i] = pother->data[i];

  Class             = pother->Class;
  IniPos            = pother->IniPos;
  Group             = pother->Group;
  wv[0]             = pother->wv[0];
  wv[1]             = pother->wv[1];
  CvTreeLabel       = pother->CvTreeLabel;
  Membership        = pother->Membership;
  Weight            = pother->Weight;

  return *this;
}
void Int32Instance::ChangeDependentColumn(int newDependentColumnIndex) {
  throw std::logic_error("Int32Instance::ShortInstance is not implemented yet.");
}
void Int32Instance::CopyDependentColumn() {
  throw std::logic_error("Int32Instance::CopyDependentColumn is not implemented yet.");
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
FloatInstance::FloatInstance(int nvars) : Instance (nvars)
{
  this->nvars = nvars;
  data = new float[nvars];
  for ( int i = 0; i < nvars; i++) {
    data[i] = 0.0;
  }
}
FloatInstance::FloatInstance(double *vals, int nvars) : Instance (nvars)
{
  this->nvars = nvars;
  data = new float[nvars];
  for ( int i = 0; i < nvars; i++) {
    data[i] = vals[i];
  }
}
FloatInstance::~FloatInstance()
{
  delete []data;
}
inline double FloatInstance::operator[](int index)
{
  return index < nvars ? data[index] : wv[index-nvars-1];
}
inline void FloatInstance::operator()(int index, double valor)
{
  if (index < nvars) data[index] = (float)valor;
  else wv[index-nvars-1] = valor;
}
FloatInstance::operator std::vector<double>()
{
  vector<double> ins;
  for (int i = 0; i < nvars; i++) {
    ins.push_back(data[i]);
  }
  ins.push_back(Class+0.5);
  ins.push_back(wv[0]);   //MultipleSplitWs1);
  ins.push_back(wv[1]); //MultipleSplitWs2);
  ins.push_back(CvTreeLabel);
  ins.push_back(Membership);
  ins.push_back(Weight);
  ins.push_back(IniPos);
  ins.push_back(Group);

  return ins;
}
void FloatInstance::DeleteColumn(int icol)
{
  nvars--;
  memmove(data+icol, data+(icol+1), sizeof(float) * (nvars - icol));
}
inline Instance& FloatInstance::operator = (const Instance & other)
{
  const FloatInstance *pother = dynamic_cast<const FloatInstance*>(&other);
  if (!pother || pother->nvars != nvars) cout << "Error assigning instance" << endl;

  for (int i = 0; i < nvars; i++) data[i] = pother->data[i];

  Class             = pother->Class;
  IniPos            = pother->IniPos;
  Group             = pother->Group;
  wv[0]             = pother->wv[0];
  wv[1]             = pother->wv[1];
  CvTreeLabel       = pother->CvTreeLabel;
  Membership        = pother->Membership;
  Weight            = pother->Weight;

  return *this;
}
void FloatInstance::ChangeDependentColumn(int newDependentColumnIndex) {
  throw std::logic_error("FloatInstance::ShortInstance is not implemented yet.");
}
void FloatInstance::CopyDependentColumn() {
  throw std::logic_error("FloatInstance::CopyDependentColumn is not implemented yet.");
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// The 11 values stored in ComboIntance are:
//    1-2: Location of the center of the detection
//    3-4: Scale or size of detection
//    5-6: Delta x & y to the center of the object from 1-2
//    7-8: Scale x & y of the object
//    9:   Binning of the deltas
//    10:  Binning of the scales
//    11:  Class of the object
ComboInstance::ComboInstance(int nvars) : Instance (nvars)
{
  this->nvars = nvars;
  data_location_bb = new float[11];
  data_description = nvars > 11 ? new unsigned char[nvars-11] : 0;
  for ( int i = 0; i < nvars-11; i++) {
    data_description[i] = 0;
  }
  for ( int i = 0; i < 11; i++) {
    data_location_bb[i] = 0;
  }
}
ComboInstance::ComboInstance(double *vals, int nvars) : Instance (nvars)
{
  this->nvars = nvars;
  data_location_bb = new float[11];
  data_description = nvars > 11 ? new unsigned char[nvars-11] : 0;
  for ( int i = 0; i < nvars; i++) {
    if (i < nvars-11) {
      data_description[i] = vals[i];
    }
    else {
      data_location_bb[i-nvars+11] = vals[i];
    }
  }
}
ComboInstance::~ComboInstance()
{
  delete []data_description;
  delete []data_location_bb;
}
inline double ComboInstance::operator[](int index)
{
  return index < nvars-11 ? data_description[index] : 
         (index < nvars ? data_location_bb[index-nvars+11] : wv[index-nvars-1]);
}
inline void ComboInstance::operator()(int index, double valor)
{
  if (index < nvars-11) data_description[index] = (unsigned char)(valor + 0.5);
  else if (index < nvars) data_location_bb[index-nvars+11] = valor;
  else wv[index-nvars-1] = valor;
}
ComboInstance::operator std::vector<double>()
{
  vector<double> ins;
  for (int i = 0; i < nvars; i++) {
    if (i < nvars-11) {
      ins.push_back(data_description[i]);
    }
    else {
      ins.push_back(data_location_bb[i-nvars+11]);
    }
  }
  ins.push_back(Class+0.5);
  ins.push_back(wv[0]);   //MultipleSplitWs1);
  ins.push_back(wv[1]); //MultipleSplitWs2);
  ins.push_back(CvTreeLabel);
  ins.push_back(Membership);
  ins.push_back(Weight);
  ins.push_back(IniPos);
  ins.push_back(Group);

  return ins;
}
void ComboInstance::DeleteColumn(int icol)
{
  nvars--;
  memmove(data_description+icol, data_description+(icol+1), 
                                  sizeof(unsigned char) * (nvars - icol));
}
inline Instance& ComboInstance::operator = (const Instance & other)
{
  const ComboInstance *pother = dynamic_cast<const ComboInstance*>(&other);
  if (!pother || pother->nvars != nvars) cout << "Error assigning instance" << endl;

  for (int i = 0; i < nvars-11; i++) data_description[i] = pother->data_description[i];
  for (int i = 0; i < 11; i++) data_location_bb[i] = pother->data_location_bb[i];

  Class             = pother->Class;
  IniPos            = pother->IniPos;
  Group             = pother->Group;
  wv[0]             = pother->wv[0];
  wv[1]             = pother->wv[1];
  CvTreeLabel       = pother->CvTreeLabel;
  Membership        = pother->Membership;
  Weight            = pother->Weight;

  return *this;
}
void ComboInstance::ChangeDependentColumn(int newDependentColumnIndex) {
  throw std::logic_error("ComboInstance::ShortInstance is not implemented yet.");
}
void ComboInstance::CopyDependentColumn() {
  throw std::logic_error("ComboInstance::CopyDependentColumn is not implemented yet.");
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
MultilabelInstance::MultilabelInstance(int nvars) : BasicInstance (nvars - num_multilabels)
{
}
MultilabelInstance::MultilabelInstance(double *vals, int nvars) : BasicInstance (vals, nvars - num_multilabels)
{
  for(int i = 0; i < num_multilabels; i++) {
    labels[i] = (bool)vals[nvars - num_multilabels + i];
  }
}
inline double MultilabelInstance::operator[](int index)
{
  return index < nvars ? BasicInstance::operator[](index) : 
         GetMultilabel(index - nvars);
}
inline void MultilabelInstance::operator()(int index, double valor)
{
  if (index < nvars) BasicInstance::operator()(index, valor);
  else SetMultilabel(index - nvars, (int)valor);
}
MultilabelInstance::operator std::vector<double>()
{
  vector<double> ins = BasicInstance::operator std::vector<double>();
  for (int i = 0; i < num_multilabels; i++) {
    ins.push_back(GetMultilabel(i));
  }

  return ins;
}
inline Instance& MultilabelInstance::operator = (const Instance & other)
{
  BasicInstance::operator=(other); 

  const MultilabelInstance *pother = dynamic_cast<const MultilabelInstance*>(&other);
  this->labels = pother->labels;

  return *this;
}
void MultilabelInstance::ChangeDependentColumn(int newDependentColumnIndex) {
  throw std::logic_error("MultilabelInstance::ShortInstance is not implemented yet.");
}
void MultilabelInstance::CopyDependentColumn() {
  throw std::logic_error("MultilabelInstance::CopyDependentColumn is not implemented yet.");
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
MultiobjectiveInstance::MultiobjectiveInstance(int nvars) : BasicInstance (nvars)
{
  objectives.reserve(num_multiobjectives > 1 ? num_multiobjectives : 1);
}
MultiobjectiveInstance::MultiobjectiveInstance(double *vals, int nvars) : BasicInstance (vals, nvars)
{
  objectives.reserve(num_multiobjectives);
  for(int i = 0; i < num_multiobjectives; i++) {
    objectives[i] = vals[i + GetNumIndependentVars()];
  }
}
inline double MultiobjectiveInstance::operator[](int index)
{
  return BasicInstance::operator[](index);
}
inline void MultiobjectiveInstance::operator()(int index, double valor)
{
  BasicInstance::operator()(index, valor);
  if (index >= GetNumIndependentVars()) {
    SetMultiobjective(index - GetNumIndependentVars(), valor);
  }
}
MultiobjectiveInstance::operator std::vector<double>()
{
  vector<double> ins = BasicInstance::operator std::vector<double>();
  for (int i = 0; i < num_multiobjectives; i++) {
    ins.push_back(GetMultiobjective(i));
  }

  return ins;
}
inline Instance& MultiobjectiveInstance::operator = (const Instance & other)
{
  BasicInstance::operator=(other); 

  const MultiobjectiveInstance *pother = dynamic_cast<const MultiobjectiveInstance*>(&other);
  this->objectives = pother->objectives;

  return *this;
}
void MultiobjectiveInstance::ChangeDependentColumn(int newDependentColumnIndex) {

  // First, we call the super method
  BasicInstance::ChangeDependentColumn(newDependentColumnIndex);

  // If the new dependent column is an objective...
  if (newDependentColumnIndex >= GetNumIndependentVars()) {
    objectives[newDependentColumnIndex-GetNumIndependentVars()] = (*this)[newDependentColumnIndex];
  }

  // Either dependent or independent...
  objectives[num_multiobjectives-1] = (*this)[nvars-1];
}
void MultiobjectiveInstance::CopyDependentColumn() {

  // First, we call the super method
  BasicInstance::CopyDependentColumn();
  
  // Insert the new objective
  // We assume the number of objectives has already been updated
  objectives.resize(num_multiobjectives, objectives[num_multiobjectives-2]);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
double GiniCriterium::Impurity(double *Pop, int NumClass)
{

  double aux = 0.0;
  double TotPop = 0.0;

  // Gini criterion
  for(int j = 0; j < NumClass; j++) {
    TotPop += Pop[j];
    aux -=  Pop[j]*Pop[j];
  }
  return (TotPop + aux/TotPop);// = n*i(t) = n*(1 - sum_j(p(j|t)^2)

}
double GiniCriterium::Initialise(int first, int last, int curSplit, double gamma, int iX, bool FuzzFlag)
{

  bool flagleft = true;
  Instances& instances = data->GetM();

  if (PopL) delete []PopL;
  if (PopR) delete []PopR;

  PopL = new double[ ((NomData*)data)->NumClass ];
  PopR = new double[ ((NomData*)data)->NumClass ];
  
  for (int j = 0; j < ((NomData*)data)->NumClass; j++) {
    PopL[j] = PopR[j] = 0.0;
  }

  if (FuzzFlag) {
    if (curSplit == last) {
      nL = R[0] = 0.0;
    }
    else {  
      flagleft = false; 
      nR = R[1] = 0.0;
    } 
  }
  else {
    R[0] = R[1] = nL = nR = 0.0;       // Initialise
  }

  for (int i = first; i <= last; i++) {
    if (instances[i].GetMembership() > 0.0) {
      if ( (gamma > 1 && i <= curSplit) // normal
                               || (gamma < 1 && instances[i][iX]+gamma < 0)) {
        PopL[instances[i].GetClass()] += instances[i].GetMembership();
        nL += instances[i].GetMembership();
      }
      else {
        PopR[instances[i].GetClass()] += instances[i].GetMembership();
        nR += instances[i].GetMembership();
      }
    }
  }
  

  if (nL > 0.0) {
    if(!FuzzFlag   || (FuzzFlag && flagleft) ) {
      R[0] = Impurity(PopL, ((NomData*)data)->NumClass);
    }
  }
  if (nR > 0.0) {
    if(!FuzzFlag  || (FuzzFlag && !flagleft) ) { 
      R[1] = Impurity(PopR, ((NomData*)data)->NumClass);
    }
  }

  if ( !FuzzFlag )
    return (nR < nL ? nR : nL);
  
  if (flagleft) return nL;
  
  return nR; 
}
double GiniCriterium::MovePoint(int curSplit, double gamma, int iX,int first,int last)
{
  Instances& instances = data->GetM();

  if (instances[curSplit].GetMembership() > 0.0) {
    if (gamma > 1 || (gamma<1 && instances[curSplit][iX]+gamma <0)) {
      PopL[instances[curSplit].GetClass()] -= instances[curSplit].GetMembership();
      PopR[instances[curSplit].GetClass()] += instances[curSplit].GetMembership();

      nL -= instances[curSplit].GetMembership();
      nR += instances[curSplit].GetMembership();
    } 
    else {
      PopL[instances[curSplit].GetClass()] += instances[curSplit].GetMembership();
      PopR[instances[curSplit].GetClass()] -= instances[curSplit].GetMembership();  

      nL += instances[curSplit].GetMembership();
      nR -= instances[curSplit].GetMembership();
    } 
    
    R[0] = R[1] = 0.0;
  
    if(nL > 0.0) R[0] = Impurity(PopL, ((NomData*)data)->NumClass);
    if(nR > 0.0) R[1] = Impurity(PopR, ((NomData*)data)->NumClass);
  }

  return (nR < nL ? nR : nL);
}
void GiniCriterium::FixLeaf(Node* n, bool add)
{
  Instances& instances = data->GetM();
  int NumClass = ((NomData*)data)->NumClass;
  double *NodePop = new double[NumClass];

  for(int i = 0; i < NumClass; i++ ) NodePop[i] = 0.0;
  ((NomData*)data)->FindPops(NodePop, n->first, n->last);

  if (add) {
    for(int i = 0; i < NumClass; i++ ) n->NodePop[i] += NodePop[i];
  }
  else {
    if (n->NodePop) delete []n->NodePop;
    n->NodePop = NodePop;
  }

  if (n->last - n->first + 1 == 0) {
;//    printf("Node w 0 evidence, class stays as %d\n", n->iClass);
  }
  else {
    n->iClass = ((NomData*)data)->WhichClass(n->NodePop);
  }

  if (((NomData*)data)->Getatdrawchosefathers() && n->parent && n->iClass != n->parent->iClass &&
      n->NodePop[n->iClass]>0 && n->NodePop[n->iClass] == n->NodePop[n->parent->iClass]) {
    n->iClass = n->parent->iClass;
  }

  if (n->d) {
    n->d->Rleaf = 0.0;
    for (int j = 0; j < NumClass; ++j) {
      if(j != n->iClass) {
        n->d->Rleaf += n->NodePop[j];
      }
    }
  }

  if ( data->GetNPrune() && n->d ) {
    n->d->RPrune = 0;
    n->d->firstTest = data->GetNumData();
    n->d->lastTest = data->GetNumData() + data->GetNPrune() - 1;
    data->FindMembershipTest(n);
    for(int i = n->d->firstTest; i <= n->d->lastTest; ++i) {
      if(instances[i].GetClass() != n->iClass)
        n->d->RPrune += instances[i].GetMembership();
    }
  }
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
double MultilabelGiniCriterium::Initialise(int first, int last, int curSplit, double gamma, int iX, bool FuzzFlag)
{
  //Select a random class label to split 


  int c = (int)((double)MultilabelInstance::GetNumMultilabels()  * TDebugRand::Rand() / (RAND_MAX + 1.0));
  MultilabelInstance::SetSplitlabelIndex(c);

//cout << "ml splitin on... " << iX << " with label number " << MultilabelInstance::GetSplitlabelIndex() << endl;

  return GiniCriterium::Initialise(first, last, curSplit, gamma, iX, FuzzFlag);
}
void MultilabelGiniCriterium::FixLeaf(Node* n, bool add)
{
  if (MultilabelInstance::GetSplitlabelIndex() < 0) {
    int c = (int)((double)MultilabelInstance::GetNumMultilabels()  * TDebugRand::Rand() / (RAND_MAX + 1.0));
    MultilabelInstance::SetSplitlabelIndex(c);
//cout << "FLml  with label number " << MultilabelInstance::GetSplitlabelIndex() << endl;
  }

  GiniCriterium::FixLeaf(n, add);

  int nmultilabels = MultilabelInstance::GetNumMultilabels();
  Matriz *clss = new Matriz(n->last - n->first + 1, nmultilabels);

  if (n->parent && n->parent->info) {
    delete (Matriz *)n->parent->info;
    n->parent->info = 0;
  }

  Instances& instances = data->GetM();
  for (int i = n->first; i <= n->last; i++) {
    for (int j = 0; j < nmultilabels; j++) {
      (*clss)[i - n->first][j] = (int)((MultilabelInstance&)instances[i]).GetMultilabel(j);
    }
  }
  n->info = clss;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void BoostedGiniCriterium::FixLeaf(Node* n, bool add)
{
  GiniCriterium::FixLeaf(n, add);

//  if (!n->parent) {root = n; return;}

  // Calculamos error
  double t = 0.0;
  for (int i = 0; i < ((NomData*)data)->NumClass; i++)
    t += n->NodePop[i];
  double leafError = 1.0 - n->NodePop[n->iClass] / t; 

//double leafError = data->Error2(root, -2, 0, data->GetNTrain()-1);

  // Si es mayor de 0.5 seguimos ejecucion normal
  //   (solo puede pasar si hay mas de 2 clases)
  if (leafError >= 0.5) return;

  // Actualizamos pesos
  Instances& instances = data->GetM();
  double KMal  = 1.0 / (2.0 * leafError);          //Mal clasificado
  double KBien = 1.0 / (2.0 * (1.0 - leafError));  //Bien clasificado
  t = 0.0;

  for(int i = n->first; i <= n->last; i++) {
  //for(int i = 0; i <= data->GetNTrain()-1; i++) {
    double newWeight =  instances[i].GetWeight();
    newWeight *= instances[i].GetClass() != n->iClass ? KMal : KBien;
    data->SetDatWeight(i, newWeight);
    t += newWeight;
  }
  for(int i = n->first; i <= n->last; i++) {
  //for(int i = 0; i <= data->GetNTrain()-1; i++) {
    double newWeight = instances[i].GetWeight();
    data->SetDatWeight(i, newWeight * (1.0 + n->last - n->first) / t);
    //data->SetDatWeight(i, newWeight * ((double)data->GetNTrain()) / t);
  }
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
double MSECriterium::Initialise(int first, int last, int curSplit, double gamma, int iX, bool FuzzFlag)
{

  bool flagleft = true;
  int Nvar = data->GetNumVar();
  Instances& instances = data->GetM();

  if (FuzzFlag) {
    if (curSplit == last) {
      nL = R[0] =  s1L = s2L  = 0.0;
    }
    else {  
      flagleft = false; 
      nR = R[1] = s1R = s2R = 0.0;
    } 
  }
  else {
    s1L = s2L = s1R = s2R = nL = nR = R[0] = R[1] = 0.0;
  }

  // Initialise
  for(int i = first; i <= last; i++) {
    if (instances[i].GetMembership() > 0.0) {

      if( (gamma > 1 && i <= curSplit) // normal
                                     || (gamma < 1 && instances[i][iX]+gamma < 0)) {
        s1L += instances[i].GetMembership()*instances[i][Nvar-1];
        s2L += instances[i].GetMembership()*instances[i][Nvar-1]*instances[i][Nvar-1];
        nL  += instances[i].GetMembership();
      } 
      else {
        s1R += instances[i].GetMembership()*instances[i][Nvar-1];
        s2R += instances[i].GetMembership()*instances[i][Nvar-1]*instances[i][Nvar-1];
        nR  += instances[i].GetMembership();
      }
    }
  }

  if (nL >0.0) {
    if (!FuzzFlag   || (FuzzFlag && flagleft)) {
      R[0] = s2L - s1L*s1L/nL;
    }
  }

  if (nR > 0.0) {
    if (!FuzzFlag  || (FuzzFlag && !flagleft)) {
      R[1] = s2R - s1R*s1R/nR;
    }
  }

  if (!FuzzFlag) return (nR < nL ? nR : nL);

  if (flagleft)  return nL;

  return nR; 
}

double MSECriterium::MovePoint(int curSplit, double gamma, int iX,int first,int last)
{
  int Nvar = data->GetNumVar();
  Instances& instances = data->GetM();

  if(instances[curSplit].GetMembership()>0.0) {
    if(gamma > 1) {
      s1L -= instances[curSplit].GetMembership() * instances[curSplit][Nvar-1]; 
      s2L -= instances[curSplit].GetMembership() * instances[curSplit][Nvar-1] * instances[curSplit][Nvar-1];
      s1R += instances[curSplit].GetMembership() * instances[curSplit][Nvar-1]; 
      s2R += instances[curSplit].GetMembership() * instances[curSplit][Nvar-1] * instances[curSplit][Nvar-1];
      nL  -= instances[curSplit].GetMembership(); 
      nR  += instances[curSplit].GetMembership();
    }
    else {
      if(instances[curSplit][iX]+gamma >=0) {
        s1R -= instances[curSplit].GetMembership() * instances[curSplit][Nvar-1];
        s2R -= instances[curSplit].GetMembership() * instances[curSplit][Nvar-1] * instances[curSplit][Nvar-1];
        nR  -= instances[curSplit].GetMembership();
        s1L += instances[curSplit].GetMembership() * instances[curSplit][Nvar-1];
        s2L += instances[curSplit].GetMembership() * instances[curSplit][Nvar-1] * instances[curSplit][Nvar-1];
        nL  += instances[curSplit].GetMembership();
      }
      else {
        s1L -= instances[curSplit].GetMembership() * instances[curSplit][Nvar-1];
        s2L -= instances[curSplit].GetMembership() * instances[curSplit][Nvar-1] * instances[curSplit][Nvar-1];
        nL  -= instances[curSplit].GetMembership();
        s1R += instances[curSplit].GetMembership() * instances[curSplit][Nvar-1];
        s2R += instances[curSplit].GetMembership() * instances[curSplit][Nvar-1] * instances[curSplit][Nvar-1];
        nR  += instances[curSplit].GetMembership();
      }
    }
  
    if (nL > 0.0) R[0] = s2L - s1L*s1L/nL;
    if (nR > 0.0) R[1] = s2R - s1R*s1R/nR;
  }

  return nR < nL ? nR : nL;
}

void MSECriterium::FixLeaf(Node* n, bool add)
{

if (add) {
printf("MSECriterium NOT READY TO ADD EVIDENCE!!!");
exit(1);
}

  int Nvar = data->GetNumVar();
  Instances& instances = data->GetM();

  n->d->Rleaf  = 0.0;
  n->fClass = 0.0;

  double MembershipTotal = 0;
  for(int i = n->first; i <= n->last; i++) {
    double v   = instances[i][Nvar-1];
    MembershipTotal += instances[i].GetMembership();
    n->d->Rleaf  += instances[i].GetMembership() * v * v;
    n->fClass += instances[i].GetMembership() * v;
  }

  n->fClass /= MembershipTotal;
  n->d->Rleaf  -= MembershipTotal * n->fClass * n->fClass;
  
  if (data->GetNPrune()) {/*CalcPruneError(n);*/}
  
}

double MultipleMSECriterium::MultipleVariance(std::vector<double> sum, std::vector<double> sumSquares, double numberInstances, int numberObjectives)
{
  std::vector<double> variances(numberObjectives, 0.0);
  for (int nObjective = 0; nObjective < numberObjectives; nObjective++) {
    variances[nObjective] = sumSquares[nObjective] - pow(sum[nObjective], 2) / numberInstances;
  } // s2R - s1R*s1R/nR
  return std::accumulate(variances.begin(), variances.end(), 0.0) / numberObjectives;
}

double MultipleMSECriterium::Initialise(int first, int last, int curSplit, double gamma, int iX, bool FuzzFlag)
{
  int numberObjectives = MultiobjectiveInstance::GetNumMultiobjectives();
  
  nL = nR = R[0] = R[1] = 0.0;
  s1L = std::vector<double>(numberObjectives, 0.0);
  s2L = std::vector<double>(numberObjectives, 0.0);
  s1R = std::vector<double>(numberObjectives, 0.0);
  s2R = std::vector<double>(numberObjectives, 0.0);

  // Initialise
  for(int i = first; i <= last; i++) {
    MultiobjectiveInstance& instance = dynamic_cast<MultiobjectiveInstance&> (data->GetInstance(i));
    if (instance.GetMembership() > 0.0) {
      if( (gamma > 1 && i <= curSplit) // normal
                                     || (gamma < 1 && instance[iX]+gamma < 0)) {
        for (int nObjective = 0; nObjective < numberObjectives; nObjective++) {
          s1L[nObjective] += instance.GetMembership() * instance.GetMultiobjective(nObjective);
          s2L[nObjective] += instance.GetMembership() * pow(instance.GetMultiobjective(nObjective), 2);
        }
        nL += instance.GetMembership();
      }
      else {
        for (int nObjective = 0; nObjective < numberObjectives; nObjective++) {
          s1R[nObjective] += instance.GetMembership() * instance.GetMultiobjective(nObjective);
          s2R[nObjective] += instance.GetMembership() * pow(instance.GetMultiobjective(nObjective), 2);
        }
        nR += instance.GetMembership();
      }
    }
  }

  if (nL > 0.0) {
    R[0] = MultipleMSECriterium::MultipleVariance(s1L, s2L, nL, numberObjectives);
  }

  if (nR > 0.0) {
    R[1] = MultipleMSECriterium::MultipleVariance(s1R, s2R, nR, numberObjectives);
  }

  return (nR < nL ? nR : nL);
}

double MultipleMSECriterium::MovePoint(int curSplit, double gamma, int iX,int first, int last)
{
  int numberObjectives = MultiobjectiveInstance::GetNumMultiobjectives();
  MultiobjectiveInstance& instance = dynamic_cast<MultiobjectiveInstance&> (data->GetInstance(curSplit));

  if(instance.GetMembership()>0.0) {
    if(gamma > 1) {
      for (int nObjective = 0; nObjective < numberObjectives; nObjective++) {
        s1L[nObjective] -= instance.GetMembership() * instance.GetMultiobjective(nObjective); 
        s2L[nObjective] -= instance.GetMembership() * pow(instance.GetMultiobjective(nObjective), 2);
        s1R[nObjective] += instance.GetMembership() * instance.GetMultiobjective(nObjective); 
        s2R[nObjective] += instance.GetMembership() * pow(instance.GetMultiobjective(nObjective), 2);
      }
      nL  -= instance.GetMembership(); 
      nR  += instance.GetMembership();
    }
    else {
      if(instance[iX]+gamma >=0) {
        for (int nObjective = 0; nObjective < numberObjectives; nObjective++) {
          s1L[nObjective] += instance.GetMembership() * instance.GetMultiobjective(nObjective);
          s2L[nObjective] += instance.GetMembership() * pow(instance.GetMultiobjective(nObjective), 2);
          s1R[nObjective] -= instance.GetMembership() * instance.GetMultiobjective(nObjective);
          s2R[nObjective] -= instance.GetMembership() * pow(instance.GetMultiobjective(nObjective), 2);
        }
        nL  += instance.GetMembership();
        nR  -= instance.GetMembership();
      }
      else {
        for (int nObjective = 0; nObjective < numberObjectives; nObjective++) {
          s1L[nObjective] -= instance.GetMembership() * instance.GetMultiobjective(nObjective);
          s2L[nObjective] -= instance.GetMembership() * pow(instance.GetMultiobjective(nObjective), 2);
          s1R[nObjective] += instance.GetMembership() * instance.GetMultiobjective(nObjective);
          s2R[nObjective] += instance.GetMembership() * pow(instance.GetMultiobjective(nObjective), 2);
        }
        nL  -= instance.GetMembership();
        nR  += instance.GetMembership();
      }
    }
  
    if (nL > 0.0) R[0] = MultipleMSECriterium::MultipleVariance(s1L, s2L, nL, numberObjectives);
    if (nR > 0.0) R[1] = MultipleMSECriterium::MultipleVariance(s1R, s2R, nR, numberObjectives);
  }

  return nR < nL ? nR : nL;
}

void MultipleMSECriterium::FixLeaf(Node* n, bool add)
{
  if (add) {
    printf("MultipleMSECriterium NOT READY TO ADD EVIDENCE!!!");
    exit(1);
  }
  
  int numberObjectives = MultiobjectiveInstance::GetNumMultiobjectives();
  Instances& instances = data->GetM();

  n->NodePop  = new double[numberObjectives];
  n->d->Rleaf = 0.0;

  double MembershipTotal = 0.0;
  for(int i = n->first; i <= n->last; i++) {
    MembershipTotal += instances[i].GetMembership();
  }
  
  for (int nObjective = 0; nObjective < numberObjectives; nObjective++) {
    n->NodePop[nObjective] = 0.0;
    for(int i = n->first; i <= n->last; i++) {
      MultiobjectiveInstance& instance = dynamic_cast<MultiobjectiveInstance&> (data->GetInstance(i));
      double currentMembership = instance.GetMembership();
      double v = instance.GetMultiobjective(nObjective);
      n->NodePop[nObjective] += currentMembership * v;
      n->d->Rleaf += currentMembership * pow(v, 2);
    }
    n->NodePop[nObjective] /= MembershipTotal;
    n->d->Rleaf -= MembershipTotal * pow(n->NodePop[nObjective], 2);
  }
  n->d->Rleaf /= numberObjectives;

  // Guardamos el valor predicho del ultimo objetivo en fClass
  n->fClass = n->NodePop[numberObjectives-1];

  if (data->GetNPrune()) {/*CalcPruneError(n);*/}
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
ComboCriterium::ComboCriterium(Data *data, std::vector<double> *alphas) 
{
  SetData(data);
  PopR = PopL = 0;
  hdeltaR = hdeltaL = 0;
  hscaleR = hscaleL = 0;
  if (alphas) {
    this->alphas = new std::vector<double>(*alphas);
  }
  else {
    this->alphas = 0;
  }
  alpha = 0;
}
ComboCriterium::~ComboCriterium()
{
  if (PopL) delete []PopL;
  if (PopR) delete []PopR;
  if (hdeltaL) delete hdeltaL;
  if (hdeltaR) delete hdeltaR;
  if (hscaleL) delete hscaleL;
  if (hscaleR) delete hscaleR;
  if (alpha) delete []alpha;
  if (alphas) delete alphas;
}

double ComboCriterium::Impurity(double *Pop, int NumClass, Matriz &delta, 
                                     Matriz &scale, double *alpha, int NumBins)
{
  double impurity;
  double TotPop = 0.0;

  // Combo criterion
  impurity = 0.0;

  for(int j = 0; j < NumClass; j++) {
    TotPop += Pop[j];
  }

#ifdef ALLINONE
    if (alpha[0] > 0.0) {
      impurity += alpha[0] * Pop[0] * (Entropy(scale[0], NumBins) + Entropy(delta[0], NumBins));
    }
#else
  for(int j = 0; j < NumClass; j++) {
    if (alpha[j] > 0.0) {
      impurity += alpha[j] * Pop[j] * (Entropy(scale[j], NumBins) + Entropy(delta[j], NumBins));
    }
  }
#endif

  impurity += alpha[NumClass] * TotPop * Entropy(Pop, NumClass);
 
  return impurity; 
}
double ComboCriterium::Entropy(double *Pop, int NumClass)
{

  double entropy = 0.0;
  double TotPop = 0.0;

  // Entropy criterion
  for(int j = 0; j < NumClass; j++) {
    TotPop += Pop[j];
  }
  for(int j = 0; j < NumClass; j++) {
    double pj = Pop[j]/TotPop;
    if (pj > 0.0) entropy += pj * log(pj) / log(2.0);
  }
 
  return -entropy; 
}
double ComboCriterium::Initialise(int first, int last, int curSplit, double gamma, int iX, bool )
{

  Instances& instances = data->GetM();
  int NumClass = ((NomData*)data)->NumClass;
  int NVar = data->GetNumVar();
  int NumBins = data->GetNumTerms(data->GetNumVarNom()-1);

  if (PopL) delete []PopL;
  if (PopR) delete []PopR;
  if (hdeltaL) delete hdeltaL;
  if (hdeltaR) delete hdeltaR;
  if (hscaleL) delete hscaleL;
  if (hscaleR) delete hscaleR;
  if (alpha) delete []alpha;

  PopL = new double[ NumClass ];
  PopR = new double[ NumClass ];
  hdeltaL = new Matriz( NumClass, NumBins );
  hdeltaR = new Matriz( NumClass, NumBins );
  hscaleL = new Matriz( NumClass, NumBins );
  hscaleR = new Matriz( NumClass, NumBins );

  alpha = new double[ NumClass + 1];
  
  for (int j = 0; j < NumClass; j++) {
    PopL[j] = PopR[j] = 0.0;
    alpha[j] = alphas ? (*alphas)[j] : log(2.0) / log(NumBins);
  }
  alpha[NumClass] = alphas ? (*alphas)[NumClass] : log(2.0) / log(NumClass);

  R[0] = R[1] = nL = nR = 0.0;       // Initialise
  for (int i = first; i <= last; i++) {
    int ic = instances[i].GetClass();
    if (instances[i].GetMembership() > 0.0) {
      if ( (gamma > 1 && i <= curSplit) // normal
                               || (gamma < 1 && instances[i][iX]+gamma < 0)) {
        PopL[ic] += instances[i].GetMembership();
#ifdef ALLINONE
        (*hdeltaL)[0][(int)instances[i][NVar-2]] += instances[i].GetMembership();
        (*hscaleL)[0][(int)instances[i][NVar-3]] += instances[i].GetMembership();
#else
        (*hdeltaL)[ic][(int)instances[i][NVar-2]] += instances[i].GetMembership();
        (*hscaleL)[ic][(int)instances[i][NVar-3]] += instances[i].GetMembership();
#endif
        nL += instances[i].GetMembership();
      }
      else {
        PopR[ic] += instances[i].GetMembership();
#ifdef ALLINONE
        (*hdeltaR)[0][(int)instances[i][NVar-2]] += instances[i].GetMembership();
        (*hscaleR)[0][(int)instances[i][NVar-3]] += instances[i].GetMembership();
#else
        (*hdeltaR)[ic][(int)instances[i][NVar-2]] += instances[i].GetMembership();
        (*hscaleR)[ic][(int)instances[i][NVar-3]] += instances[i].GetMembership();
#endif
        nR += instances[i].GetMembership();
      }
    }
  }
  

  if (nL > 0.0) {
    R[0] = Impurity(PopL, NumClass, *hscaleL, *hdeltaL, alpha, NumBins);
  }
  if (nR > 0.0) {
    R[1] = Impurity(PopR, NumClass, *hscaleL, *hdeltaL, alpha, NumBins);
  }

  return (nR < nL ? nR : nL);
}
double ComboCriterium::MovePoint(int curSplit, double gamma, int iX,int first,int last)
{
  Instances& instances = data->GetM();
  int NVar = data->GetNumVar();
  int NumBins = data->GetNumTerms(data->GetNumVarNom()-1);

  if (instances[curSplit].GetMembership() > 0.0) {
    int ic = instances[curSplit].GetClass();
    if (gamma > 1 || (gamma<1 && instances[curSplit][iX]+gamma <0)) {
      PopL[ic] -= instances[curSplit].GetMembership();
      PopR[ic] += instances[curSplit].GetMembership();

#ifdef ALLINONE
      (*hdeltaL)[0][(int)instances[curSplit][NVar-2]] -= instances[curSplit].GetMembership();
      (*hdeltaR)[0][(int)instances[curSplit][NVar-2]] += instances[curSplit].GetMembership();

      (*hscaleL)[0][(int)instances[curSplit][NVar-3]] -= instances[curSplit].GetMembership();
      (*hscaleR)[0][(int)instances[curSplit][NVar-3]] += instances[curSplit].GetMembership();
#else
      (*hdeltaL)[ic][(int)instances[curSplit][NVar-2]] -= instances[curSplit].GetMembership();
      (*hdeltaR)[ic][(int)instances[curSplit][NVar-2]] += instances[curSplit].GetMembership();

      (*hscaleL)[ic][(int)instances[curSplit][NVar-3]] -= instances[curSplit].GetMembership();
      (*hscaleR)[ic][(int)instances[curSplit][NVar-3]] += instances[curSplit].GetMembership();
#endif

      nL -= instances[curSplit].GetMembership();
      nR += instances[curSplit].GetMembership();
    } 
    else {
      PopL[ic] += instances[curSplit].GetMembership();
      PopR[ic] -= instances[curSplit].GetMembership();  

#ifdef ALLINONE
      (*hdeltaL)[0][(int)instances[curSplit][NVar-2]] += instances[curSplit].GetMembership();
      (*hdeltaR)[0][(int)instances[curSplit][NVar-2]] -= instances[curSplit].GetMembership();

      (*hscaleL)[0][(int)instances[curSplit][NVar-3]] += instances[curSplit].GetMembership();
      (*hscaleR)[0][(int)instances[curSplit][NVar-3]] -= instances[curSplit].GetMembership();
#else
      (*hdeltaL)[ic][(int)instances[curSplit][NVar-2]] += instances[curSplit].GetMembership();
      (*hdeltaR)[ic][(int)instances[curSplit][NVar-2]] -= instances[curSplit].GetMembership();

      (*hscaleL)[ic][(int)instances[curSplit][NVar-3]] += instances[curSplit].GetMembership();
      (*hscaleR)[ic][(int)instances[curSplit][NVar-3]] -= instances[curSplit].GetMembership();
#endif

      nL += instances[curSplit].GetMembership();
      nR -= instances[curSplit].GetMembership();
    } 
    
    R[0] = R[1] = 0.0;
  
    if(nL > 0.0) R[0] = Impurity(PopL, ((NomData*)data)->NumClass, *hscaleL, *hdeltaL, alpha, NumBins);
    if(nR > 0.0) R[1] = Impurity(PopR, ((NomData*)data)->NumClass, *hscaleL, *hdeltaL, alpha, NumBins);
  }

  return (nR < nL ? nR : nL);
}

void ComboCriterium::FixLeaf(Node* n, bool add)
{
if (add) {
printf("ComboCriterium NOT READY TO ADD EVIDENCE!!!");
exit(1);
}

  Instances& instances = data->GetM();
  int NumClass = ((NomData*)data)->NumClass;
  int NVar = data->GetNumVar();
  int NumBins = data->GetNumTerms(data->GetNumVarNom()-1);

  if (NumClass) {
    if (n->NodePop) delete []n->NodePop;
    n->NodePop = new double[NumClass];
    for(int i = 0; i < NumClass; i++ ) n->NodePop[i] = 0.0;
  }
  n->d->Rleaf = 0.0;
  n->iClass = 0;

  ((NomData*)data)->FindPops(n->NodePop, n->first, n->last);
  n->iClass = ((NomData*)data)->WhichClass(n->NodePop);

  int count_positives = 0;
  for (int j = 0; j < NumClass; ++j) {
    if(j != n->iClass) {
      n->d->Rleaf += n->NodePop[j] /3.0;
    }

    if (j > 0) {
      count_positives += n->NodePop[j];
    }
  }

  /*ComboSplitNodeInfo *csni = new ComboSplitNodeInfo();
  if (n->parent==0) { //root node
    csni->indexes = new int[count_positives];
  }
  else if (n->parent->child == n) {
    csni->indexes = ((ComboSplitNodeInfo*)n->parent->info)->indexes;
  }
  else {
    int despl = ((ComboSplitNodeInfo*)n->parent->child->info)->count;
    csni->indexes = ((ComboSplitNodeInfo*)n->parent->child->info)->indexes + despl;
  }

  csni->count = count_positives;
  for (int i = n->first, k = 0; i <= n->last; i++) {
    if (instances[i].GetClass() > 0) {
      csni->indexes[k] = instances[i].GetIniPos();
      k++;
    }
  }
  if (n->parent && n->parent->info) {
    delete (ComboSplitNodeInfo *)n->parent->info;
    n->parent->info = 0;
  }

  n->info = csni;
    */ 

  ComboSplitNodeInfo *csni = new ComboSplitNodeInfo();

  if (n->parent && n->parent->info) {
    delete ((ComboSplitNodeInfo *)n->parent->info)->instances_data;
    delete (ComboSplitNodeInfo *)n->parent->info;
    n->parent->info = 0;
  }


//  csni->hist_scale = new Matriz( NumClass, NumBins );
//  csni->hist_delta = new Matriz( NumClass, NumBins );
  csni->instances_data = new Matriz( count_positives, 7);

  Matriz pop_scale(1, NumBins);
  Matriz pop_delta(1, NumBins);
  for (int i = n->first, ip = 0; i <= n->last; i++) {
    int ic = instances[i].GetClass();
  //  (*csni->hist_delta)[ic][(int)instances[i][NVar-2]] += instances[i].GetMembership();
    //(*csni->hist_scale)[ic][(int)instances[i][NVar-3]] += instances[i].GetMembership();
    if ( ic > 0 ) {
      pop_delta[0][(int)instances[i][NVar-2]]++;     //Binning deltas
      pop_scale[0][(int)instances[i][NVar-3]]++;     //Binning scales
      (*csni->instances_data)[ip][0] = instances[i][NVar-7]; //Delta x to obj
      (*csni->instances_data)[ip][1] = instances[i][NVar-6]; //Delta y to obj
      (*csni->instances_data)[ip][2] = instances[i][NVar-9]; //Detection x scale
      (*csni->instances_data)[ip][3] = instances[i][NVar-8]; //Detection y scale
      (*csni->instances_data)[ip][4] = ic;
      //(*csni->instances_data)[ip][5] = instances[i][NVar-3];
      //(*csni->instances_data)[ip][6] = instances[i][NVar-2];
      (*csni->instances_data)[ip][5] = instances[i].GetIniPos();
      (*csni->instances_data)[ip][6] = instances[i][NVar-2];
      ip++;
    }
  }


  int i_class_delta = ((NomData*)data)->WhichClass(pop_delta[0], NumBins);
  int i_class_scale = ((NomData*)data)->WhichClass(pop_scale[0], NumBins);

  for (int i = 0; i < NumBins; i++) {
    if (i_class_delta != i) {
      n->d->Rleaf += pop_delta[0][i]/3.0;
    }
    if (i_class_scale != i) {
      n->d->Rleaf += pop_scale[0][i]/3.0;
    }
  }

  if (i_class_delta==0) {
    delete csni->instances_data;
    csni->instances_data = new Matriz( 0, 7);
  }
  n->info = csni;
//for (int i=0;i<n->depth;i++) cout << "====";
//cout << "D" << endl;
//csni->hist_scale->saveToStream(cout);
//csni->hist_delta->saveToStream(cout);
}

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


NomData::NomData(char *filename, Instance *builder)
  : Data(filename, builder)
{
  init();
}
NomData::NomData(int NTotalData, int NVar, Instance *builder)
  : Data(NTotalData, NVar, builder)
{
  init();
}
NomData::NomData(int NTotalData, int NVar, double**data, Instance *builder)
  : Data(NTotalData, NVar, data, builder)
{
  init();
}

NomData::~NomData()
{
  delete[] Pop;
//  delete[] PopR;
//  delete[] PopL;
  delete splitCriterium;
}

void NomData::init()
{
  atdrawchosefathers = false;
  NumClass = NomTerms[Nnom].size();

//  printf("Class variables to %d\n", NumClass);
  Pop = new double[NumClass];
  if (DepVarType==NOM) {
    splitCriterium = new GiniCriterium(this);
  }
  else if (DepVarType==ORD) {
    splitCriterium = new MSECriterium(this);
  }
  else {
    splitCriterium = 0;
    cout << "No split criterium set." << endl;
  }
//  PopR = new double[NumClass];
//  PopL = new double[NumClass];
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
    instances[i].Membership = /*instances[i].Weight =*/1.0;//???????????????


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
      instances[i].Weight = weight[instances[i].Class];
      if(i<Num)
        total += instances[i].Weight;
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
    instances[i].Membership = 1.0;


  FindPops(class_weights, 0, NTotal-1);

  for(int i = 0; i< NumClass; ++i)
    class_weights[i] = NTotal/(NumClass*class_weights[i]);

  for(int i = 0; i < NTotal ; ++i) {
    instances[i].Weight = class_weights[ instances[i].Class ];
  }

  delete[] class_weights;

}

void NomData::FindPops(double* HoldPop, int first, int last)
{
  for (int j=0; j < NumClass; j++)
    HoldPop[j]=0.0;

  for(int i=first; i<= last; i++)
    HoldPop[instances[i].GetClass()] += instances[i].Membership;
}

int NomData::WhichClass(double* pt, int NumClass)
{

  double max=0.0;
  int maxj = -1;

  if (NumClass <= 0) NumClass = this->NumClass;

  for(int j = 0; j < NumClass; j++, pt++) {
    if (j == 0 || *pt>max) {
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
      double aux2 = ((i==instances[i].Class) - Pop[i]/TotPop);
      aux1 += aux2*aux2;

    }
    aux += aux1*instances[i].Membership;
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
      double aux2 = ((i==instances[i].Class) - Pop[i]/TotPop);
      aux1 += aux2*aux2;

    }
    aux += aux1*instances[i].Membership;
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


//double NomData::Initialise(int first, int last, int curSplit, double gamma, int iX, bool FuzzFlag)
//{

//  cout << "Shouldn't come here" << endl;
  //double min_in_node = splitCriterium->Initialise(first, last, curSplit, gamma, iX, FuzzFlag);
//  R[0] = 0.0;
//  R[1] = splitCriterium->GetLastSplitImpurity();
//  return 0.0; // min_in_node;

/*  bool flagleft = true;

  if(FuzzFlag) {
    if(curSplit == last) {
      nL = R[0] = 0.0;
    }
    else {  
      flagleft = false; 
      nR = R[1] = 0.0;
    } 
  }
  else {
    R[0] = R[1] = nL = nR = 0.0;                 // Initialise
  }
  
  for (int j=0; j < NumClass; j++)
    PopL[j] = PopR[j] = 0.0;

  for(int i=first;i<=last;i++) {
    if(instances[i].Membership > 0.0) {
      if( (gamma > 1 && i <= curSplit) // normal
        || (gamma < 1 && instances[i][iX]+gamma < 0))
      {
        PopL[instances[i].Class] += instances[i].Membership;
        nL += instances[i].Membership;
      } else {
        PopR[instances[i].Class] += instances[i].Membership;
        nR += instances[i].Membership;
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
  
  return nR;*/ 
//}

//double NomData::MovePoint(int curSplit, double gamma, int iX,int first,int last)
//{
//  cout << "Shouldn't come here" << endl;
//  double min_in_node = splitCriterium->MovePoint(curSplit, gamma, iX, first, last);
//  R[0] = 0.0;
//  R[1] = splitCriterium->GetLastSplitImpurity();
//  return 0.0;//min_in_node;

/*
  if(instances[curSplit].Membership > 0.0)
  {
    if(gamma > 1 || (gamma<1 & instances[curSplit][iX]+gamma <0))
    {
      PopL[instances[curSplit].Class] -= instances[curSplit].Membership;
      PopR[instances[curSplit].Class] += instances[curSplit].Membership;

      nL -= instances[curSplit].Membership;
      nR += instances[curSplit].Membership;

    } else {
      PopL[instances[curSplit].Class] += instances[curSplit].Membership;
      PopR[instances[curSplit].Class] -= instances[curSplit].Membership;  

      nL += instances[curSplit].Membership;
      nR -= instances[curSplit].Membership;

    } 

    
    R[0] = R[1] = 0.0;
  
    if(nL > 0.0) R[0] = FindImpurity(PopL,first,last);
    if(nR > 0.0) R[1] = FindImpurity(PopR,first,last);
  }

  return (nR < nL ? nR : nL);*/
//}

void NomData::FixLeaf(Node* n, bool add)
{
  splitCriterium->FixLeaf(n, add);

/*  if(NumClass) {
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
*//*    char puf[128];
    sprintf(puf, "JAJAJAJAJALLL%f", (float)n->parent->fSplit);
    if (console) (*console)(puf);
  *//*}

  for (int j=0; j<NumClass; ++j)
  {
    if(j != n->iClass)
      n->d->Rleaf += n->NodePop[j];
  }

*/  /*// ASG: Gain-ratio//////////
  n->entropy = FindImpurity(n->NodePop);
  // Gain-ratio //////////////*/

/*
  if(NPrune)
  {
    n->d->RPrune = 0;
    n->d->firstTest = Ndata; n->d->lastTest = Ndata+ NPrune - 1;
    FindMembershipTest(n);
    for(int i = n->d->firstTest; i <= n->d->lastTest; ++i)
    {
      if(instances[i].Class != n->iClass)
        n->d->RPrune += instances[i].Membership;
    }
  }*/
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

        if(instances[i].GetClass() != cur->iClass)
          cur->d->RleafCV += instances[i].Weight;  

        if(cur->child && cur->SplitType == ORD)
        {
          instances[i](Nvar, -cur->fSplit);
          if(cur->coeff)
          {
            for(int j=0;j<NVarsOrdSplit;j++) 
              instances[i](Nvar, instances[i][Nvar] + cur->coeff[j]*instances[i][j]);
          } else 
            instances[i](Nvar, instances[i][Nvar] + instances[i][cur->att]);
          
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
            SortOn(Nvar, cur->d->firstTest, cur->d->lastTest);
            cur->child->d->firstTest = cur->child->d->lastTest = cur->d->firstTest;
            cur->child->d->lastTest--;
            while(cur->child->d->lastTest < cur->d->lastTest && instances[cur->child->d->lastTest+1][Nvar] < 0.0)
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
        if(instances[i].GetClass() != cur->iClass && instances[i].Membership > 0.0)
          cur->d->RleafCV += instances[i].Membership;  
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
  if (n->IsLeaf(K)) {
    if (n->parent && n->parent->child == n) {
      if (MultiobjectiveInstance::GetNumMultiobjectives() == -1) {
        sprintf(buf,"%s is ", VarNames[Nvar-1].c_str());
        if (DepVarType==NOM) {
          sprintf(buf+strlen(buf),"%s ", NomTerms[Nnom].at(n->iClass).c_str());
        } else {
          sprintf(buf+strlen(buf),"%f ", n->fClass);
        }
      } else {
        if (DepVarType==NOM) {
          throw std::invalid_argument("Can't describe multi-objective nominal nodes");
        } else {
          int first_objective = Nvar - MultiobjectiveInstance::GetNumMultiobjectives();
          for (int i=first_objective; i<Nvar; i++) {
            sprintf(buf+strlen(buf),"%s is ", VarNames[i].c_str());
            sprintf(buf+strlen(buf),"%f ", n->NodePop[i-first_objective]);
            if (i != Nvar - 1) {
              sprintf(buf+strlen(buf),"and ");
            }
          }
        }
      }
      sprintf(buf+strlen(buf),"(K=%d) (%g", n->K, n->d ? n->d->Rleaf : 0);
      sprintf(buf+strlen(buf),"(%g)/", n->d ? n->d->RsubTree : 0);
      sprintf(buf+strlen(buf),"%d) ", n->last - n->first + 1);
      sprintf(buf+strlen(buf),"(a=%g) T=%d\n", n->alpha, n->T);
    }
    else {
      if (MultiobjectiveInstance::GetNumMultiobjectives() == -1) {
        if (DepVarType==NOM) {
          sprintf(buf,"ELSE var %s is %s (K=%d) (%g(%g)/%d) (a=%g) T=%d\n",
                          VarNames[Nvar-1].c_str(),
                          NomTerms[Nnom].at(n->iClass).c_str(), 
                          n->K, 
                          n->d ? n->d->Rleaf : 0, 
                          n->d ? n->d->RsubTree : 0, 
                          (n->last - n->first + 1),
                          n->alpha, n->T);
        } else {
          sprintf(buf,"ELSE var %s is %f (K=%d) (%g(%g)/%d) (a=%g) T=%d\n",
                          VarNames[Nvar-1].c_str(),
                          n->fClass, 
                          n->K, 
                          n->d ? n->d->Rleaf : 0, 
                          n->d ? n->d->RsubTree : 0, 
                          (n->last - n->first + 1),
                          n->alpha, n->T);
        }
      } else {
        if (DepVarType==NOM) {
          throw std::invalid_argument("Can't describe multi-objective nominal nodes");
        } else {
          int first_objective = Nvar - MultiobjectiveInstance::GetNumMultiobjectives();
          sprintf(buf,"ELSE ");
          for (int i=Nvar-MultiobjectiveInstance::GetNumMultiobjectives(); i<Nvar; i++) {
            sprintf(buf+strlen(buf),"var %s is %f ",
                          VarNames[i].c_str(),
                          n->NodePop[i-first_objective]);
            if (i != Nvar - 1) {
              sprintf(buf+strlen(buf),"and ");
            }
          }
          sprintf(buf+strlen(buf),"(K=%d) (%g(%g)/%d) (a=%g) T=%d\n",
                      n->K, 
                      n->d ? n->d->Rleaf : 0, 
                      n->d ? n->d->RsubTree : 0, 
                      (n->last - n->first + 1),
                      n->alpha, n->T);
        }
      }
    }
  } 
  else if (n->SplitType==NOM)  {
    if (n->parent) {
      if(n->parent->child == n) {
        sprintf(buf,"IF the value of (%s) is in ",VarNames[n->att].c_str());
      }
      else {
        sprintf(buf,"ELSE IF the value of (%s) is in ",
                                                     VarNames[n->att].c_str());
      }
    }
    else {
      sprintf(buf,"IF the value of %s is in ",VarNames[n->att].c_str());
    }

    int flag =0;
    for(unsigned i=0;i<NomTerms[n->att-Nordfuzz].size(); ++i)  {
      if(1<<i & n->NomSplit) {
        if(flag==0) {
          sprintf(buf+strlen(buf),"{ %s",
                                      NomTerms[n->att-Nordfuzz].at(i).c_str());
          flag=1;
        }
        else {
          sprintf(buf+strlen(buf),", %s", 
                                      NomTerms[n->att-Nordfuzz].at(i).c_str());
        }
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
    else { 
      sprintf(buf,"IF (%g) ",-(n->fSplit));
    }
    if (n->coeff) {
      for(int i=0;i<NVarsOrdSplit;i++) {
        if(fabs(n->coeff[i]) > 0.0) {
          sprintf(buf+strlen(buf)," + (%g) %s",n->coeff[i],VarNames[i].c_str());
        }
      }
    }
    else {
      sprintf(buf+strlen(buf)," + %s", VarNames[n->att].c_str());
    }
    sprintf(buf+strlen(buf)," < 0 (K=%d) (%g(%g)/%d) (a=%g) T=%d\n",
        n->K, n->d ? n->d->Rleaf : 0, n->d ? n->d->RsubTree : 0, (n->last - n->first + 1), n->alpha, n->T);
  }
  else if(n->SplitType == FUZZ) {
    if(n->d && n->d->fuzzindex) {
      if(n->parent)  {
        if(n->parent->child == n) {
          sprintf(buf,"IF fuzzification(%d) of (%g) ", n->d->fuzzindex,
                                                                 -(n->fSplit));
        }
        else {
          sprintf(buf,"ELSE IF fuzzification(%d) of (%g) ", n->d->fuzzindex,
                                                                 -(n->fSplit));
        }
      }
    else {
      sprintf(buf,"IF fuzzification(%d) of (%g) ",n->d->fuzzindex,-(n->fSplit));
    }

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
        if(n->parent->child == n) {
          sprintf(buf,"IF the value of (%s) is %s", VarNames[n->att].c_str(), 
                                              n->d->FuzzSplit->GetName().c_str());
        }
        else {
          sprintf(buf,"ELSE IF the value of (%s) is  %s", 
                    VarNames[n->att].c_str(), n->d->FuzzSplit->GetName().c_str());
        }
      }
      else {
        sprintf(buf,"IF the value of %s is %s", VarNames[n->att].c_str(), 
                                              n->d->FuzzSplit->GetName().c_str());
      }
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
    int y= instances[i].GetClass();
    int ypred; 
    
    while(cur)
    {
     //  Node* temp = cur;

      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {
        cur->d->firstTest = cur->d->lastTest = i;
        FindMembershipTest(cur);

        unity += instances[i].Membership;

        ypred = cur->iClass;
        if(ypred != y)
          ress += instances[i].Membership;
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
        double aux = cur->d->tempMU*instances[i].Weight;
        if(cur->parent) aux *= cur->parent->d->cumMU;
        int Class = instances[i].Class;
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

        double mu =  cur->ParamMem(instances[i],NVarsOrdSplit);
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
    int Class = instances[i].Class;
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

        double mu = cur->ParamMem(instances[i],NVarsOrdSplit);
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
      cur->FixDers(instances[i],K,NVarsOrdSplit);
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
    int Class = instances[i].GetClass();
    while(cur)  {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf) {
        if(Class != cur->iClass) {
          double aux = instances[i].Weight*cur->d->tempMU;
          if(cur->parent) aux *= cur->parent->d->cumMU;
          ss += aux;
        }
      }
      else {
        double mu =  cur->ParamMem(instances[i],NVarsOrdSplit);
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
    total += instances[i].Weight;
  }
  return ss/total;  // if all data points have weight (instances[i].Weight =1),
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
    int Class = instances[i].Class;    
    instances[i].CvTreeLabel = 0.0;
    while(cur)
    {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {    
        if(Class != cur->iClass)
        {  
          
          double aux = instances[i].Weight*cur->d->tempMU;
          if(cur->parent) aux *= cur->parent->d->cumMU;
          ss += aux;
          if(aux>0.99) instances[i].CvTreeLabel = 1.0;
        }
      } else {

        double mu =  cur->ParamMem(instances[i],NVarsOrdSplit);
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
    total += instances[i].Weight;
  }  
  return ss/total;  // if all data points have weight (instances[i].Weight =1),
            // total = end - begin +1
}

double NomData::Fuzziness(int i, double fuzz)
{

  int N = NTrain + NSel;

  if(i<N) 
    return fuzz;
  if(instances[i][Nvar] > instances[N-1][Nvar])
    return instances[N-1][Nvar+1];
  if( instances[i][Nvar] <  instances[0][Nvar])
    return 0;

  double val = N/2.0;
  double delta = val;

  while((int)(delta /= 2.0) >= 1)
  {
    if(instances[i][Nvar] > instances[(int)val][Nvar])
      val += delta;
    else val -= delta;
  }
  
  return instances[(int)val][Nvar+1];
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
        double aux = cur->d->tempMU*instances[i].Weight;
        if(cur->parent) aux *= cur->parent->d->cumMU;
      
        total += aux;
        fuzz += aux*(1.0-aux);

        if(i==begin)
        {
          for(int j=0; j<NumClass; ++j)
            cur->NodePop[j]  = 0.0;
        }
        
        cur->NodePop[instances[i].Class] += aux;
        
        if (i == end)   
          cur->iClass = WhichClass(cur->NodePop);
          
        
        
      } else {

        double mu =  cur->ParamMem(instances[i],NVarsOrdSplit);
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
    int Class = instances[i].Class;    
    while(cur)
    {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {
        double aux = cur->d->tempMU*instances[i].Weight;
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

        double mu =  cur->ParamMem(instances[i],NVarsOrdSplit);
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
    int Class = instances[i].Class;    
    while(cur)
    {
      bool curIsLeaf = cur->IsLeaf(K);
      if(curIsLeaf)
      {
        double aux = cur->d->tempMU*instances[i].Weight;
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
        double mu = cur->ParamMem(instances[i],NVarsOrdSplit);
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
  SortOn(0, 0, NTotal);
  if(index <0)
  {
    for(int i=0; i<NTotal; ++i)
      Class[i] = instances[i].Class;
    return;
  }

  for(int i=0; i < NTotal; ++i)
    instances[i].Class = (int)(0.5 + ((int)Class[i] == index));
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

        double mu =  cur->ParamMem(instances[i],NVarsOrdSplit);
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

  while(cur) {

    if (cur->IsLeaf(K)) {
      if (finalnode) *finalnode = cur;
      break;
    }

    left = GoLeft(cur, idat);

    cur = left ? cur->child : cur->child->sib;
  }

  return cur->iClass;
}
//Calculate the prune errors and pops forthetree, it alse rearranges the data
void NomData::CalculateAndArrangeTestData(Node * root)
{
/*  n->RPrune = 0;
  n->d->firstTest = Ndata; n->d->lastTest = Ndata+ NPrune - 1;
  FindMembershipTest(n);
  for(int i = n->d->firstTest; i <= n->d->lastTest; ++i)
  {
    if(instances[i].Class != n->iClass)
      n->RPrune += instances[i].Membership;
  }*/
}
double NomData::NextDifferentValue(int var, int index)
{
  double val = instances[index][var];
  for(int i=index+1;i<NTotal;i++) {
    if (val!=instances[i][var]) return instances[i][var];
  }
  return val;
}
NomData *NomData::GenerateSinteticData(int magnitude, int ini, int fin)
{
  int ntot = (fin-ini+1)*magnitude;
  NomData *nd = new NomData(ntot, Nvar, instanceBuilder);
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
  for( int i = ini; i <= fin; i++ ) {
    for( int j = 0; j < magnitude; j++ ) {
      nd->instances[i*magnitude+j] = instances[i];
    }
  }
  //Desvio los datos entre su valor superior e inferior
  double desde, rango, ant;
//  nd->SortOn(0, 0, ntot-1);
//  rango = nd->NextDifferentValue(0, magnitude) - nd->m[0][0];
//  desde = nd->m[0][0] - rango;
//  rango *= 2;
  int k;
  for(k = 0; k < nd->Nordfuzz; k++) {
    nd->SortOn(k, 0, ntot-1);
    SortOn(k, ini, fin);
    rango = NextDifferentValue(k, ini) - instances[ini][k];
    desde = instances[ini][k] - rango;
    ant = instances[ini][k];
    rango *= 2.0;
    for(int i=ini;i<=fin;i++) {
      for(int j=0;j<magnitude;j++) {
        nd->instances[i*magnitude+j](k, drand()*rango + desde);
      }
      desde = ant;
      if (i<fin-1) {
        rango = NextDifferentValue(k, i+1) - ant;
        ant = instances[i+1][k];
      }
      else {
        rango = (instances[fin][k] - ant)*2.0;
      }
    }
  }
  for(;k<Nvar-1;k++) {
  }


  nd->Pop  = new double[NumClass];
//  nd->PopR = new double[NumClass];
//  nd->PopL = new double[NumClass];
  nd->ponder = false;
  nd->Ponderate(Ndata);
  nd->iClass = WhichClass(nd->Pop);
  return nd;
}
NomData *NomData::GenerateSinteticData2(int ntot, int ini, int fin)
{
  NomData *nd = new NomData(ntot, Nvar, instanceBuilder);
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
  for(k = 0; k < nd->Nordfuzz; k++) {
    SortOn(k, ini, fin);
    for(int j=0;j<ntot;j++) {
      int idx = ini + (int)(drand()*(fin-ini-1));
      desde = instances[idx][k];
      rango = instances[idx+1][k] - instances[idx][k];
      nd->instances[j](k, drand()*rango + desde);
    }
  }
  for(;k<Nvar-1;k++) {
  }


  nd->Pop  = new double[NumClass];
//  nd->PopR = new double[NumClass];
//  nd->PopL = new double[NumClass];
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

void NomData::RemoveInstancesOfClass(int class_index)
{
  int instances_to_remove = 0;

  for(int i = 0; i < NTotal ; ++i) {
    if (instances[i].GetClass() == class_index) {
      instances[i](Nvar, 1);
      instances_to_remove++;
    }
    else {
      instances[i](Nvar, 0);
    }
  }
  SortOn(Nvar);

  redim(NTotal - instances_to_remove);

  OriginalOrder();
  for(int i = 0; i < NTotal; i++) {
    instances[i].IniPos = i;
  }
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
      int c = (int)( nc * TDebugRand::Rand() / (RAND_MAX + 1.0));
      if (d->GetDatClass(i) <= c) c++;
//cout<<d->m[i][Nvar-1]<<"->" << c << endl;
      d->instances[i].SetClass(c); //Aqui se cambia la clase
//k++;
    }
  }
//cout<<"tot " << k << "/" << d->GetNTotal() << endl;
  return d;
}

//------------------------------------  NomData::SubstituteClassLabels  -----
//---------------------------------------------------------------------------
void NomData::SubstituteClassLabels(Data *data, int attribute_idx)
{
  if (data->GetNTotal() != GetNTotal()) {
    cout << "Different Number of instances: Nothing done." << endl;
    return;
  }

  if (attribute_idx < 0) {
    attribute_idx = data->GetNumVar()-1;
  }
  bool IsNom = attribute_idx < data->GetNordfuzz() ? false : true;

  NomTerms[GetNumVarNom()].clear();
  if (IsNom) {
    NomTerms[GetNumVarNom()] = ((NomData*)data)->NomTerms[attribute_idx - data->GetNordfuzz()];
    for(int i=0; i < GetNTotal(); i++) {
      cout << "Previous class: " << instances[i].GetClass();
      instances[i].SetClass((int)data->GetValueVar(i, attribute_idx));
      cout << " New: " << instances[i].GetClass() << endl;
    }
  }
  else {
    for(int i=0; i < GetNTotal(); i++) {
      ostringstream os;
      os << data->GetValueVar(i, attribute_idx);
      instances[i].SetClass((int)DameRepresentacion(os.str(), GetNumVarNom()));
    }
  }
  NumClass = NomTerms[GetNumVarNom()].size();
  delete[] Pop;
//  delete[] PopR;
//  delete[] PopL;
  delete splitCriterium;
  init();
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//--------------------------------------------------------------  Data  -----
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


Data::Data(char *filename, Instance *builder)
{
  instanceBuilder =  builder ? builder->NewInstance(0) : new BasicInstance(0);

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
Data::Data(int NTotalData, int _NVar, Instance *builder)
{
  instanceBuilder =  builder ? builder->NewInstance(0) : new BasicInstance(0);

  DefaultInitValues();
  SetDefaultNames(_NVar);
  redim(NTotalData);
//  init(NTotalData, _NVar);
}

Data::Data(int NTotalData, int _NVar, double **data, Instance *builder)
{
  instanceBuilder =  builder ? builder->NewInstance(0) : new BasicInstance(0);

  DefaultInitValues();
  SetDefaultNames(_NVar);
  redim(NTotalData);

  for(int i=0;i<NTotal;i++) {
    for (int j=0;j<Nvar-1;j++) instances[i](j, data[i][j]);
    //Depvar
    char txt[256];
    sprintf(txt, "%g", data[i][Nvar-1]);
    instances[i].SetClass((int)DameRepresentacion(txt, 0));
  }
}


Data::~Data()
{
  delete instanceBuilder;

  for(unsigned i=0;i<reg_eventos.size();i++) reg_eventos[i]->OnDelete(this);

  instances.clear();
//  for(int i = 0; i < NTotal; i++) {
//    delete[] m[i];
//  }

//  delete[] m; m = NULL;
  delete[] vshift; vshift = NULL;
  delete[] vscale; vscale = NULL;
  delete[] Gvshift; Gvshift = NULL;
  delete[] Gvscale; Gvscale = NULL;
//  if(VarNames) delete[] VarNames; VarNames =NULL;
//  if(NomTerms) delete[] NomTerms;  NomTerms =NULL;
  if (GroupCount) delete []GroupCount;

  delete []NomTerms;
  delete []VarNames;
}
double Data::GetValueVar(int i, int j)
{
  if (j < Nvar-1)  return instances[i][j]; 
  if (j == Nvar-1) return instances[i].GetClass();
  if (j == Nvar)   return instances[i][Nvar];
  if (j == Nvar+1) return instances[i][Nvar+1];
  if (j == Nvar+2) return instances[i].CvTreeLabel;
  if (j == Nvar+3) return instances[i].Membership;
  if (j == Nvar+4) return instances[i].Weight;
  if (j == Nvar+5) return instances[i].IniPos;
  /*if (j == Nvar+6)*/ return instances[i].Group;
}
void Data::SetValueVar(int i, int j, double value)
{
  if (j < Nvar-1)  instances[i](j, value); 
  if (j == Nvar-1) instances[i].SetClass((int)value);
  if (j == Nvar)   instances[i](Nvar, value);
  if (j == Nvar+1) instances[i](Nvar+1, value);
  if (j == Nvar+2) instances[i].CvTreeLabel = value;
  if (j == Nvar+3) instances[i].Membership = value;
  if (j == Nvar+4) instances[i].Weight = value;
  if (j == Nvar+5) instances[i].IniPos = (int)value;
  if (j == Nvar+6) instances[i].Group = (int)value;
}
void Data::SetDefaultNames(int nvar)
{
  int numberObjectives = MultiobjectiveInstance::GetNumMultiobjectives();
  bool multiObjective = numberObjectives > 0;

  //Variables "horizontales" (o las que tienen que ver con atributos...)
  Nvar = nvar;
  Nord = multiObjective ? Nvar-numberObjectives : Nvar-1;
  Nfuzz=Nnom=0;
  Nordfuzz = Nord+Nfuzz;
  NVarsOrdSplit = Nord;
  DepVarType = NOM;
  OriginalVarType.resize(Nvar);
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
  instances.clear();
  NTotal = 0;
  //Variables "verticales" (o las que tienen que ver con no. ejemplos de train, test...)
  NSel = 0;
  SetNTrain(NTotal/2);

  //Iterative growing and pruning mathod variables
  GroupsCreated = false;
  GroupCount = 0;

  splitCriterium = 0;
}
void Data::redim(int ndatos)
{

  if (NTotal == ndatos) return;

  if (NTotal > ndatos) {//Reducimos
    for(int i = ndatos; i < NTotal; i++) {
      instances.pop_back();
    }
  }
  else {
    //double **mm;

    /*if (n_reservados > ndatos) {//Ampliando 
      mm = m;
    }
    else {//Ampliando y cambiando de puntero
      n_reservados = ndatos;//n_reservados * 2 > ndatos ? n_reservados * 2 : ndatos;
      mm = new double*[n_reservados];
      //Los datos extras los marcamos a cero para que casque si alguien se pasa
      for(int i = ndatos; i < n_reservados; i++) mm[i] = 0;
      for(int i = 0; i < NTotal; i++) mm[i] = m[i];
    }*/
    //instances.resize(NTotal);
    for(int i = NTotal; i < ndatos; i++) {
      // Extra elements for working space later
      instances.push_back(instanceBuilder->NewInstance(Nvar-1));
      // (0->Nvar-2) Independent variables ; (Nvar-1) Dependent variable
      for (int j = 0; j < Nvar-1; j++) instances[i](j, 0.0);
      instances[i].Class = -1;
      // (Nvar) and (Nvar+1) Multiple split working space
      instances[i].SetWorkingVar(0, drand());
      instances[i].SetWorkingVar(1, instances[i].GetWorkingVar(0));
      // (Nvar+2) Label for CV group
      instances[i].CvTreeLabel = (int) (1+(Ncv)*drand());
      // (Nvar+3) Membership
      instances[i].Membership = 1.0;
      // (Nvar+4) Weight (only for NomData with non-homogeneous distribution of classes)
      instances[i].Weight = 1.0;
      // (Nvar+5) Initial position
      instances[i].IniPos = i;
      // (Nvar+6) Group number for Iterative growing pruning method
      instances[i].Group = 1;
    }
   // if (m && m != instances) delete []m;
   //  m = instances;
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
//    instances[i]->CvTreeLabel = (int) (1+(Ncv)*drand()); // Label for CV group
}

Data* Data::DataFromFile(char *filename, Instance* builder)
{
  string nombre = filename;
  int pos = nombre.find_last_of('.');
  string ext = nombre.substr(pos);
  if (ext.compare(".cre")==0) {
    return LoadFromCre(filename, builder);
  }
  else if (ext.compare(".asc")==0) {
    return LoadFromAsc(filename, builder);
  }
  else {
    return LoadFromCre(filename, builder);
  }
  return 0;
}
Data* Data::LoadFromCre(char *filename, Instance* builder)
{
//  ifstream ifs(filename);
//  return new NomData(ifs, 0);
  return new NomData(filename, builder);
}
Data* Data::LoadFromAsc(char *filename, Instance* builder)
{
  string sfilename =filename;
  NomData *dat = new NomData(builder);

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
  int numberObjectives = MultiobjectiveInstance::GetNumMultiobjectives();
  bool multiObjective = numberObjectives > 0;
  
  fscanf(f, "%d", &Nvar);
//  init(NTotal?NTotal:1, Nvar);
//  printf("Nomtermsi, ... to %d\n", Nvar);
  NomTerms = new std::vector<std::string>[Nvar];
  VarNames = new string[Nvar];
  vector<string> varnames;
  OriginalVarType.resize(Nvar);
  
  Nord = Nfuzz = Nnom = iNom = 0;
  for(int i=0;i<Nvar;i++) {
    fscanf(f, "%s %d", varname, &tipo);
    if (tipo==0) {      // ORD
      varnames.insert(varnames.begin()+Nord+Nfuzz, varname);
      Nord++; 
      OriginalVarType[i] = ORD;
    }
    else if (tipo==1) { // FUZZ
      varnames.insert(varnames.begin()+Nord+Nfuzz, varname);
      Nfuzz++; 
      OriginalVarType[i] = FUZZ;
    }
    else if (tipo==2) { // NOM
      varnames.push_back(varname);
      OriginalVarType[i] = NOM;
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
    else if (tipo==3) { // IGN
      OriginalVarType[i] = IGN;
    }
    else if (tipo==4) { // GROUP
      OriginalVarType[i] = GROUP;
    }
    else if (tipo==5) { // WEIGHT
      OriginalVarType[i] = WEIGHT;
    }
  }

  //Variable dependiente
  DepVarType = OriginalVarType[Nvar-1];
  if (DepVarType==NOM) {
    Nnom -= multiObjective ? numberObjectives : 1;
  }
  else if (DepVarType==ORD) {
    Nord -= multiObjective ? numberObjectives : 1;
    fgets(linea, 1024, f);
  }

  if (Nvar != Nord+Nfuzz+Nnom+(multiObjective ? numberObjectives : 1)) {
    Nvar = Nord + Nfuzz + Nnom + (multiObjective ? numberObjectives : 1);
  }
  Nordfuzz = Nord+Nfuzz;
  NVarsOrdSplit = Nord;
  for (int i = 0; i < Nvar; i++) VarNames[i] = varnames[i];

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
  char *tok, *ok;
  int ilinea = 0;
  int iOrd, iNom, /*iIgn,*/ iCol;
  int iOrdIni=0, iNomIni=0, iIgnIni=0;

  //Initialising counters
  if (colIni==-1) colIni = 0;
  if (colFin==-1 || colFin>=Nvar) colFin = OriginalVarType.size()-1;//Nvar-1;
  for (int i = 0; i < colIni; i++) {
    if (OriginalVarType[i]==ORD) iOrdIni++;
    else if (OriginalVarType[i]==NOM) iNomIni++;
    else if (OriginalVarType[i]==IGN) iIgnIni++;
  }

  //loop for reading the data
  bool NotSaid = true;
  while(fgets(linea, 0x100000, f)) {
    int ll = strlen(linea);
    if (linea[ll-1] != '\n') {
      cout << "Warning, line " << (ilinea+1) << " is too long" << endl;
    }
    else {
      linea[ll-1] = ' ';
    }
    if (ll < Nvar) break;
    if (NTotal==ilinea) {//Reservamos más espacio
      redim( (int) (NTotal > 10 ? NTotal * 1.5 : 20) );
    }
//    if (NTotal>=70 && 0==ilinea%(NTotal/70)) printf("=");
    iOrd = iOrdIni;
    iNom = iNomIni;
/*    iIgn = iIgnIni;*/
    iCol = colIni;
    int natt_read = 0;
    for(tok = strtok(linea, " \r\t"); tok && iCol<=colFin; 
                                            tok = strtok(0, " \r\t"),iCol++) {
      if (OriginalVarType[iCol]==ORD) {
        double val;

        if (tok[0]=='?') val = DBL_MAX;
        else {
          val = strtod(tok, &ok);
          if (tok==ok) cout << "Error reading line " << (ilinea+1) << " var " << (iOrd+iNom+1) << endl;
        }

        if (iCol == Nvar-1) { //Dependent Var is ORD
          iOrd = Nvar-1;
        }
        instances[ilinea](iOrd, val);
        if (NotSaid && instances[ilinea][iOrd] != val) {
          cout << "*****************************************************" << endl;
          cout << "*****************************************************" << endl;
          cout << "Losing precision in instance " << (ilinea+1) << " var " << (iOrd+iNom+1);
          cout << ". From " << val << " to " << instances[ilinea][iOrd] << endl;
          cout << "*****************************************************" << endl;
          cout << "*****************************************************" << endl;
          NotSaid = false;
        }
        iOrd++;
      }
      else if (OriginalVarType[iCol]==NOM) {
        if (Nordfuzz + iNom == Nvar - 1) instances[ilinea].SetClass((int)DameRepresentacion(tok, iNom));
        else instances[ilinea](Nordfuzz + iNom, DameRepresentacion(tok, iNom));
        iNom++;
      }
      else if (OriginalVarType[iCol]==WEIGHT) {
        double val = strtod(tok, &ok);
        if (tok==ok) cout << "Error reading line " << (ilinea+1) << " var " << (iOrd+iNom+1) << endl;
        instances[ilinea].Weight = val;
      }
      else if (OriginalVarType[iCol]==GROUP) {
        int val = strtol(tok, &ok, 10);
        if (tok==ok) cout << "Error reading line " << (ilinea+1) << " var " << (iOrd+iNom+1) << endl;
        instances[ilinea].Group = val;
      }
      else if (OriginalVarType[iCol]==IGN) {
        //Do nothing
      }
      natt_read++;
    }
    if (natt_read != colFin + 1 ) {
      cout << "Warning, line " << (ilinea+1) << " read " << natt_read;
      cout << " attributes instead of " << (colFin+1) << endl;
    }
    else if (tok) {
      cout << "Warning, line " << (ilinea+1) << " seems to be longer than should be" << endl;
    }
    ilinea++;
  }

  //Ajustamos el tamaño
  if (NTotal != ilinea) redim(ilinea);
}
double Data::DameRepresentacion(string val, int ivarnom, bool force)
{
  vector<string>::iterator k;
  int ic;
  bool noesta = true;
  vector<string> &noms = NomTerms[ivarnom];
  double ret;
  
  //Busco 'val' en noms (teniendo en cuenta que la lista 
  //  de etiquetas esta ordenada por nombre)
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
      if ( imcol == Nvar-1 ) {
        for (int j = 0; j < NTotal && instances[j].Class >= 0; j++) {
          if ( instances[j].Class >= ic ) instances[j].Class++;
        }
      }
      else {
        for (int j = 0; j < NTotal && instances[j][imcol]>0.0; j++) {
          if ( instances[j][imcol] > ic ) instances[j](imcol, instances[j][imcol] + 1);
        }
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
    double hold = instances[i][iatt];
    instances[i](iatt, instances[elem][iatt]);
    instances[elem](iatt, hold);
  }
}
void Data::ChangeDependentColumn(int newDependentColumnIndex)
{
  for(int i = 0; i < NTotal; i++) {
    instances[i].ChangeDependentColumn(newDependentColumnIndex);
  }
}
void Data::CopyDependentColumn(int new_num_multi)
{
  // Multi-objetivo aumentado
  if (new_num_multi > 0) {
    if (MultilabelInstance::GetNumMultilabels() > 0) {
      MultilabelInstance::SetNumMultilabels(new_num_multi);
    }
    if (MultiobjectiveInstance::GetNumMultiobjectives() > 0) {
      MultiobjectiveInstance::SetNumMultiobjectives(new_num_multi);
    }
  }
  
  for(int i = 0; i < NTotal; i++) {
    instances[i].CopyDependentColumn();
  }

  // Nuevo nombre de clase
  string* newArray = new string[Nvar + 1];
  std::copy(VarNames, VarNames + Nvar, newArray);
  newArray[Nvar] = VarNames[Nvar-1] + "_b";
  delete []VarNames;
  VarNames = newArray;

  this->Nvar += 1;
}
void Data::Scramble(int begin,int end)
{

  for(int i = begin; i <= end ; ++i)
  {
    instances[i](Nvar, drand());
  }
  SortOn(Nvar,begin, end);
}


void Data::ScaleVariableIQ(int iVar, int scfirst,int first, int last)
{
  SortOn(iVar,scfirst,last);
  int N = last-scfirst+1;
  double val1 = instances[scfirst+N/4][iVar];
  double val2 = instances[scfirst+(3*N)/4][iVar];

  if(val2-val1 > 1e-6) {
    vscale[iVar] = 1.0/(val2-val1);
    vshift[iVar] = 0.5-vscale[iVar]*val2;
    for(int i = first; i <= last; i++) {
      instances[i](iVar, vscale[iVar] * instances[i][iVar] + vshift[iVar]);
    }
  } 
  else {
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
      int Class = instances[i].Class;
      if(Class != cur->iClass) nerror+=instances[i].Weight;
    }
  }
  else 
  {
    leftLast = GetLeftLast(cur, begin, end, K);
    nerror = Error2(cur->child, K, begin, leftLast) +
                                    Error2(cur->child->sib, K, leftLast+1, end);
  }

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
    instances[i](iVar, (instances[i][iVar] -vshift[ilabel])/vscale[ilabel]);

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
  if (nGroups == last-first+1) {
    for(j = first; j <= last; j++) {
      SetValueVar(j, icol, j + igrp);
      GroupCount[j] = 1;
    }
    return;
  }

  //Sort by class
  SortByClass(first, last);

//  srand(time(0) + rand());//randomize();
  if (DepVarType == NOM) {
    iclass=-1;
    for(i = first;i <= (unsigned)last; i++) {
      //A random number is thrown
      instances[i](Nvar, TDebugRand::Rand());//%10000;
      //The number of classes and number of elements of each class are obtained
      if (iclass!=GetDatClass(i)){ //(instances[i].Class)) {
        iclass=GetDatClass(i);     //instances[i].Class;
        nnclass.resize(nnclass.size()+1);
      }
      nnclass.at(nnclass.size()-1)++;
    }
  }
  else {
    nnclass.resize(1);
    nnclass[0] = last - first + 1;
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
  for(i = 0; i < nnclass.size(); i++) {
    SortOn(Nvar, ClassBegin, ClassBegin+nnclass.at(i)-1);
    nelem = nnclass.at(i)/nGroups;
    nelem = nelem*nGroups + ClassBegin;
    //A group number is assigned
    for(j = ClassBegin; j < nelem; j++) {
      SetValueVar(j, icol, igrp + j % nGroups);
      GroupCount[ (int)GetValueVar(j, icol) - igrp ]++;
    }
    //The elements from nGroups*nelem to nnclass.at(i) are added at random to
    // the groups
    for(; j < ClassBegin+nnclass.at(i); j++) {
      SetValueVar(j, icol, igrp + gr % nGroups);
      gr++;
      GroupCount[(int)GetValueVar(j, icol) - igrp]++;
    }
    ClassBegin += nnclass.at(i);
  }

  //Sort by group
  SortOn(icol, first, last);
}

static int theSortAtt;
static bool InvertSort=false;
static int CV;

//  -------------------------------------------------------------------
//    Comparing functions
//  -------------------------------------------------------------------
//  -------------------------------------------------------------------
//  -------------------------------------------------------------------
//
inline bool compare_variables(Instance* i1, Instance* i2)
{
  return (*i1)[theSortAtt] < (*i2)[theSortAtt];
}
inline bool compare_variables_inv(Instance* i1, Instance* i2)
{
  return (*i1)[theSortAtt] > (*i2)[theSortAtt];
}
inline bool compare_class(Instance* i1, Instance* i2)
{
  return i1->GetClass() < i2->GetClass();
}
inline bool compare_class_inv(Instance* i1, Instance* i2)
{
  return i1->GetClass() > i2->GetClass();
}
inline bool compare_ws1(Instance* i1, Instance* i2)
{
  return i1->GetWorkingVar(0) < i2->GetWorkingVar(0);
}
inline bool compare_ws1_inv(Instance* i1, Instance* i2)
{
  return i1->GetWorkingVar(0) > i2->GetWorkingVar(0);
}
inline bool compare_ws2(Instance* i1, Instance* i2)
{
  return i1->GetWorkingVar(1) < i2->GetWorkingVar(1);
}
inline bool compare_ws2_inv(Instance* i1, Instance* i2)
{
  return i1->GetWorkingVar(1) > i2->GetWorkingVar(1);
}
bool compare_cv(Instance* i1, Instance* i2)
{
  int A = (int) i1->GetCvTreeLabel(); 
  int B = (int) i2->GetCvTreeLabel();
  if(A == CV && B != CV) return false;
  if(A != CV && B == CV) return true;
  return false;
}
bool compare_cv_inv(Instance* i1, Instance* i2)
{
  int A = (int) i1->GetCvTreeLabel(); 
  int B = (int) i2->GetCvTreeLabel();
  if(A == CV && B != CV) return true;
  if(A != CV && B == CV) return false;
  return false;
}
inline bool compare_membership(Instance* i1, Instance* i2)
{
  return i1->GetMembership() < i2->GetMembership();
}
inline bool compare_membership_inv(Instance* i1, Instance* i2)
{
  return i1->GetMembership() > i2->GetMembership();
}
inline bool compare_weight(Instance* i1, Instance* i2)
{
  return i1->GetWeight() < i2->GetWeight();
}
inline bool compare_weight_inv(Instance* i1, Instance* i2)
{
  return i1->GetWeight() > i2->GetWeight();
}
inline bool compare_inipos(Instance* i1, Instance* i2)
{
  return i1->GetIniPos() < i2->GetIniPos();
}
inline bool compare_inipos_inv(Instance* i1, Instance* i2)
{
  return i1->GetIniPos() > i2->GetIniPos();
}
inline bool compare_group(Instance* i1, Instance* i2)
{
  return i1->GetGroup() < i2->GetGroup();
}
inline bool compare_group_inv(Instance* i1, Instance* i2)
{
  return i1->GetGroup() > i2->GetGroup();
}
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
  //instances.sort(int iVar, int first, int last, bool _InvertSort);

  if (iVar < Nvar-1) {
    stable_sort(instances.begin() + first, instances.begin() + last + 1, 
                  InvertSort ? compare_variables_inv : compare_variables);
  }
  else if (iVar == Nvar-1) {
    stable_sort(instances.begin() + first, instances.begin() + last + 1, 
                  InvertSort ? compare_class_inv : compare_class);
  }
  else if (iVar == Nvar) {
    stable_sort(instances.begin() + first, instances.begin() + last + 1, 
                  InvertSort ? compare_ws1_inv : compare_ws1);
  }
  else if (iVar == Nvar+1) {
    stable_sort(instances.begin() + first, instances.begin() + last + 1, 
                  InvertSort ? compare_ws2_inv : compare_ws2);
  }
  else if (iVar == Nvar+2) {
    stable_sort(instances.begin() + first, instances.begin() + last + 1, 
                  InvertSort ? compare_cv_inv : compare_cv);
  }
  else if (iVar == Nvar+3) {
    stable_sort(instances.begin() + first, instances.begin() + last + 1, 
                  InvertSort ? compare_membership_inv : compare_membership);
  }
  else if (iVar == Nvar+4) {
    stable_sort(instances.begin() + first, instances.begin() + last + 1, 
                  InvertSort ? compare_weight_inv : compare_weight);
  }
  else if (iVar == Nvar+5) {
    stable_sort(instances.begin() + first, instances.begin() + last + 1, 
                  InvertSort ? compare_inipos_inv : compare_inipos);
  }
  else if (iVar == Nvar+6) {
    stable_sort(instances.begin() + first, instances.begin() + last + 1, 
                  InvertSort ? compare_group_inv : compare_group);
  }
}

int Data::PreSort(int first, int last)
{
  int i;
  SortOn(Nvar+3,first,last);
  for(i=first; i<=last; ++i)
  {
    if (instances[i].Membership > MIN_BELONG_LOW)
      break;
  }
  return i;
}


void Data::SortCV(int cv, int begin, int end)
{
  CV = cv;
  theSortAtt = Nvar+2;
  if(end < 0) end = NGrow-1;
  stable_sort(instances.begin() + begin, instances.begin() + end + 1, compare_cv);
}

int Data::FindFirstTest(int cv, int first, int last)
{
  int ret = last+1;
  while(ret > first && fabs(instances[ret-1].CvTreeLabel - cv ) < 0.5)
    ret--;
  return ret;
}

int Data::SetCV(int cv)
{
  Ndata = NGrow;

  if(cv != 0)
  {
    SortCV(cv);
    while(Ndata && fabs(instances[Ndata-1].CvTreeLabel - cv) < 0.1) Ndata--;
    if(DepVarType == NOM) Ponderate(Ndata);
  }


  return Ndata;
}

void Data::MarkOrder()
{ //It wont work if someone touches the values in column c
  // before reseting the order
  for(int i = 0; i < GetNTotal(); i++) 
    instances[i].SetWorkingVar(1, i);
}
void Data::ResetOrder()
{
  int c = Data::MultipleSplitWs2 + Nvar;
  SortOn(c);
}
void Data::SetNTrain(int first, int last)
{
  int c = Data::MultipleSplitWs2 + Nvar;
  for(int i = 0; i < first; i++)            instances[i](Nvar+1, 1);
  for(int i = first; i <= last; i++)        instances[i](Nvar+1, 0);
  for(int i = last+1; i < GetNTotal(); i++) instances[i](Nvar+1, 2);
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
  for(i = 0; i < NVarsOrdSplit; i++) {
    if (attributes_to_use && !((*attributes_to_use)[i])) continue;
    iSplit=-1;
    chi = FindOrdSplit(scfirst, n->last, i, Nmin,n->d->MembershipTotal, &iSplit);
    if(iSplit >=0 && (*chiBest < 0.0 || chi < *chiBest)) {
      n->att = i;
      int iS2 = iSplit+1;
      while(instances[iS2].Membership < MIN_BELONG_LOW || instances[iSplit][i] == instances[iS2][i])
        iS2++;
      n->fSplit = 0.5 * ( instances[iSplit][i] + instances[iS2][i] );
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
        instances[i](Nvar, 0.0);
        for(int j=0;j<NVarsOrdSplit;j++) {
          instances[i](Nvar, instances[i][Nvar] + instances[i][j]*n->coeff[j]);
        }
      }
      iSplit =-1;

      FindOrdSplit(scfirst, n->last, Nvar, Nmin, n->d->MembershipTotal, &iSplit);
      if(iSplit >= 0) {
        int iS2 = iSplit+1;
        while(instances[iS2].Membership < MIN_BELONG_LOW)
          iS2++;
        n->fSplit = 0.5 * ( instances[iSplit][Nvar] + instances[iS2][Nvar]);
        tempmin = min_in_child;
        chi = FindMultiSplit(scfirst, n->last, Nmin,n->d->MembershipTotal, use,
                                             n->coeff, &(n->fSplit), &tempmin);
        EliminateVariables(scfirst, n->last, Nmin, n->d->MembershipTotal, 
                                   n->coeff,  use, chi,&(n->fSplit), &tempmin);
        for(i=n->first;i<=n->last;i++) instances[i](Nvar, instances[i][Nvar] - n->fSplit);
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
    else  { // DE if(n->coeff) 
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
      //instances[i](Nvar, instances[i][n->att] - (n->fSplit));
      instances[i].SetWorkingVar(0, instances[i][n->att] - (n->fSplit));
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
      //while((*leftLast) < n->last && instances[(*leftLast)+1][Nvar] < 0.0) 
      while((*leftLast) < n->last && instances[(*leftLast)+1].GetWorkingVar(0) < 0.0) 
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


double Data::FindFuzzSplit(int first, int last, int iX, int Nmin, double NMembers, FuzzySet** pfz)
{
     return 0.0;
}


double Data::FindNomSplit(int first, int last, int iX, int Nmin, double NMembers, int* pNomSplit)
{
  *pNomSplit = -1; // there is an error
  
  bool MinInParent = Nmin < 0 ? true : false;
  bool MinInLeaf   = !MinInParent; 

  Nmin = Nmin < 0 ? -Nmin : Nmin;

  if ( iX < Nordfuzz ) return 0.0; // variable not nominal
  if ( MinInParent && NMembers <   Nmin ) return 0.0; // Not enough data to split
  if ( MinInLeaf   && NMembers < 2*Nmin ) return 0.0; // Not enough data to have Nmin in each child 

  int dim = NomTerms[iX-Nordfuzz].size();

  double chi =-1.0;
  
  if(DepVarType == NOM && NomTerms[Nnom].size() == 2) {
    double** Pop = new double*[dim];
    
    double* Norm = new double[dim];
    int i;
    for(i = 0; i < dim; i++ ) {
      Pop [i] = new double[2];
      Pop[i][0] = i + 0.5;
      Pop[i][1] = 0.0;
      Norm[i] =0.0;
    }

    for (i = first; i <= last; i++ ) {    
      if(instances[i].GetClass() == 1)
        Pop[ (int)instances[i][iX] ][1] += instances[i].Membership;
      Norm[ (int)instances[i][iX] ] += instances[i].Membership;
    }
    


    for (i = 0; i< dim; ++i) {
      if(Norm[i] >= 1e-10) Pop[i][1] /= Norm[i];
      else                 Pop[i][1] = 1.0; 
    }

    theSortAtt = 1;
    qsort(Pop,dim,sizeof(double**),compAtt);
    
    int k = 0;
    for (i = 0; i< dim ; ++i) {
    
      k += (1 << (int)Pop[i][0]);

      int count = SortNomOn(k,iX, first,last);
      double min_in_node = splitCriterium->Initialise(first,last, first+count-1);
      if( (!MinInLeaf || min_in_node >= Nmin)  && 
          (chi < 0.0 || (splitCriterium->GetLastSplitImpurity() <= chi)) ) {
        chi = splitCriterium->GetLastSplitImpurity(); 
        *pNomSplit= k;
        min_in_child = min_in_node;
      }

    }

    for(i=0; i<dim; ++i)
      delete[] Pop[i];
  
    delete[] Norm;
    delete[] Pop;
  }
  else {

    for (int i=1; i < pow(2.0,dim-1); ++i) {
    //  membership of sets = i in binary form
      //R[0] = R[1] = 0.0;          
      int count = SortNomOn(i,iX, first,last);
      double min_in_node = splitCriterium->Initialise(first,last, first+count-1);
      if( (!MinInLeaf || min_in_node >= Nmin)  && 
          (chi<0.0 || (splitCriterium->GetLastSplitImpurity() <= chi)) ) {
        chi = splitCriterium->GetLastSplitImpurity();
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
    if((1<<(int)instances[i][iX]) & Membership)
    {
      instances[i](Nvar, 0);
      ++count;
    } else {
      instances[i](Nvar, 1);
    }
  }
  theSortAtt = Nvar;
  stable_sort(instances.begin() + first, instances.begin() + last + 1, compare_ws1);
  return count;
}



double Data::FindOrdSplit(int first, int last, int iX, int Nmin,double NMembers, int *iSplit, double gamma)
{
  *iSplit = -1; // Indicates an error

  bool MinInParent = Nmin < 0 ? true : false;
  bool MinInLeaf   = !MinInParent; 

  Nmin = Nmin < 0 ? -Nmin : Nmin;

  if ( MinInParent && NMembers <   Nmin ) return 0.0; // Not enough data to split
  if ( MinInLeaf   && NMembers < 2*Nmin ) return 0.0; // Not enough data to have Nmin in each child 
                                          
  int iVar = iX;
  if(gamma < 1) iVar = Nvar;

  SortOn(iVar,first,last);
    // We need a split between different values
  int curSplit = last;
  if(gamma > 1) {                                  // For single ordinal split
    curSplit--;//GMM
    while( curSplit >= (first) && 
           fabs(instances[curSplit][iVar] - instances[curSplit+1][iVar]) < 1e-10 &&  
           (!MinInLeaf || NMembers >= Nmin) ) {

      curSplit--;
      NMembers -= instances[curSplit+1].Membership;
    }
//    if(curSplit == (first-1) || NMembers < Nmin) return 0.0;  // Error: no good split
  }
  else {
    curSplit = last;
  }

  // Begin main algorithm
  
  double chi= 0.0;
  double min_in_node = splitCriterium->Initialise(first, last, curSplit,gamma, iX);  // ASG: ask James why iX and not iVar?
      
  while(curSplit >= first) {
    if( curSplit < last &&
        fabs(instances[curSplit][iVar] - instances[curSplit+1][iVar]) >= 1e-10 && //check for new best
        (!MinInLeaf || min_in_node >= Nmin) &&
        ((splitCriterium->GetLastSplitImpurity() < chi) || *iSplit <0) &&
        instances[curSplit].Membership>= MIN_BELONG_LOW ) {

      *iSplit = curSplit; 
      chi = splitCriterium->GetLastSplitImpurity();
      min_in_child = min_in_node;
    }
    
    // Move current point from left to right
    if(instances[curSplit].Membership > 0.0) { 
      min_in_node = splitCriterium->MovePoint(curSplit,gamma,iX,first,last);
    }

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
    instances[i](Nvar, 0.0);  
    for(int j=0;j<NVarsOrdSplit;j++) instances[i](Nvar, instances[i][Nvar] + instances[i][j]*a[j]);
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
        for(i = first; i <= last; i++) {
          instances[i](Nvar+1, ( instances[i][Nvar] - (*c) ) / ( instances[i][iX] + gamma ));
        }
        chi = FindOrdSplit(first, last, iX, Nmin, NMembers, &iSplit, gamma);

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
        
          

          deltaBest = 0.5 * ( instances[iSplit][Nvar+1] + instances[iS2][Nvar+1] );
          
          gammaBest = gamma;
        }
      } // end gamma-loop
      a[iX] -= deltaBest;
      *c += deltaBest*gammaBest;
      for(i = first; i <= last; i++) instances[i](Nvar, instances[i][Nvar] - instances[i][iX] * deltaBest);
    } // end iX loop

    int iSplit;
    chiPrev = chiLoop;
    // Normalise
    double sum = 0.0;

    for(i = 0; i < NVarsOrdSplit; i++) if(fabs(a[i]) > sum) sum = fabs(a[i]);
    for(i = 0; i < NVarsOrdSplit; i++) a[i] /= sum;
    for(i = first; i <= last; i++) instances[i](Nvar, instances[i][Nvar] / sum);
    *c /= sum;
    chiLoop = FindOrdSplit(first, last, Nvar, Nmin  ,NMembers, &iSplit); // ASG: if problems change to Nmin - 1 avoids problems related to not finding a split
    if(iSplit >=0 /*&& (*chiBest < 0.0 || chi < *chiBest)*/)
    {
      int iS2 = iSplit+1;
      while(instances[iS2].Membership < MIN_BELONG_LOW )
        iS2++;
      *c = 0.5 * ( instances[iSplit][Nvar] + instances[iS2][Nvar] );
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
    instances[i](Nvar, 0.0);
    for(int j=0;j<NVarsOrdSplit;j++)
      instances[i](Nvar, instances[i][Nvar] + a[j]*instances[i][j]);
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
    for(int iX = 0; iX <NVarsOrdSplit; iX++)
    {
      if(use[iX] == 0) continue;
      arem = a[iX]; a[iX] = 0.0;

      for(int i = first; i <= last; i++) {
        instances[i](Nvar, instances[i][Nvar] - arem * instances[i][iX]);
      }
      
      chi = FindOrdSplit(first, last, Nvar, Nmin, NMembers, &iSplit);

    //  if(iSplit >= 0) // To Review ASG
      {
        double deltaTemp = I-chi;
        if(flag == 0 || deltaTemp < deltaMin) {deltaMin = deltaTemp;}
        if(flag == 0 || deltaTemp > deltaMax) {deltaMax = deltaTemp;imax = iX;}
      }
      for(int i = first; i <= last; i++) {
        instances[i].SetWorkingVar(0, instances[i].GetWorkingVar(0) + arem * instances[i][iX]);
      }
      a[iX] = arem; arem = 0.0; flag++;
    }
    if( flag > 1 ) {
      if( delta - deltaMax < cartBeta * (delta - deltaMin) ) {
        stop = 0;
        for(int i = first;i <= last; i++) {
          instances[i].SetWorkingVar(0, instances[i].GetWorkingVar(0) - a[imax] * instances[i][imax]);
        }
        a[imax] = 0.0;
        use[imax] = 0;
        chimax = I - deltaMax;
      } 
      else {
        stop = 1;
      }
    }
  } while(!stop);
  chimax = FindOrdSplit(first, last, Nvar, Nmin, NMembers, &iSplit);
  if(iSplit < 0) 
  {
    for(i = 0; i<NVarsOrdSplit; ++i)
    {
      a[i] = old_a[i];
    }
    delete []old_a;
    return -1;   // ASG: Recover old coeffs
  }

  delete []old_a;
  int iS2 = iSplit+1;
  while(instances[iS2].Membership < MIN_BELONG_LOW)
    iS2++;
  *fSplit = 0.5 * ( instances[iSplit][Nvar] + instances[iS2][Nvar] );
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
    instances[i].Membership = 1.0*instances[i].Weight;



  while(cur->parent)
  {
    for(int i= n->first; i<= n->last; ++i)
    {
      switch(cur->parent->SplitType)
      {
      case NOM:

        if((1<<(int)instances[i][cur->parent->att]) & cur->parent->NomSplit)
        {
          if(!cur->sib)
          {
            instances[i].Membership = 0.0;
          }
        }
        else if(cur->sib)
        {
          instances[i].Membership = 0.0;
        }
        break;
      case ORD:
        instances[i].SetWorkingVar(0, -cur->parent->fSplit);

        if (cur->parent->coeff && !cur->parent->IsUnivariateSplit()) // !cur->parent->univariate_split 
        {
          for(int j=0;j<NVarsOrdSplit;j++) {
            instances[i].SetWorkingVar(0, instances[i].GetWorkingVar(0) + 
                                     cur->parent->coeff[j] * instances[i][j]);
          }
//if ( (instances[i][cur->parent->att] - cur->parent->fSplit) * instances[i][Nvar]<=0.0) 
//printf("MS<-------------%g--------%g-------------------------\n", 
//instances[i][Nvar], instances[i][cur->parent->att] - cur->parent->fSplit);
        }
        else {
          instances[i].SetWorkingVar(0, instances[i].GetWorkingVar(0) + 
                                              instances[i][cur->parent->att]);
/*if (instances[i][Nvar]==0.0) 
printf("MS<-------------%g--------%g-------------------------\n", instances[i][Nvar], m[i][cur->parent->att]-cur->parent->fSplit);*/
        }


        if(instances[i].GetWorkingVar(0) <=0.0) // ASG to Verify
        {
          if(!cur->sib)  // R Node
          {
            instances[i].Membership = 0.0;
          }
        }
        else if(cur->sib)  // L Node
        {
          instances[i].Membership = 0.0;
        }
        break;

      case FUZZ:

  /*      int iX = cur->parent->att;
        if(iX==Nvar)
        {
          instances[i][Nvar] = -cur->parent->fSplit;
          for(int j=0;j<NVarsOrdSplit;j++)
            instances[i][Nvar] += cur->parent->coeff[j]*instances[i][j];
        }

        if(cur->sib)
        {
          instances[i].Membership *= cur->parent->d->FuzzSplit->Membership((m[i][iX]- Gvshift[iX])/Gvscale[iX]);
        }
        else
        {
          instances[i].Membership *= (1.0 - cur->parent->d->FuzzSplit->Membership((m[i][iX]- Gvshift[iX])/Gvscale[iX]));
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
    total += instances[i].Membership;

  return total;
}
bool Data::GoLeft(Node* n, int idat)
{
  bool goleft = true;

  if (n->SplitType == NOM) {
    if( ! ((1 << (int)instances[idat][n->att] ) & n->NomSplit ) ) {
      goleft = false;
    }
  }
  else if (n->SplitType == ORD) {
    instances[idat].SetWorkingVar(0, -n->fSplit);

    if (n->coeff && !n->IsUnivariateSplit()) {
      for(int j = 0; j < NVarsOrdSplit; j++) {
        instances[idat].SetWorkingVar(0, instances[idat].GetWorkingVar(0) + 
                                         n->coeff[j] * instances[idat][j]);
      }
    }
    else {
      instances[idat].SetWorkingVar(0, instances[idat].GetWorkingVar(0) + 
                                                  instances[idat][n->att]);
    }

    if(instances[idat].GetWorkingVar(0) > 0.0) {
      goleft = false;
    }
  }

  return goleft;
}

void Data::FindMembershipTest(const Node* n)
{
  const Node* cur = n;

  for(int i= n->d->firstTest; i<= n->d->lastTest; ++i)
    instances[i].Membership = 1.0*instances[i].Weight;

  
  while(cur->parent)
  {  
    for(int i= n->d->firstTest; i<= n->d->lastTest; ++i)
    {    
      switch(cur->parent->SplitType)
      {
      case NOM:
        
        if((1<<(int)instances[i][cur->parent->att]) & cur->parent->NomSplit)
        {  
          if(!cur->sib)
            instances[i].Membership = 0.0;
        }
        else if(cur->sib)
          instances[i].Membership = 0.0;
        break;

      case ORD:
        instances[i].SetWorkingVar(0, -cur->parent->fSplit);
        if(cur->parent->coeff) {
          for(int j=0;j<NVarsOrdSplit;j++) 
            instances[i].SetWorkingVar(0, instances[i].GetWorkingVar(0) + cur->parent->coeff[j] * instances[i][j]);
        }
        else {
          instances[i].SetWorkingVar(0, instances[i].GetWorkingVar(0) + instances[i][cur->parent->att]);
        }
      
        
                
        if(instances[i][Nvar] <=0.0) // ASG to Verify 
        {
          if(!cur->sib)  // R Node
            instances[i].Membership = 0.0;
        }
        else if(cur->sib)  // L Node
          instances[i].Membership = 0.0;
        break;
        
      case FUZZ:
/*        int iX = cur->parent->att;

        if(iX == Nvar)
        {
          instances[i][Nvar] = -cur->parent->fSplit;
          for(int j=0;j<NVarsOrdSplit;j++)
            instances[i][Nvar] += cur->parent->coeff[j]*instances[i][j];
        }

        if(cur->sib)
        {
          instances[i].Membership *= cur->parent->d->FuzzSplit->Membership((m[i][iX]- Gvshift[iX])/Gvscale[iX]);
        }
        else
        {
          instances[i].Membership *= (1.0 - cur->parent->d->FuzzSplit->Membership((m[i][iX]- Gvshift[iX])/Gvscale[iX]));
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
    for(int i = n->first; i <= n->last; ++i) {
      if(instances[i].Membership > MIN_BELONG_HIGH) {
        double temp;
        if(n->att == Nvar) {
          instances[i].SetWorkingVar(0, -n->fSplit);
          for(int j=0;j<NVarsOrdSplit;j++) {
            instances[i].SetWorkingVar(0, instances[i].GetWorkingVar(0) + 
                                                n->coeff[j]*instances[i][j]);
          }
          temp = instances[i][Nvar];
        }
        else {
          temp = (instances[i][n->att] - Gvshift[n->att]) / Gvscale[n->att];
        }

        if(i == n->first || temp > d2) {
          d2 = temp;
        }
        if(i==n->first || temp < d1) {
          d1 = temp;
        }
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
          val += coeff[j] * instances[i][j];
    else val = instances[i][att];
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
        double aux = cur->d->tempMU*instances[i].Weight;
        if(cur->parent) aux *= cur->parent->d->cumMU;
        
        total += aux;
        fuzz += aux*(1.0-aux);

      } else {

        double mu =  cur->ParamMem(instances[i],NVarsOrdSplit);
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
    instances[i].Weight = 1.0;
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
      if (instances[i].Membership==0.0) {break;}
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
    for(i=begin;i<=end;i++) instances[i][Nvar] = m[i][n->att] - (n->fSplit);
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
    instances[i].SetWorkingVar(0, TDebugRand::Rand());
    instances[i].CvTreeLabel = (int) (1+(Ncv)*drand()); // Label for CV group
    //The number of classes and number of elements of each class are obtained
    if (iclass!=(instances[i].Class)) {
      iclass=instances[i].Class;
      nnclass.resize(nnclass.size()+1);
    }
    nnclass.at(nnclass.size()-1)++;
  }

  //Each class is ordered by the random number just thrown
  int ClassBegin=0;
  for(i=0;i<nnclass.size();i++) {
    SortOn(Nvar, ClassBegin, ClassBegin+nnclass.at(i)-1);
    //The elements in the class are numbered
    for(j=0;j<nnclass[i];j++) {
      instances[ j+ClassBegin ].SetWorkingVar(0, j);
    }
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
    instances[i].SetWorkingVar(0, iclass + 1 + TDebugRand::Rand());
    instances[i].CvTreeLabel = (int) (1+(Ncv)*drand()); // Label for CV group
    //The number of classes and number of elements of each class are obtained
    if (iclass!=(instances[i].Class)) {
      iclass=instances[i].Class;
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
      instances[ j+ClassBegin ].SetWorkingVar(0, -1);
    ClassBegin += nnclass.at(i);
  }

  //Sort by group
  SortOn(Nvar, 0, NTotal-1);
}
void Data::SaveToFile(char *FileName, int ini, int fin)
{
  int lon = strlen(FileName);

  if (lon >= 6 && strcmpci(FileName+lon-5, ".arff") == 0 ) {
    SaveAsArff(FileName, ini, fin);
  }
  else if (lon >= 7 && strcmpci(FileName+lon-6, ".names") == 0 ) {
    SaveAsC45(FileName);
  }
  else {
    ofstream ofs(FileName);
    SaveToStream(ofs, ini, fin);
    ofs.close();
  }
}
string Data::GetExampleAsString(int iEx, string sep)
{
  int iNom, iOrd;
  std::ostringstream buf;
 
  if (sep.size()==0) sep="\t";

  buf.precision(14);
  iNom = iOrd = 0;
  for(int i = 0; i < (int)OriginalVarType.size() - 1; i++) {
    if (OriginalVarType[i]==ORD) {
      if (instances[iEx][iOrd]==DBL_MAX) buf << "?" << sep;
      else                               buf << instances[iEx][iOrd] << sep;
      iOrd++;
    }
    else if (OriginalVarType[i]==NOM) {
      if (NomTerms[iNom].size() == 0) {
        buf << (int)instances[iEx][Nordfuzz + iNom] << sep;
      }
      else {
        buf << NomTerms[iNom][(int)instances[iEx][Nordfuzz + iNom]] << sep;
      }
      iNom++;
    }
    else if (OriginalVarType[i]==GROUP) {
      buf << GetDatGroup(iEx) << sep;
    }
    else if (OriginalVarType[i]==WEIGHT) {
      buf << instances[iEx].Weight << sep;
    }
  }


  if ( DepVarType == NOM ) {
    if (NomTerms[iNom].size() == 0) {
      buf << instances[iEx].Class;
    }
    else {
      buf << NomTerms[iNom][instances[iEx].Class];
    }
  }
  else {
    buf << instances[iEx][Nvar-1];
  }

  return buf.str();
}
string Data::GetNamesAsString()
{
  std::ostringstream buf;
  int iNom, iOrd;
  
  buf << OriginalVarType.size();//Nvar;
  iNom = iOrd = 0;
  for (int i = 0; i < (int)OriginalVarType.size(); i++) {
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
    else if (OriginalVarType[i]==IGN) {
      buf << "IGNORE\t3";
    }
    else if (OriginalVarType[i]==GROUP) {
      buf << "GROUP\t4";
    }
    else if (OriginalVarType[i]==WEIGHT) {
      buf << "WEIGHT\t5";
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
      out << instances[i][j] << "\t";
    //Nominal
    for(; j < Nvar-1; ++j)
      out << NomTerms[j-Nordfuzz][(int)instances[i][j]] << "\t";

    //Depvar
    if (DepVarType==NOM)
      out << NomTerms[Nnom][(int)instances[i][j]] << "\n";
    else //ordinal
      out << instances[i][j] << "\n";
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
      if (instances[i][j]==DBL_MAX) fprintf(f, "?,");
      else fprintf(f, "%f," , instances[i][j]);

    //Nominal
    for(;j<Nvar-1;j++)
      fprintf(f, "%s,", NomTerms[j-Nordfuzz][(int)instances[i][j]].c_str());

    //Depvar
    if (DepVarType==NOM)
      fprintf(f, "%s\n", NomTerms[Nnom][instances[i].Class].c_str());
    else //ordinal
      fprintf(f, "%f\n", instances[i][j]);
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
void Data::SaveAsArff(char *fichero, int ini, int fin)
{
  unsigned iNom, iOrd;
  char fn[512], n[256];

  if (ini < 0) ini = 0;
  if (fin < 0) fin = GetNTotal() - 1;

  //Se anade la extension .arff si no la trae ya puesta
  int lon = strlen(fichero);
  strcpy(fn, fichero);
  if (lon<6 || strcmpci(fichero+lon-5, ".arff")!=0)
    strcat(fn, ".arff");
  
  strcpy(n, fn);
  lon = strlen(fichero);
  n[lon-5] = '\0';

  //
  //Se abre el fichero 
  //
  FILE *f = fopen(fn, "wt");
  if (f==NULL) return;

  
  fprintf(f, "@RELATION %s\n", n);
  fprintf(f, "\n");

  iNom = iOrd = 0; 
  for (int i = 0; i < (int)OriginalVarType.size(); i++) {
    if (OriginalVarType[i]==ORD) {
      fprintf(f, "@ATTRIBUTE\t%s\tREAL", VarNames[iOrd++].c_str());
    }
    else if (OriginalVarType[i]==NOM) {
      fprintf(f, "@ATTRIBUTE\t%s\t{%s", VarNames[Nordfuzz + iNom].c_str(), NomTerms[iNom][0].c_str());
      for (unsigned iv = 1; iv < NomTerms[iNom].size(); iv++) {
        fprintf(f, ",%s", NomTerms[iNom][iv].c_str());
      }
      fprintf(f, "}");
      iNom++;
    }
    fprintf(f, "\n");
  }

  fprintf(f, "@DATA\n");
  for(int i = ini; i <= fin; i++) {
    fprintf(f, "%s\n", GetExampleAsString(i, ",").c_str());
  }

  fclose(f);

}
//------------------------------------------------------  DeleteColumn  -----
//---------------------------------------------------------------------------
void Data::DeleteColumn(int icol)
{
  for(int i = 0; i < NTotal; i++) {
    instances[i].DeleteColumn(icol);
  }

  string *vn = new string[Nvar - 1];
  for(int i = 0; i < icol; i++) vn[i] = VarNames[i];
  for(int i = icol + 1; i < Nvar; i++) vn[i - 1] = VarNames[i];
  delete []VarNames;
  VarNames = vn;

  int ovt_icol = -1;
  for(int i = 0; i < (int)OriginalVarType.size(); i++) {
    if (OriginalVarType[i]==ORD || OriginalVarType[i]==NOM || 
                                        OriginalVarType[i]==FUZZ) {
      ovt_icol++;
      if (ovt_icol == icol) {
        OriginalVarType.erase(OriginalVarType.begin()+i);
        break;
      }
    }
  }

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
  for(int i = ini; i <= fin; i++) {
    int inipos =  instances[lastNTotal + i - ini].IniPos;
    instances[lastNTotal + i - ini] = d2->instances[ i ];
    instances[lastNTotal + i - ini].IniPos = inipos;
//    for(int j = 0; j < Nvar + 5; j++) { // dont copy the original order as it makes no sense
//      m[ lastNTotal + i - ini ][ j ] = d2->m[ i ][ j ];
//    }
//    m[ lastNTotal + i - ini ][ Nvar+6 ] = d2->m[ i ][ Nvar+6 ];
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
  fprintf(f, "%s", GetNamesAsString().c_str());
  rewind(f);
  clon->GetNames(f);
  fclose(f);

  //Copiamos los datos
  clon->init();
  int num_dif=fin-ini+1;
  clon->redim( num_dif * factor );

  for(int i = ini; i <= fin; i++) {
    for(int k = 0; k < factor; k++) {
      clon->instances[ k * num_dif + i - ini ] = instances[i];
    }
  }
  
  clon->SetSplitCriterium(this->splitCriterium);

  return clon;
}
//------------------------------------------------  GetBootstrapSample  -----
//---------------------------------------------------------------------------
Data* Data::GetBootstrapSample(int N, bool Weighted, bool IncludeOobAsTest, std::vector<int> *idxs)
{
  Data *bootstrap = Clone(0,0);//NewData();

  //Por omision se extraen igual numero de ejemplos que datos de entrenamiento
  if (N<=0) N=NTrain;

  // Cargamos los nombres y redimensionamos los datos
//  FILE *f=tmpfile();
 // fprintf(f, GetNamesAsString().c_str());
//  rewind(f);
 // bootstrap->GetNames(f);
  //fclose(f);
  //bootstrap->init();

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
    bootstrap->instances[i] = instances[elem];
  }

  if (Weighted)
    ResetWeights();

  int nNoUsados=0;
  if (IncludeOobAsTest) {
    //Los ejemplos no utilizados en el bootstrap los pongo como test del dataset
    if ( ! idxs ) {
      for(int i=0; i<NTrain;i++) {
        if (incluido[i]) continue;
        bootstrap->instances[N+nNoUsados] =  instances[i];
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
      instances[i].SetWorkingVar(0, 0);
      nClass++;
    }
    else {
      instances[i].SetWorkingVar(0, 1);
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

  int n = (int) ((N + 0.5) / NumClass);

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
  MINomData *clon = new MINomData(data->instanceBuilder);

  // Cargamos los nombres
  FILE *f=tmpfile();
  fprintf(f, "%s", data->GetNamesAsString().c_str());
  rewind(f);
  clon->GetNames(f);
  fclose(f);

  //Copiamos los datos
  clon->init();
  int n_total = data->GetNTotal();
  clon->redim(n_total);
  for(int i = 0; i < n_total; i++) {
    clon->instances[ i ] = data->instances[ i ];
    clon->instances[ i ].Group = (int) data->instances[ i ][ i_mi_column ];
  }

  // Marcamos la columna de mi
  clon->MarkOrder();
  clon->SortByGroup();
  int i_mi = -1;
  double last_value = clon->instances[ 0 ].Group - 1.0;
  for(int i = 0; i < n_total; i++) {
    if (last_value !=  clon->instances[ i ].Group ) i_mi++;
    last_value = clon->instances[ i ].Group;
    clon->instances[ i ].Group = i_mi;
  }
  clon->ResetOrder();

  clon->DeleteColumn( i_mi_column );
  clon->SetNTrain( data->GetNTrain() );
  clon->ResetWeights();

  return clon;
}

//-------------------------------------------------  SetSplitCriterium  -----
//---------------------------------------------------------------------------
void Data::SetSplitCriterium(SplitCriterium *sc)
{
  SplitCriterium *newSplit = 0;
  GiniCriterium *newSplit1  = dynamic_cast<GiniCriterium*>(sc);
  MSECriterium *newSplit2   = dynamic_cast<MSECriterium*>(sc);
  ComboCriterium *newSplit3 = dynamic_cast<ComboCriterium*>(sc);
  MultilabelGiniCriterium *newSplit4 = dynamic_cast<MultilabelGiniCriterium*>(sc);
  BoostedGiniCriterium *newSplit5  = dynamic_cast<BoostedGiniCriterium*>(sc);
  MultipleMSECriterium *newSplit6   = dynamic_cast<MultipleMSECriterium*>(sc);

  if (newSplit4) {
    newSplit = new MultilabelGiniCriterium(this);
  }
  else if (newSplit1) {
    newSplit = new GiniCriterium(this);
  }
  else if (newSplit2) {
    newSplit = new MSECriterium(this);
  }
  else if (newSplit3) {
    newSplit = new ComboCriterium(this);
  }
  else if (newSplit5) {
    newSplit = new BoostedGiniCriterium(this);
  }
  else if (newSplit6) {
    newSplit = new MultipleMSECriterium(this);
  }

  if (newSplit == 0) {
    cout << "Warning: split type not found" << endl;
  }
  else {
    if (splitCriterium) delete splitCriterium;
    splitCriterium = newSplit;
  }

}
void Data::SetSplitCriterium(string splitClassName)
{
  SplitCriterium *newSplit = 0;

  if ( splitClassName.compare("GiniCriterium")==0 ) {
    newSplit = new GiniCriterium(this);
  }
  else if ( splitClassName.compare("MSECriterium")==0 ) {
    newSplit = new MSECriterium(this);
  }
  else if ( splitClassName.compare("MultipleMSECriterium")==0 ) {
    newSplit = new MultipleMSECriterium(this);
  }
  else if ( splitClassName.compare("ComboCriterium")==0 ) {
    newSplit = new ComboCriterium(this);
  }
  else if ( splitClassName.compare("MultilabelGiniCriterium")==0 ) {
    newSplit = new MultilabelGiniCriterium(this);
  }
  else if ( splitClassName.compare("BoostedGiniCriterium")==0 ) {
    newSplit = new BoostedGiniCriterium(this);
  }

  if (newSplit == 0) {
    cout << "Warning: split " << splitClassName << " not found" << endl;
    exit(1);
  }
  else {
    if (splitCriterium) delete splitCriterium;
    splitCriterium = newSplit;
  }
}
//---------------------------------------------------------------------------
double Data::DistanciaEjem(int d1, int d2)
{
  double dist = 0.0;

  for (int i = 0; i < Nvar-1; i++)
    dist += pow( instances[d1][i] - instances[d2][i], 2.0);

  return sqrt(dist);
}
//---------------------------------------------------------------------------
double Data::Min(int icol)
{
  double min = numeric_limits<double>::infinity(); 

  for(int i=0; i<NTotal;i++)
    if (instances[i][icol]<min) 
      min = instances[i][icol];

  return min;
}
double Data::Max(int icol)
{
  double max = -numeric_limits<double>::infinity(); 

  for(int i=0; i<NTotal;i++)
    if (instances[i][icol]>max) 
      max = instances[i][icol];

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
  for(unsigned k = 0; k < ini_c_0.size(); k++) {
    for(int i = 0;i < tam_c_1[k]; i++) {
      for(int j = 0; j < Nvar-1; j++) {
        int pos = ini_c_0[k] + 
              (int) ((double)tam_c_0[k]*TDebugRand::Rand()/(RAND_MAX+1.0));
        d->instances[tam_acum+i](j, instances[pos][j]);
      }
      int pos = ini_c_0[k] + 
              (int) ((double)tam_c_0[k]*TDebugRand::Rand()/(RAND_MAX+1.0));
      d->instances[tam_acum+i].Class = instances[pos].Class;
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
        for(int iAtt = 0; iAtt < Nvar-1; iAtt++) {
          double dif = instances[pos][iAtt] - instances[i][iAtt];
          d->instances[iSal](iAtt, instances[i][iAtt] + gap*dif);
        }
        int iclass = gap>drand() ? pos : i;     //Clase
        d->instances[iSal].Class = instances[iclass].Class; //Clase
        if ( instances[pos].Class !=instances[i].Class ) iSal++; //Generate only frontier examples
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
    d->instances[i].SetWorkingVar(0, 10);
    for(int j=0;j<knn;j++) {
      int pos = (int)distancias[j][1];
      if ( instances[pos].Class != instances[i].Class ) {//Distinta clase -> se aniade
        d->instances[i].SetWorkingVar(0, 1);
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
MINomData::MINomData(char *filename, int mi_attrib, Instance *builder) : 
                             NomData(filename, builder)
{
  NMITrain = 0;
  NMITotal = 0;
  int n_total = GetNTotal();

  if (OriginalVarType[mi_attrib]!=GROUP) {
    for(int i = 0; i < n_total; i++) {
      instances[ i ].Group = (int) instances[ i ][ mi_attrib ];
    }
  }

  // Marcamos la columna de mi
  MarkOrder();
  SortByGroup();
  int i_mi = -1;
  double last_value = instances[ 0 ].Group - 1.0;
  for(int i = 0; i < n_total; i++) {
    if ( last_value !=  instances[ i ].Group ) i_mi++;
    last_value = instances[ i ].Group;
    instances[ i ].Group = i_mi;
  }
  ResetOrder();

  if (OriginalVarType[mi_attrib]!=GROUP) {
    cout << "Deleting column " << mi_attrib << endl;
    DeleteColumn( mi_attrib );
  }
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
      ((MINomData*)bootstrap)->instances[imi] = instances[k];
      bootstrap->SetDatGroup(imi, igrp); //Renumero las instancias
      //bootstrap->SetDatWeight(imi, (double)NTrain/NMITrain/(posiciones[elem+1]-posiciones[elem])); //Las doy un peso inv proporcional al no de instancias en el ejemplo
      imi++;
      if ( imi == size ) {
       // if (i < N/2 ) {
          size = (int) (0.5 + (double)(N+5.0)*size/(i+1));
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
          ((MINomData*)bootstrap)->instances[imi] = instances[k];
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
  Weights.resize(NMITotal, 1.0);
  cout << endl; 
}
//---------------------------------------------------------------------------
void MINomData::SetNTrain(int n)
{
  NomData::SetNTrain(n); 
  if (n == GetNTotal()) {
    int i_mi = -1;
    double last_value = instances[ 0 ].Group - 1.0;
    for(int i = 0; i < n; i++) {
      if ( last_value !=  instances[ i ].Group ) i_mi++;
      last_value = instances[ i ].Group;
      instances[ i ].Group = i_mi;
    }
  }
  NMITrain = CountMIInstances(0, n - 1);
  NMITotal = CountMIInstances(0, GetNTotal() - 1 );// No debe estar aqui...
//cout << NMITrain << " " << NMITotal << " " << n << " " << GetNTotal() << endl;
  if ( Weights.size() != (unsigned)NMITotal ) {
    Weights.clear();
    Weights.resize(NMITotal, 1.0);
    //ResetWeights();
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

  SetNTrain(GetNTotal());

  ResetWeights();
  
}

//----------------------------------  MINomData::SubstituteClassLabels  -----
//---------------------------------------------------------------------------
void MINomData::SubstituteClassLabels(Data *data, int attribute_idx)
{
  MINomData *mind = dynamic_cast<MINomData*>(data);

  if (mind) { //Se marca ejemplo a ejemplo
    NomData::SubstituteClassLabels(data, attribute_idx); 
    return;
  }
  else if (!mind && data->GetNTotal() != NMITotal ) {
    cout << "Different Number of mi / instances: Nothing done." << endl;
    return;
  }

  if (attribute_idx < 0) {
    attribute_idx = data->GetNumVar()-1;
  }
  bool IsNom = attribute_idx < data->GetNordfuzz() ? false : true;

  cout << "Old classes: ";
  for(unsigned i = 0; i < NomTerms[GetNumVarNom()].size(); i++) {
    cout << NomTerms[GetNumVarNom()][i] << ", ";
  }
  cout << endl; 
  NomTerms[GetNumVarNom()].clear();
  if (IsNom) {
    NomTerms[GetNumVarNom()] = data->GetTermVector(attribute_idx - data->GetNordfuzz());
  }

  SortByGroup();
  data->OriginalOrder();
  int last_value =  GetDatGroup(0);
  vector<int> new_hist(NomTerms[GetNumVarNom()].size(), 0);
  for(int i = 0, imi = 0; i < data->GetNTotal(); i++) {
    while ( imi < GetNTotal() && last_value ==  GetDatGroup(imi) ) {
      if (IsNom) {
        instances[imi].SetClass((int) data->GetValueVar(i, attribute_idx));
      }
      else {
        ostringstream os;
        os << data->GetValueVar(i, attribute_idx);
        instances[imi].SetClass((int)DameRepresentacion(os.str(), GetNumVarNom()));
      }
      imi++;
    }
    new_hist[instances[imi-1].GetClass()]++;
    if (imi < GetNTotal()) last_value =  GetDatGroup(imi);
  }

  NumClass = NomTerms[GetNumVarNom()].size();
  delete[] Pop;
//  delete[] PopR;
//  delete[] PopL;
  delete splitCriterium;
  init();
  cout << "New classes: ";
  for(unsigned i = 0; i < NomTerms[GetNumVarNom()].size(); i++) {
    cout << NomTerms[GetNumVarNom()][i] << "(" << new_hist[i] << "), ";
  }
  cout << endl; 
}
//-------------------------------------------------------------  Clone  -----
//---------------------------------------------------------------------------
Data* MINomData::Clone(int ini, int fin, int factor)
{
  MINomData *mind =(MINomData*) Data::Clone(ini, fin, factor);

//  mind->SetNTrain(GetNTrain());
  //mind->ResetWeights();

  return mind;
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