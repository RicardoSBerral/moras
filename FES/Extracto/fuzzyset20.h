//---------------------------------------------------------------------------

#ifndef Fuzzyset20H
#define Fuzzyset20H
//---------------------------------------------------------------------------
#include <stdio.h>
#include "fuzzydef20.h"
#include <string> 
#include "definitions20.h"


class FuzzySet
{
public:
	FuzzySet(const char* name, int type);
	FuzzySet(FuzzySet*);
        virtual ~FuzzySet(){}

	static void Defaults(int id,double& p1,double& p2,double& p3,double Domain[2]);

	double ValueFromTruth(double TruthVal);
	double Membership(double Val);
	void Normalise();
	double GetHeight();

	bool IsEmpty()
		{return _IsEmpty;}
	bool IsOK()
		{return (_Status == OK);}
	bool IsNormal();

	int SetDomain(double dom_min, double dom_max);
	void GetDomain(double& d1, double& d2){d1 = _Domain[0];d2=_Domain[1];}
	int SetAlpha(double a)
		{if(a < 0 || a > 1) return 1; _Alpha = a; return 0;}
	double GetAlpha(){return _Alpha;}
	const std::string GetDesc(){return _Description;}
	FuzzySet* Next()	{return _Next;}
	void SetNext(FuzzySet* n){_Next = n;}
	void Draw();
	virtual void PrintMe();
	double* GetValues(){return _MemVec;}
	double GetParam(int i){return _Params[i];}
  std::string& GetName(){return _Name;}
//2_  char* GetName(){return _Name;}
	int GetType(){return _CurveType;}
	void SetType(int t){_CurveType = t;}
	void SetName(const char* name){_Name = name;}
//2_	void SetName(const char* name){strcpy(_Name, name);} //2
	void Generate(double p1, double p2, double p3=0);
	int CheckParams(double p1, double p2, double p3);
	void Write(FILE* os);
	void Read(FILE* is);
protected:
	std::string _Name,_Description;
//	char _Name[256],_Description[256];
	int _CurveType;
	bool _IsEmpty;
	int _Order;
	double _Domain[2];
	double _Params[4];
	double _Alpha;
	double _MemVec[VECMAX];
	FuzzySet* _Next;
	int _Status;
};
/*
class LinearSet : public FuzzySet
{
public:
	LinearSet(const char* name);
	void Generate(double low, double high, int type);
	void PrintMe();
};

class S_Set : public FuzzySet
{
public:
	S_Set(const char* name);
	void Generate(double left, double flex, double right, int type);
	void PrintMe();
};

class PiSet : public FuzzySet
{
public:
	PiSet(const char* name);
	void Generate(double Center, double Width);
	void PrintMe();
};

class BetaSet : public FuzzySet
{
public:
	BetaSet(const char* name);
	void Generate(double Center, double Width, double Weight = 1.0);
	void PrintMe();
};

class GaussianSet : public FuzzySet
{
public:
	GaussianSet(const char* name);
	void Generate(double Center, double Width);
	void PrintMe();
};
class InterpSet : public FuzzySet
{
public:
	InterpSet(const char* name);
	void Generate(double pts[], int n);
	void PrintMe();
};
class TriangleSet : public FuzzySet
{
public:
	TriangleSet(const char* name);
	void Generate(double Left, double Center, double Right);
	void PrintMe();
};

class ShoulderSet : public FuzzySet
{
public:
	ShoulderSet(const char* name);
	void Generate(double Edge, double Floor, int type);
	void PrintMe();
};

FuzzySet* CreateFuzzySet(const char* name, int Type, double Dom[2],
	double p1=0, double p2=0, double p3=0, double p4=0);
*/

#endif
