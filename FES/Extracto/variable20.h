//---------------------------------------------------------------------------

#ifndef Variable20H
#define Variable20H
//---------------------------------------------------------------------------
#include "fuzzyset20.h"
#include <string> 

#define fzBOOL 0
#define fzINT 1
#define fzFLOAT 2

class Property
{
public:
	Property(const char* name, int type, double dom[],Property* child)
		: Name(name)
		{Types=NULL;Type = type; next = child;Domain[0] = dom[0];Domain[1] = dom[1];value = 0.5*(dom[0]+dom[1]);}
	Property(Property* prop);
	~Property();
	std::string Name;
    int num;
	FuzzySet*  AddType(const char* name,FuzzySet* FS)
		{FS->SetNext(Types);Types = FS; return Types;}
	void Decorate(int num, bool FuzzFlag=true);
	int Type;
	Property* next;
	FuzzySet* Types;
	double Domain[2];
	void Read(FILE* is);
	void Write(FILE* os);
	double value;
	double GetVal(){return value;}
};
#endif
