//---------------------------------------------------------------------------


#include "variable20.h"
#include<string>

//---------------------------------------------------------------------------

Property::Property(Property* prop)
{	

	Domain[0] = prop->Domain[0];
	Domain[1] = prop->Domain[1];
	Name  = prop->Name;
	num= prop->num;
	Type = prop->Type;
	Types = NULL;
	value = prop->value;
}
Property::~Property()
{
	FuzzySet * cur = Types;
	while(cur)
	{
		Types=cur->Next();
		delete cur;
		cur = Types;
	}
}

void Property::Read(FILE* is)
{
	char buf[255];
	fscanf(is,"\tVariable: %s\n",buf); Name = buf;
	fscanf(is,"\tType = %d Domain = (%lf,%lf)\n",&Type,&(Domain[0]),&(Domain[1]));

	FuzzySet *cur = NULL,*old=NULL;
	int count;
	fscanf(is,"\tNumber of terms: %d\n",&count);

	for(int i=0;i<count;i++)
	{
		cur = new FuzzySet("dummy",0);
		cur->Read(is);
		if(!old) Types = cur;
		else old->SetNext(cur);
		old = cur;
	}
}

void Property::Write(FILE* os)
{
	fprintf(os,"\tVariable: %s\n",Name.c_str());
	fprintf(os,"\tType = %d Domain = (%f,%f)\n",Type,Domain[0],Domain[1]);

	int count = 0;
	FuzzySet *cur = Types;
	while(cur){count++; cur = cur->Next();}

	fprintf(os,"\tNumber of terms: %d\n",count);

	cur = Types;
	for(int i=0;i<count;i++)
	{
		cur->Write(os);
		cur = cur->Next();
	}
}

void Property::Decorate(int num,bool FuzzFlag)
{
	double d0 = Domain[0];
	double d1 = Domain[1];

   char ret[64];

	if(FuzzFlag)
	{
		double pos = d0+(d1-d0)/(num+1);
		double halfw = (d1-d0)/(num+1);
		FuzzySet* temp = new FuzzySet("Low",SHOULDER_DEC);
	
		temp->SetDomain(d0,d1);
		double p0=pos,p1=pos+halfw,p2=pos+2*halfw;
		temp->Generate(pos,pos+halfw,pos+2*halfw);
		//	temp->SetNext(Types); ASG: otherwise, list of fuzzysets becomes larger and larger 
		temp->SetNext(NULL);	

		Types = temp;	

		for(int i=1;i<num-1;i++)
		{
			if(num-1 == 2) sprintf(ret,"Medium");
			else if(num-1 == 4)
			{
				if(i == 1) sprintf(ret,"Medium_Low");
				if(i == 2) sprintf(ret,"Medium");
				if(i == 3) sprintf(ret,"Medium_High");
			} else sprintf(ret,"Set_%d",i);

		//	temp = new FuzzySet(ret,TRIANGLE);
			temp = new FuzzySet(ret,SHOULDER_DEC);  // ASG: Different fuzz strategy
			temp->SetDomain(d0,d1);
			temp->Generate(p0,p1,p2);
			temp->SetNext(Types);
			Types = temp;
			p0 += halfw; p1 += halfw; p2 += halfw;
		}

	//	temp = new FuzzySet("High",SHOULDER_INC);
		temp = new FuzzySet("High", SHOULDER_DEC);
		temp->SetDomain(d0,d1);
		temp->Generate(p0,p1,p2);
		temp->SetNext(Types);
		Types = temp;
	}
	else
	{
		double segment = (d1-d0)/num;
		double init = d0;
		FuzzySet* temp = new FuzzySet("Low",SQUARE);	
		temp->SetDomain(d0,d1);
		temp->Generate(init,init+segment,0);
		//	temp->SetNext(Types); ASG: otherwise, list of fuzzysets becomes larger and larger 
		temp->SetNext(NULL);	

		Types = temp;	

		for(int i=1;i<num-1;i++)
		{
			if(num-1 == 2) sprintf(ret,"Medium");
			else if(num-1 == 4)
			{
				if(i == 1) sprintf(ret,"Medium_Low");
				if(i == 2) sprintf(ret,"Medium");
				if(i == 3) sprintf(ret,"Medium_High");
			} else sprintf(ret,"Set_%d",i);
			temp = new FuzzySet(ret,SQUARE);
			temp->SetDomain(d0,d1);
			init +=segment;
			temp->Generate(init,init+segment,0);
			temp->SetNext(Types);
			Types = temp;
		}
		temp = new FuzzySet("High",SQUARE);
		init +=segment;
		temp->SetDomain(d0,d1);
		temp->Generate(init,init+segment,0);
		temp->SetNext(Types);
		Types = temp;
	}
}
