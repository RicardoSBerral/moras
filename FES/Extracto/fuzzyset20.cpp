 #include<string>
//---------------------------------------------------------------------------

#include "fuzzyset20.h"

//---------------------------------------------------------------------------

#include <math.h>

//extern void Display(const char* s);
double SCurveValue(double val, double left,double flex, double right);

FuzzySet::FuzzySet(FuzzySet* fs)
{
  _Name = fs->_Name;
  _Description = fs->_Description;;
  _CurveType = fs->_CurveType;
  _IsEmpty = fs->_IsEmpty;
  _Order = fs->_Order;
  _Domain[0] = fs->_Domain[0];
  _Domain[1] = fs->_Domain[1];
  _Params[0] = fs->_Params[0];
  _Params[1] = fs->_Params[1];
  _Params[2] = fs->_Params[2];
  _Params[3] = fs->_Params[3];
  _Alpha = fs->_Alpha;
  _Status = fs->_Status;
  _Next = NULL;
  for(int i=0;i<VECMAX;i++){_MemVec[i] = fs->_MemVec[i];}

}



void FuzzySet::Generate(double p0, double p1, double p2)
{
  if(_Domain[1] <= _Domain[0]) {_Status = BadDomain; return;}
  if(CheckParams(p0,p1,p2)) return;

  _Params[0] = p0;
  _Params[1] = p1;
  _Params[2] = p2;

  double dwidth = _Domain[1]-_Domain[0];
  for(int i=0;i<VECMAX;i++)
  {
    double val = _Domain[0]+((double)i/VECMAX)*dwidth;
    double temp;
    _MemVec[i] = 0.0;
    switch(_CurveType)
    {
      case LIN_INC:    // p0 = high, p1 = low
      case LIN_DEC:
        if(val > p1) _MemVec[i] = 1.0;
        else if(val > p0) _MemVec[i] = (val-p0)/(p1-p0);
        break;
      case SIG_INC:    //p0 = left, p1 = flex, p2 = right
      case SIG_DEC:
        _MemVec[i] = SCurveValue(val,p0,p1,p2);
        break;
      case PISET:      // p0 = Center, p1 = Width
        _MemVec[i] = 1.0-SCurveValue(val,p0,p0+0.5*p1,p0+p1);
        temp = SCurveValue(val,p0-p1,p0-0.5*p1,p0);
        if(temp < _MemVec[i]) _MemVec[i] = temp;
        break;
      case BETASET:    // p0 = Center, p1 = Flex, p2 = Weight
        _MemVec[i] = 1.0/(1+p2*((val-p0)/p1)*((val-p0)/p1));
        break;
      case GAUSSIAN:    // p0 = Center, p1 = Width
        _MemVec[i] = exp(-p1*(val-p0)*(val-p0));
        break;
      case TRIANGLE:    // p0 = Left, p1 = Center, p2 = Right
        if(val <= p0 || val > p2) _MemVec[i] = 0.0;
        else if(val < p1) _MemVec[i] = (val-p0)/(p1-p0);
        else if(val <= p2) _MemVec[i] = 1-(val-p1)/(p2-p1);
        break;
      case SQUARE:
        if(val <= p0 || val >= p1) _MemVec[i] =0.0;
        else _MemVec[i] = 1.0;
        break;
      case SHOULDER_INC:    // p0 = Edge, p1 = Floor,
      case SHOULDER_DEC:
        if(val < p0) _MemVec[i] = 0.0;
        else if(val < p1) _MemVec[i] = (val-p0)/(p1-p0);
        else _MemVec[i] = 1.0;
        break;
    }
    if(_CurveType == LIN_DEC || _CurveType == SIG_DEC || _CurveType == SHOULDER_DEC)
      _MemVec[i] = 1.0-_MemVec[i];
  }

}

int FuzzySet::CheckParams(double p0, double p1, double p2)
{

  switch(_CurveType)
  {
    case LIN_INC:    // p0 = high, p1 = low
    case LIN_DEC:
      if(p0 >= p1) {_Status = BadDomain;return 1;}
      break;
    case SIG_INC:    //p0 = left, p1 = flex, p2 = right
    case SIG_DEC:
      if(p0>=p1 || p1>=p2){_Status = BadDomain;return 1;}
      break;
    case PISET:      // p0 = Center, p1 = Width
      break;
    case BETASET:    // p0 = Center, p1 = Flex, p2 = Weight
      break;
    case GAUSSIAN:    // p0 = Center, p1 = Width
      break;
    case TRIANGLE:    // p0 = Left, p1 = Center, p2 = Right
      break;
    case SQUARE:
      break;
    case SHOULDER_INC:    // p0 = Edge, p1 = Floor,
    case SHOULDER_DEC:
      break;
  }
  return 0;
}

void FuzzySet::PrintMe()
{
/*  char buf[255];
  switch(_CurveType)
  {
    case LIN_INC:    // p0 = high, p1 = low
    case LIN_DEC:
      sprintf(buf,"        LowEdge  = %g",_Params[0]);Display(buf);
      sprintf(buf,"   HighEdge = %g\r\n",_Params[1]);Display(buf);
      break;
    case SIG_INC:    //p0 = left, p1 = flex, p2 = right
    case SIG_DEC:
      sprintf(buf,"        LeftEdge  = %g",_Params[0]);Display(buf);
      sprintf(buf,"   Flex      = %g",_Params[1]);Display(buf);
      sprintf(buf,"   RightEdge = %g\r\n",_Params[2]);Display(buf);
      break;
    case PISET:      // p0 = Center, p1 = Width
        sprintf(buf,"        Center  = %g",_Params[0]);Display(buf);
      sprintf(buf,"   Width      = %g\r\n",_Params[1]);Display(buf);
      break;
    case BETASET:    // p0 = Center, p1 = Flex, p2 = Weight
      sprintf(buf,"        Center  = %g",_Params[0]);Display(buf);
      sprintf(buf,"   Flex    = %g",_Params[1]);Display(buf);
      sprintf(buf,"   Weight  = %g\r\n",_Params[2]);Display(buf);
      break;
    case GAUSSIAN:    // p0 = Center, p1 = Width
      sprintf(buf,"        Center  = %g",_Params[0]);Display(buf);
      sprintf(buf,"   Width    = %g\r\n",_Params[1]);Display(buf);
      break;
    case TRIANGLE:    // p0 = Left, p1 = Center, p2 = Right
        sprintf(buf,"        Left  = %g",_Params[0]);Display(buf);
      sprintf(buf,"   Center    = %g",_Params[1]);Display(buf);
      sprintf(buf,"   Right    = %g\r\n",_Params[2]);Display(buf);
      break;
    case SHOULDER_INC:    // p0 = Edge, p1 = Floor,
    case SHOULDER_DEC:
        sprintf(buf,"        Edge  = %g",_Params[0]);Display(buf);
      sprintf(buf,"   Floor    = %g",_Params[1]);Display(buf);
      break;
  }
  return;
  */
}



FuzzySet::FuzzySet(const char* name, int type)
{
  _Name = name;
  _Next = NULL;
  _IsEmpty = true;
  _Alpha=0;
  _Status = OK;
  _Domain[0] = _Domain[1] = 0.0;
  _CurveType = type;
   _Description = "FuzzySet";
}
#define PLOTROWS   26
#define PLOTCOLS   70
#define NUDGE 0.00001

void FuzzySet::Draw()
{
  /*
  char        WkArea[PLOTROWS][PLOTCOLS+1];
  int         i,j,k,Sidx,PltHgt,N,HorzPos,VertPos;
  int         ScalingIdx=1,ScaleCtl  =4;
  double      x;
  char buf[1024];

  const char  Symbol='.';
  const int   PlotHeight[]    = {0,11,21,26,51};
  const double ScalingFactor[] = {0,10.0,20.0,25.0,50.0};

//--Blank out the plot area and then put a std::string terminator
//--at the end of each row. This is used in the fprint output.
  for(i=0;i<PLOTROWS;i++)
   {
     for(j=0;j<PLOTCOLS;j++)  WkArea[i][j]=' ';
    WkArea[i][PLOTCOLS]='\0';
   }
//--Fill out the plot MxN area.
  Sidx=ScalingFactor[ScalingIdx];
  for(k=0;k<VECMAX;k+=ScaleCtl)
    {
     VertPos=(int)(Sidx*(_MemVec[k]+NUDGE));
     HorzPos= k / ScaleCtl;
     if((VertPos+1>=0)&&(VertPos+1<PLOTROWS))
      WkArea[VertPos+1][HorzPos]=Symbol;
    }
//--Now we write out the plot area
  buf[0] = '\0';
  Display("\r\n");
  Display("\r\n");
  sprintf(buf,"%s%s\r\n","        FuzzySet:    ",_Name);
  Display(buf);
  sprintf(buf,"%s%s\r\n","        Description: ",_Description);
  Display(buf);
  PrintMe();
  PltHgt=PlotHeight[ScalingIdx];
  for(N=PltHgt,i=0;i<PltHgt;i++,--N)
    {
    j=PltHgt-i;
     x=(j-1)/ScalingFactor[ScalingIdx];
     sprintf(buf,"%10.2f%s\r\n",x,&WkArea[N][0]);
     Display(buf);
    }
//  sprintf(buf,"%s%s",
//   "          0___|___|___|___|___|___|___|___|___|___|___|___|",
//   "___|___|___|___0\r\n");
  sprintf(buf,"%s%s",
   "          0________________________________________________",
   "_______________0\r\n");

   Display(buf);
  x=(_Domain[1]-_Domain[0])/8.0;
  Display("     ");
  for(i=0;i<9;i++) {sprintf(buf,"%8.2f",_Domain[0]+(x*i));Display(buf);}
  Display("\r\n");
  sprintf(buf,
    "%s%10.2f%s%10.2f\r\n", "        Domained:   ",
     _Domain[0]," to ",_Domain[1]);
    Display(buf);
  if(_Alpha>0)
  {
    sprintf(buf,
    "%s%10.2f\r\n",        "         AlphaCut:   ",_Alpha);
    Display(buf);
  }
  Display("\r\n");
  */
}

int FuzzySet::SetDomain(double dom_min, double dom_max)
{
  if(dom_min >= dom_max) return 1;
  _Domain[0] = dom_min;
  _Domain[1] = dom_max;
  return 0;
}

bool FuzzySet::IsNormal()
{
  for(int i=0;i<VECMAX;i++)
    if(_MemVec[i] == 1.0) return true;
  return false;
}

double FuzzySet::GetHeight()
{
  double maxval = _MemVec[0];
  for(int i=1;i<VECMAX;i++)
    if(_MemVec[i] > maxval) maxval = _MemVec[i];
  return maxval;
}

void FuzzySet::Normalise()
{
  double maxval = GetHeight();
  for(int i=0;i<VECMAX;i++)
    _MemVec[i] /= maxval;
}

double FuzzySet::Membership(double Val)
{
   if(Val<_Domain[0])
    return( _MemVec[0] < _Alpha ? 0.0 : _MemVec[0]);
   if(Val>_Domain[1])
  return( _MemVec[VECMAX-1] < _Alpha ? 0.0 : _MemVec[VECMAX-1]);
  double Range=_Domain[1] - _Domain[0];
   if(Range==0)
   {
    _Status = BadRange;
    return(0);
   }
  int pos=(int)(((Val-_Domain[0])/Range)*(VECMAX-1));
  return(_MemVec[pos]<_Alpha ? 0.0 : _MemVec[pos]);
}

double FuzzySet::ValueFromTruth(double TruthVal)
{
  _Status = OK;
  for(int i=0;i<VECMAX;i++)
    if(TruthVal >= _MemVec[i])
    {
      double Range=_Domain[1]-_Domain[0];
      return (_Domain[0]+(i*Range)/(VECMAX-1));
    }
  _Status = NotFound;
  return 0;
}
void FuzzySet::Read(FILE*  is)
{

  char buf[255];


  fscanf(is,"\tFuzzy Set: %[^'\n']\n",buf);
  _Name = buf;
  fscanf(is,"\tFuzzy Description: %[^'\n']\n",buf);
  _Description=buf;
  fscanf(is,"\t_CurveType: %d\n",&_CurveType);
fscanf(is,"\t_IsEmpty: %d\n",&_IsEmpty);
fscanf(is,"\t_Order: %d\n",&_Order);
fscanf(is,"\t_Domain: %lf %lf\n",&_Domain[0],&_Domain[1]);
fscanf(is,"\t_Params: %lf %lf %lf %lf\n",&_Params[0],&_Params[1],&_Params[2],&_Params[3]);
fscanf(is,"\t_Alpha: %lf\n",&_Alpha);
fscanf(is,"\t_Status: %d\n",&_Status);
for(int i=0;i<VECMAX;i++) fscanf(is,"%lf ",&_MemVec[i]);
}

void FuzzySet::Write(FILE*  os)
{
  fprintf(os,"\tFuzzy Set: %s\n",_Name.c_str());
  fprintf(os,"\tFuzzy Description: %s\n",_Description.c_str());
  fprintf(os,"\t_CurveType: %d\n",_CurveType);
fprintf(os,"\t_IsEmpty: %d\n",_IsEmpty);
fprintf(os,"\t_Order: %d\n",_Order);
fprintf(os,"\t_Domain: %f %f\n",_Domain[0],_Domain[1]);
fprintf(os,"\t_Params: %f %f %f %f\n",_Params[0],_Params[1],_Params[2],_Params[3]);
fprintf(os,"\t_Alpha: %f\n",_Alpha);
fprintf(os,"\t_Status: %d\n",_Status);
for(int i=0;i<VECMAX;i++) fprintf(os,"%f ",_MemVec[i]);
fprintf(os,"\n");

}

double SCurveValue(double val, double left,double flex, double right)
{
  if(val > right) return 1.0;
  else if(val > flex) return 1.0-0.5*(val-right)*(val-right)/((right-flex)*(right-flex));
  else if(val > left) return 0.5*(val-left)*(val-left)/((flex-left)*(flex-left));
  else return 0.0;
}


void FuzzySet::Defaults(int id,double& p1,double& p2,double& p3,double Domain[2])
{
  double w = Domain[1] - Domain[0];
  p1 = p2 = p3 = 0.0;
  switch(id)
  {
    case LIN_INC:    // p0 = high, p1 = low
    case LIN_DEC:
      p1=0.33*Domain[0];p2=0.66*Domain[1];
      break;
    case SIG_INC:    //p0 = left, p1 = flex, p2 = right
    case SIG_DEC:
      p1=Domain[0]+0.33*w;p2=Domain[0]+0.5*w; p3=Domain[0]+0.66*w ;
      break;
    case PISET:      // p0 = Center, p1 = Width
      p1=0.5*(Domain[0]+Domain[1]),p2=0.25*(Domain[1]-Domain[0]);
      break;
    case BETASET:    // p0 = Center, p1 = Flex, p2 = Weight
      p1=Domain[0]+0.5*w,p2=0.25*w, p3=1.0;
      break;
    case GAUSSIAN:    // p0 = Center, p1 = Width
      p1=0.5*(Domain[0]+Domain[1]),p2=0.1/(Domain[1]-Domain[0]);
      break;
    case TRIANGLE:    // p0 = Left, p1 = Center, p2 = Right
      p1=Domain[0]+0.33*w,p2=Domain[0]+0.5*w, p3=Domain[0]+0.66*w ;
      break;
    case SQUARE:
      p1 = Domain[0] + 0.25*w; 
      p2 = Domain[0] + 0.75*w; 
      break;
    case SHOULDER_INC:    // p0 = Edge, p1 = Floor,
    case SHOULDER_DEC:
      p1=0.33*Domain[0],p2=0.66*Domain[1];
      break;
  }
}
