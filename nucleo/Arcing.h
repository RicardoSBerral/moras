//---------------------------------------------------------------------------
// ClassExplorer Pro generated header file
// Created by  on 04/12/00, 18:48:28
//---------------------------------------------------------------------------
#ifndef BaggingH
#define BaggingH
//#include "Classifier.h"
#include "Ensemble.h"
#include <vector>

//---------------------------------------------------------------------------
class GPTree;
//---------------------------------------------------------------------------
//----------------------------------------------------  ArcingTemplate  -----
//---------------------------------------------------------------------------

class ArcingTemplate : public Ensemble
{

  protected:
    virtual void StartArc() {return;}
    virtual Data* Resample(int iClas) {return data;}
    virtual void AfterBuilding(int ic, bool &Continue, bool &ValidClassifer) {}
    virtual void EndArc() {return;}

  public:
    ArcingTemplate();
    virtual ~ArcingTemplate();

    virtual void Build(Data *data, FuncionDeProgreso *fp=0);

};
//---------------------------------------------------------------------------
//------------------------------------------------------  GroupBagging  -----
//---------------------------------------------------------------------------

class GroupBagging : public ArcingTemplate
{

  protected:
    virtual Data* Resample(int iClas);

  public:
    GroupBagging();
    virtual ~GroupBagging();

};
//---------------------------------------------------------------------------
//-----------------------------------------------------------  Bagging  -----
//---------------------------------------------------------------------------

class Bagging : public ArcingTemplate
{

  protected:
    int ResampleSize;
    bool WithReplacement;
    bool EqualizeClasses;

    virtual void EndArc();
    virtual Data* Resample(int iClas);
    virtual void AfterBuilding(int ic, bool &Continue, bool &ValidClassifer);

  public:
    Bagging(int _ResampleSize=-1, bool _WithReplacement=true, 
                                                bool _EqualizeClasses=false);
    virtual ~Bagging();

  public:
    double *EvidenceWeights;
};
//---------------------------------------------------------------------------
//-----------------------------------------------------------  Wagging  -----
//---------------------------------------------------------------------------

class Wagging : public ArcingTemplate
{

  protected:
    virtual Data* Resample(int iClas);

  public:
    Wagging();
    virtual ~Wagging();

};
//---------------------------------------------------------------------------
//----------------------------------------------------------  Boosting  -----
//---------------------------------------------------------------------------

class Boosting : public ArcingTemplate
{
  protected:
    virtual void StartArc();
    virtual void EndArc();
    virtual Data* Resample(int iClas);
    virtual void AfterBuilding(int ic, bool &Continue, bool &ValidClassifer);
    virtual void RandomReweight(Data *data, int first, int last);

    //Funciones para guardar y leer de fichero
    bool ForceContinuationByReweighting;
    bool NewBootstrap;
    std::vector<int> ReweightingPoints;
    bool AsInWebb00multiboosting;
    bool AsInBauer99empirical;
    double MinWeight;

  public:
    Boosting();
    virtual ~Boosting();

  public:
    //
    virtual std::string Info(int value=0);

    //Funciones para guardar y leer de fichero
    virtual void Guardar(std::ostream &salida, int version=0);
    virtual void Leer(std::istream &in, int version=0);
  
    bool GetForceContinuationByReweighting() {
      return ForceContinuationByReweighting;
    }
    void SetForceContinuationByReweighting(bool value) {
      ForceContinuationByReweighting = value;
    }
    double GetMinWeight() { return MinWeight;}
    void SetMinWeight(double value) { MinWeight = value;}
    void SetAsInWebb00multiboosting();
    void SetAsInBauer99empirical();
    void SetAdaBoost();
    void SetDefaultClassifiers();
};
//---------------------------------------------------------------------------
//--------------------------------------------------------  MIBoosting  -----
//---------------------------------------------------------------------------

class MIBoosting : public Boosting
{

  protected:
    virtual void StartArc();
    virtual void EndArc();
    virtual Data* Resample(int iClas);
    virtual void AfterBuilding(int ic, bool &Continue, bool &ValidClassifer);

    void ResetWeights();

  protected:
    double* Weights;
    bool Resampling;
  
  public:
    MIBoosting();
    virtual ~MIBoosting();

    virtual double Error(int first, int last);
  
    virtual void SecuencialError(int first, int last, std::vector<double> *errors,
                std::vector<double> *inderrs=0, std::vector<int> *final_class=0);
    virtual void SecuencialClassify(int ElementIndex,
                         std::vector<int> *classes, std::vector<int> *indclss=0);
};
//---------------------------------------------------------------------------
//---------------------------------------------------------  MIBagging  -----
//---------------------------------------------------------------------------
class MIBagging : public MIBoosting   // Mal heredado, hay que pasarlo a 
// herencia multiple y que herede de MIClassifier y de MIEnsemble
{
  protected:
    int ResampleSize;
    std::vector<std::vector<double> > DistributionWeights;
    std::vector<int> idxs_oob;
    std::vector<std::vector<int> > idxs_oobs;
    double **dout;
    Data *ndout;
    Data *ndout_all;
    int **ndout_all_votes;
    bool MaxSingleInstanceMargin; 
    int PositiveClass; 
    bool UseOOBForNodePops; 
    int OutputLevel;
    Matriz *q;

  protected:
    virtual void CalcDiffs(Data* d, int locualo);
    virtual void StartArc();
    virtual void EndArc();
    virtual Data* Resample(int iClas);
    virtual void AfterBuilding(int ic, bool &Continue, bool &ValidClassifer);
  
  public:
    MIBagging(int _ResampleSize=-1);
    virtual ~MIBagging();

    virtual double Error(int first, int last);
  
    virtual void SecuencialClassify(int ElementIndex,
                         std::vector<int> *classes, std::vector<int> *indclss=0);

    void SetOutputLevel(int level) { OutputLevel = level; }
};
//---------------------------------------------------------------------------
//-------------------------------------------  BoostingWithBaggingInfo  -----
//---------------------------------------------------------------------------

class BoostingWithBaggingInfo : public Boosting
{
protected:
  double *percents;
  double Umbral;
  Ensemble *ens;
  std::vector<int> excluidos;
  char auxinfo[256];

  virtual void StartArc();
  virtual void AfterBuilding(int ic, bool &Continue, bool &ValidClassifer);
  virtual void EndArc();

public:
  BoostingWithBaggingInfo(double Umbral=0.90, Ensemble *ens=0);
  virtual ~BoostingWithBaggingInfo();

//  void Build(Data *data);
  std::string Info(int value=0);
};
//---------------------------------------------------------------------------
//---------------------------------------------------  BoostingGPTrees  -----
//---------------------------------------------------------------------------
class BoostingGPTrees : public Boosting
{
protected:
  int NMin;
  int nMaxClassifiers;
  void (*salida)(char *);

public:
  BoostingGPTrees(void (*_salida)(char *), int nClassifiers, int _NMin=0);

  void Init(Data *data);

};
//---------------------------------------------------------------------------
//------------------------------------------------  InvBoostingGPTrees  -----
//---------------------------------------------------------------------------
class InvBoostingGPTrees : public BoostingGPTrees
{
protected:
  virtual void AfterBuilding(int ic, bool &Continue, bool &ValidClassifer);

public:
  InvBoostingGPTrees(void (*_salida)(char *), int nClassifiers, int _NMin=0);
};
//---------------------------------------------------------------------------
//--------------------------------------------------  IndepClassifiers  -----
//---------------------------------------------------------------------------

class IndepClassifiers : public ArcingTemplate
{
protected:
  virtual void StartArc();
  virtual Data* Resample(int iClas);
  virtual void AfterBuilding(int ic, bool &Continue, bool &ValidClassifer);
  virtual void EndArc();

  double factor;
  double factor_incr;
  double e_ok;
  double e_nook;
  double e_train;
  bool AsBreiman;
  bool UseDummyClass;
  int dummy_class;
  bool GenerateOutputsAsKuncheva;
  int **outputs;

public:
  IndepClassifiers(double factor=0.386, bool UseDummyClass=false, 
                                                       double factor_incr=0.0);
  virtual ~IndepClassifiers();

  virtual void Init(Data *data);
  virtual void SetData(Data *data);
  virtual int Classify(std::vector<double> distrb);
  using ArcingTemplate::Classify;//para no ocultar Classify(int)

  void SetAsBreiman(bool Value) {AsBreiman = Value;}
  void SetGenerateOutputsAsKuncheva(bool Value) {GenerateOutputsAsKuncheva = Value;}
  std::string Info(int value=0);
};
//---------------------------------------------------------------------------
//-----------------------------------------------------  PastingIVotes  -----
//---------------------------------------------------------------------------

class PastingIVotes : public ArcingTemplate
{
protected:
  int ResampleSize;
  double error_estimate;
  double p;

protected:
  virtual Data* Resample(int iClas);
  virtual void AfterBuilding(int ic, bool &Continue, bool &ValidClassifer);

public:
  PastingIVotes(int );
  virtual ~PastingIVotes();

};
//---------------------------------------------------------------------------
//---------------------------------------------  CompetentModelBagging  -----
//---------------------------------------------------------------------------
#define COMPBAG
class CompetentModelBagging : 
#ifdef COMPBAG
                              public Bagging
{
protected:
  virtual void CalcWeights(int iClas);
#else
                              public Boosting
{
protected:
  virtual void CalcWeights(Classifier *Clasf, bool &Continue);
#endif
  std::vector<Classifier*> Referees;
  int clasevotes;
  int clasevotesC;
  int clasevotesW;
  int mconf[2][2];
  bool UseReferees;
  
 public:
  CompetentModelBagging(){UseReferees = true;}
  virtual int Classify(int ElementIndex);
  virtual double Error(int first, int last);
  
  void SetUseReferees(bool Value) {UseReferees = Value;}
  bool GetUseReferees() {return UseReferees;}
  
  void GetVotes(int ElementIndex, std::vector<int> &WithRefs, 
                                                      std::vector<int> &NoRefs);
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

class BaggingReg : public ClsfReg
{
  private:
    static BaggingReg *reg;
    static BaggingReg *autoreg(){return new BaggingReg();}
    BaggingReg() : ClsfReg("Bagging", "Baggnig", true){}

  public:
    Classifier *CreateClassifier(){ return new Bagging();}

    static BaggingReg *Reg(){return reg;}
};
//---------------------------------------------------------------------------
class BoostingWithBaggingInfoReg : public ClsfReg
{
    static BoostingWithBaggingInfoReg *reg;
    static BoostingWithBaggingInfoReg *autoreg() {
      return new BoostingWithBaggingInfoReg();
    }
    BoostingWithBaggingInfoReg() : ClsfReg("BoostingWithBaggingInfo",
                                           "BoostingWithBaggingInfo", true){}
  public:
    Classifier *CreateClassifier(){ return new BoostingWithBaggingInfo();}

    static BoostingWithBaggingInfoReg *Reg(){return reg;}
};
//---------------------------------------------------------------------------
class BoostingReg : public ClsfReg
{
  private:
    static BoostingReg *reg;
    static BoostingReg *autoreg(){return new BoostingReg();}
    BoostingReg() : ClsfReg("Boosting", "Boosting", true){}

  public:
    Classifier *CreateClassifier(){ return new Boosting();}

    static BoostingReg *Reg(){return reg;}
};
//---------------------------------------------------------------------------
#endif
