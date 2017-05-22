//---------------------------------------------------------------------------
// 
// 
//---------------------------------------------------------------------------
#ifndef EvaluationsH
#define EvaluationsH

#include <vector>
#include <string>

class Classifier;
class Data;
class Matriz;
//---------------------------------------------------------------------------
//--------------------------------------------------------  Evaluation  -----
//---------------------------------------------------------------------------
class Evaluation
{
  public:
    Evaluation(Classifier *clf, Data *data, std::string params=std::string()){}
    virtual ~Evaluation(){}

  public:
    virtual double MainResult()=0;
    virtual std::string ToString();
    virtual int NumberOfOutputs()=0;
    virtual std::string OutputName(int iOutput)=0;
    virtual std::string Output(int i)=0;

  public:
    static Evaluation* EvaluationByClassName(std::string EvaluationClassName, 
                           Classifier *c, Data *data, std::string params=std::string());
};
class ErrorEvaluation : public Evaluation
{
  protected:
    double error;
    Matriz *matriz;

  public:
    ErrorEvaluation(Classifier *clf, Data *data, std::string params="");
    virtual ~ErrorEvaluation();

  public:
    virtual double MainResult() {return error;};
    virtual int NumberOfOutputs();
    virtual std::string OutputName(int iOutput);
    virtual std::string Output(int i);
};
class MSEEvaluation : public Evaluation
{
  protected:
    double rms;

  public:
    MSEEvaluation(Classifier *clf, Data *data, std::string params="");
    virtual ~MSEEvaluation(){}

  public:
    virtual double MainResult() {return rms;};
    virtual int NumberOfOutputs();
    virtual std::string OutputName(int iOutput);
    virtual std::string Output(int i);
};
class MultipleMSEEvaluation : public Evaluation
{
  protected:
    double rms;

  public:
    MultipleMSEEvaluation(Classifier *clf, Data *data, std::string params="");
    virtual ~MultipleMSEEvaluation(){}

  public:
    virtual double MainResult() {return rms;};
    virtual int NumberOfOutputs();
    virtual std::string OutputName(int iOutput);
    virtual std::string Output(int i);
};
class ComboEvaluation : public Evaluation
{
  protected:
    double error;
    Matriz *matriz;

  public:
    ComboEvaluation(Classifier *clf, Data *data, std::string params="");
    virtual ~ComboEvaluation();

  public:
    virtual double MainResult() {return error;};
    virtual int NumberOfOutputs();
    virtual std::string OutputName(int iOutput);
    virtual std::string Output(int i);
};
class MultilabelHammingLossEvaluation : public Evaluation
{
  protected:
    double error;
    Matriz *matrizz;

  public:
    MultilabelHammingLossEvaluation(Classifier *clf, Data *data, std::string params="");
    virtual ~MultilabelHammingLossEvaluation();

  public:
    virtual double MainResult() {return error;};
    virtual int NumberOfOutputs();
    virtual std::string OutputName(int iOutput);
    virtual std::string Output(int i);
};
class EvidenceEvaluation : public Evaluation
{
  protected:
    std::vector< std::vector<double> > seq_errors;

  public:
    EvidenceEvaluation(Classifier *clf, Data *data, std::string params="");
    virtual ~EvidenceEvaluation();

  public:
    virtual double MainResult() {return 0.0;};
    virtual int NumberOfOutputs();
    virtual std::string OutputName(int iOutput);
    virtual std::string Output(int i);
};
class StackingDataEvaluation : public Evaluation
{
  protected:
    std::string base_file_name;

  public:
    StackingDataEvaluation(Classifier *clf, Data *data, std::string params="");
    virtual ~StackingDataEvaluation();

  public:
    virtual double MainResult() {return 0.0;};
    virtual int NumberOfOutputs();
    virtual std::string OutputName(int iOutput);
    virtual std::string Output(int i);
};
//---------------------------------------------------------------------------
#endif
