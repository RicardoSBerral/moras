//---------------------------------------------------------------------------
// 
//
//---------------------------------------------------------------------------
#ifndef NaiveBayesH
#define NaiveBayesH
#include "Classifier.h"

//---------------------------------------------------------------------------
class NaiveBayes : public Classifier
{
protected:

private:

public:
  NaiveBayes();
  virtual ~NaiveBayes();

  //Funcies virtuales sobreescritas
  virtual void Build(Data *data, FuncionDeProgreso *fp);
  virtual std::string Info(int value=0);

  virtual std::vector<double> UnnormalizedDistribution(int ElementIndex);
  std::vector<double> Margen(int ini, int fin, int div=-1);

};

#endif

