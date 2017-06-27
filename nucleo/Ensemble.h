//---------------------------------------------------------------------------
// ClassExplorer Pro generated header file
// Created by  on 04/12/00, 18:48:28
//---------------------------------------------------------------------------
#ifndef EnsembleH
#define EnsembleH
#include "Classifier.h"

//---------------------------------------------------------------------------
class Ensemble : public Classifier
{
protected:
  void QuickSort(int const AHigh, int iLo, int iHi, double *W=0);
  void DoQuickSort(int const AHigh, int iLo, int iHi, double *A);

  std::vector<Classifier*> Classifiers;
  std::vector<double> ClassifierWeights;

  bool UseWeights;
  int ClassWithWeigths;
  int ClassWithoutWeigths;
  int ClassifiersToUse;

private:
  std::vector<int> OrdenOriginal;
  std::vector<double> PesosOriginales;

  void ResetOrdenOriginal();

public:
  Ensemble();
  virtual ~Ensemble();

  //Funcies virtuales sobreescritas
  virtual void Init(Data * data);
  virtual void Build(Data *data, FuncionDeProgreso *fp=0)=0;
  virtual void SetData(Data *data);
  virtual std::string Info(int value=0);
  //virtual int Classify(int ElementIndex);

  virtual double Average(int ElementIndex);
  virtual double Average(int ElementIndex, int AverageIndex);
  virtual std::vector<double> MultipleAverage(int ElementIndex);
  virtual std::vector<double> UnnormalizedDistribution(int ElementIndex);
//  virtual double ClassificationCertainty(int ElementIndex, int &Class/*, 
//		  bool ExcluirTrain=false*/);
  std::vector<double> Margen(int ini, int fin, int div=-1);
  virtual void SecuencialError(int first, int last, std::vector<double> *errors,
              std::vector<double> *inderrs=0, std::vector<int> *final_class=0);
  virtual void RegressionError(int first, int last, std::vector<double> *errors,
              int dVar=-1);
  virtual void SecuencialClassify(int ElementIndex,
                       std::vector<int> *classes, std::vector<int> *indclss=0);
  //

  //Gestion de la lista de clasificadores
  void AddClassifier(Classifier *Clasf);
  void RemoveClassifier(int i);
  int Count();
  Classifier *GetClassifier(int Index);
  void SetClassifier(int Index, Classifier *c);
  void Exchange(int Index1, int Index2);

  //Orden de los clasificadores
   std::vector<int> DameOrdenOriginal() {return OrdenOriginal;}
  void OrdenarClasificadoresPorPeso(double *W=0);
  void OrdenarClasificadoresPorOrdenOriginal();
  void PonerEsteOrdenComoOrdenOriginal() { ResetOrdenOriginal(); }

  //Gets y Sets varios
  double GetClassifierWeight(int Index);
  double GetWeightInUse(int Index);
  std::vector<double>& GetUsingWeights();
  std::vector<double>& GetWeights(){return GetUsingWeights();}
  void SetClassifierWeight(int Index, double Weight);
  bool GetUseWeights() {return UseWeights;}
  void SetUseWeights(bool Value) {UseWeights=Value;}
  int GetClassifiersToUse() {
    return ClassifiersToUse<=0 ? Count() : ClassifiersToUse;
  }
  void SetClassifiersToUse(int Value) {
    Value = (Value>Count()) ? Count() : Value;
    ClassifiersToUse = Value;
  }
  int GetLastClassWithWeigths(){return ClassWithWeigths;}
  int GetLastClassWithoutWeigths(){return ClassWithoutWeigths;}

  //Funciones para guardar y leer de fichero
  virtual void Guardar(std::ostream &salida, int version=0);
  virtual void Leer(std::istream &in, int version=0);

  //Otras
  int **MatClasif0(Data *data);

  // Computes a matrix accounting for the strength/diversity
  // as defined in Zhang06 Ensemble Pruning Via SDP
  double** computeGMatrix(Data *data);

  // Computes the G matrix and multiplies it by a vector 
  // of ones of size of the ensemble. 
  // This is equivalent to compute a strenght/diversity 
  // measure of a (pruned) ensemble.
  double computeGValue(Data *data);

  // Computes a diversity-strenght measure for the selected classifiers
  // in the array.
  static double evaluateGMatrix(double** G, int size, int *selectedClassifiers, int n);
 
  // Computes isSelected' * G * isSelected
  // isSelected[i] must be 1 if the i classifiers is selected
  //                       0 otherwise
  static double evaluateGMatrix(double** G, int size, int *isSelected);
};

//---------------------------------------------------------------------------
/*class SimpleStatistics
{
  protected:
    std::vector<double> pop;
    double mn;
    double stdv;

    void Calculate();

  public:
    SimpleStatistics(std::vector<double> pop) {
      this->pop = pop;
      Calculate();
    }

    double mean() {return mn;}
    double stdev() {return stdv;}
};*/
//---------------------------------------------------------------------------
/*class EnsembleStatistics
{
  protected:
    Ensemble *ens;
    int **Mcd;
    Data *data;

  public:
    EnsembleStatistics(Ensemble *ens, int **Mcd);
    EnsembleStatistics(Ensemble *ens, Data *data);

    std::vector<double> Distribs();
    std::vector<double> Certezas();
    std::vector<double> Margen();
    std::vector<double> IndivError();
    std::vector<double> SecuencialError();
    double Error();

    void ThStat(double &mean, double &stdev);
    void QStat(double &mean, double &stdev);
    double PairQStat(int i, int j);
    void KStat(double &mean, double &stdev);
    double PairKStat(int i, int j);
};*/
//---------------------------------------------------------------------------
#endif

