//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
using namespace std;
//---------------------------------------------------------------------------
//----------------------------------------------------------  NaiveBayes  -----
//---------------------------------------------------------------------------

NaiveBayes::NaiveBayes():Classifier()
{
}

NaiveBayes::~NaiveBayes()
{
}
//---------------------------------------------------------------------------
virtual void NaiveBayes::Build(Data *data, FuncionDeProgreso *fp)
{
  
}
//---------------------------------------------------------------------------
vector<double> NaiveBayes::UnnormalizedDistribution(int ElementIndex)
{
  vector<double> distrb;
  vector<double>& weights = GetUsingWeights();

  for (int i=0;i<((NomData*)data)->NumClass;i++) distrb.push_back(0.0);

  for(int i=0;i<GetClassifiersToUse();i++) {
    Classifier *c = Classifiers.at(i);
    int clase = c->Classify(ElementIndex);

    distrb[clase] += weights[i];
  }


  return distrb;
}
