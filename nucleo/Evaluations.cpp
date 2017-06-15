//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#include "Evaluations.h"

#include "Arcing.h"
#include "GPTree.h"
#include "data20.h"
#include "node20.h"
#include "Matriz.h"
#include "Utils.h"

#include <stdlib.h>

#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
string Evaluation::ToString() 
{
  string s;
  ostringstream out;

  for(int i = 0; i < NumberOfOutputs(); i++) {
    if (i > 0) out << endl;
    out << Output(i); 
  }

  return out.str();
}
Evaluation* Evaluation::EvaluationByClassName(std::string EvaluationClassName, 
                                                Classifier *c, Data *data, string params) 
{ //This should be done with a class register and the classes should be able to clone themselves
  if (EvaluationClassName.compare("ErrorEvaluation")==0) {
    return new ErrorEvaluation(c, data, params);
  }
  else if (EvaluationClassName.compare("MSEEvaluation")==0) {
    return new MSEEvaluation(c, data, params);
  }
  else if (EvaluationClassName.compare("MultipleMSEEvaluation")==0) {
    return new MultipleMSEEvaluation(c, data, params);
  }
  else if (EvaluationClassName.compare("ComboEvaluation")==0) {
    return new ComboEvaluation(c, data, params);
  }
  else if (EvaluationClassName.compare("MultilabelHammingLossEvaluation")==0) {
    return new MultilabelHammingLossEvaluation(c, data, params);
  }
  else if (EvaluationClassName.compare("EvidenceEvaluation")==0) {
    return new EvidenceEvaluation(c, data, params);
  }
  else if (EvaluationClassName.compare("StackingDataEvaluation")==0) {
    return new StackingDataEvaluation(c, data, params);
  }
  return 0;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
ErrorEvaluation::ErrorEvaluation(Classifier *clf, Data *data, string params) : Evaluation(clf, data)
{
  matriz = new Matriz(((NomData*)data)->NumClass);
  error = clf->Error(data);
  for(int i = 0; i < (int)clf->rc.size(); i++) {
    (*matriz)[clf->rc[i]][clf->pc[i]] += data->GetInstance(i).GetWeight();
  }
}
ErrorEvaluation::~ErrorEvaluation()
{
  delete matriz;
}
int ErrorEvaluation::NumberOfOutputs()
{
  return 2;
}
std::string ErrorEvaluation::OutputName(int iOutput)
{
  return iOutput == 0 ? "error" : "confusion matrix";
}
std::string ErrorEvaluation::Output(int iOutput)
{
  ostringstream out;

  out.precision(14);
  if ( iOutput == 0 ) {
    out << error;
  }
  else if ( iOutput == 1 ) {
    for (int i = 0; i < matriz->filas(); i++) {
      for (int j = 0; j < matriz->columnas()-1; j++) {
        out << (*matriz)[i][j] << "\t";
      }
      out << (*matriz)[i][matriz->columnas()-1] << ";\t";
    }
  }
  else {
  }

  return out.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
MSEEvaluation::MSEEvaluation(Classifier *clf, Data *data, string params) : Evaluation(clf, data)
{
  Data *holdData = clf->GetData();

  clf->SetData(data);

  rms = 0.0;
  double ww = 0.0;
  int Nvar = data->GetNumVar();
  for (int i = 0 ; i < data->GetNTotal(); i++) {
    double d = clf->Average(i) - data->GetInstance(i)[Nvar-1];
    // printf("Average[%d] = %f, Nvar-1 = %f, d = %f\n", i, clf->Average(i), data->GetInstance(i)[Nvar-1], d);
    rms += data->GetInstance(i).GetWeight() * d * d;
    ww += data->GetInstance(i).GetWeight();
  }

  rms /= ww;

  clf->SetData(holdData);
}
int MSEEvaluation::NumberOfOutputs()
{
  return 1;
}
string MSEEvaluation::OutputName(int iOutput)
{
  return "MSE";
}
std::string MSEEvaluation::Output(int iOutput)
{
  ostringstream out;

  out.precision(14);
  if ( iOutput == 0 ) {
    out << rms;
  }
  else {
  }

  return out.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
MultipleMSEEvaluation::MultipleMSEEvaluation(Classifier *clf, Data *data, string params) : Evaluation(clf, data)
{
  Data *holdData = clf->GetData();
  int numberObjectives = MultiobjectiveInstance::GetNumMultiobjectives();

  clf->SetData(data);

  rms = 0.0;
  double ww = 0.0;
  for (int i = 0 ; i < data->GetNTotal(); i++) {
    std::vector<double> averages = clf->MultipleAverage(i);
    MultiobjectiveInstance& instance = dynamic_cast<MultiobjectiveInstance&> (data->GetInstance(i));
    for (int nObjective = 0; nObjective < numberObjectives; nObjective++) {
      double d = averages[nObjective] - instance.GetMultiobjective(nObjective);
      // printf("Average[%d][%d] = %f, obj = %f, d = %f\n", i, nObjective, averages[nObjective], instance.GetMultiobjective(nObjective), d);
      rms += instance.GetWeight() * d * d / numberObjectives;
      ww += instance.GetWeight() / numberObjectives;
    }
  }
  rms /= ww;

  clf->SetData(holdData);
}
int MultipleMSEEvaluation::NumberOfOutputs()
{
  return 1;
}
string MultipleMSEEvaluation::OutputName(int iOutput)
{
  return "MultipleMSE";
}
std::string MultipleMSEEvaluation::Output(int iOutput)
{
  ostringstream out;

  out.precision(14);
  if ( iOutput == 0 ) {
    out << rms;
  }
  else {
  }

  return out.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
ComboEvaluation::ComboEvaluation(Classifier *clf, Data *data, string parmas) : Evaluation(clf, data)
{
  Instances& instances = data->GetM();
  int NVar = data->GetNumVar();
  Ensemble *ens = (Ensemble*)clf;
  ens->SetData(data);

  string fn = GenerateNewFileName("locations_", ".txt");
  FILE *f = fopen(fn.c_str(),"w");
  for (int i = 0; i < ens->GetClassifiersToUse(); i++) {
    DecisionTree *dt = (DecisionTree*)ens->GetClassifier(i);
    for (int j = 0; j < data->GetNTotal(); j++) {
      dt->Classify(j);
      ComboSplitNodeInfo *csi = (ComboSplitNodeInfo*)dt->GetLastClassifyingNode()->info;
      Matriz &md = *csi->instances_data;
      if (dt->GetLastClassifyingNode()->iClass==0) continue; //Si el nodo es en mayoria background, nada
      //for (int k = 0; k < md.filas(); k++) {
      //}
      for (int k = 0; k < md.filas(); k++) {
        if (md[k][6]==0) continue;
  //    for (int k = 0; k < csi->count; k++) {
     //   fprintf(f, "%d %d\n", csi->indexes[k], instances[j].GetIniPos());
        //fprintf(f, "%g %g %g %g %g %g %d %d %d\n", md[k][0] * instances[j][NVar-9] +  instances[j][NVar-11], 
        /*double det_scale_x = instances[j][NVar-9];
        double det_scale_y = instances[j][NVar-8];
        double det_loc_x = instances[j][NVar-11];
        double det_loc_y = instances[j][NVar-10];*/
        fprintf(f, "%g %g %g %g %g %g %g %g %d\n", 
          /* Delta_x_evidence * scale_X_instance / scale_x_evidence + pos_x_instance */
          md[k][0], // * instances[j][NVar-9] / md[k][2] +  instances[j][NVar-11],
          /* Delta_y_evidence * scale_y_instance / scale_y_evidence + pos_y_instance */
          md[k][1], // * instances[j][NVar-8] / md[k][3] +  instances[j][NVar-10],
          /* scale_X_instance */
          instances[j][NVar-9],
          /* scale_y_instance */
          instances[j][NVar-8],
          /* Index sift evidence*/
          md[k][5],
          //md[k][2] * instances[j][NVar-9],
          /* Index sift instance */
          (double)instances[j].GetIniPos(),
          //md[k][3] * instances[j][NVar-8],
          /* pos_x_instance */
          instances[j][NVar-11],
          /* pos_y_instance */
          instances[j][NVar-10],
          /* class_evidence */
          (int)md[k][4]);
          //, (int)md[k][5],(int)md[k][6]);
      }
    }
  }

  fclose(f);
}
ComboEvaluation::~ComboEvaluation()
{
}
int ComboEvaluation::NumberOfOutputs()
{
  return 1;
}
std::string ComboEvaluation::OutputName(int iOutput)
{
  return "stats";
}
std::string ComboEvaluation::Output(int iOutput)
{
  ostringstream out;

  out.precision(14);

  return out.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
MultilabelHammingLossEvaluation::MultilabelHammingLossEvaluation(Classifier *clf, Data *data, 
                                                          string parmas) : Evaluation(clf, data)
{
  Instances& instances = data->GetM();
  Ensemble *ens = (Ensemble*)clf;
  ens->SetData(data);
  error = 0;

  for (int j = 0; j < data->GetNTotal(); j++) {
    vector<double> acum(MultilabelInstance::GetNumMultilabels(), 0.0);
    int ninst = 0;
    for (int i = 0; i < ens->GetClassifiersToUse(); i++) {
      DecisionTree *dt = (DecisionTree*)ens->GetClassifier(i);
      dt->Classify(j);
      Matriz *md = (Matriz*)dt->GetLastClassifyingNode()->info;
      for (int ik = 0; ik < md->filas(); ik++) {
        ninst++;
        for (int jk = 0; jk < MultilabelInstance::GetNumMultilabels(); jk++) {
          acum[jk] += (*md)[ik][jk];
        }
      }
    }
    for (int jk = 0; jk < MultilabelInstance::GetNumMultilabels(); jk++) {
      bool on = ((MultilabelInstance&)instances[j]).GetMultilabel(jk);
      bool predicion_on = acum[jk]/ninst > 0.5;
      if ( (~predicion_on && on) || (predicion_on && ~on) ) error++;
    }
  }
cout << error << endl;

  error /= MultilabelInstance::GetNumMultilabels()*data->GetNTotal();

}
MultilabelHammingLossEvaluation::~MultilabelHammingLossEvaluation()
{
}
int MultilabelHammingLossEvaluation::NumberOfOutputs()
{
  return 1;
}
std::string MultilabelHammingLossEvaluation::OutputName(int iOutput)
{
  return "Hamming loss";
}
std::string MultilabelHammingLossEvaluation::Output(int iOutput)
{
  ostringstream out;

  out << error;

  out.precision(14);

  return out.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
EvidenceEvaluation::EvidenceEvaluation(Classifier *clf, Data *data, string params) : Evaluation(clf, data)
{
  Instances& instances = data->GetM();
  int NClass =((NomData*) data)->NumClass;
  Ensemble *ens = (Ensemble*)clf;
  ens->SetData(data);
  vector<int> min_inst_per_leaf;

  if (params.empty()) {
    min_inst_per_leaf.push_back(0);
  }
  else {
    istringstream iss(params);
    while (!iss.eof()) {
      int val;
      iss >> val;
      if ( (iss.rdstate() & istringstream::failbit ) != 0 ) {
        cerr << "Error processing string: " << params << endl;
        exit(1);
      }
      min_inst_per_leaf.push_back(val);
    }
    sort(min_inst_per_leaf.begin(), min_inst_per_leaf.end());
  }

  vector< vector<double> > errors(min_inst_per_leaf.size()*2, 
                             vector<double>(ens->GetClassifiersToUse(), 0.0));

double reme=0.0;
double mm, Mm;

  for (int j = 0; j < data->GetNTotal(); j++) {
    vector< vector<double> > evidences(min_inst_per_leaf.size(), vector<double>(NClass, 0.0));
mm = 10000;
Mm = 0;
vector<double> vv(NClass, 0.0);
vector<double> Vv(NClass, 0.0);
    for (int i = 0; i < ens->GetClassifiersToUse(); i++) {
      double *v;
      DecisionTree *dt = (DecisionTree*)ens->GetClassifier(i);
      dt->Classify(j);
      Node *node = dt->GetLastClassifyingNode();
      for(unsigned int ie = 0; ie < min_inst_per_leaf.size(); ie++ ) {
        int node_evi = node->last - node->first + 1;
        int parent_evi = node->parent ? node->parent->last - node->parent->first + 1 : node_evi;
        //bool EnoughEvidence = node_evi >= min_inst_per_leaf[ie] && parent_evi-node_evi >= min_inst_per_leaf[ie];
        bool EnoughEvidence = parent_evi >= min_inst_per_leaf[ie];
        while (!EnoughEvidence && node->parent) {
          node = node->parent;
          node_evi = node->last - node->first + 1;
          parent_evi = node->parent ? node->parent->last - node->parent->first + 1 : node_evi;
          //EnoughEvidence = node_evi >= min_inst_per_leaf[ie] && parent_evi-node_evi >= min_inst_per_leaf[ie];
          EnoughEvidence = parent_evi >= min_inst_per_leaf[ie];
        }
        v = node->NodePop;
        for (int k = 0; k < NClass; k++) {
          evidences[ie][k] += node->NodePop[k];
        }
        errors[ie][i] += Classifier::WhichClass(evidences[ie]) == instances[j].GetClass() ? 
                                                   0 : instances[j].GetWeight();
        vector<double> we = evidences[ie]; 
        for (int k = 0; k < NClass; k++) {
          we[k] *= ((Bagging*)ens)->EvidenceWeights[k];
        }
        errors[min_inst_per_leaf.size()+ ie][i] += Classifier::WhichClass(we) == instances[j].GetClass() ? 
                                                   0 : instances[j].GetWeight();
      }
double mmm = 0.0;
for(int k = 0; k < NClass; k++) 
  mmm+=v[k];
reme += mmm;
if (mmm<mm) {mm=mmm;for(int k = 0; k < NClass; k++)vv[k]=v[k];} 
if (mmm>Mm) {Mm=mmm;for(int k = 0; k < NClass; k++)Vv[k]=v[k];}
    }
//printf("%g(", mm);
//for(int k = 0; k < NClass-1; k++)printf("%g,", vv[k]);
//printf("%g) %g(", vv[NClass-1],  Mm);
//for(int k = 0; k < NClass-1; k++)printf("%g,", Vv[k]);
//printf("%g)\n", Vv[NClass-1]);
  }

  for(unsigned int ie = 0; ie < errors.size(); ie++ ) {
    for (int i = 0; i < ens->GetClassifiersToUse(); i++) {
      errors[ie][i] /= data->GetNTotal();
    }
  }

  reme /= (ens->GetClassifiersToUse()*data->GetNTotal());
  errors.push_back(vector<double>());
  (*(errors.end()-1)).push_back(reme);
  (*(errors.end()-1)).push_back(mm);
  (*(errors.end()-1)).push_back(Mm);

  seq_errors = errors;
}
EvidenceEvaluation::~EvidenceEvaluation()
{
}
int EvidenceEvaluation::NumberOfOutputs()
{
  return 1;
}
std::string EvidenceEvaluation::OutputName(int iOutput)
{
  return "Evidence";
}
std::string EvidenceEvaluation::Output(int iOutput)
{
  ostringstream out;

  out.precision(14);
  for(unsigned int i = 0; i < seq_errors.size(); i++) {
    for(unsigned int j = 0; j < seq_errors[i].size(); j++) {
      out << seq_errors[i][j] << " ";
    }
    if ( i < seq_errors.size() - 1 ) out << endl;
  }

  return out.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//----------------i-----------------------------------------------------------
StackingDataEvaluation::StackingDataEvaluation(Classifier *clf, Data *data, 
                                         string params) : Evaluation(clf, data)
{
  FILE *fe, *fv;

  int nclases = ((NomData*)data)->NumClass;
  Ensemble *ens = (Ensemble*)clf;
  ens->SetData(data);

  MINomData *midata = dynamic_cast<MINomData*>(data);

  string set = params.empty() ? "train" : params;
  bool train = set.compare("train")==0 ? true : false;

  string n = GenerateNewFileName("evidence_" + set + "_",".cre");
  fe = fopen(n.c_str(), "w");
  n = GenerateNewFileName("vote_" + set + "_",".cre");
  fv = fopen(n.c_str(), "w");

  base_file_name = n.substr(7);

  if (midata) {
    fprintf(fe, "%d\n%d\n", midata->GetMITotal(), nclases+1);
    fprintf(fv, "%d\n%d\n", midata->GetMITotal(), nclases+1);
  }
  else {
    fprintf(fe, "%d\n%d\n", data->GetNTotal(), nclases+1);
    fprintf(fv, "%d\n%d\n", data->GetNTotal(), nclases+1);
  }

  for(int k = 0; k < nclases; k++) {
    ostringstream out;
    out << "Att" << k;
    string aux = out.str();
    fprintf(fe, "%s 0\n", aux.c_str());
    fprintf(fv, "%s 0\n", aux.c_str());
  }
  fprintf(fe, "Class 2\n");
  fprintf(fv, "Class 2\n");

  data->MarkOrder();
  data->OriginalOrder();

  for (int i = 0; i < data->GetNTotal(); i++) { // loop through single instance
    vector<double> evidence(nclases, 0.0);
    vector<double> votes(nclases, 0.0);
    int n_classifiers_voting;
    do {
      n_classifiers_voting = 0;
      for (int j = 0; j < ens->GetClassifiersToUse(); j++) { //For each classifier in the ensemble
        int iori = data->GetDatIniPos(i);
        if (train && ens->GetClassifier(j)->UsedInOriginalTrainingData(iori)) continue;

        vector<double> d = ens->GetClassifier(j)->UnnormalizedDistribution(i);
        int clase = ens->GetClassifier(j)->Classify(i);

        n_classifiers_voting++;

        votes[clase] += 1.0;
        for(int k = 0; k < nclases; k++) evidence[k] += d[k];
      }
     
      i++;
    //Loop to group multiple instances output
    } while (midata && i<data->GetNTotal() && data->GetDatGroup(i-1)==data->GetDatGroup(i));  
    i--;

    double norm_cons = 0.0;
    for(int k = 0; k < nclases; k++) {
       norm_cons += evidence[k];
    }
    for(int k = 0; k < nclases; k++) {
      fprintf(fe, "%.8g ", evidence[k]/norm_cons);
      fprintf(fv, "%.8g ", votes[k]/n_classifiers_voting);
    }
    fprintf(fe, "%d\n", data->GetDatClass(i));
    fprintf(fv, "%d\n", data->GetDatClass(i));
  }

  data->ResetOrder();

  fclose(fe);
  fclose(fv);
}
StackingDataEvaluation::~StackingDataEvaluation()
{
}
int StackingDataEvaluation::NumberOfOutputs()
{
  return 1;
}
std::string StackingDataEvaluation::OutputName(int iOutput)
{
  return "Evidence";
}
std::string StackingDataEvaluation::Output(int iOutput)
{
  return base_file_name;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//----------------i-----------------------------------------------------------

