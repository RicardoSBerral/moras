#include "Language.h"
#include <stdio.h>

StringRepository* StringRepository::sr = 0;

StringRepository* StringRepository::GetStringRepository()
{
  if (sr==0) {
#ifdef IN_ENGLISH
    sr = new StringRepositoryEnglish();
#else
    sr = new StringRepositorySpanish();
#endif
 } 

  return sr;
}

std::string StringRepository::GetString(std::string Key)
{
  std::string value = GetStringRepository()->strings[Key];

  if (value.length()==0) {
    printf("Key %s with no associated value\n", Key.c_str());
  }

  return value;
}

StringRepositoryEnglish::StringRepositoryEnglish():StringRepository()
{
   strings["LoadClassifier"] 	= "LoadClassifier";
   strings["LoadDataset"] 	= "LoadDataset";
   strings["ChangeDependentColumn"] 	= "ChangeDependentColumn";
   strings["CopyDependentColumn"] 	= "CopyDependentColumn";
   strings["GenerateWrongLabels"] 	= "GenerateWrongLabels";
   strings["Corrected"] 	= "Corrected";
   strings["LoadLabels"] 	= "LoadLabels";
   strings["FilterData"] 	= "FilterData";
   strings["SetSplitCriterium"] 	= "SetSplitCriterium";
   strings["Classify"] 	= "Classify";
   strings["Error"] 	= "Error";
   strings["ResubstituteNodeStats"] 	= "ResubstituteNodeStats";
   strings["EvaluateClassifier"] 	= "EvaluateClassifier";
   strings["Mrse"] 	= "Mrse";
   strings["SequentialError"] 	= "SequentialError";
   strings["ClassificationMatrix"] 	= "ClassificationMatrix";
   strings["ConfusionMatrix"] 	= "ConfusionMatrix";
   strings["Margin"] 	= "Margin";
   strings["MarginOfInstances"] 	= "MarginOfInstances";
   strings["Certainty"] 	= "Certainty";
   strings["MarginSum"] 	= "MarginSum";
   strings["EstimErrorVal"] 	= "EstimErrorVal";
   strings["EstimErrorTest"] 	= "EstimErrorTest";
   strings["DiversityMeasures"] 	= "DiversityMeasures";
   strings["EnsembleMeasures"] 	= "EnsembleMeasures";
   strings["SetPropertyValue"] 	= "SetPropertyValue";
   strings["SetNTrain"] 	= "SetNTrain";
   strings["GeneratePartitions"] 	= "GeneratePartitions";
   strings["SaveDataset"] 	= "SaveDataset";
   strings["SaveClassifier"] 	= "SaveClassifier";
   strings["ClassifierInfo"] 	= "ClassifierInfo";
   strings["BuildIGPEnsemble"] 	= "BuildIGPEnsemble";
   strings["BuildBaggingCART"] 	= "BuildBaggingCART";
   strings["BuildBaggingC45"] 	= "BuildBaggingC45";
   strings["BuildBaggingNNet"] 	= "BuildBaggingNNet";
   strings["BuildRandomForest"] 	= "BuildRandomForest";
   strings["BuildBaggingComMod"] 	= "BuildBaggingComMod";
   strings["BuildBoosting"] 	= "BuildBoosting";
   strings["BuildBoostingC45"] 	= "BuildBoostingC45";
   strings["BuildClassSwitching"] 	= "BuildClassSwitching";
   strings["BuildBoostingWithBaggingInfo"] 	= "BuildBoostingWithBaggingInfo";
   strings["BuildCART"] 	= "BuildCART";
   strings["BuildC45"] 	= "BuildC45";
   strings["BuildNNet"] 	= "BuildNNet";
   strings["EnsembleReport"] 	= "EnsembleReport";
   strings["ClassificationMap2D"] 	= "ClassificationMap2D";
   strings["TestAll"] 	= "TestAll";
   strings["Clear"] 	= "Clear";
}

StringRepositorySpanish::StringRepositorySpanish():StringRepository()
{
   strings["LoadClassifier"] 	= "CargaClasificador";
   strings["LoadDataset"] 	= "CargaDatos";
   strings["ChangeDependentColumn"] 	= "CambiaColumnaDependiente";
   strings["CopyDependentColumn"] 	= "CopiaColumnaDependiente";
   strings["GenerateWrongLabels"] 	= "GeneraEtiquetasErroneas";
   strings["Corrected"] 	= "Corregidos";
   strings["LoadLabels"] 	= "CargaEtiquetas";
   strings["FilterData"] 	= "FiltraDatos";
   strings["SetSplitCriterium"] 	= "SetSplitCriterium";
   strings["Classify"] 	= "Clasificar";
   strings["Error"] 	= "ClasErr";
   strings["ResubstituteNodeStats"] 	= "ResubstituyeStatsNodos";
   strings["EvaluateClassifier"] 	= "EvaluaClasificador";
   strings["Mrse"] 	= "Mrse";
   strings["SequentialError"] 	= "ErrSecClas";
   strings["ClassificationMatrix"] 	= "MatrizClasif";
   strings["ConfusionMatrix"] 	= "MatrizConfusion";
   strings["Margin"] 	= "Margen";
   strings["MarginOfInstances"] 	= "MargenEjemplos";
   strings["Certainty"] 	= "Certeza";
   strings["MarginSum"] 	= "SumaMargen";
   strings["EstimErrorVal"] 	= "EstimErrorVal";
   strings["EstimErrorTest"] 	= "EstimErrorTest";
   strings["DiversityMeasures"] 	= "DiversityMeasures";
   strings["EnsembleMeasures"] 	= "EnsembleMeasures";
   strings["SetPropertyValue"] 	= "SetPropertyValue";
   strings["SetNTrain"] 	= "SetNTrain";
   strings["GeneratePartitions"] 	= "CreaConjuntosTrainAleat";
   strings["SaveDataset"] 	= "GuardaDatos";
   strings["SaveClassifier"] 	= "GuardaClasificador";
   strings["ClassifierInfo"] 	= "InfoClasificador";
   strings["BuildIGPEnsemble"] 	= "ConstruyeConjuntoIGP";
   strings["BuildBaggingCART"] 	= "ConstruyeBaggingCART";
   strings["BuildBaggingC45"] 	= "ConstruyeBaggingC45";
   strings["BuildBaggingNNet"] 	= "ConstruyeBaggingNNet";
   strings["BuildRandomForest"] 	= "ConstruyeRandomForest";
   strings["BuildBaggingComMod"] 	= "ConstruyeBaggingComMod";
   strings["BuildBoosting"] 	= "ConstruyeBoosting";
   strings["BuildBoostingC45"] 	= "ConstruyeBoostingC45";
   strings["BuildClassSwitching"] 	= "ConstruyeClassSwitching";
   strings["BuildBoostingWithBaggingInfo"] 	= "ConstruyeBoostingWithBaggingInfo";
   strings["BuildCART"] 	= "ConstruyeCART";
   strings["BuildC45"] 	= "ConstruyeC45";
   strings["BuildNNet"] 	= "ConstruyeNNet";
   strings["EnsembleReport"] 	= "EnsembleReport";
   strings["ClassificationMap2D"] 	= "MapaDeClasif2D";
   strings["TestAll"] 	= "TestAll";
   strings["Clear"] 	= "Clear";
}
