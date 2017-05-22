/**
 * Fichero de cabecera del objeto que implementa un clasificador basado en redes neuronales
 * 
 * Autor:  Daniel Hernández Lobato, Aitor Sanchez Martinez
 * Fecha: 4/06/2005
 */

#include "Classifier.h"
#include "data20.h"
#include <doublefann.h>
#include <vector>
#include <map>



//variables a modificar para el estudio

#define MAX_EPOCAS 1000


// Heredaremos de Classifier e implementaremos los metodos necesarios para definir el calsificador

class NNet : public Classifier, DataEvents {

  // Atributos de la clase NNet
 protected: 

    static std::map<Data *, fann_train_data*> rep_datos;

    virtual void OnDelete(Data *);


 private:
  
 
  
  double a_momento;          // Momento utilizado en el aprendizaje on line

  // Entero que indica el numero de veces seguidas que debe
  // incrementarse el error de validacion antes de para el entrenamiento.
  int a_reglaValidacion;

  // Indica el tamaño del conjunto de validacion en tanto por 1.
  double a_tamanioValidacion; 

  // Indica el tamaño del conjunto de validacion en tanto por 1.
  double a_coefAprendizaje; 


  // Vector de medias usado en la normalizacion.
  std::vector<double> a_medias;
  // Vector de desviaciones usado en la normalizacion.
  std::vector<double> a_desviaciones;

  struct fann_train_data *a_patrones; // Datos con los que trabajar.
  struct fann *a_nn;                  // Red Neruonal.


  unsigned int a_neuronas_capa_oculta_ini;    //numero de neuronas en capa oculta
  unsigned int a_numeroEpocas;     // Número de epocas que entrenaremos la red
  unsigned int a_neuronas_capa_oculta_fijo;    //numero fijo de neuronas inicial (0 si hay que buscarlo)
  unsigned int a_max_neuronas;     // Número maximo de neuronas que entrenaremos la red


 public:

  unsigned a_neuronas_capa_oculta_fin;    //numero final de neuronas ocultas
 

  // Metodos privados de la clase NNet
  
 private:
  
  void normalizarDatos(struct fann_train_data *data);
  void calcularVectoresNormalizacion(struct fann_train_data *data);

  struct fann_train_data * transformarDatosFANN(Data *datos, bool training = false);
  struct fann_train_data * extraerSubConjuntoFANN(struct fann_train_data *data, int inicio, int fin);


 /**
   *  Metodo que calcula el maximo de neuronas permitido para la capa oculta de la red
   *  Se calcula comoLA raiz cuadrada del numero que se le pasa como parametro (el numero de datos de train) 
   *  Return : el numero potencia de 2 inmediatamente superior a la la raiz cuadrada del numero pasado
   */

  unsigned int sacarMaximoNeuronas(unsigned int num_data);


  /**
   *  Metodo de calculo del error actual de una red neuronal
   *  struct fann *ann: red actual de la red
   *  struct fann_train_data *data. datos de entrenamiento actuales
   *  Return : error actual de la red.
   */

  unsigned int checkError(struct fann *ann,struct fann_train_data *train);


  /*  Funcion que imita la funcion de la libreria fann "fann_train_on_data" con la excepcion de que este metodo
   *  utiliza el error estandar sobre la salida esperada en lugar del MSE 
   *  
   *  struct fann *ann: red neuronal sobre la que trabajar
   *  struct fann_train_data *data datos de entrenamiento
   *  unsigned int max_epochs: maximo numero de epocas que entrenar
   *  unsigned int epochs_between_check_error: cada x numero de epocas se chequea el error
   *  desired_error: error deseado
 
   *  Return: el error obtenido )  
   */
  unsigned int own_fann_train_on_data(struct fann *ann, struct fann_train_data *train, unsigned int max_epochs, unsigned int epochs_between_check_error, unsigned int desired_error);


  // Metodos publicos de la clase NNet
  
 public:

  NNet(unsigned int num_neuronas=0,unsigned int num_epocas=MAX_EPOCAS);

  virtual ~NNet();

  virtual void Build(Data *data, FuncionDeProgreso *fp = 0);
  //  virtual double Error(int first, int last);
  virtual void SetData(Data *datos);
  //    virtual int Classify(int index);
  virtual std::vector<double> UnnormalizedDistribution(int ElementIndex);

};

