/**
 * Fichero que contiene la implementacion de la clase NNet
 *
 * Autor: Daniel Hernandez Lobato, Aitor Sanchez Martinez
 * Fecha 04/06/2005
 */


#include "NNet.h"
#include <map>
#include <iostream>
#include <fann_internal.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>


//EL DEFINE PARA EL NUMERO MAXIMO DE EPOCAS SE ENCUENTRA EN NNET.H
//#define MAX_EPOCAS 1000

//variables a modificar para el estudio
#define EPOCAS_ENTRE_CHEQUEOS 250
#define NUM_NEURONAS_CAPA_OCULTA_INI  2
#define COEF_APRENDIZAJE 0.7 
#define NUM_CAPAS_RED 3
#define ERROR_DESEADO 0
#define FULLY_CONECTED 1


//CON ESTA VARIABLE SE ACTIVAN O DESACTIVAN LOS COMENTARIOS POR PANTALLA 
#define DEBUG  0


//CON ESTA VARIABLE SE ACTIVA (1) O DESACTIVA (0) EL REORDENAMIENTO DE DATOS PARA EL ENTRENAMIENTO 
#define SHUFFLE  1


/*se define eta variable
 * 0 si se quiere como maximo de neuronas la potencia de 2 inmediatamente superior o igual a srt(num_data)
 * Otro numero (potencia de 2 para optimizacion del codigo)
 */
#define MAX_NEURONAS 0


std::map<Data *, fann_train_data*> NNet::rep_datos;


using namespace std;

/**
 * Constructor de la clase NNet
 * unsigned int num_neuronas: numero de neuronas en la capa oculta con que debe crear la red
 */

NNet::NNet(unsigned int num_neuronas,unsigned int num_epocas) {


  data = NULL;
  a_patrones = NULL;
  a_nn = NULL;

  a_numeroEpocas = num_epocas;

  // No usado en esta implementacion aunque la red se crea con ella no se usa aqui, la red lo controla
  a_coefAprendizaje=COEF_APRENDIZAJE;   
  

  a_max_neuronas=MAX_NEURONAS;
  
    

  //sirven para ir modificando las neuronas de la capa oculta
  //mediante busqueda binaria (2(x)-4(x)-8(x)-16(v)-12(x)-14(v)) o (64(v)-32(v)-16(x)-24(v)...)
  a_neuronas_capa_oculta_ini=NUM_NEURONAS_CAPA_OCULTA_INI;
  

  //sirve para saber el numero de neuronas fijo (0 si hay que calcularlo) 
  a_neuronas_capa_oculta_fijo=num_neuronas;

  //sirve para guardar el numero optimo de neuronas ocultas
  a_neuronas_capa_oculta_fin=1;

}



/**
 * Destructor de la clase NNet
 */

NNet::~NNet() {

//  if (a_patrones != NULL) {
//    fann_destroy_train(a_patrones);
//  }

  if (a_nn != NULL)
    fann_destroy(a_nn);

  a_medias.clear();
  a_desviaciones.clear();
}



/**
 *  Metodo que entrena una red neuronal para clasificar unos datos
 *
 *  data : datos con los que construir la red.
 */

void NNet::Build(Data *datos, FuncionDeProgreso *fp)
{

  /*NOT USED
    int nIncrementosValidacion = a_reglaValidacion;
    double coefAprendizaje = a_coefAprendizaje, mseTrainning, mseValidacion;
    double mseTrainningAnterior = 10000 , mseValidacionAnterior = 10000;
    struct fann_train_data *validacion;
  */


  struct fann_train_data *entrenamiento;
  struct fann_train_data *datos_fann;


  unsigned int seguir=1;
  

  // variable para el calculo del error
  unsigned int error=0;



  //variables a utilizar para la busqueda binaria del numero minimo de neuronas
  unsigned neuronas_capa_oculta=a_neuronas_capa_oculta_ini;
  unsigned neuronas_capa_oculta_previo=0;
  unsigned minimo_neuronas=1; 

  time_t timer;


  if (datos == NULL) {
    std::cerr << "Error intento de crear clasificador NNet sin datos" << std::endl;
    return;
  }

  Init(datos);
  
  // Transformamos los datos

  datos_fann = transformarDatosFANN(datos, true);

  // Primero Calculamos los vectores de Normalizacion
  
  calcularVectoresNormalizacion(datos_fann);

  // se normalizan los datos
  normalizarDatos(datos_fann);


  // se llama al metodo que guarda los datos (que transforma y normaliza) 
  // SE UTILIZARAN LOS VECTOREES DE NORMALIZACION CALCULADOS CON DATOS_FANN -> ES CORRECTO?? ->>>> si
  
  SetData(datos); 



  // Obtenemos el conjunto de entrenamiento (en nuestro caso todo el conjunto de datos)
  
  if (SHUFFLE==1) fann_shuffle_train_data(datos_fann);
  

  entrenamiento = extraerSubConjuntoFANN(datos_fann, 0, datos_fann->num_data - 1);
 
  

  //ALGORITMO SEGUIDO PARA EL APRENDIZAJE

  //se fija el numero de neuronas de capa oculta inicial
  neuronas_capa_oculta=a_neuronas_capa_oculta_ini;
  neuronas_capa_oculta_previo=0;
  minimo_neuronas=1;
  unsigned neuronas_previo_temp=0;
  unsigned error_min=0;
  struct fann *a_nn_aux=NULL;


  //se determina el numero maximo de neuronas para la red dependiendo de la variable MAX_NEURONAS
  if (MAX_NEURONAS!=0)
    a_max_neuronas=MAX_NEURONAS;
  else
   a_max_neuronas=sacarMaximoNeuronas(entrenamiento->num_data);
   
  if (DEBUG) cout<<"Numero maximo de neuronas: "<<a_max_neuronas<<" (numero datos = "<<entrenamiento->num_data<<")\n";


  //si no hay numero de neuronas definido se crea por busqueda de la solucion optima
  if (a_neuronas_capa_oculta_fijo==0){


    //se mide el tiempo de origen
    timer=time(NULL);


    //bucle hasta le maximo de nuronas para comprobar (numero razonable~100), 
    //o hasta que no se quiera seguir (min encontrado)
    while (neuronas_capa_oculta<=a_max_neuronas &&seguir==1){

    

      // Construimos la red


      a_nn = fann_create(FULLY_CONECTED,COEF_APRENDIZAJE, NUM_CAPAS_RED,
			 entrenamiento->num_input ,neuronas_capa_oculta, //NUMERO DE NEURONAS EN CAPA OCULTA (PARAM) 
			 entrenamiento->num_output);


      if (DEBUG) cout<<"Numero actual de neuronas: "<<neuronas_capa_oculta<<"\n";

  
      //lo primero que hay que hacer es fijar el metodo de aprendizaje de la red
      //se utilizara FANN_TRAIN_RPROP, que aunque es el metodo por defecto, se especifica de todas formas
      //para tenerlo presente. Este metodo es adaptativo, luego no utiliza el coeficiente de aprendizaje o "learning rate"

      fann_set_training_algorithm(a_nn,FANN_TRAIN_RPROP);


      //se entrena la red el numero de epocas que le digamos, teniendo en cuenta que se chequea segun el intervalo deseado.
      error=own_fann_train_on_data(a_nn,entrenamiento,a_numeroEpocas, EPOCAS_ENTRE_CHEQUEOS, ERROR_DESEADO);

      //la funcion callback pone el error global del entrenameinto de la neurona, asi que en este punto podemos
      //trabajar con dicho error

      neuronas_previo_temp=neuronas_capa_oculta;

      //se  mira si se tiene que actualizar las neuronas de la capa oculta
    
      if (error>ERROR_DESEADO){
    
	//si se llega a numero de nuronas impar el bueno es el guardado par anterior
	if (neuronas_capa_oculta%2==1) {
	  seguir=0;
	  error=error_min;
	}else{
      
	  //se mira a ver si se tiene que subir o no al doble de neuronas en funcion de lo anterior
	  if (neuronas_capa_oculta_previo>neuronas_capa_oculta){
	    neuronas_capa_oculta=(neuronas_capa_oculta+neuronas_capa_oculta_previo)/2;
      
	  }else{
	    if (a_nn_aux==NULL) 
	      neuronas_capa_oculta= neuronas_capa_oculta*2;
	    else
	      neuronas_capa_oculta=(neuronas_capa_oculta+minimo_neuronas)/2;
	  }
	}

	//en este punto se puede destruir la red, puesto que sabemos que no es ganadora (salvo que se acabe el bucle)
	if (!(a_nn_aux==NULL && neuronas_capa_oculta>a_max_neuronas)) {
	  fann_destroy(a_nn);
	}else{
	  a_nn_aux=a_nn;
	  minimo_neuronas=neuronas_previo_temp;
	  seguir=0;
	}

      }else if (error<=ERROR_DESEADO){

 
	minimo_neuronas=neuronas_capa_oculta;
      
	//se guarda el error minimo
	error_min=error;
      
	//guarda la red neuronal que por ahora va ganando      
	if (a_nn_aux!=NULL) fann_destroy(a_nn_aux);
	a_nn_aux=a_nn; 
 
	//si hay error cero se mira si es el minimo alcanzado (si es un numero impar estamos seguros de que es el minimo)
	if (neuronas_capa_oculta%2==1) seguir=0;
	else{
      
	  if (neuronas_capa_oculta_previo>neuronas_capa_oculta){
	    neuronas_capa_oculta=neuronas_capa_oculta-((neuronas_capa_oculta_previo-neuronas_capa_oculta)/2);
	  }else{  
	
	    //se situa el numero de neuronas en la mitad inferior del intervalo entre el numero anterior y el actual
	    neuronas_capa_oculta=(neuronas_capa_oculta+neuronas_capa_oculta_previo)/2;
	  }
	}

      }

    
      
      //actualiza el numero de neuronas previo (sera el previo de la siguiente iteracion)
      neuronas_capa_oculta_previo=neuronas_previo_temp;

    }

    //se guarda la red final
    a_nn=a_nn_aux;

    //se actualiza la variable para que sea accesible desde el exterior
    a_neuronas_capa_oculta_fin=minimo_neuronas;


  }else{
    

    //se mide el tiempo de origen
    timer=time(NULL);

    //hay ya numero de neuronas definido!!, se utiliza ese
//cout << "INPUTS=" << entrenamiento->num_input << endl;
    a_nn = fann_create(FULLY_CONECTED,COEF_APRENDIZAJE, NUM_CAPAS_RED,
			 entrenamiento->num_input ,a_neuronas_capa_oculta_fijo, //NUMERO DE NEURONAS EN CAPA OCULTA 
			 entrenamiento->num_output);
    

    //se actualiza la variable para que sea accesible desde el exterior
    a_neuronas_capa_oculta_fin=a_neuronas_capa_oculta_fijo;


    if (DEBUG) cout<<"Numero de neuronas fijo: "<<a_neuronas_capa_oculta_fin<<"\n";
    
    fann_set_training_algorithm(a_nn,FANN_TRAIN_RPROP);

    //se entrena la red el numero de epocas que le digamos, teniendo en cuenta que se chequea segun el intervalo desead
    error=own_fann_train_on_data(a_nn,entrenamiento,a_numeroEpocas, 0, 0);

  }



  if (DEBUG) {
    cout<<"Tiempo: "<< difftime(time(NULL),timer)<<"sec\tNeuronas: "<<a_neuronas_capa_oculta_fin<<"\tError: "<<error<<"\n";
    //cout<<"Datos de la red resultante del entrenamiento:\n\n";
    //fann_print_parameters(a_nn);
 
  }
 

  fann_destroy_train(entrenamiento);
  fann_destroy_train(datos_fann);
}



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
unsigned int NNet::own_fann_train_on_data(struct fann *ann, struct fann_train_data *train, unsigned int max_epochs, unsigned int epochs_between_check_error, unsigned int desired_error)
{
  unsigned int error=100000;
  unsigned int i;

	
  if(epochs_between_check_error && epochs_between_check_error < 1){
    cout<<"Epochs_between_reports < 1!!!";
    return error;
  }
 

  /* some training algorithms need stuff to be cleared etc. before training starts.
   */
  
  /* NO HACE FALTA PORQUE ES UNA RED NUEVA CADA VEZ QUE SE LLAMA AQUI..Y SE LLAMA INTERNAMENTE AL RESETEO DE ARRAYS
    if(ann->training_algorithm == FANN_TRAIN_RPROP ||
     ann->training_algorithm == FANN_TRAIN_QUICKPROP){
  }*/
   
  for(i = 1; i <= max_epochs; i++){
		
    /* train */
    fann_train_epoch(ann, train);
	  
		
    /* print current output */
    if((epochs_between_check_error && (i % epochs_between_check_error == 0)) || (i == max_epochs) ){

      //se coge el error habitual standar!!!
      error=checkError(ann,train);
      
      if (DEBUG) cout<<"\tEpocas actuales: "<<i<<" Error actual: "<<error<<"\n";
	    
      if (desired_error && error <= desired_error) return error;
    } 

	 
  }
  
  //si se ha llegado hasta aqui no se ha conseguido error cero y se ha llegado al numero de epocas
  //se deberia tener el error calculado arriba por llegar
  return error;
}



/**
 *  Metodo de calculo del error actual de una red neuronal
 *  struct fann *ann: red actual de la red
 *  struct fann_train_data *data. datos de entrenamiento actuales
 *  Return : error actual de la red.
 */

unsigned int NNet::checkError(struct fann *ann,struct fann_train_data *train){
    
  // la variable para almacenar datos resultantes de clasificacion, es decir, el porcentaje 
  // de acierto de un patron sobre una clase
  fann_type *result_class;

  fann_type aux=-1.0,max= -1.0; //sirve para ir almacenando el maximo obtenido de porcentaje a una clase  
  fann_type o_aux=-1.0,o_max= -1.0; //sirve para ir almacenando el maximo obtenido de porcentaje a una clase (salidas)  
    
  unsigned class_winner=0; //sirve para guardar la clase ganadora
  unsigned o_class_winner=0; //sirve para guardar la clase ganadora (salidas esperadas)

  unsigned int j=0,k=0;
  unsigned int error=0;


  for (j=0;j<train->num_data;j++){
    
    result_class = fann_run(a_nn, train->input[j]);

    k=0;
    aux=-1.0,max= -1.0;
    o_aux=-1.0,o_max= -1.0;
    class_winner=0;
    o_class_winner=0;

    for(k=0;k<a_nn->num_output;k++){
	 
      aux=result_class[k];
      o_aux= train->output[j][k];


      //si hay dos clases con el mismo porcentaje de acierto, se queda la primera, de ahi que quitemos el = en
      //la asignacion

      if (aux>max){
        max=aux;
        class_winner=k; //se guarda quien es la neurona de salida ganadora
 
      }
      if (o_aux>o_max){
        o_max=o_aux;
        o_class_winner=k; //se guarda quien es la neurona de salida ganadora (de las esperadas)
 
      }
    }

    //ahora se comprueba si la neurona ganadora es la misma que con el conjunto de entrenamiento
    if (o_class_winner!=class_winner)
      error++;
  }

  return error;
}



/**
 *  Metodo que calcula el maximo de neuronas permitido para la capa oculta de la red
 *  Se calcula comoLA raiz cuadrada del numero que se le pasa como parametro (el numero de datos de train) 
 *  Return : el numero potencia de 2 inmediatamente superior a la la raiz cuadrada del numero pasado
 */

unsigned int NNet::sacarMaximoNeuronas(unsigned int num_data){
  
  double pro=1;
  double sqr=sqrt(num_data);
  
  while(1){
    pro*=2;
    if (pro>=sqr) return (unsigned)pro;
  }


}



/**
 *  Metodo que calcula los vectores de normalizacion de los datos
 *
 *  data : datos con los que construir la red.
 */

void NNet::calcularVectoresNormalizacion(struct fann_train_data *datos)
{
  a_medias.clear();
  a_desviaciones.clear();

  // Primero normalizamos los datos
  
  for (unsigned i = 0 ; i < datos->num_input ; i++) {
    double media = 0;

    for (unsigned j = 0 ; j < datos->num_data ; j++) 
      media += datos->input[ j ][ i ];

    a_medias.push_back((double) media / datos->num_data);
  }
  
  for (unsigned i = 0 ; i < datos->num_input ; i++) {
    double desv = 0;

    for (unsigned j = 0 ; j < datos->num_data ; j++) 
      desv += (datos->input[ j ][ i ] - a_medias[ i ]) * (datos->input[ j ][ i ] - a_medias[ i ]);

    a_desviaciones.push_back(sqrt((double) desv / datos->num_data));
  }
}

/**
 *  Metodo que normaliza los datos
 *
 *  data : datos con los que construir la red.
 */

void NNet::normalizarDatos(struct fann_train_data *datos)
{
  for (unsigned j = 0 ; j < datos->num_data ; j++) 
    for (unsigned i = 0 ; i < datos->num_input ; i++) 
      datos->input[ j ][ i ] = (double) (datos->input[ j ][ i ] - a_medias[ i ]) / a_desviaciones[ i ];
}


/**
 *  Metodo que fija los nuevos datos de la red 
 *
 *  data : datos con los que calcular el error de clasificacion.
 */

void NNet::SetData(Data *datos) {

  Classifier::SetData(datos);

  if (datos != NULL) {
    
//    if (a_patrones != NULL)
//      fann_destroy_train(a_patrones);

    a_patrones = NNet::rep_datos[datos];

    if (!a_patrones) {
      this->a_patrones = transformarDatosFANN(datos);
      normalizarDatos(this->a_patrones);
      NNet::rep_datos[data] = a_patrones;
      datos->AddListener(this);
    }
  }

  //this->data = datos;
}

//---------------------------------------------------------------------------
void NNet::OnDelete(Data *data)
{
  fann_train_data *ftd = NNet::rep_datos[data];
  if (ftd) {
    if (ftd==a_patrones) a_patrones = 0;
    fann_destroy_train(ftd);
    NNet::rep_datos[data] = 0;
  } 
}


/**
 *  Metodo que transforma los datos al tipo fann 
 *
 *  data : datos con los que calcular el error de clasificacion.
 */

struct fann_train_data * NNet::transformarDatosFANN(Data *datos, bool training) 
{
  unsigned int num_input, num_output, num_data, i, j, num_input_sin_expandir;
  unsigned int line = 1;
  fann_type *data_input, *data_output;
  struct fann_train_data* data = (struct fann_train_data *)malloc(sizeof(struct fann_train_data));

  if(data == NULL){
    fprintf(stderr, "Error al pedir memoria para datos fann.\n");
    return NULL;
  }

  if (training == true)
	  num_data = datos->GetNTrain();
  else
	  num_data = datos->GetNTotal();

  // Calculamos el numero de entradas

  num_input = datos->GetNumVarOrd() + datos->GetNumVarFuz();

  for (int u = 0 ; u < datos->GetNumVarNom(); u++)
	num_input += (datos->GetTermVector( u )).size();

 
  num_input_sin_expandir = datos->GetNumVar() - 1;

  num_output =  ((NomData*)datos)->NumClass;
  
  
  struct fann_error *errdat = (struct fann_error *)data;
  errdat->errstr = NULL;
  errdat->errno_f = 0;
  errdat->error_log = stderr;


  data->num_data = num_data;
  data->num_input = num_input;
  data->num_output = num_output;
  data->input = (fann_type **)calloc(num_data, sizeof(fann_type *));
  if(data->input == NULL){
    fprintf(stderr, "Error al pedir memoria para datos fann.\n");
    fann_destroy_train(data);
    return NULL;
  }
  
  data->output = (fann_type **)calloc(num_data, sizeof(fann_type *));
  if(data->output == NULL){
    fprintf(stderr, "Error al pedir memoria para datos fann.\n");
    fann_destroy_train(data);
    return NULL;
  }
  
  data_input = (fann_type *)calloc(num_input*num_data, sizeof(fann_type));
  if(data_input == NULL){
    fprintf(stderr, "Error al pedir memoria para datos fann.\n");
    fann_destroy_train(data);
    return NULL;
  }

  data_output = (fann_type *)calloc(num_output*num_data, sizeof(fann_type));
  if(data_output == NULL){
    fprintf(stderr, "Error al pedir memoria para datos fann.\n");
    fann_destroy_train(data);
    return NULL;
  }
  
  for(i = 0; i != num_data; i++){
    data->input[i] = data_input;
    data_input += num_input;
    
    for(j = 0; j != num_input; j++)
      data->input[ i ][ j ] = 0.0;

    for(unsigned u = 0, v = 0, j = 0; j != num_input_sin_expandir; j++) {
	if (j < (unsigned)datos->GetNumVarOrd()){
	      data->input[ i ][ u ] = datos->GetValueVar(i, j)  ==  1.7976931348623157e+308 ? 
		      0 : datos->GetValueVar(i, j);
	      u++;
	} else {
		for (unsigned k = 0 ; k < (datos->GetTermVector( v )).size() ; k++, u++) {
			if ( ((unsigned)datos->GetValueVar(i, j)) == k)
				data->input[ i ][ u ] = 1;
			else
				data->input[ i ][ u ] = 0;
		}

		v++;
	}


    }


    line++;
    
    data->output[ i ] = data_output;
    data_output += num_output;
    

    unsigned clase = datos->GetDatClass( i );

    if (clase < 0) {
      fprintf(stderr, "Error al construir datos libreria fann. Clase no registrada.\n");
      fann_destroy_train(data);
      return NULL;
    }

    if (false) {
	data->output[ i ][ 0 ] = clase == 0 ? 1 : -1;
    } else {	    
    for (unsigned k = 0 ; k != num_output ; k++)
      if (k == clase)
        data->output[ i ][ k ] = 1;
      else
        data->output[ i ][ k ] = 0;
    }
  }
  return data;
}




/**
 *  Metodo que devuelve la clase de un dato
 *
 *  data : datos con los que calcular el error de clasificacion.
 */

std::vector<double> NNet::UnnormalizedDistribution(int index)
{
  fann_type *pred;
  std::vector<double> ret; 

  // Calculamos la salida de la red
  pred = fann_run(a_nn, a_patrones->input[ index ]);

  //double tot = 0.0;
  for (unsigned j = 0 ; j < a_patrones->num_output ; j++) {
    //tot += pred[j];
    ret.push_back(pred[j]); 
  }
  //for (unsigned j = 0 ; j < a_patrones->num_output ; j++) {
    //ret[j] /= tot;
  //} 

  return ret;
}



/**
 * Funcion que copia una parte de los datos 
 * dada por inicio y fin.
 */


struct fann_train_data * NNet::extraerSubConjuntoFANN(struct fann_train_data *data, int inicio, int fin)
{
  struct fann_train_data * dest;
  unsigned int x, x2;
  int len = fin - inicio + 1;

  if ( (dest = (struct fann_train_data *)malloc(sizeof(struct fann_train_data))) == NULL ) {
    fprintf(stderr,"Error no hay memoria para crear estructura de datos fann.\n");
    return NULL;
  }

  //fann_init_error_data((struct fann_error *)dest);

  struct fann_error *errdat = (struct fann_error *)data;
  errdat->errstr = NULL;
  errdat->errno_f = 0;
  errdat->error_log = stderr;

  dest->num_data = len;
  dest->num_input = data->num_input;
  dest->num_output = data->num_output;

  if ( ((dest->input  = (fann_type **)calloc(len, sizeof(fann_type *))) == NULL) ||
       ((dest->output = (fann_type **)calloc(len, sizeof(fann_type *))) == NULL) ) {
    fprintf(stderr,"Error no hay memoria para crear estructura de datos fann.\n");
    fann_destroy_train(dest);
    return NULL;
  }


  fann_type *data_input, *data_output;

  if ( ((data_input  = (fann_type *)calloc(dest->num_input * len,  sizeof(fann_type))) == NULL) ||
       ((data_output = (fann_type *)calloc(dest->num_output * len, sizeof(fann_type))) == NULL) ) {
    fprintf(stderr,"Error no hay memoria para crear estructura de datos fann.\n");
    fann_destroy_train(dest);
    return NULL;
  }

  for (x = inicio, x2 = 0 ; x <= (unsigned)fin ; x++, x2++ ) {

    dest->input[x2] = &(data_input[ x2 * dest->num_input ]);
    dest->output[x2] = &(data_output[ x2 * dest->num_output ]);

    memcpy(dest->input[x2],  data->input[x],  dest->num_input  * sizeof(fann_type));
    memcpy(dest->output[x2], data->output[x], dest->num_output * sizeof(fann_type));
  }
  return dest;
}




 /*
   int NNet::Classify(int index) 
   {
   fann_type *pred;
   double max = 10000;

   // Calculamos la salida de la red

   pred = fann_run(a_nn, a_patrones->input[ index ]);

   unsigned j, pk = 0;

   for (j = 0 ; j < a_patrones->num_output ; j++) {
   if (fabs(pred[ j ] - 1) < max) {
   pk = j;
   max = fabs(pred[ j ] - 1);
   }
   }

   return pk;
   }*/
