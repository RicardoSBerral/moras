//---------------------------------------------------------------------------

#ifndef LectorH
#define LectorH
#include <string>
#include <vector>
#include <map>
#include <stdlib.h>
#include <climits>
//---------------------------------------------------------------------------
typedef void ( *EventoSalida)(char *);
//typedef enum TTipo { tpDesconocido, tpCadena, tpNumero, tpEntero, tpBooleano, tpDatos } Tipo;
//---------------------------------------------------------------------------
class Tipo;
class Variable;
class Matriz;
//---------------------------------------------------------------------------
//-------------------------------------------------------  Excepciones  ---
//---------------------------------------------------------------------------
class MandatoError
{
  std::string msj_error;

  public:
    MandatoError(std::string err) : msj_error(err){;}

  public:
    std::string que(){return msj_error;}
};
class ParseError : public MandatoError
{
  int en;

  public:
    ParseError(std::string err, int pos) : MandatoError(err), en(pos){;}

  public:
    int donde(){return en;}
};
class ParametroError : public MandatoError
{
  public:
    ParametroError(std::string err): MandatoError(err){;}
};
class TipoIncorrecto : public ParametroError
{
  public:
    TipoIncorrecto(std::string err): ParametroError(err){;}
};
class ParametroNoOpcional : public ParametroError
{
  public:
    ParametroNoOpcional(std::string err): ParametroError(err){;}
};
class RangoNoValido : public ParametroError
{
  public:
    RangoNoValido(std::string err): ParametroError(err){;}
};
class NombreFicheroNoValido : public ParametroError
{
  public:
    NombreFicheroNoValido(std::string err): ParametroError(err){;}
};
//---------------------------------------------------------------------------
// Clase principal del parser
//---------------------------------------------------------------------------
class Mandato
{
  protected:
    void salida(char* sal);

  public:
    static std::string VERDADERO;
    static std::string FALSO;
    static int NivelSalida;
    static EventoSalida Salida;
    static std::vector<Mandato*> Mandatos;

  public:
    Mandato();
    virtual ~Mandato();

    virtual std::string ComoCadena()=0;
    virtual double ComoNumero()=0;
    virtual int ComoEntero()=0;
    virtual bool ComoBooleano()=0;
    virtual void* ComoDatos()=0;
//    virtual std::vector<double> ComoVectorNumeros(){;}
  //  virtual std::vector<std::vector<double>*> ComoMatrizNumeros(){;}

    virtual Tipo* DameTipo();
    virtual bool EsTipo(Tipo *tipo); 

    virtual bool EsCadena() {return false;}
    virtual bool EsNumero() {return false;}
    virtual bool EsEntero() {return false;}
    virtual bool EsBooleano() {return false;}
    virtual bool EsDatos() {return false;}
    virtual bool EsVectorNumeros() {return false;}
    virtual bool EsMatrizNumeros() {return false;}

    virtual Mandato* Ejecutar()=0;

    static int BuscaFueraDeParentesis(std::string exp, std::string cod);
    static int BuscaFueraDeParentesisInv(std::string exp, std::string cod);
    static void TrimExtendida(std::string *exp);
    static void Trim(std::string *exp);
    static Mandato* valueOf(std::string exp);
    static std::string ACad(double num);
    static std::string ACad(int num);

    static void EliminarMandatos();
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class Tipo
{
  private:
    static Tipo* _TCadena;
    static Tipo* _TNumerico;
    static Tipo* _TEntero;
    static Tipo* _TBooleano;
    static Tipo* _TDatos;
    static Tipo* _TMatriz;
    static Tipo* _TNinguno;
    static Tipo* _TCualquiera;
    static Tipo* _TFicheroEntrada;
    static Tipo* _TFicheroSalida;

  public:
    static Tipo* TCadena()    {return _TCadena;}
    static Tipo* TNumerico()  {return _TNumerico;}
    static Tipo* TEntero()    {return _TEntero;}
    static Tipo* TBooleano()   {return _TBooleano;}
    static Tipo* TDatos()       {return _TDatos;}
    static Tipo* TMatriz()       {return _TMatriz;}
    static Tipo* TNinguno()       {return _TNinguno;}
    static Tipo* TCualquiera()     {return _TCualquiera;}
    static Tipo* TFicheroEntrada() {return _TFicheroEntrada;}
    static Tipo* TFicheroSalida()  {return _TFicheroSalida;}

  protected:
    virtual bool ConvertibleA(Tipo *tipo)=0;

  public:
    virtual ~Tipo(){}

    virtual bool EsTipo(Tipo *tipo);

    static bool EsTipoBasico(Tipo *tipo);

    virtual std::string nombre()=0;

    virtual Variable *NuevaVariable(std::string nombre, Mandato *inic=0)=0;
};
class TipoCadena : public Tipo
{
  protected:
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoCadena*>(tipo);
    }

  public:
    virtual Variable *NuevaVariable(std::string nombre, Mandato *inic=0);
    virtual std::string nombre(){return "cadena";}

};
class TipoCadenaFicheroLectura : public TipoCadena
{
  protected:
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoCadenaFicheroLectura*>(tipo);
    }

  public:
    virtual std::string nombre(){return "file_in";}
};
class TipoCadenaFicheroEscritura : public TipoCadena
{
  protected:
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoCadenaFicheroEscritura*>(tipo);
    }

  public:
    virtual std::string nombre(){return "file_out";}
};
class TipoNumerico : public Tipo
{
  protected: 
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoNumerico*>(tipo);
    }

  public:
    virtual Variable *NuevaVariable(std::string nombre, Mandato *inic=0);
    virtual std::string nombre(){return "numero";}

};
class TipoEntero : public TipoNumerico
{
  protected:
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoNumerico*>(tipo);
    }

  public:
    virtual Variable *NuevaVariable(std::string nombre, Mandato *inic=0);
    virtual std::string nombre(){return "entero";}
};
class TipoBooleano : public Tipo
{
  protected: 
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoBooleano*>(tipo);
    }

  public:
    virtual Variable *NuevaVariable(std::string nombre, Mandato *inic=0);
    virtual std::string nombre(){return "booleano";}
    
};
class TipoDatos : public Tipo
{
  protected: 
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoDatos*>(tipo);
    }

  public:
    virtual Variable *NuevaVariable(std::string nombre, Mandato *inic=0);
    virtual std::string nombre(){return "datos";}
    
};
class TipoMatriz : public TipoDatos
{
  protected:
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoMatriz*>(tipo);
    }

  public:
    virtual Variable *NuevaVariable(std::string nombre, Mandato *inic=0);
    virtual std::string nombre(){return "matriz";}

};
/*class TipoCompuesto : public Tipo
{
  protected:
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoCompuesto*>(tipo);
    }

  public:
    virtual Variable *NuevaVariable(std::string nombre, Mandato *inic=0);
    std::string nombre(){return "compuesto";}

};*/
class TipoNinguno : public Tipo
{
  protected: 
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<TipoNinguno*>(tipo);
    }

  public:
    virtual Variable *NuevaVariable(std::string nombre, Mandato *inic=0){
      return 0;
    }
    virtual std::string nombre(){return "ninguno";}
    
};
class TipoCualquiera : public Tipo
{
  protected: 
    virtual bool ConvertibleA(Tipo *tipo){
      return dynamic_cast<Tipo*>(tipo);
    }

  public:
    virtual Variable *NuevaVariable(std::string nombre, Mandato *inic=0){
      return 0;
    }
    virtual std::string nombre(){return "cualquiera";}
    
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class Variable : public Mandato
{
  private:
    std::string FNombre;
    std::string FEspacioDeNombres;

  public:
    Variable(std::string nombre);
    Variable();
    virtual ~Variable();

  public:
    virtual void Set(Mandato *)=0;

    virtual void SetCadena(std::string)=0;
    virtual void SetNumero(double)=0;
    virtual void SetEntero(int)=0;
    virtual void SetBooleano(bool)=0;
    virtual void SetDatos(void*)=0;

    std::string Nombre() {return NombreCorto();}
    std::string NombreCompleto() {
      return UneNombre(FNombre, FEspacioDeNombres);
    }
    std::string NombreCorto() {return FNombre;}
    std::string EspacioDeNombres() {return FEspacioDeNombres;}

    Mandato* Ejecutar() {return this;}

  //Parte estatica
  private:
    static int VarsAnon;
    static std::map<std::string, Variable*> Variables;
    static std::map<std::string, std::vector<Variable *> * > esps;

    static bool NombreValido0(std::string nom);
    static std::vector<Variable *> * edn(std::string Namespace);

  public:
    static Variable* VariablePorNombre(std::string nombre,
                  std::string Namespace="", bool force=false, Mandato *inic=0);
    static int NumVariables(std::string Namespace="");
    static Variable* VariablePorNumero(int i, std::string Namespace="");
    static void EliminarVarDeListas(std::string nombre);
    static void EliminarVars(std::string Namespace="");
    static bool NombreValido(std::string nom, std::string Namespace="");
    static std::vector<Variable*> VariablesPorTipo(Tipo *tipo,
                                                     std::string Namespace="");
    static std::string UneNombre(const std::string &Nombre, 
                                          const std::string &EspacioDeNombres);
    static void ParteNombre(const std::string NombreLargo, 
                      std::string &EspacioDeNombres, std::string &NombreCorto);
    static void print();
};
//---------------------------------------------------------------------------
class IVariableDesglosable
{
  public:
    virtual ~IVariableDesglosable() {}

  public:
    virtual int NumElementos() = 0;
    virtual Variable* Elemento(int i) = 0;
};
//---------------------------------------------------------------------------
class VariableMutante : public Variable
{
  protected:
    Variable *FVar;

  public:
    VariableMutante(std::string nombre, Variable *valor) : Variable(nombre) 
      { FVar = valor; }
    VariableMutante(std::string nombre) : Variable(nombre) { FVar = 0; }
    VariableMutante() : Variable() { FVar = 0; }

    Variable *GetVar(){ return FVar;}
    void SetVar(Variable *v){ FVar = v;}

    virtual Tipo* DameTipo(); 
    virtual bool EsCadena();
    virtual bool EsNumero();
    virtual bool EsEntero();
    virtual bool EsBooleano();
    virtual bool EsDatos();
    virtual bool EsVectorNumeros();
    virtual bool EsMatrizNumeros();

    virtual std::string ComoCadena(){return FVar->ComoCadena();}
    virtual double ComoNumero()     {return FVar->ComoNumero();}
    virtual int ComoEntero()        {return FVar->ComoEntero();}
    virtual bool ComoBooleano()     {return FVar->ComoBooleano();}
    virtual void* ComoDatos()       {return FVar->ComoDatos();}

    virtual void Set(Mandato *m){FVar->Set(m);}

    virtual void SetCadena(std::string v){FVar->SetCadena(v);}
    virtual void SetNumero(double v)     {FVar->SetNumero(v);}
    virtual void SetEntero(int v)        {FVar->SetEntero(v);}
    virtual void SetBooleano(bool v)     {FVar->SetBooleano(v);}
    virtual void SetDatos(void* v)       {FVar->SetDatos(v);}

    Mandato* Ejecutar() {if(FVar) FVar->Ejecutar(); return this;}
};
//---------------------------------------------------------------------------
class VariableCadena : public Variable
{
  std::string FValor;
  public:
    VariableCadena(std::string nombre, std::string valor) : Variable(nombre) { FValor = valor; }
    VariableCadena(std::string nombre) : Variable(nombre) { FValor = ""; }
    VariableCadena() : Variable() { FValor = ""; }

    virtual Tipo* DameTipo(){return Tipo::TCadena();} 

    virtual bool EsCadena() {return true;}
    virtual std::string ComoCadena();
    virtual double ComoNumero();
    virtual int ComoEntero();
    virtual bool ComoBooleano();
    virtual void* ComoDatos();

    virtual void Set(Mandato *m){SetCadena(m->ComoCadena());}

    virtual void SetCadena(std::string);
    virtual void SetNumero(double);
    virtual void SetEntero(int);
    virtual void SetBooleano(bool);
    virtual void SetDatos(void*);
};
//---------------------------------------------------------------------------
class VariableNumero : public Variable
{
  protected:
    double FValor;

  public:
    VariableNumero(std::string nombre, double valor) : Variable(nombre) { FValor = valor; }
    VariableNumero(double valor) : Variable() { FValor = valor; }
    VariableNumero() : Variable() { FValor = 0.0; }

    virtual Tipo* DameTipo(){return Tipo::TNumerico();} 

    virtual bool EsNumero() {return true;}
    virtual std::string ComoCadena(){ return ACad(FValor);}
    virtual double ComoNumero(){ return FValor; }
    virtual int ComoEntero(){ return (int)FValor; }
    virtual bool ComoBooleano() {return (0!=ComoEntero()); }
    virtual void* ComoDatos() { return this; }

    virtual void Set(Mandato *m){SetNumero(m->ComoNumero());}

    virtual void SetCadena(std::string cad){FValor = atof(cad.c_str());}
    virtual void SetNumero(double valor) {FValor = valor;}
    virtual void SetEntero(int valor) {FValor = valor;}
    virtual void SetBooleano(bool valor) {FValor = valor;}
    virtual void SetDatos(void* valor) {}
};
//---------------------------------------------------------------------------
class VariableEntero : public VariableNumero
{
  public:
    VariableEntero(std::string nombre, int valor) : 
                                            VariableNumero(nombre, valor) { }
    VariableEntero(int valor) : VariableNumero(valor) {  }
    VariableEntero() : VariableNumero(0.0) { }

    virtual Tipo* DameTipo(){return Tipo::TEntero();} 

    virtual bool EsEntero() {return true;}

    virtual void Set(Mandato *m){SetEntero(m->ComoEntero());}
};
//---------------------------------------------------------------------------
class VariableBooleano : public VariableNumero
{
  private:
    static VariableBooleano *GenVarTrue();
    static VariableBooleano *GenVarFalse();
  public:
    static VariableBooleano *VarTrue;
    static VariableBooleano *VarFalse;

  public:
    VariableBooleano(std::string nombre, bool valor) : 
                                            VariableNumero(nombre, valor) { }
    VariableBooleano(bool valor) : VariableNumero(valor) {  }
    VariableBooleano() : VariableNumero(0.0) {  }

    virtual Tipo* DameTipo(){return Tipo::TBooleano();} 

    virtual bool EsBooleano() {return true;}

    virtual void Set(Mandato *m){SetBooleano(m->ComoBooleano());}

    virtual std::string ComoCadena();
};
//---------------------------------------------------------------------------
class VariableNumeroHistorico : public VariableNumero
{
  protected:
    std::vector<double> FValores;
  public:
    VariableNumeroHistorico(std::string nombre, double valor) :
                   VariableNumero(nombre, valor) { FValores.push_back(FValor); }
    VariableNumeroHistorico(double valor) : VariableNumero(valor)
                                                 { FValores.push_back(FValor); }
    VariableNumeroHistorico() : VariableNumero() { FValores.push_back(FValor); }

    virtual void SetCadena(std::string var) {
      VariableNumero::SetCadena(var);
    }
    virtual std::string ComoCadena() {
    std::string cad = "";
      for (unsigned i=0;i<FValores.size();i++) {
        cad += ACad(FValores[i]) + "\t";
      }

      return cad;
    }
    virtual double ComoNumero() {
      return VariableNumero::ComoNumero();
    }
    virtual void SetNumero(double valor) {
      VariableNumero::SetNumero(valor);
      FValores.push_back(FValor);
    }
    std::vector<double> DameHistoricos(){ return FValores; }
    void Inicializar() { FValores.clear(); }
};
//---------------------------------------------------------------------------
class VariableDatos : public Variable
{
  protected:
    void* FValor;
    Tipo *tipo;

  public:
    VariableDatos(std::string nombre, void* valor, Tipo *tipo=0) :
                        Variable(nombre) { FValor = valor; this->tipo = tipo; }
    VariableDatos(void* valor, Tipo *tipo=0) : Variable()
                                         { FValor = valor; this->tipo = tipo; }
    VariableDatos(Tipo *tipo=0) : Variable() { FValor = 0; this->tipo = tipo; }

    virtual Tipo* DameTipo(){return tipo ? tipo : Tipo::TDatos();}

    virtual bool EsDatos() {return true;}
    virtual std::string ComoCadena(){return (FValor) ? "DATOS" : "(VACIO)";}
    virtual double ComoNumero(){return (int)FValor;}
    virtual int ComoEntero(){return (int)FValor;}
    virtual bool ComoBooleano(){ return (FValor!=0); }
    virtual void* ComoDatos() {return FValor; }

    virtual void Set(Mandato *m){SetDatos(m->ComoDatos());}

    virtual void SetCadena(std::string valor) {FValor = (void*)valor.c_str();}
    virtual void SetNumero(double ) {}
    virtual void SetEntero(int) {}
    virtual void SetBooleano(bool){}
    virtual void SetDatos(void* valor) {FValor = valor;}
};
//---------------------------------------------------------------------------
class VariableMatriz : public VariableDatos
{

  public:
    VariableMatriz(std::string nombre, Matriz* valor=0) :
                                               VariableDatos(nombre, valor) { }
    VariableMatriz(Matriz* valor) : VariableDatos(valor) { }

    virtual Tipo* DameTipo(){return Tipo::TMatriz();}
    virtual bool EsDatos() {return true;}

    virtual std::string ComoCadena();
    virtual void SetCadena(std::string valor);

    Matriz *GetMatriz(){ return (Matriz*)ComoDatos(); }
};
//---------------------------------------------------------------------------
class MandatoAsignacion : public Mandato
{
  Mandato *FMandato;
  Variable *FVar;
  public:
    MandatoAsignacion(Mandato *mand, std::string NomVar);
    MandatoAsignacion(Mandato *mand, Variable *Var);

    virtual Tipo* DameTipo(){return FVar->DameTipo();} 

    virtual bool EsCadena() {return FVar->EsCadena();}
    virtual bool EsNumero() {return FVar->EsNumero();}
    virtual bool EsEntero() {return FVar->EsEntero();}
    virtual bool EsBooleano() {return FVar->EsBooleano();}
    virtual bool EsDatos() {return FVar->EsDatos();}

    virtual std::string ComoCadena() {return FVar->ComoCadena();}
    virtual double ComoNumero() {return FVar->ComoNumero();}
    virtual int ComoEntero() {return FVar->ComoEntero();}
    virtual bool ComoBooleano() {return FVar->ComoBooleano();}
    virtual void* ComoDatos() {return FVar->ComoDatos();}

    Mandato* Ejecutar();
};
//---------------------------------------------------------------------------
class Cadena : public Mandato
{
  protected:
  std::string cad;
  public:
    Cadena(std::string cad) : Mandato() { this->cad = cad; }
    Cadena(char* cad) : Mandato() { this->cad = cad; }
    Cadena() : Mandato() { this->cad = ""; }

    virtual Tipo* DameTipo(){return Tipo::TCadena();} 

    virtual bool EsCadena() {return true;}
    virtual std::string ComoCadena();
    virtual double ComoNumero();
    virtual int ComoEntero();
    virtual bool ComoBooleano();
    virtual void* ComoDatos();

    Mandato* Ejecutar() {return this;}
};
//---------------------------------------------------------------------------
class Concatenacion : public Cadena
{
  Mandato *sum1;
  Mandato *sum2;
  public:
    Concatenacion(Mandato *sum1, Mandato *sum2) : Cadena() {
      this->sum1 = sum1;
      this->sum2 = sum2;
    }

    virtual Mandato* Ejecutar() {
      sum1->Ejecutar();
      sum2->Ejecutar();
      cad = sum1->ComoCadena() + sum2->ComoCadena();
      return this;
    }

};
//---------------------------------------------------------------------------
class Numero : public Mandato
{
  protected:
  double num;
  bool ejec;
  public:
    Numero(double num) { this->num = num; ejec = true;}
    Numero(char* cad) { this->num = atof(cad); ejec = true;}
    Numero() {
      num = 0;
      ejec = false;
    }

    virtual Tipo* DameTipo(){return Tipo::TNumerico();} 

    virtual bool EsNumero() {return true;}
    virtual std::string ComoCadena();
    virtual double ComoNumero();
    virtual int ComoEntero();
    virtual bool ComoBooleano();
    virtual void* ComoDatos();

    Mandato* Ejecutar() {return this;}
};
//---------------------------------------------------------------------------
class Suma : public Numero
{
  Mandato *sum1;
  Mandato *sum2;
  public:
    Suma(Mandato *sum1, Mandato *sum2) : Numero() {
      this->sum1 = sum1;
      this->sum2 = sum2;
    }

    virtual Mandato* Ejecutar() {
      sum1->Ejecutar();
      sum2->Ejecutar();
      num = sum1->ComoNumero() + sum2->ComoNumero();
      ejec = true;
      return this;
    }

};
//---------------------------------------------------------------------------
class Resta : public Numero
{
  Mandato *sum1;
  Mandato *sum2;
  public:
    Resta(Mandato *sum1, Mandato *sum2) : Numero() {
      this->sum1 = sum1;
      this->sum2 = sum2;
    }

    virtual Mandato* Ejecutar() {
      sum1->Ejecutar();
      sum2->Ejecutar();
      num = sum1->ComoNumero() - sum2->ComoNumero();
      ejec = true;
      return this;
    }

};
//---------------------------------------------------------------------------
class Producto : public Numero
{
  Mandato *sum1;
  Mandato *sum2;
  public:
    Producto(Mandato *sum1, Mandato *sum2) : Numero() {
      this->sum1 = sum1;
      this->sum2 = sum2;
    }

    virtual Mandato* Ejecutar() {
      sum1->Ejecutar();
      sum2->Ejecutar();
      num = sum1->ComoNumero() * sum2->ComoNumero();
      ejec = true;
      return this;
    }

};
//---------------------------------------------------------------------------
class Division : public Numero
{
  Mandato *div1;
  Mandato *div2;
  public:
    Division(Mandato *div1, Mandato *div2) : Numero() {
      this->div1 = div1;
      this->div2 = div2;
    }

    virtual Mandato* Ejecutar() {
      div1->Ejecutar();
      div2->Ejecutar();
      num = div1->ComoNumero() / div2->ComoNumero();
      ejec = true;
      return this;
    }

};
//---------------------------------------------------------------------------
class Parametro : public Mandato
{
  protected:
    Mandato *mnd;
    bool opcional;
    Tipo *tipo;
    std::string nombre;
    std::string descripcion;
    std::string poromision;

  public:
    Parametro(bool o=false, Tipo *t=Tipo::TNinguno(), std::string n="",
                                                             std::string d="") {
      mnd = 0;
      PonPropiedades(o, t, n.size()==0 && t ? t->nombre() : n, d);
    }

  public:
    virtual bool Opcional() {return opcional;}
    virtual void checkValue(Mandato *mnd);

    virtual Tipo* DameTipo(){return tipo;}

    void SetMandato(Mandato *mnd){ this->mnd = mnd;}
    bool Asignado(){return mnd;}

    void PonNombre(std::string nombre) {this->nombre = nombre;}
    std::string Nombre() {return nombre;}
    void PonDescripcion(std::string descrip) {this->descripcion = descrip;}
    std::string Descripcion() {return descripcion;}

    virtual std::string ComoCadena() {return mnd->ComoCadena();}
    virtual double ComoNumero() {return mnd->ComoNumero();}
    virtual int ComoEntero() {return mnd->ComoEntero();}
    virtual bool ComoBooleano() {return mnd->ComoBooleano();}
    virtual void* ComoDatos() {return mnd->ComoDatos();}

    virtual Mandato* Ejecutar();

    void PonPropiedades(bool o=false, Tipo *t=Tipo::TNinguno(),
                                           std::string n="", std::string d="");

    std::string ValorPorOmision() {return poromision;}
    void SetPorOmision(std::string val) {poromision = val;}

};
class ParametroNumero : public Parametro
{
  protected:
    double menor;
    double mayor;

  public:
    ParametroNumero(double min=INT_MIN, double max=INT_MAX) : menor(min), mayor(max){ ; }

    virtual void checkValue(Mandato *mnd);

    virtual Tipo* DameTipo(){return Tipo::TNumerico();}

    bool Rango() {return mayor>=menor;}
    double Min(){return menor;}
    double Max(){return mayor;}
};
class ParametroEntero : public Parametro
{
  protected:
    int menor;
    int mayor;

  public:
    ParametroEntero(int min=1, int max=0) : menor(min), mayor(max){ ; }

    virtual void checkValue(Mandato *mnd);

    virtual Tipo* DameTipo(){return Tipo::TEntero();}

    bool Rango() {return mayor>=menor;}
    int Min(){return menor;}
    int Max(){return mayor;}
};
class ParametroCadena : public Parametro
{
  protected:
    std::vector<std::string> CadenasValidas;

  public:
    ParametroCadena(std::vector<std::string> *CadenasValidas=0){
      if (CadenasValidas)
        this->CadenasValidas = *CadenasValidas;
    }

    virtual void checkValue(Mandato *mnd);

    virtual Tipo* DameTipo(){return Tipo::TCadena();}

    virtual std::string ComoCadena();
    virtual int ComoEntero();

    int NumCadenasValidas(){ return CadenasValidas.size();}
    std::string CadenaValida(int i){ return CadenasValidas[i];}
};
class ParametroFicheroEntrada : public Parametro
{
  protected:
    std::string Filtros;

  public:
    ParametroFicheroEntrada(std::string filtros="") {Filtros = filtros;}

    virtual void checkValue(Mandato *mnd);
    virtual Tipo* DameTipo(){return Tipo::TFicheroEntrada();}

    std::string DameFiltros() {return Filtros;}
};
class ParametroFicheroSalida : public Parametro
{
  protected:
    std::string Filtros;

  public:
    ParametroFicheroSalida(std::string filtros=""){ Filtros = filtros;}

    virtual void checkValue(Mandato *mnd);
    virtual Tipo* DameTipo(){return Tipo::TFicheroSalida();}

    std::string DameFiltros() {return Filtros;}
};
class ParametroBooleano : public ParametroCadena
{
  public:
    ParametroBooleano(bool o=false, std::string n="", std::string d="")
                                                       : ParametroCadena() {
      this->CadenasValidas.push_back(Mandato::VERDADERO);
      this->CadenasValidas.push_back(Mandato::FALSO);
      PonPropiedades(o, 0, n, d);
    }
};
//---------------------------------------------------------------------------
class Funcion : public Mandato
{
  protected:
   std::vector<Mandato*> Parametros;

  public:
    Funcion(std::vector<Mandato *> params);
    Funcion();
    ~Funcion();

  //C++ does not allow you to call pure
  //virtual functions from the base constructor (while Java does) 
  virtual Tipo* TipoRetorno(){ return DameTipo();}
  virtual unsigned NumParams()=0;
  virtual Parametro* ParamInfo(int i);

  int NumParamsAsignados();
//  std::vector<Mandato *> *G*etParametros();
  //void Param(int i, Mandato *param){(*GetParametros())[i] = param;}
  Mandato* Param(int i){return ParamInfo(i);}
};
//---------------------------------------------------------------------------
class FuncionDatos : public Funcion
{
  protected:
   void* datos;
  public:
    FuncionDatos(std::vector<Mandato *> k):Funcion(k) { datos=0; }
    FuncionDatos():Funcion() { datos=0; }
    ~FuncionDatos() { }

    virtual bool EsDatos() {return true;}

    virtual double ComoNumero() {return 0;}
    virtual std::string ComoCadena() {return "";}
    virtual int ComoEntero()    {return (long int)datos;}
    virtual bool ComoBooleano() {return datos;}
    virtual void* ComoDatos()   {return datos;}
};
//---------------------------------------------------------------------------
class IntefazFuncion
{
  public:
    virtual ~IntefazFuncion() {}

  public:
    virtual Funcion* CrearFuncion()=0;
    virtual std::string NombreFuncion()=0;
};
class ListaInterfazFunciones
{
  protected:
    std::vector<IntefazFuncion*> iFunciones;
    std::vector<std::string> FunsCat;
    std::vector<Funcion*> InfoFunciones;
    std::vector<std::string> Categorias;

  public:
    ListaInterfazFunciones(){ Categorias.push_back("");}

  public:
    void AnadirInterfazFuncion(IntefazFuncion *inf, std::string cat="");
    Funcion *CrearFuncionDeCadena(std::string *cadena);
    void CrearParametrosDeCadena(std::string *cadena, Funcion *f);

    int NumFunciones(){return size();}
    int size(){return iFunciones.size();}
    IntefazFuncion *iFuncion(int i){return iFunciones[i];}
    Funcion *InfoFuncion(int i){return InfoFunciones[i];}
    Funcion *InfoFuncionPorNombre(std::string nom);

    std::string DameCategoria(int i);
    int NumCategorias(){return Categorias.size();}

    std::vector<int> IndicesFuncionesPorCat(std::string cat);
    std::vector<IntefazFuncion*> InterfacesPorCat(std::string cat);
    std::vector<IntefazFuncion*> InterfacesPorTipo(Tipo *tipo);

  private:
    static ListaInterfazFunciones _ListaInterfazFunciones;
  public:
    static ListaInterfazFunciones &lif(){
      return _ListaInterfazFunciones;
    }
};
ListaInterfazFunciones &ifns();
//---------------------------------------------------------------------------
class Seno : public IntefazFuncion, Funcion
{
 //  std::vector<Mandato*> Parametros;
   double res;
  public:
    Seno(std::vector<Mandato *> mands) : Funcion(mands) { res = 0.0; }
    Seno() : Funcion() { res = 0.0; }

  virtual unsigned NumParams(){return 1;}
  virtual Funcion* CrearFuncion() {return new Seno();}
  virtual std::string NombreFuncion() {return "sen";}

    virtual Mandato* Ejecutar();
    virtual bool EsNumero() {return true;}
    virtual std::string ComoCadena() {return ACad(res);}
    virtual double ComoNumero(){return res;}
    virtual int ComoEntero(){return (int)(res);}
    virtual bool ComoBooleano(){return res!=0.0;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
class Coseno : public IntefazFuncion, Funcion
{
  // std::vector<Mandato*> Parametros;
   double res;
  public:
    Coseno(std::vector<Mandato *> mands) : Funcion(mands) { res = 0.0; }
    Coseno() : Funcion() { res = 0.0; }

  virtual unsigned NumParams(){return 1;}
  virtual Funcion* CrearFuncion() {return new Coseno();}
  virtual std::string NombreFuncion() {return "cos";}

    virtual Mandato* Ejecutar();
    virtual bool EsNumero() {return true;}
    virtual std::string ComoCadena() {return ACad(res);}
    virtual double ComoNumero(){return res;}
    virtual int ComoEntero(){return (int)(res);}
    virtual bool ComoBooleano(){return res!=0.0;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};

#endif
