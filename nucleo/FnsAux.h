//---------------------------------------------------------------------------

#ifndef FnsAuxH
#define FnsAuxH
#include "Lector.h"
//---------------------------------------------------------------------------
void DarDeAltaFnsAux();
//---------------------------------------------------------------------------
std::string &replace_all(std::string s_find, std::string s_rpl, std::string &s);
//---------------------------------------------------------------------------
class Matriz;
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class DameNombreFichero : public IntefazFuncion, Funcion
{
  protected:
    std::string txt;

  public:
    DameNombreFichero(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    DameNombreFichero() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 1;}
    virtual Funcion* CrearFuncion() {return new DameNombreFichero();}
    virtual std::string NombreFuncion() {return "DameNombreFichero";}

    virtual Mandato* Ejecutar();

    virtual bool EsCadena(){return true;}
    virtual std::string ComoCadena() {return txt;}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class NombreProc : public IntefazFuncion, Funcion
{
  public:
    NombreProc(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    NombreProc() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 1;}
    virtual Funcion* CrearFuncion() {return new NombreProc();}
    virtual std::string NombreFuncion() {return "NombreProc";}

    virtual Mandato* Ejecutar();

    virtual bool EsCadena(){return true;}
    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class Append : public IntefazFuncion, Funcion
{
//   std::vector<Mandato*> Parametros;
  public:
    Append(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    Append() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new Append();}
    virtual std::string NombreFuncion() {return "Append";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return Param(0)->ComoCadena();}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class Mean : public IntefazFuncion, Funcion
{
  protected:
   Matriz *means;

  public:
    Mean(std::vector<Mandato *> mands) : Funcion(mands) {
      CreateParams();
      means=0;
    }
    Mean() : Funcion() { CreateParams(); means=0; }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 1;}
    virtual Funcion* CrearFuncion() {return new Mean();}
    virtual std::string NombreFuncion() {return "Mean";}

    virtual Mandato* Ejecutar();

    virtual bool EsDatos(){return true;}
    virtual Tipo* DameTipo(){return Tipo::TMatriz();}

    virtual std::string ComoCadena();
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return (long int)means;}
    virtual bool ComoBooleano(){return means;}
    virtual void* ComoDatos(){return means;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class StartRand : public IntefazFuncion, Funcion
{
  protected:
   std::vector<double> means;

  public:
    StartRand(std::vector<Mandato *> mands) : Funcion(mands) { CreateParams(); }
    StartRand() : Funcion() { CreateParams(); }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new StartRand();}
    virtual std::string NombreFuncion() {return "StartRand";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena(){return "";};
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
//    static boolean DadaDeAlta;
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class Print : public IntefazFuncion, Funcion
{
  protected:
    std::string txt;

  public:
    Print(std::vector<Mandato *> mands) : Funcion(mands) {CreateParams();  }
    Print() : Funcion() {CreateParams();  }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 1;}
    virtual Funcion* CrearFuncion() {return new Print();}
    virtual std::string NombreFuncion() {return "Print";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return txt;}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class Call : public IntefazFuncion, Funcion
{
  protected:
    std::vector<Mandato *> script;
    std::vector<std::string> txt;

    void LoadScript();
    bool ScriptLoaded();

  public:
    Call(std::vector<Mandato *> mands) : Funcion(mands) {CreateParams();  }
    Call() : Funcion() {CreateParams();  }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 1;}
    virtual Funcion* CrearFuncion() {return new Call();}
    virtual std::string NombreFuncion() {return "Call";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class Time : public IntefazFuncion, Funcion
{
  double secs;

  public:
    Time(std::vector<Mandato *> mands) : Funcion(mands) {CreateParams();  }
    Time() : Funcion() {CreateParams();  }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 2;}
    virtual Funcion* CrearFuncion() {return new Time();}
    virtual std::string NombreFuncion() {return "Time";}

    virtual bool EsNumero() {return true;}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return ACad(secs);}
    virtual double ComoNumero(){return secs;}
    virtual int ComoEntero(){return (int)(0.5 + secs);}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class Repeat : public IntefazFuncion, Funcion
{
  public:
    Repeat(std::vector<Mandato *> mands) : Funcion(mands) {CreateParams();  }
    Repeat() : Funcion() {CreateParams();  }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new Repeat();}
    virtual std::string NombreFuncion() {return "Repeat";}

    virtual Mandato* Ejecutar();
    virtual std::string ComoCadena() {return "";}
    virtual double ComoNumero(){return 0;}
    virtual int ComoEntero(){return 0;}
    virtual bool ComoBooleano(){return false;}
    virtual void* ComoDatos(){return 0;}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class If : public IntefazFuncion, Funcion
{
  protected:
    Mandato *mand;

  public:
    If(std::vector<Mandato *> mands) : Funcion(mands) {CreateParams();  }
    If() : Funcion() {CreateParams();  }

  public:
    void CreateParams();
    virtual unsigned NumParams(){return 3;}
    virtual Funcion* CrearFuncion() {return new If();}
    virtual std::string NombreFuncion() {return "If";}

    virtual Mandato* Ejecutar();

    virtual std::string ComoCadena() {return mand->ComoCadena();}
    virtual double ComoNumero(){return mand->ComoNumero();}
    virtual int ComoEntero(){return mand->ComoEntero();}
    virtual bool ComoBooleano(){return mand->ComoBooleano();}
    virtual void* ComoDatos(){return mand->ComoDatos();}

    static void DarDeAlta();
};
//---------------------------------------------------------------------------
#endif

