//---------------------------------------------------------------------------

#ifndef Util20H
#define Util20H
//---------------------------------------------------------------------------
#include <string>
#include <stdio.h>

#ifndef WIN32
        #include <termios.h>
        #include <unistd.h> // read()
#else
        #include <conio.h> // kbhit()
#endif

//---------------------------------------------------------------------------
typedef enum ModoObtencionNumerosAleat {moNormal, moDeRandYGuardaEnFichero, 
                                        moDeFichero} TModoObtencionNumerosAleat;
//---------------------------------------------------------------------------
class TDebugRand
{
    static std::string NombreFichero;
  public:
    static int Peticiones;
    static TModoObtencionNumerosAleat mona;
    static FILE *f;
    static void ResetFile();
    static void Start();
    static void End();
    static int Rand();
    static void SetFichero(std::string algo);
    static void SetModo(TModoObtencionNumerosAleat modo);
};

class keyboard
{
public:

      keyboard();
    ~keyboard();
    int kbhit();
    int getch();

private:

#ifndef WIN32
    struct termios initial_settings, new_settings;
#endif
    int peek_character;

};
#endif
