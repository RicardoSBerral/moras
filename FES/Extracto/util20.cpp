//---------------------------------------------------------------------------
#include "util20.h"
#include <stdlib.h>
#ifndef WIN32
        #include <unistd.h> // read()
#else
        #include <conio.h> // kbhit()
#endif

///---------------------------------------------------------------------------

std::string TDebugRand::NombreFichero = "temp.alea";
TModoObtencionNumerosAleat TDebugRand::mona = moDeRandYGuardaEnFichero;
FILE* TDebugRand::f = 0;
int TDebugRand::Peticiones = 0;

void TDebugRand::ResetFile()
{
  if (f) fclose(f);
  f = fopen(NombreFichero.c_str(),"w");
  fclose(f);
  f=0;
  Peticiones = 0;
}

void TDebugRand::Start()
{
  if (f) {fclose(f); f=0;}
  if (mona==moDeRandYGuardaEnFichero) {
    ResetFile();
    f = fopen(NombreFichero.c_str(),"wb");
  }
  else if (mona==moDeFichero) {
    f = fopen(NombreFichero.c_str(),"rb");
  }
}
void TDebugRand::End()
{
  if (f) fclose(f);
  f=0;
  Peticiones = 0;
}

int TDebugRand::Rand()
{
  int res;
  Peticiones++;
  if (!f || mona==moNormal) {
    res=rand();
  }
  else if (mona==moDeRandYGuardaEnFichero) {
    res=rand();
    fwrite(&res,1,sizeof(int),f);
  }
  else {
    fread(&res,1,sizeof(int),f);
  }
  return res;
}
void TDebugRand::SetFichero(std::string algo)
{
  NombreFichero = algo;
  if (f) {
    Start();
  }
}
void TDebugRand::SetModo(TModoObtencionNumerosAleat modo)
{
  mona = modo;
}
//----------------------------
keyboard::keyboard()
{
#ifndef WIN32
    tcgetattr(0,&initial_settings);
    new_settings = initial_settings;
    new_settings.c_lflag &= ~ICANON;
    new_settings.c_lflag &= ~ECHO;
    new_settings.c_lflag &= ~ISIG;
    new_settings.c_cc[VMIN] = 1;
    new_settings.c_cc[VTIME] = 0;
    tcsetattr(0, TCSANOW, &new_settings);
#endif
    peek_character=-1;
}

keyboard::~keyboard()
{
#ifndef WIN32
    tcsetattr(0, TCSANOW, &initial_settings);
#endif
}

int keyboard::kbhit()
{
#ifndef WIN32
    unsigned char ch;
    int nread;

    if (peek_character != -1) return 1;
    new_settings.c_cc[VMIN]=0;
    tcsetattr(0, TCSANOW, &new_settings);
    nread = read(0,&ch,1);
    new_settings.c_cc[VMIN]=1;
    tcsetattr(0, TCSANOW, &new_settings);

    if (nread == 1)
    {
        peek_character = ch;
        return 1;
    }
    return 0;
#else
    return ::kbhit();
#endif
}

int keyboard::getch()
{
#ifndef WIN32
char ch;

    if (peek_character != -1)
    {
        ch = peek_character;
        peek_character = -1;
    } else read(0,&ch,1);

    return ch;
#else
    return ::getch(); 
#endif
}

////////////////////////////////////////////////////////////////
