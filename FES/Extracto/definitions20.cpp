#include "definitions20.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int strcmpci(const char *s1, const char *s2)
{
  static int S1cap = 256;
  static char *S1 = (char*)malloc(S1cap + 1);
  static int S2cap = 256;
  static char *S2 = (char*)malloc(S2cap + 1);

  int s1lon = strlen(s1);
  if (s1lon>S1cap) {
    S1cap = s1lon;
    S1 = (char*)realloc(S1, S1cap+1);
  }
  int s2lon = strlen(s2);
  if (s2lon>S2cap) {
    S2cap = s2lon;
    S2 = (char*)realloc(S2, S2cap+1);
  }

  //Pasamoms las cadenas a mayusculas
  for(int i=0;i<s1lon;i++)
    S1[i] = toupper(s1[i]);
  S1[s1lon] = '\0';
  for(int i=0;i<s2lon;i++)
    S2[i] = toupper(s2[i]);
  S2[s2lon] = '\0';

  return strcmp(S1, S2);
}
