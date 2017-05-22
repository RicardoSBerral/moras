#include<stdio.h>
#include<string.h>

int main(int argn, char**argv)
{
  char cadena[ 100024 ];
  int i = 0;
  FILE *f = stdin;

  while ( !feof(f) ) {
    fgets(cadena, 100024, f);
    if (feof(f)) break;

    if ( argn == 1 ) { 
       i = strlen(cadena);
      while(cadena[ i ] != '\t' && cadena[ i ] != ' ' && i >= 0) i--;

      printf("%s", cadena + i + 1);
    }
    else {
      int nt = 0;
      char *tok = strtok(cadena, "\t ");
      while (tok != NULL) {
        nt++;
        tok = strtok(0, "\t ");
      }
      printf("%d\n", nt);
    }
  }

  return 0;
}
