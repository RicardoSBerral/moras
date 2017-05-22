#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <getopt.h>

int main(int argc,char**argv)
{
  FILE *f;
  extern char *optarg;
  char *format=0, *tok, o;
  char leido[0x100001];
  double val;
  bool mult = false;
  bool remove_text = false;

  /*  Process options */ 
  while ( (o = getopt(argc, argv, "rf:x:")) != EOF ) {
    switch (o) {
      case 'f':   
        format = optarg;
        break;
      case 'x':
        char *fin;
        val = strtod(optarg, &fin);
        if (fin==optarg) exit(1);
        mult = true;
        break;
      case 'r':
        remove_text = true;
        break;
      case '?':   printf("unrecognised option\n");
        exit(1);
    }
  }


  f = stdin;
  if (!f) {
    exit(1);
  }

  leido[0x10000] = '\0';
  while(fgets(leido, 0x100000, f)) {
    char *cad = leido;
    char *end;
    double num;
    while (*cad != '\0') {
      while (isspace(*cad)) {
        if ( !remove_text ) printf("%c", *cad);
        if ( remove_text && *cad == '\n' ) printf("%c", *cad);
        cad++;
      }
      num = strtod(cad, &end);
      if (end != cad) { /* Algo ha leido */
        printf(format, mult ? num*val : num);
        cad = end;
      }
      else if (*cad != '\0') {
        if ( !remove_text ) printf("%c", *cad);
        if ( remove_text && *cad == '\n' ) printf("%c", *cad);
        cad++;
      }
    }
  }

  fclose(f);

  return 0;
}

