#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

int main(int argc,char**argv)
{
  FILE *f;
  extern char *optarg;
  char *filename=0, o;
  char leido[0x100000];

  /*  Process options */ 
  while ( (o = getopt(argc, argv, "f:")) != EOF ) {
    switch (o) {
      case 'f':   
        filename = optarg;
        break;
      case '?':   printf("unrecognised option\n");
        exit(1);
    }
  }

  f= filename ? fopen(filename, "r") : stdin;
  if (!f) {
    printf("Error openening file: %s\n", filename);
    exit(1);
  }

  bool evenline = false;
  while(fgets(leido, 0x100000, f)) {
    if (evenline) printf("%s", leido);
    evenline = !evenline;
  }

  fclose(f);

  return 0;
}

