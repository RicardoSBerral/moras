#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

int main(int argc,char**argv)
{
  FILE *f;
  extern char *optarg;
  char *filename=0, o;
  char leido[0x100000];
  bool sec[1024];
  int n_sec = 0;
  int cada = 1;
  int prim = 0;
  int ulti = -1;
  bool inverse = false;

  /*  Process options  */
  while ( (o = getopt(argc, argv, "f:c:p:u:s:i")) != EOF ) {
    switch (o) {
      case 'f':
        filename = optarg;
        break;
      case 'c':
        cada = atoi(optarg);
        break;
      case 's':
        n_sec = strlen(optarg);
        for(int i=0;i<n_sec;i++) {
          sec[i] = (optarg[i] == '1');
        }
        break;
      case 'p':
        prim = atoi(optarg);
        break;
      case 'i':
        inverse = true;
        break;
      case 'u':
        ulti = atoi(optarg);
        break;
      case '?':
        printf("unrecognised option\n");
        exit(1);
    }
  }

  f= filename ? fopen(filename, "r") : stdin;
  if (!f) {
    printf("Error openening file: %s\n", filename);
    exit(1);
  }

  int i=-1, k=0, i_sec=0;
  while(fgets(leido, 0x100000, f)) {
    i++;
    if (i<prim) continue;
    if (i>ulti && ulti>0) continue;
    if (inverse) {
      if (k%cada!=0) {
        if (n_sec==0) printf("%s", leido);
        else {
          if (sec[i_sec]) printf("%s", leido);
          i_sec++;
          if (i_sec==n_sec) i_sec=0;
        }
      }
    }
    else {
      if (k%cada==0) {
        if (n_sec==0) printf("%s", leido);
        else {
          if (sec[i_sec]) printf("%s", leido);
          i_sec++;
          if (i_sec==n_sec) i_sec=0;
        }
      }
    }
    k++;
  }

  return 0; 
}
