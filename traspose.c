#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <vector>
#include <string>
#include <string.h>

using namespace std;

vector< vector<string>* >* leermatriz(FILE *f);

int main(int argc,char**argv)
{
  FILE *f;
  extern char *optarg;
  char *filename=0, o;

  /*  Process options  */
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
 
  leermatriz(f);  

  return 0;    
}

vector< vector<string>* >* leermatriz(FILE *f)
{
  if (!f) return 0;

  char leido[1<<20];
  char *next;
  unsigned maxlon, fils;
  vector< vector<string>* >*m = new vector< vector<string>* >;


  fils = 0;
  maxlon = 0;
  while (!feof(f)) {
    bool nosevayantodaviaaunhaymas = false;
    fgets(leido, sizeof(leido), f); 
    if (strlen(leido)==sizeof(leido)-1 && leido[sizeof(leido)-2]!='\n')
      nosevayantodaviaaunhaymas = true;
    if (feof(f)) break;
    m->push_back(new vector<string>);
    next = strtok(leido, " \t\r\n");
    while (next) {
      (*m)[fils]->push_back((string)next);
      if ((*m)[fils]->size()>maxlon) maxlon = (*m)[fils]->size();
      next = strtok(0, " \t\r\n");
    };
    if (!nosevayantodaviaaunhaymas) fils++;
  };

  for(unsigned i=0;i<maxlon;i++) {
    for(unsigned j=0;j<fils;j++) {
      if ((*m)[j]->size()>i) printf("%s\t", (*(*m)[j])[i].c_str());
      else printf("\t");
    }
    printf("\n");
  }

  return m;
}
